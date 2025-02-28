#include <Interpreters/SuccinctRangeFilter.h>
#include <city.h>
#include <Columns/ColumnArray.h>
#include <Columns/ColumnNullable.h>
#include <Columns/ColumnLowCardinality.h>
#include <DataTypes/DataTypeArray.h>
#include <DataTypes/DataTypeNullable.h>
#include <DataTypes/DataTypeLowCardinality.h>
#include <queue>

#include <Common/logger_useful.h>

namespace DB
{
namespace ErrorCodes
{
    extern const int BAD_ARGUMENTS;
}

// Builds the LOUDS-DS structure from a trie.
// void buildLOUDSDS(std::unique_ptr<TrieNode> & root, size_t l_depth) {
//     // Here, you would traverse the trie and fill in surf.dense and surf.sparse.
//     // Also, compute l_depth using the ratio R (default R=64) as described.
//     l_depth = l_depth;
//     // For illustration, we leave this unimplemented.
// }


// TODO: replace the loop with a specialized rank data structure
// Compute the child position in the dense part using a rank/select operation.
size_t SuccinctRangeFilter::childPositionDense(size_t currentPos, uint8_t label, size_t level) const
{
    // Ensure the level is valid.
    if (level >= surf.dense.d_hasChild.size())
    {
        throw std::out_of_range("childPositionDense: level out of range");
    }
    
    // Get the bitset representing the "hasChild" flags at this level.
    const std::bitset<256>& hasChild = surf.dense.d_hasChild[level];
    
    // Compute the rank: count the number of bits set (i.e. children present) in positions [0, label)
    size_t rank = 0;
    for (size_t i = 0; i < label; ++i)
    {
        if (hasChild.test(i))
            ++rank;
    }
    
    // Check if the target label itself exists as a child.
    if (!hasChild.test(label))
    {
        // In a well-formed trie, this function should only be called when the label exists.
        throw std::invalid_argument("childPositionDense: label not found in hasChild bitset");
    }
    
    // The new child node's position in the next level is computed by adding the rank to the current position.
    return currentPos + rank;
}


// For the sparse part, compute the node boundaries.
size_t SuccinctRangeFilter::getSparseNodeStart(size_t pos) const
{
    // If pos is 0 or out-of-bounds, return pos as it must be a node start.
    if (pos == 0 || pos >= surf.sparse.s_LOUDS.size())
        return pos;

    // Walk backward until a node boundary is found.
    // A node boundary is indicated by a true value in s_LOUDS.
    while (pos > 0 && !surf.sparse.s_LOUDS[pos])
    {
        pos--;
    }
    return pos;
}


size_t SuccinctRangeFilter::getSparseNodeEnd(size_t pos) const
{
    // If pos is out of bounds, return pos.
    if (pos >= surf.sparse.s_labels.size())
        return pos;
    
    // Starting from the next position, search for the next node boundary.
    size_t i = pos + 1;
    while (i < surf.sparse.s_labels.size() && !surf.sparse.s_LOUDS[i])
    {
        ++i;
    }
    // 'i' is now the index of the next node's start (or the end of s_labels).
    return i;
}

size_t SuccinctRangeFilter::rankTrue(const std::vector<bool>& vec, size_t pos) const {
    size_t count = 0;
    // Ensure we don't go past the vector size.
    for (size_t i = 0; i < pos && i < vec.size(); ++i) {
        if (vec[i]) {
            ++count;
        }
    }
    return count;
}

// Helper: find the index of the count-th false in 'vec'.
// count is 1-indexed; i.e. count == 1 returns the index of the first false.
size_t SuccinctRangeFilter::selectFalse(const std::vector<bool>& vec, size_t count) const {
    size_t zeroCount = 0;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (!vec[i]) {
            ++zeroCount;
            if (zeroCount == count)
                return i;
        }
    }
    // If not found, return vec.size() (could also throw an error in production code).
    return vec.size();
}

// TODO: replace the loop with a specialized rank data structure
// Compute the child position in the sparse part.
size_t SuccinctRangeFilter::childPositionSparse(size_t currentPos) const
{
    // Compute the number of true bits in s_hasChild from index 0 to currentPos (inclusive).
    // This gives the "node index" among those nodes that have children.
    size_t num = rankTrue(surf.sparse.s_hasChild, currentPos + 1);

    // In our LOUDS encoding for the sparse part, a false bit in s_LOUDS marks the end of a node.
    // The children of the current node start immediately after the (num+1)-th false.
    size_t boundaryIndex = selectFalse(surf.sparse.s_LOUDS, num + 1);

    // The child node's starting position in s_labels is one position after the boundary.
    return boundaryIndex + 1;
}

// Find the smallest label in the dense node at [currentPos, ...) that is >= target.
size_t SuccinctRangeFilter::findLowerBoundInDense(/*size_t currentPos, */size_t level, uint8_t target) const {
    const std::bitset<256>& labels = surf.dense.d_labels[level];
    // Search for the next set bit starting at 'target'.
    for (size_t b = target; b < 256; ++b) {
        if (labels.test(b))
            return b;  // b is valid (0..255)
    }
    // Return 256 as a sentinel value indicating no label was found.
    return 256;
}
// Check whether there is a next label in the node at the given level.
bool SuccinctRangeFilter::canMoveNext(const Iterator & iter, size_t level) const
{
    // Check for dense node.
    if (level < l_depth)
    {
        // In the dense representation, the current position is the label value.
        uint16_t currentLabel = iter.levelPositions[level];
        // Scan for the next set bit (label) in the bitset.
        for (uint16_t next = currentLabel + 1; next < 256; ++next)
        {
            if (surf.dense.d_labels[level].test(next))
                return true;
        }
        return false;
    }
    else
    {
        // For sparse nodes, iter.levelPositions[level] is an index into s_labels.
        size_t currentPos = iter.levelPositions[level];
        // Get the start and end boundaries of the current sparse node.
        // size_t nodeStart = getSparseNodeStart(currentPos);
        size_t nodeEnd = getSparseNodeEnd(currentPos);
        // If there's another label in the node after currentPos, we can move next.
        return (currentPos + 1 < nodeEnd);
    }
}


// Given a current position in a node, return the next position.
size_t SuccinctRangeFilter::nextPosition(size_t pos, size_t level) const
{
    // Dense levels: use the bitset representing labels at this level.
    if (level < l_depth)
    {
        const std::bitset<256>& labels = surf.dense.d_labels[level];
        // Search for the next set bit after 'pos'.
        for (size_t i = pos + 1; i < 256; ++i)
        {
            if (labels.test(i))
            {
                return i;
            }
        }
        // No next label found in this dense node: return an invalid index.
        return std::numeric_limits<size_t>::max();
    }
    else
    {
        // Sparse levels: assume that labels for a node are stored contiguously in s_labels.
        // The LOUDS bitvector (s_LOUDS) marks the beginning of a new node (true indicates a new node).
        size_t next = pos + 1;
        if (next < surf.sparse.s_labels.size() && !surf.sparse.s_LOUDS[next])
        {
            // Next position is within the current node.
            return next;
        }
        // No next label exists in the current node.
        return surf.sparse.s_labels.size(); // using vector size as an invalid marker.
    }
}

// Check whether the node at the given position and level has a child.
bool SuccinctRangeFilter::hasChild(size_t pos, size_t level) const
{
    if (level < l_depth) {
        // For dense nodes, assume pos maps directly to a label index.
        return surf.dense.d_hasChild[level].test(pos);
    }
    else {
        return surf.sparse.s_hasChild[pos];
    }
}

// Get the left-most child position for a given node.
size_t SuccinctRangeFilter::leftMostChild(size_t pos, size_t level) const
{
    // In a dense node, find the first label that is set.
    if (level < l_depth) {
        for (uint16_t b = 0; b < 256; b++) {
            if (surf.dense.d_labels[level].test(b))
                return b;  // simplified: returning label value as position
        }
    }
    else {
        // For sparse, assume the left-most child is the next element.
        return pos + 1;
    }
    return pos;
}

// SuccinctRangeFilter::SuccinctRangeFilter(size_t ds_ratio_)
// {
//     ds_ratio = ds_ratio_;
// }

// bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b)
// {
//     if (a.ds_ratio != b.ds_ratio)
//         return false;
//     else
//         return true;
    
// }

// void SuccinctRangeFilter::add(const char * data, size_t len)
// {
//     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "add {} {}", data, len);

// }

// bool SuccinctRangeFilter::find(const char * data, size_t len)
// {
//     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "find {} {}", data, len);
//     if (len > 3)
//     {
//         return true;
//     }
//     else
//     {
//         return false;
//     }

// }

size_t SuccinctRangeFilter::getWriteSize()
{
    size_t d_labels_size = surf.dense.d_labels.size() * sizeof(std::bitset<256>);
    size_t d_has_child_size = surf.dense.d_hasChild.size() * sizeof(std::bitset<256>);
    size_t d_is_prefix_tree_size = surf.dense.d_isPrefixKey.size() * sizeof(bool);
    size_t d_values_size = surf.dense.d_values.size() * sizeof(uint64_t);

    size_t s_labels_size = surf.sparse.s_labels.size() * sizeof(uint8_t);
    size_t s_has_child_size = surf.sparse.s_hasChild.size() * sizeof(bool);
    size_t s_louds_size = surf.sparse.s_LOUDS.size() * sizeof(bool);
    size_t s_values_size = surf.sparse.s_values.size() * sizeof(uint64_t);

    size_t louds_dense_size = d_labels_size + d_has_child_size + d_is_prefix_tree_size + d_values_size;
    size_t louds_sparse_size = s_labels_size + s_has_child_size + s_louds_size + s_values_size;
    
    return louds_dense_size + louds_sparse_size;
}

SuccinctRangeFilter::SuccinctRangeFilter(std::unique_ptr<TrieNode> root, size_t l_depth_) // TODO: change this so that it make full dense and sparse and then finds depth
{
    l_depth = l_depth_;
    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "SuccinctRangeFilter");
    // Clear existing data
    surf.dense = LOUDSDenseTrie();
    surf.sparse = LOUDSSparseTrie();

    if (!root)
    {
        return;
    }

    // ------------------------------------------------------------
    // PHASE 0: Construct d_values and s_values
    // ------------------------------------------------------------
    
    std::queue<KeyItem> key_queue;
    key_queue.push({ root.get(), std::vector<char>() });
    while (!key_queue.empty())
    {
        auto [node_ptr, key] = key_queue.front();
        key_queue.pop();

        if (!node_ptr->children.empty())
        {
            for (auto & kv : node_ptr->children)
            {
                std::vector<char> temp_key = key;
                char ch = kv.first;
                TrieNode * child_ptr = kv.second.get();
                temp_key.emplace_back(ch);
                key_queue.push({ child_ptr, temp_key });
            }
            if (node_ptr->is_terminal)
            {
                std::vector<char> temp_key = key;
                temp_key.emplace_back(char(0xFF));
                key_queue.push({ new TrieNode(), temp_key });
            }
        }
        else
        {
            // char * str_ptr = key.data();
            if (key.size() <= l_depth)
            {
                if (key.back() == char(0xFF))
                {
                    key.pop_back();
                }
                // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_value: {}", key);
                std::vector<char>* buffer_ptr = new std::vector<char>(key);
                // buffer_ptr = &key;
                // auto vec_ptr = std::make_unique<std::vector<char>>(key);
                surf.dense.d_values.emplace_back(reinterpret_cast<uint64_t>(buffer_ptr));
            }
            else
            {
                if (key.back() == char(0xFF))
                {
                    key.pop_back();
                }
                // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_value: {}", key);
                std::vector<char>* buffer_ptr = new std::vector<char>(key);
                // buffer_ptr = &key;
                surf.sparse.s_values.emplace_back(reinterpret_cast<uint64_t>(buffer_ptr));
            }
        }
    }
    

    // ------------------------------------------------------------
    // PHASE 1: BFS up to depth < l_depth for the DENSE portion
    // ------------------------------------------------------------


    std::queue<BFSItem> bfs_queue;
    bfs_queue.push({ root.get(), 0 });

    // We'll keep the BFS order of nodes that belong in the dense part
    std::vector<TrieNode*> dense_nodes;
    // We also record which nodes belong in the sparse part
    std::vector<TrieNode*> sparse_roots;  // all nodes at depth == l_depth

    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "enqueueing");

    while (!bfs_queue.empty())
    {
        auto [node_ptr, depth] = bfs_queue.front();
        bfs_queue.pop();

        if (!node_ptr->children.empty())
        {
            if (depth < l_depth)
            {
                // This node goes into the DENSE structure
                dense_nodes.push_back(node_ptr);

                // Enqueue children if they remain below l_depth
                for (auto & kv : node_ptr->children)
                {
                    // char c = kv.first;
                    TrieNode * child_ptr = kv.second.get();
                    // If child depth is still < l_depth, we keep BFSing in dense
                    if (depth + 1 < l_depth)
                    {
                        bfs_queue.push({ child_ptr, depth + 1 });
                    }
                    else
                    {
                        // The child is exactly at depth l_depth => goes to the SPARSE part
                        // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "child pushed to sparse {}", kv.first);
                        // sparse_roots.push_back(node_ptr);
                        sparse_roots.push_back(child_ptr);
                    }
                }
            }
        }
        // else
        // {
        //     char * str_ptr = key.data();
        //     if (depth > l_depth)
        //     {
        //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_value: {}", key);
        //         surf.dense.d_values.emplace_back(reinterpret_cast<uint64_t>(str_ptr));
        //     }
        //     else
        //     {
        //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_value: {}", key);
        //         surf.sparse.s_values.emplace_back(reinterpret_cast<uint64_t>(str_ptr));
        //     }
        // }
    }

    // ------------------------------------------------------------
    // Build the LOUDS-DENSE portion from 'dense_nodes'
    // (A minimal example that sets bits according to child edges.)
    // ------------------------------------------------------------


    // We need one entry per dense node in BFS order
    const size_t dense_count = dense_nodes.size();
    surf.dense.d_labels.resize(dense_count);
    surf.dense.d_hasChild.resize(dense_count);
    surf.dense.d_isPrefixKey.resize(dense_count, false);

    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "building LOUDS_DENSE {}", dense_count);

    char c = ' ';
    std::vector<char> key;

    // We might store a value in d_values for each terminal node (example).
    // The real logic depends on your design.
    // Let's store a dummy ID = i if node i is terminal.
    for (size_t i = 0; i < dense_count; ++i)
    {
        TrieNode * node = dense_nodes[i];
        surf.dense.d_isPrefixKey[i] = node->is_terminal;

        // For each child (with depth < l_depth) we encountered,
        // set the bits in d_labels and d_hasChild
        for (auto & kv : node->children)
        {
            c = kv.first;
            TrieNode * child_ptr = kv.second.get();

            // // We only set bits if the child's BFS node is also in dense_nodes
            // // (meaning it is depth < l_depth).
            // // We check if child_ptr is in the BFS portion for dense:
            // // (In a real design, you'd need a map from child_ptr->dense_index or something.)
            // // For simplicity, let's do a quick check:
            // auto it = std::find(dense_nodes.begin(), dense_nodes.end(), child_ptr);
            // if (it != dense_nodes.end())
            // {
            // Mark label c
            surf.dense.d_labels[i].set(static_cast<unsigned char>(c));
            // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "set d_labels: {} {}", i, c);

            // If child has any children => hasChild = true
            if (!child_ptr->children.empty())
            {
                surf.dense.d_hasChild[i].set(static_cast<unsigned char>(c));
            }
            // else
            // {

            //     surf.dense.d_values.push_back(static_cast<uint64_t>(c));
            // }
            
            // }
            // else
            // {
            //     // This child is in the sparse part => from the perspective of the dense node,
            //     // we treat it as if it doesn't continue with children (hasChild=0).
            //     // Optionally, you might treat it as a "terminal" here or store a pointer offset.
            //     surf.dense.d_labels[i].set(static_cast<unsigned char>(c));
            //     // no set on d_hasChild => remains false
            // }
        }

        // If you want to store a d_value for terminal nodes:
        // if (node->is_terminal)
        // {
        //     // Just push back an ID
        //     surf.dense.d_values.push_back(static_cast<uint64_t>(c));
        // }
    }

    // ------------------------------------------------------------
    // PHASE 2: BFS for the SPARSE portion
    // gather all nodes at depth >= l_depth
    // ------------------------------------------------------------
    // We'll do a BFS over the sub-tries whose roots we found.
    // To avoid duplicates, we might need a visited set if a node can appear multiple times.
    // For simplicity, we'll assume no overlap.

    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "enqueueing sparse");

    std::queue<TrieNode*> sparse_queue;
    for (TrieNode * sp_root : sparse_roots)
        sparse_queue.push(sp_root);

    // We'll store them in BFS order
    std::vector<TrieNode*> sparse_nodes;
    while (!sparse_queue.empty())
    {
        auto * node = sparse_queue.front();
        sparse_queue.pop();
        sparse_nodes.push_back(node);

        // Enqueue children
        for (auto & kv : node->children)
        {
            sparse_queue.push(kv.second.get());
        }
    }

    // ------------------------------------------------------------
    // Build the LOUDS-SPARSE portion from 'sparse_nodes'
    // ------------------------------------------------------------
    // For each node in BFS order, we add its labels to s_labels.
    // s_LOUDS[i] = true if s_labels[i] is the first label of a node;
    //             false if it continues the same node.
    // If node->is_terminal => we add 0xFF as a special label at the start.

    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "building LOUDS_SPARSE");

    for (size_t i = 0; i < sparse_nodes.size(); ++i)
    {
        TrieNode * node = sparse_nodes[i];

        // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "sparse node is_terminal: {} children: {}", node->is_terminal, node->children.size());

        // Indicate that the *next* label we write belongs to a new node
        bool first_label_of_node = true;

        // If this node is terminal => put 0xFF as a special label
        if (node->is_terminal && !node->children.empty())
        {
            // surf.sparse.s_labels.push_back(0xFF);
            // surf.sparse.s_hasChild.push_back(false);     // 0xFF doesn't lead to a real child
            // surf.sparse.s_LOUDS.push_back(first_label_of_node); 
            // first_label_of_node = false;

            // // Optionally store a corresponding value
            // surf.sparse.s_values.push_back(static_cast<uint64_t>(i));
            // node->children.push_back(std::make_pair(0xFF, nullptr));
            node->children[0xFF] = std::make_unique<TrieNode>();
        }

        // For each real child label, push it onto s_labels
        for (auto & kv : node->children)
        {
            uint8_t val = static_cast<uint8_t>(kv.first);
            surf.sparse.s_labels.push_back(val);
            // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_label {}", kv.first);

            TrieNode * child_ptr = kv.second.get();
            bool has_subtrie = !child_ptr->children.empty(); 
            surf.sparse.s_hasChild.push_back(has_subtrie);

            // s_LOUDS: set to true if this is the first label for this node;
            // once set, we reset first_label_of_node so subsequent labels are false.
            surf.sparse.s_LOUDS.push_back(first_label_of_node);
            first_label_of_node = false;
        }
    }

    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "LOUDS:");

    // for (size_t i = 0; i < surf.dense.d_labels.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_labels: {} {}", i, surf.dense.d_labels[i].to_string());
    // }
    // for (size_t i = 0; i < surf.dense.d_hasChild.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_hasChild: {} {}", i, surf.dense.d_hasChild[i].to_string());
    // }
    // for (size_t i = 0; i < surf.dense.d_isPrefixKey.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_isPrefixKey: {} {}", i, surf.dense.d_isPrefixKey[i]);
    // }
    // for (size_t i = 0; i < surf.dense.d_values.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} {}", i, surf.dense.d_values[i]);
    // }
    // for (size_t i = 0; i < surf.dense.d_values.size(); ++i)
    // {
    //     std::vector<char> * buffer_ptr = reinterpret_cast<std::vector<char>*>(surf.dense.d_values[i]);
    //     if (buffer_ptr)
    //     {
    //         std::vector<char> buffer = *buffer_ptr;
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} {}", i, buffer);
    //     }
    //     else
    //     {
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} nullptr", i);
    //     }

    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} {}", i, surf.dense.d_values[i]);
    // }

    // for (size_t i = 0; i < surf.sparse.s_labels.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_labels: {} {}", i, surf.sparse.s_labels[i]);
    // }
    // for (size_t i = 0; i < surf.sparse.s_hasChild.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_hasChild: {} {}", i, surf.sparse.s_hasChild[i]);
    // }
    // for (size_t i = 0; i < surf.sparse.s_LOUDS.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_LOUDS: {} {}", i, surf.sparse.s_LOUDS[i]);
    // }
    // for (size_t i = 0; i < surf.sparse.s_values.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_values: {} {}", i, surf.sparse.s_values[i]);
    // }
    // for (size_t i = 0; i < surf.sparse.s_values.size(); ++i)
    // {
    //     std::vector<char> * buffer_ptr = reinterpret_cast<std::vector<char>*>(surf.sparse.s_values[i]);
    //     if (buffer_ptr)
    //     {
    //         std::vector<char> buffer = *buffer_ptr;
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_values: {} {}", i, buffer);
    //     }
    //     else
    //     {
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_values: {} nullptr", i);
    //     }
    // }
}

std::optional<uint64_t> SuccinctRangeFilter::ExactKeySearch(const std::string & key) const {
    size_t level = 0;
    // Starting positions in the dense or sparse arrays.
    size_t posDense = 0;   // pointer into the current dense level (conceptually)
    size_t posSparse = 0;  // pointer into the sparse part (once the cutoff is reached)
    
    // First, search the dense portion (levels 0 .. cutoff-1).
    for (size_t i = 0; i < key.size() && level < l_depth; i++, level++) {
        uint8_t target = static_cast<uint8_t>(key[i]);
        // In the current dense level, check if the label exists.
        if (!surf.dense.d_labels[level].test(target))
            return std::nullopt; // key byte not found
        // Check if this label has a child.
        if (!surf.dense.d_hasChild[level].test(target)) {
            // Leaf reached in dense portion; if the key length exactly matches, return its value.
            if (i == key.size() - 1 && surf.dense.d_isPrefixKey[posDense])
                return surf.dense.d_values[posDense];
            else
                return std::nullopt;
        }
        // Otherwise, update posDense using a (hypothetical) rank/select on d_hasChild.
        posDense = childPositionDense(posDense, target, level);
    }
    
    // If key remains and we are in the sparse portion:
    for (size_t i = level; i < key.size(); i++, level++) {
        uint8_t target = static_cast<uint8_t>(key[i]);
        // Determine the boundaries of the current node in the sparse arrays.
        size_t nodeStart = getSparseNodeStart(posSparse);
        size_t nodeEnd = getSparseNodeEnd(posSparse);
        bool found = false;
        for (size_t j = nodeStart; j < nodeEnd; j++) {
            if (surf.sparse.s_labels[j] == target) {
                found = true;
                // If no child exists, then we have reached a leaf.
                if (!surf.sparse.s_hasChild[j])
                    return surf.sparse.s_values[j];
                // Otherwise, update posSparse to the child node position.
                posSparse = childPositionSparse(j);
                break;
            }
        }
        if (!found)
            return std::nullopt;
    }
    
    // If the loop completes exactly, then the key exactly matched a prefix.
    // Check the final node’s flag and return its value if present.
    if (level < l_depth) {
        if (surf.dense.d_isPrefixKey[posDense])
            return surf.dense.d_values[posDense];
    } else {
        // (Assume similar handling in the sparse portion.)
        // For brevity, we return nullopt if the node is not a valid key.
        return std::nullopt;
    }
    return std::nullopt;
}

Iterator SuccinctRangeFilter::LowerBound(const std::string & key) const {
    Iterator iter;
    size_t level = 0;
    size_t posDense = 0;   // current pointer in dense part
    size_t posSparse = 0;  // current pointer in sparse part
    
    // Process dense levels.
    for (; level < key.size() && level < l_depth; level++) {
        uint8_t target = static_cast<uint8_t>(key[level]);
        // In a dense node, find the smallest label >= target.
        size_t candidate = findLowerBoundInDense(/*posDense, */level, target);
        if (candidate == 256) { // not found in current node; need to backtrack
            // For simplicity, mark iterator invalid.
            iter.valid = false;
            return iter;
        }
        iter.levelPositions.push_back(posDense); // store current position (for illustration)
        // If candidate does not have a child, we have reached a leaf.
        if (!surf.dense.d_hasChild[level].test(candidate)) {
            iter.currentLevel = level;
            return iter;
        }
        posDense = childPositionDense(posDense, candidate, level);
    }
    
    // Process sparse levels.
    for (; level < key.size(); level++) {
        uint8_t target = static_cast<uint8_t>(key[level]);
        size_t nodeStart = getSparseNodeStart(posSparse);
        size_t nodeEnd = getSparseNodeEnd(posSparse);
        // uint8_t candidate = 255;
        size_t candidatePos = 0;
        bool found = false;
        for (size_t j = nodeStart; j < nodeEnd; j++) {
            uint8_t label = surf.sparse.s_labels[j];
            if (label >= target) {
                // candidate = label;
                candidatePos = j;
                found = true;
                break;
            }
        }
        if (!found) {
            iter.valid = false;
            return iter;
        }
        iter.levelPositions.push_back(candidatePos);
        if (!surf.sparse.s_hasChild[candidatePos]) {
            iter.currentLevel = level;
            return iter;
        }
        posSparse = childPositionSparse(candidatePos);
    }
    
    iter.currentLevel = level;
    return iter;
}

bool SuccinctRangeFilter::MoveToNext(Iterator & iter) const {
    // Start from the current leaf position and try to move forward in the node.
    size_t level = iter.currentLevel;
    // Attempt to find a next label in the current node.
    if (canMoveNext(iter, level)) {
        // Move within the same node.
        iter.levelPositions[level] = nextPosition(iter.levelPositions[level], level);
        // Descend to the left-most key in the subtree.
        while (hasChild(iter.levelPositions[level], level)) {
            level++;
            size_t childPos = leftMostChild(iter.levelPositions[level - 1], level);
            iter.levelPositions.push_back(childPos);
        }
        iter.currentLevel = level;
        return true;
    }
    else {
        // Backtrack: move up until we find a node where we can advance.
        while (level > 0) {
            level--;
            if (canMoveNext(iter, level)) {
                iter.levelPositions[level] = nextPosition(iter.levelPositions[level], level);
                // Remove deeper levels and descend again.
                iter.levelPositions.resize(level + 1);
                while (hasChild(iter.levelPositions[level], level)) {
                    level++;
                    size_t childPos = leftMostChild(iter.levelPositions[level - 1], level);
                    iter.levelPositions.push_back(childPos);
                }
                iter.currentLevel = level;
                return true;
            }
        }
    }
    // If we cannot backtrack further, we have reached the end.
    iter.valid = false;
    return false;
}

}
