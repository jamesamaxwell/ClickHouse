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

SuccinctRangeFilter::SuccinctRangeFilter(std::unique_ptr<TrieNode> root, size_t l_depth) // TODO: change this so that it make full dense and sparse and then finds depth
{
    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "SuccinctRangeFilter");
    // Clear existing data
    surf.dense = LOUDSDenseTrie();
    surf.sparse = LOUDSSparseTrie();

    if (!root)
        return; // empty trie => nothing to do

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
            // else
            // {
            //     // Node is at or beyond l_depth => goes to sparse
            //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "parent pushed to sparse");
            //     for (auto & kv : node_ptr->children)
            //     {
            //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "{}", kv.first);
            //     }
            //     sparse_roots.push_back(node_ptr);
            // }
        }
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
                surf.dense.d_hasChild[i].set(static_cast<unsigned char>(c));
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
        if (node->is_terminal)
        {
            // Just push back an ID
            surf.dense.d_values.push_back(static_cast<uint64_t>(c));
        }
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

    // LOG_DEBUG(getLogger("SuccinctRangeFilter"), "built LOUDS_SPARSE {}", surf.sparse.s_labels.size());

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
}

}
