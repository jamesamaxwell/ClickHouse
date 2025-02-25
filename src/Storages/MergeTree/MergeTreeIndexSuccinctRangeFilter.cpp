#include <Storages/MergeTree/MergeTreeIndexSuccinctRangeFilter.h>

#include <Columns/ColumnArray.h>
#include <Columns/ColumnConst.h>
#include <Columns/ColumnFixedString.h>
#include <Columns/ColumnNullable.h>
#include <Columns/ColumnString.h>
#include <Columns/ColumnTuple.h>
#include <Columns/ColumnsNumber.h>
#include <Common/FieldVisitorsAccurateComparison.h>
#include <Common/HashTable/ClearableHashMap.h>
#include <Common/HashTable/Hash.h>
#include <DataTypes/DataTypeArray.h>
#include <DataTypes/DataTypeMap.h>
#include <DataTypes/DataTypeNullable.h>
#include <DataTypes/DataTypeTuple.h>
#include <DataTypes/DataTypesNumber.h>
#include <IO/WriteHelpers.h>
#include <Interpreters/ExpressionAnalyzer.h>
#include <Interpreters/TreeRewriter.h>
#include <Interpreters/castColumn.h>
#include <Interpreters/convertFieldToType.h>
#include <Interpreters/misc.h>
#include <Parsers/ASTFunction.h>
#include <Parsers/ASTIdentifier.h>
#include <Parsers/ASTLiteral.h>
#include <Parsers/ASTSelectQuery.h>
#include <Parsers/ASTSubquery.h>
#include <Storages/MergeTree/MergeTreeData.h>
#include <Storages/MergeTree/MergeTreeIndexUtils.h>
#include <Storages/MergeTree/RPNBuilder.h>
#include <base/types.h>

// Delete these eventually
#include <Interpreters/BloomFilterHash.h>
#include <Interpreters/BloomFilter.h>

#include <Common/logger_useful.h>

namespace DB
{


namespace ErrorCodes
{
    extern const int BAD_ARGUMENTS;
    extern const int ILLEGAL_COLUMN;
    extern const int ILLEGAL_TYPE_OF_ARGUMENT;
    extern const int INCORRECT_QUERY;
    extern const int NUMBER_OF_ARGUMENTS_DOESNT_MATCH;
    extern const int LOGICAL_ERROR;
}

ColumnWithTypeAndName getPreparedSetInfo(const ConstSetPtr & prepared_set)
{
    if (prepared_set->getDataTypes().size() == 1)
        return {prepared_set->getSetElements()[0], prepared_set->getElementsTypes()[0], "dummy"};

    Columns set_elements;
    for (auto & set_element : prepared_set->getSetElements())

        set_elements.emplace_back(set_element->convertToFullColumnIfConst());

    return {ColumnTuple::create(set_elements), std::make_shared<DataTypeTuple>(prepared_set->getElementsTypes()), "dummy"};
}

MergeTreeIndexGranuleSuccinctRangeFilter::MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, const TrieNode & root, size_t total_rows_)
    : ds_ratio(ds_ratio_), total_rows(total_rows_)
{
    num_columns = 1;
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter {}", root.is_terminal);

    // if (root.children.size() == 0)
    // {
    //     LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "root has no children");
    // }
    // else
    // {
    //     for (auto & child : root.children)
    //     {
    //         LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter root child {}", child.first);
    //     }
    //     for (auto & child : root.children)
    //     {
    //         for (auto & grandchild : child.second->children)
    //         {
    //             LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter root grandchild {}", grandchild.first);
    //         }
    //     }
    // }

    // Filter superfluous nodes in the trie
    auto [pruned_root, total_terminals] = pruneSubtree(root);
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "trie pruned {}", total_terminals);

    // if (pruned_root == nullptr)
    // {
    //     LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "pruned_root is null");
    // }
    // else
    // {
    //     for (auto & child : pruned_root->children)
    //     {
    //         LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter pruned child {}", child.first);

    //     }
    //     for (auto & child : pruned_root->children)
    //     {
    //         for (auto & grandchild : child.second->children)
    //         {
    //             LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter pruned grandchild {}", grandchild.first);
    //         }
    //     }
    // }

    size_t l_depth = findLargestDepth(pruned_root, ds_ratio);
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "LOUDS_DENSE depth {}", l_depth);

    if (pruned_root == nullptr)
    {
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "pruned_root is null");
    }
    else
    {
        SuccinctRangeFilterPtr surf = std::make_shared<SuccinctRangeFilter>(std::move(pruned_root), l_depth);

        dense_nodes = surf->getFilter().dense.d_labels.size();
        d_values = surf->getFilter().dense.d_values.size();

        sparse_nodes = surf->getFilter().sparse.s_labels.size();
        s_values = surf->getFilter().sparse.s_values.size();

        surfs.push_back(surf);
    }
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "Constructor surfs size: {}", surfs.size());
}

/**
 * Recursively prunes a subtree. Returns:
 *   ( pruned_subtree_root, number_of_terminals_in_subtree ).
 *
 * If the subtree contains exactly 1 terminal, all children are dropped
 * and the returned root is marked terminal.
 */
std::pair<std::unique_ptr<TrieNode>, size_t> MergeTreeIndexGranuleSuccinctRangeFilter::pruneSubtree(const TrieNode & old_node)
{
    // if (old_node == nullptr)
    // {
    //     LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "old_node is null");
    // }
    // Create a new node that will hold the pruned version of old_node.
    auto new_node = std::make_unique<TrieNode>();
    new_node->is_terminal = old_node.is_terminal;

    // Count how many terminals are in this entire subtree (including itself).
    size_t total_terminals = old_node.is_terminal ? 1 : 0;

    // Recursively process each child.
    for (const auto & kv : old_node.children)
    {
        // kv.first  = character label
        // kv.second = unique_ptr<TrieNode>
        const TrieNode * child_ptr = kv.second.get();
        auto [pruned_child, child_terminals] = pruneSubtree(*child_ptr);

        // If child subtree has at least 1 terminal, keep it in the new node
        // (unless we decide to prune further below).
        if (child_terminals > 0)
        {
            new_node->children[kv.first] = std::move(pruned_child);
            total_terminals += child_terminals;
        }
    }

    // If exactly 1 terminal was found below this node, we can prune all children
    // and just mark this node as terminal.
    if (total_terminals == 1)
    {
        new_node->children.clear();  // drop all children
        new_node->is_terminal = true;
    }

    return {std::move(new_node), total_terminals};
}

/**
 * @brief  Finds the largest depth l_depth satisfying:
 *         (# of nodes with depth < l_depth) * ds_ratio <= (# of nodes with depth >= l_depth)
 * @param  root     Reference to the root of the trie (depth 0).
 * @param  ratio The ratio used in the comparison.
 * @return The largest depth l_depth that satisfies the condition.
 */
size_t MergeTreeIndexGranuleSuccinctRangeFilter::findLargestDepth(const std::unique_ptr<TrieNode> & root, size_t ratio)
{
    // If there's no trie at all, largest depth is 0.
    if (!root)
        return 0;

    // 1) Traverse the trie (BFS) to count how many nodes exist at each depth.
    std::queue<std::pair<const TrieNode*, size_t>> node_queue;
    node_queue.push({ root.get(), 0 });

    std::vector<size_t> count_at_depth;  // count_at_depth[d] = # of nodes at depth d
    size_t max_depth = 0;

    while (!node_queue.empty())
    {
        auto [node_ptr, depth] = node_queue.front();
        node_queue.pop();

        if (depth >= count_at_depth.size())
            count_at_depth.resize(depth + 1, 0);

        count_at_depth[depth]++;
        max_depth = std::max(max_depth, depth);

        // Enqueue children with depth+1
        for (const auto & kv : node_ptr->children)
        {
            node_queue.push({ kv.second.get(), depth + 1 });
        }
    }

    // 2) Build a prefix sum over count_at_depth for quick lookups.
    //    prefix[d] = total nodes up to (and including) depth d.
    std::vector<size_t> prefix(count_at_depth.size());
    prefix[0] = count_at_depth[0];
    for (size_t d = 1; d <= max_depth; ++d)
        prefix[d] = prefix[d - 1] + count_at_depth[d];

    const size_t total_nodes = prefix[max_depth];

    // 3) Find the largest depth d s.t. (# nodes with depth < d)*ratio <= (# nodes with depth >= d).
    for (int d = static_cast<int>(max_depth); d >= 0; --d)
    {
        size_t n_less = (d == 0) ? 0 : prefix[d - 1];
        size_t n_at_least = total_nodes - n_less;

        if (static_cast<double>(n_less) * ratio <= static_cast<double>(n_at_least))
        {
            return static_cast<size_t>(d);
        }
    }

    return 0;
}


MergeTreeIndexGranuleSuccinctRangeFilter::MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, size_t num_columns_)
    : ds_ratio(ds_ratio_), num_columns(num_columns_)
{
    total_rows = 0;
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter second constructor {}", num_columns_);
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "ds_ratio {}", ds_ratio);
    for (size_t column = 0; column < num_columns_; ++column)
        surfs.push_back(std::make_shared<SuccinctRangeFilter>(nullptr, ds_ratio));
}

bool MergeTreeIndexGranuleSuccinctRangeFilter::empty() const
{
    bool empty = surfs.empty();
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "empty granule {}", empty);
    /// One possible definition:
    /// Return true if we haven't added any rows yet, or if no filters exist.
    return empty;
}

void MergeTreeIndexGranuleSuccinctRangeFilter::serializeBinary(WriteBuffer & ostr) const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "serializeBinary");
    writeVarUInt(total_rows, ostr);
    writeVarUInt(dense_nodes, ostr);
    writeVarUInt(d_values, ostr);
    writeVarUInt(sparse_nodes, ostr);
    writeVarUInt(s_values, ostr);

    SuccinctRangeFilterPtr surf = surfs[0];

    size_t d_labels_size = dense_nodes * sizeof(std::bitset<256>);
    size_t d_has_child_size = dense_nodes * sizeof(std::bitset<256>);
    size_t d_is_prefix_key_size = dense_nodes * sizeof(bool);
    size_t d_values_size = d_values * sizeof(uint64_t);

    size_t s_labels_size = sparse_nodes * sizeof(uint8_t);
    size_t s_has_child_size = sparse_nodes * sizeof(bool);
    size_t s_louds_size = sparse_nodes * sizeof(bool);
    size_t s_values_size = s_values * sizeof(uint64_t);

    size_t louds_dense_size = d_labels_size + d_has_child_size + d_is_prefix_key_size + d_values_size;
    size_t louds_sparse_size = s_labels_size + s_has_child_size + s_louds_size + s_values_size;
    
    size_t write_size = louds_dense_size + louds_sparse_size;

    // for (const auto & surf : surfs)
    // {
    // const char * write_buffer = reinterpret_cast<const char *>(surf->getFilter().dense.d_labels.data());

    write_size = surf->getWriteSize();
    ostr.write(serializeTrie(surf->getFilter(),
                    d_labels_size,
                    d_has_child_size,
                    d_is_prefix_key_size,
                    d_values_size,
                    s_labels_size,
                    s_has_child_size,
                    s_louds_size,
                    s_values_size).data(), write_size);
    //}
}

std::vector<char> MergeTreeIndexGranuleSuccinctRangeFilter::serializeTrie(const LOUDSdsTrie &trie,
    size_t denseLabelsBytes,
    size_t denseHasChildBytes,
    size_t denseIsPrefixKeyBytes,
    size_t denseValuesBytes,
    size_t sparseLabelsBytes,
    size_t sparseHasChildBytes,
    size_t sparseLOUDSBytes,
    size_t sparseValuesBytes) const
{
    // Compute total number of bytes in the final buffer.
    size_t totalBytes = denseLabelsBytes + denseHasChildBytes + denseIsPrefixKeyBytes +
                        denseValuesBytes + sparseLabelsBytes + sparseHasChildBytes +
                        sparseLOUDSBytes + sparseValuesBytes;
    
    std::vector<char> buffer(totalBytes);
    size_t offset = 0;
    
    // Copy the dense.d_labels vector.
    // Each element is a std::bitset<256> (typically 32 bytes, but this may depend on the implementation).
    memcpy(buffer.data() + offset,
           reinterpret_cast<const char*>(trie.dense.d_labels.data()),
           denseLabelsBytes);
    offset += denseLabelsBytes;
    
    // Copy the dense.d_hasChild vector.
    memcpy(buffer.data() + offset,
           reinterpret_cast<const char*>(trie.dense.d_hasChild.data()),
           denseHasChildBytes);
    offset += denseHasChildBytes;
    
    std::vector<char> packedBoolVector = packBoolVector(trie.dense.d_isPrefixKey);
    memcpy(buffer.data() + offset,
        packedBoolVector.data(),
        packedBoolVector.size());
    offset += packedBoolVector.size();
    
    // Copy the dense.d_values vector.
    memcpy(buffer.data() + offset,
           reinterpret_cast<const char*>(trie.dense.d_values.data()),
           denseValuesBytes);
    offset += denseValuesBytes;
    
    // Copy the sparse.s_labels vector.
    memcpy(buffer.data() + offset,
           reinterpret_cast<const char*>(trie.sparse.s_labels.data()),
           sparseLabelsBytes);
    offset += sparseLabelsBytes;
    
    packedBoolVector = packBoolVector(trie.sparse.s_hasChild);
    memcpy(buffer.data() + offset,
        packedBoolVector.data(),
        packedBoolVector.size());
    offset += packedBoolVector.size();

    packedBoolVector = packBoolVector(trie.sparse.s_LOUDS);
    memcpy(buffer.data() + offset,
        packedBoolVector.data(),
        packedBoolVector.size());
    offset += packedBoolVector.size();
    
    // Copy the sparse.s_values vector.
    memcpy(buffer.data() + offset,
           reinterpret_cast<const char*>(trie.sparse.s_values.data()),
           sparseValuesBytes);
    // offset += sparseValuesBytes; // (Not really needed at the end.)
    
    return buffer;
}

std::vector<char> MergeTreeIndexGranuleSuccinctRangeFilter::packBoolVector(const std::vector<bool>& boolVec) const {
    size_t numBytes = (boolVec.size() + 7) / 8; // 8 bits per byte
    std::vector<char> packed(numBytes, 0);
    for (size_t i = 0; i < boolVec.size(); ++i) {
        if (boolVec[i]) {
            packed[i / 8] |= (1 << (i % 8));
        }
    }
    return packed;
}

void MergeTreeIndexGranuleSuccinctRangeFilter::deserializeBinary(ReadBuffer & istr, MergeTreeIndexVersion version)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "deserializeBinary");
    if (version != 1)
        throw Exception(ErrorCodes::LOGICAL_ERROR, "Unknown index version {}.", version);

    readVarUInt(total_rows, istr);
    readVarUInt(dense_nodes, istr);
    readVarUInt(d_values, istr);
    readVarUInt(sparse_nodes, istr);
    readVarUInt(s_values, istr);

    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "total_rows: {} dense_nodes: {} d_values: {} sparse_nodes: {} s_values: {}", total_rows, dense_nodes, d_values, sparse_nodes, s_values);

    size_t d_labels_size = dense_nodes * sizeof(std::bitset<256>);
    size_t d_has_child_size = dense_nodes * sizeof(std::bitset<256>);
    size_t d_is_prefix_key_size = dense_nodes * sizeof(bool);
    size_t d_values_size = d_values * sizeof(uint64_t);

    size_t s_labels_size = sparse_nodes * sizeof(uint8_t);
    size_t s_has_child_size = sparse_nodes * sizeof(bool);
    size_t s_louds_size = sparse_nodes * sizeof(bool);
    size_t s_values_size = s_values * sizeof(uint64_t);

    size_t louds_dense_size = d_labels_size + d_has_child_size + d_is_prefix_key_size + d_values_size;
    size_t louds_sparse_size = s_labels_size + s_has_child_size + s_louds_size + s_values_size;
    
    size_t write_size = louds_dense_size + louds_sparse_size;

    // std::vector<char> buffer (1024,0);
    std::vector<char> buffer(write_size);
    istr.readStrict(buffer.data(), write_size);

    SuccinctRangeFilterPtr surf = surfs[0];
    // LOUDSdsTrie &trie = surf->getFilter();
    size_t offset = 0;

    surf->getFilter().dense.d_labels.resize(d_labels_size / sizeof(std::bitset<256>));
    memcpy(surf->getFilter().dense.d_labels.data(),
            buffer.data() + offset,
            d_labels_size);
    offset += d_labels_size;

    surf->getFilter().dense.d_hasChild.resize(d_has_child_size / sizeof(std::bitset<256>));
    memcpy(surf->getFilter().dense.d_hasChild.data(),
            buffer.data() + offset,
            d_has_child_size);
    offset += d_has_child_size;

    surf->getFilter().dense.d_isPrefixKey.resize(d_is_prefix_key_size);
    for (size_t i = 0; i < d_is_prefix_key_size && i < surf->getFilter().dense.d_isPrefixKey.size(); i++) {
        surf->getFilter().dense.d_isPrefixKey[i] = (buffer[offset + (i / 8)] >> (i % 8)) & 1;
    }
    offset += (d_is_prefix_key_size + 7) / 8;

    surf->getFilter().dense.d_values.resize(d_values_size / sizeof(uint64_t));
    memcpy(surf->getFilter().dense.d_values.data(),
           buffer.data() + offset,
           d_values_size);
    offset += d_values_size;

    surf->getFilter().sparse.s_labels.resize(s_labels_size / sizeof(uint8_t));
    memcpy(surf->getFilter().sparse.s_labels.data(),
           buffer.data() + offset,
           s_labels_size);
    offset += s_labels_size;

    surf->getFilter().sparse.s_hasChild.resize(s_has_child_size);
    for (size_t i = 0; i < s_has_child_size * 8 && i < surf->getFilter().sparse.s_hasChild.size(); i++) {
        surf->getFilter().sparse.s_hasChild[i] = (buffer[offset + (i / 8)] >> (i % 8)) & 1;
    }
    offset += (s_has_child_size + 7) / 8;

    surf->getFilter().sparse.s_LOUDS.resize(s_louds_size);
    for (size_t i = 0; i < s_louds_size * 8 && i < surf->getFilter().sparse.s_LOUDS.size(); i++) {
        surf->getFilter().sparse.s_LOUDS[i] = (buffer[offset + (i / 8)] >> (i % 8)) & 1;
    }
    offset += (s_louds_size + 7) / 8;

    surf->getFilter().sparse.s_values.resize(s_values_size / sizeof(uint64_t));
    memcpy(surf->getFilter().sparse.s_values.data(),
           buffer.data() + offset,
           s_values_size);
    offset += s_values_size;
    
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "deserialized LOUDS:");
    // LOUDSdsTrie &trie = surf->getFilter();

    // for (size_t i = 0; i < trie.dense.d_labels.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_labels: {} {}", i, trie.dense.d_labels[i].to_string());
    // }
    // for (size_t i = 0; i < trie.dense.d_hasChild.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_hasChild: {} {}", i, trie.dense.d_hasChild[i].to_string());
    // }
    // for (size_t i = 0; i < trie.dense.d_isPrefixKey.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_isPrefixKey: {} {}", i, trie.dense.d_isPrefixKey[i]);
    // }
    // for (size_t i = 0; i < trie.dense.d_values.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} {}", i, trie.dense.d_values[i]);
    // }
    // for (size_t i = 0; i < trie.dense.d_values.size(); ++i)
    // {
    //     std::vector<char> * buffer_ptr = reinterpret_cast<std::vector<char>*>(trie.dense.d_values[i]);
    //     if (buffer_ptr)
    //     {
    //         std::vector<char> buffer0 = *buffer_ptr;
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} {}", i, buffer0);
    //     }
    //     else
    //     {
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} nullptr", i);
    //     }

    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "d_values: {} {}", i, trie.dense.d_values[i]);
    // }

    // for (size_t i = 0; i < trie.sparse.s_labels.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_labels: {} {}", i, trie.sparse.s_labels[i]);
    // }
    // for (size_t i = 0; i < trie.sparse.s_hasChild.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_hasChild: {} {}", i, trie.sparse.s_hasChild[i]);
    // }
    // for (size_t i = 0; i < trie.sparse.s_LOUDS.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_LOUDS: {} {}", i, trie.sparse.s_LOUDS[i]);
    // }
    // for (size_t i = 0; i < trie.sparse.s_values.size(); ++i)
    // {
    //     LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_values: {} {}", i, trie.sparse.s_values[i]);
    // }
    // for (size_t i = 0; i < trie.sparse.s_values.size(); ++i)
    // {
    //     std::vector<char> * buffer_ptr = reinterpret_cast<std::vector<char>*>(trie.sparse.s_values[i]);
    //     if (buffer_ptr)
    //     {
    //         std::vector<char> buffer0 = *buffer_ptr;
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_values: {} {}", i, buffer0);
    //     }
    //     else
    //     {
    //         LOG_DEBUG(getLogger("SuccinctRangeFilter"), "s_values: {} nullptr", i);
    //     }
    // }
}

void fillingSuccinctRangeFilter(SuccinctRangeFilterPtr & surf)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "fillingSuccinctRangeFilter");
    surf->getFilter();
}

MergeTreeIndexConditionSuccinctRangeFilter::MergeTreeIndexConditionSuccinctRangeFilter(
    const ActionsDAG * filter_actions_dag, 
    ContextPtr context_, 
    const Block & header_)
    : WithContext(context_)
    , header(header_)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexConditionSuccinctRangeFilter");
    if (!filter_actions_dag)
    {
        rpn.push_back(RPNElement::FUNCTION_UNKNOWN);
        return;
    }

    RPNBuilder<RPNElement> builder(
        filter_actions_dag->getOutputs().at(0),
        context_,
        [&](const RPNBuilderTreeNode & node, RPNElement & out) { return extractAtomFromTree(node, out); });
    rpn = std::move(builder).extractRPN();
}


bool MergeTreeIndexConditionSuccinctRangeFilter::alwaysUnknownOrTrue() const
{
    std::vector<bool> rpn_stack;

    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "alwaysUnknownOrTrue");

    for (const auto & element : rpn)
    {
        if (element.function == RPNElement::FUNCTION_UNKNOWN)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_UNKNOWN");
        }
        else if (element.function == RPNElement::ALWAYS_TRUE)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "ALWAYS_TRUE");
        }
        else if (element.function == RPNElement::FUNCTION_EQUALS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_EQUALS");
        }
        else if (element.function == RPNElement::FUNCTION_NOT_EQUALS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_NOT_EQUALS");
        }
        else if (element.function == RPNElement::FUNCTION_HAS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_HAS");
        }
        else if (element.function == RPNElement::FUNCTION_HAS_ANY)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_HAS_ANY");
        }
        else if (element.function == RPNElement::FUNCTION_HAS_ALL)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_HAS_ALL");
        }
        else if (element.function == RPNElement::FUNCTION_IN)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_IN");
        }
        else if (element.function == RPNElement::FUNCTION_NOT_IN)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_NOT_IN");
        }
        else if (element.function == RPNElement::ALWAYS_FALSE)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "ALWAYS_FALSE");
        }
        else if (element.function == RPNElement::FUNCTION_NOT)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_NOT");
        }
        else if (element.function == RPNElement::FUNCTION_AND)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_AND");
        }
        else if (element.function == RPNElement::FUNCTION_OR)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_OR");
        }
        else if (element.function == RPNElement::FUNCTION_GREATER)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_GREATER");
        }
        else if (element.function == RPNElement::FUNCTION_LESS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_LESS");
        }
        else if (element.function == RPNElement::FUNCTION_GREATER_OR_EQUALS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_GREATER_OR_EQUALS");
        }
        else if (element.function == RPNElement::FUNCTION_LESS_OR_EQUALS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), "FUNCTION_LESS_OR_EQUALS");
        }
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Unexpected function type in KeyCondition::RPNElement");

        if (element.function == RPNElement::FUNCTION_UNKNOWN
            || element.function == RPNElement::ALWAYS_TRUE)
        {
            rpn_stack.push_back(true);
        }
        else if (element.function == RPNElement::FUNCTION_EQUALS
            || element.function == RPNElement::FUNCTION_NOT_EQUALS
            || element.function == RPNElement::FUNCTION_HAS
            || element.function == RPNElement::FUNCTION_HAS_ANY
            || element.function == RPNElement::FUNCTION_HAS_ALL
            || element.function == RPNElement::FUNCTION_IN
            || element.function == RPNElement::FUNCTION_NOT_IN
            || element.function == RPNElement::ALWAYS_FALSE
            || element.function == RPNElement::FUNCTION_GREATER
            || element.function == RPNElement::FUNCTION_LESS
            || element.function == RPNElement::FUNCTION_LESS_OR_EQUALS
            || element.function == RPNElement::FUNCTION_GREATER_OR_EQUALS)
        {
            rpn_stack.push_back(false);
        }
        else if (element.function == RPNElement::FUNCTION_NOT)
        {
            // do nothing
        }
        else if (element.function == RPNElement::FUNCTION_AND)
        {
            auto arg1 = rpn_stack.back();
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "arg1 {}", arg1);
            rpn_stack.pop_back();
            auto arg2 = rpn_stack.back();
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "arg2 {}", arg2);
            rpn_stack.back() = arg1 && arg2;
        }
        else if (element.function == RPNElement::FUNCTION_OR)
        {
            auto arg1 = rpn_stack.back();
            rpn_stack.pop_back();
            auto arg2 = rpn_stack.back();
            rpn_stack.back() = arg1 || arg2;
        }
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Unexpected function type in KeyCondition::RPNElement");
    }

    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "allwaysUnknownOrTrue {}", rpn_stack[0]);
    return rpn_stack[0];
}

bool MergeTreeIndexConditionSuccinctRangeFilter::mayBeTrueOnGranule(const MergeTreeIndexGranuleSuccinctRangeFilter * granule) const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "mayBeTrueOnGranule");
    std::vector<BoolMask> rpn_stack;
    const auto & filters = granule->getFilters();
    if (!filters.empty())
    {
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "filters not empty");
    }

    for (const auto & element : rpn)
    {
        rpn_stack.emplace_back(true, true);
        if (element.function == RPNElement::FUNCTION_UNKNOWN)
        {
            rpn_stack.emplace_back(true, true);
        }
        else if (element.function == RPNElement::FUNCTION_IN
            || element.function == RPNElement::FUNCTION_NOT_IN
            || element.function == RPNElement::FUNCTION_EQUALS
            || element.function == RPNElement::FUNCTION_NOT_EQUALS
            || element.function == RPNElement::FUNCTION_HAS
            || element.function == RPNElement::FUNCTION_HAS_ANY
            || element.function == RPNElement::FUNCTION_HAS_ALL
            || element.function == RPNElement::FUNCTION_GREATER
            || element.function == RPNElement::FUNCTION_LESS
            || element.function == RPNElement::FUNCTION_LESS_OR_EQUALS
            || element.function == RPNElement::FUNCTION_GREATER_OR_EQUALS)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ first case ------------------------------------------------------ ");
            // bool match_rows = true;
            // bool match_all = element.function == RPNElement::FUNCTION_HAS_ALL;
            // const auto & predicate = element.predicate;
            // for (size_t index = 0; match_rows && index < predicate.size(); ++index)
            // {
            //     const auto & query_index_hash = predicate[index];
            //     const auto & filter = filters[query_index_hash.first];
            //     const ColumnPtr & hash_column = query_index_hash.second;

            //     match_rows = maybeTrueOnSuccinctRangeFilter(&*hash_column,
            //                                         filter,
            //                                         hash_functions,
            //                                         match_all);
            // }

            // rpn_stack.emplace_back(match_rows, true);
            // if (element.function == RPNElement::FUNCTION_NOT_EQUALS || element.function == RPNElement::FUNCTION_NOT_IN)
            //     rpn_stack.back() = !rpn_stack.back();

        }
        else if (element.function == RPNElement::FUNCTION_NOT)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ FUNCTION_NOT ------------------------------------------------------ ");

//            rpn_stack.back() = !rpn_stack.back();
        }
        else if (element.function == RPNElement::FUNCTION_OR)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ FUNCTION_OR ------------------------------------------------------ ");

            // auto arg1 = rpn_stack.back();
            // rpn_stack.pop_back();
            // auto arg2 = rpn_stack.back();
            // rpn_stack.back() = arg1 | arg2;
        }
        else if (element.function == RPNElement::FUNCTION_AND)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ FUNCTION_AND ------------------------------------------------------ ");

            // auto arg1 = rpn_stack.back();
            // rpn_stack.pop_back();
            // auto arg2 = rpn_stack.back();
            // rpn_stack.back() = arg1 & arg2;
        }
        else if (element.function == RPNElement::ALWAYS_TRUE)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ ALWAYS_TRUE ------------------------------------------------------ ");

            // rpn_stack.emplace_back(true, false);
        }
        else if (element.function == RPNElement::ALWAYS_FALSE)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ ALWAYS_FALSE ------------------------------------------------------ ");

            //rpn_stack.emplace_back(false, true);
        }
        else
            throw Exception(ErrorCodes::LOGICAL_ERROR, "Unexpected function type in KeyCondition::RPNElement");
    }

    if (rpn_stack.size() != 1)
        throw Exception(ErrorCodes::LOGICAL_ERROR, "Unexpected stack size in KeyCondition::mayBeTrueInRange");

    return rpn_stack[0].can_be_true;
}

bool MergeTreeIndexConditionSuccinctRangeFilter::extractAtomFromTree(const RPNBuilderTreeNode & node, RPNElement & out)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "extractAtomFromTree");

    {
        Field const_value;
        DataTypePtr const_type;

        if (node.tryGetConstant(const_value, const_type))
        {
            if (const_value.getType() == Field::Types::UInt64)
            {
                out.function = const_value.safeGet<UInt64>() ? RPNElement::ALWAYS_TRUE : RPNElement::ALWAYS_FALSE;
                LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " Type: UInt64 ");
                return true;
            }

            if (const_value.getType() == Field::Types::Int64)
            {
                out.function = const_value.safeGet<Int64>() ? RPNElement::ALWAYS_TRUE : RPNElement::ALWAYS_FALSE;
                LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " Type: Int64 ");
                return true;
            }

            if (const_value.getType() == Field::Types::Float64)
            {
                out.function = const_value.safeGet<Float64>() != 0.0 ? RPNElement::ALWAYS_TRUE : RPNElement::ALWAYS_FALSE;
                LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " Type: Float64 ");
                return true;
            }
        }
    }

    auto tf = traverseFunction(node, out, nullptr /*parent*/);
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "extractAtomFromTree {}", tf);
    return tf;
}

bool MergeTreeIndexConditionSuccinctRangeFilter::traverseFunction(const RPNBuilderTreeNode & node, RPNElement & out, const RPNBuilderTreeNode * parent)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "  traverseFunction  ");

    if (!node.isFunction())
        return false;

    const auto function = node.toFunctionNode();
    auto arguments_size = function.getArgumentsSize();
    auto function_name = function.getFunctionName();

    // if (function_name == "less" || function_name == "greater") // This is just for testing
    // {
    //     return true;
    // }

    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "  function name: {} {} ", function_name, arguments_size);


    if (parent == nullptr)
    {
        /// Recurse a little bit for indexOf().
        for (size_t i = 0; i < arguments_size; ++i)
        {
            auto argument = function.getArgumentAt(i);
            if (traverseFunction(argument, out, &node))
                return true;
        }
    }

    if (arguments_size != 2)
        return false;

    /// indexOf() should be inside comparison function, e.g. greater(indexOf(key, 42), 0).
    /// Other conditions should be at top level, e.g. equals(key, 42), not equals(equals(key, 42), 1).
    if ((function_name == "indexOf") != (parent != nullptr))
        return false;

    auto lhs_argument = function.getArgumentAt(0);
    auto rhs_argument = function.getArgumentAt(1);

    if (functionIsInOrGlobalInOperator(function_name))
    {
        if (auto future_set = rhs_argument.tryGetPreparedSet(); future_set)
        {
            if (auto prepared_set = future_set->buildOrderedSetInplace(rhs_argument.getTreeContext().getQueryContext()); prepared_set)
            {
                if (prepared_set->hasExplicitSetElements())
                {
                    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "traverseTreeIn");
                    const auto prepared_info = getPreparedSetInfo(prepared_set);
                    if (traverseTreeIn(function_name, lhs_argument, prepared_set, prepared_info.type, prepared_info.column, out))
                        return true;
                }
            }
        }
        return false;
    }

    if (function_name == "equals" ||
        function_name == "notEquals" ||
        function_name == "has" ||
        function_name == "mapContains" ||
        function_name == "indexOf" ||
        function_name == "hasAny" ||
        function_name == "hasAll" ||
        function_name == "less" ||
        function_name == "greater" ||
        function_name == "greaterOrEquals" ||
        function_name == "lessOrEquals")
    {
        Field const_value;
        DataTypePtr const_type;

        if (rhs_argument.tryGetConstant(const_value, const_type))
        {
            if (traverseTreeEquals(function_name, lhs_argument, const_type, const_value, out/*, parent*/))
                return true;
        }
        else if (lhs_argument.tryGetConstant(const_value, const_type) && (function_name == "equals" || function_name == "notEquals"))
        {
            if (traverseTreeEquals(function_name, rhs_argument, const_type, const_value, out/*, parent*/))
                return true;
        }

        return false;
    }

    return false;
}

bool MergeTreeIndexConditionSuccinctRangeFilter::traverseTreeIn(
    const String & function_name,
    const RPNBuilderTreeNode & key_node,
    const ConstSetPtr & prepared_set,
    const DataTypePtr & type,
    const ColumnPtr & column,
    RPNElement & out)
{
    auto key_node_column_name = key_node.getColumnName();

    if (header.has(key_node_column_name))
    {
        size_t row_size = column->size();
        size_t position = header.getPositionByName(key_node_column_name);
        const DataTypePtr & index_type = header.getByPosition(position).type;
        const auto & converted_column = castColumn(ColumnWithTypeAndName{column, type, ""}, index_type);
        out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithColumn(index_type, converted_column, 0, row_size)));

        if (function_name == "in"  || function_name == "globalIn")
            out.function = RPNElement::FUNCTION_IN;

        if (function_name == "notIn"  || function_name == "globalNotIn")
            out.function = RPNElement::FUNCTION_NOT_IN;

        return true;
    }

    if (key_node.isFunction())
    {
        auto key_node_function = key_node.toFunctionNode();
        auto key_node_function_name = key_node_function.getFunctionName();
        size_t key_node_function_arguments_size = key_node_function.getArgumentsSize();

        WhichDataType which(type);

        if (which.isTuple() && key_node_function_name == "tuple")
        {
            const auto & tuple_column = typeid_cast<const ColumnTuple *>(column.get());
            const auto & tuple_data_type = typeid_cast<const DataTypeTuple *>(type.get());

            if (tuple_data_type->getElements().size() != key_node_function_arguments_size || tuple_column->getColumns().size() != key_node_function_arguments_size)
                throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "Illegal types of arguments of function {}", function_name);

            bool match_with_subtype = false;
            const auto & sub_columns = tuple_column->getColumns();
            const auto & sub_data_types = tuple_data_type->getElements();

            for (size_t index = 0; index < key_node_function_arguments_size; ++index)
                match_with_subtype |= traverseTreeIn(function_name, key_node_function.getArgumentAt(index), nullptr, sub_data_types[index], sub_columns[index], out);

            return match_with_subtype;
        }

        if (key_node_function_name == "arrayElement")
        {
            /** Try to parse arrayElement for mapKeys index.
              * It is important to ignore keys like column_map['Key'] IN ('') because if the key does not exist in the map
              * we return the default value for arrayElement.
              *
              * We cannot skip keys that does not exist in map if comparison is with default type value because
              * that way we skip necessary granules where the map key does not exist.
              */
            if (!prepared_set)
                return false;

            auto default_column_to_check = type->createColumnConstWithDefaultValue(1)->convertToFullColumnIfConst();
            ColumnWithTypeAndName default_column_with_type_to_check { default_column_to_check, type, "" };
            ColumnsWithTypeAndName default_columns_with_type_to_check = {default_column_with_type_to_check};
            auto set_contains_default_value_predicate_column = prepared_set->execute(default_columns_with_type_to_check, false /*negative*/);
            const auto & set_contains_default_value_predicate_column_typed = assert_cast<const ColumnUInt8 &>(*set_contains_default_value_predicate_column);
            bool set_contain_default_value = set_contains_default_value_predicate_column_typed.getData()[0];
            if (set_contain_default_value)
                return false;

            auto first_argument = key_node_function.getArgumentAt(0);
            const auto column_name = first_argument.getColumnName();
            auto map_keys_index_column_name = fmt::format("mapKeys({})", column_name);
            auto map_values_index_column_name = fmt::format("mapValues({})", column_name);

            if (header.has(map_keys_index_column_name))
            {
                /// For mapKeys we serialize key argument with bloom filter

                auto second_argument = key_node_function.getArgumentAt(1);

                Field constant_value;
                DataTypePtr constant_type;

                if (second_argument.tryGetConstant(constant_value, constant_type))
                {
                    size_t position = header.getPositionByName(map_keys_index_column_name);
                    const DataTypePtr & index_type = header.getByPosition(position).type;
                    const DataTypePtr actual_type = BloomFilter::getPrimitiveType(index_type);
                    out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithField(actual_type.get(), constant_value)));
                }
                else
                {
                    return false;
                }
            }
            else if (header.has(map_values_index_column_name))
            {
                /// For mapValues we serialize set with bloom filter

                size_t row_size = column->size();
                size_t position = header.getPositionByName(map_values_index_column_name);
                const DataTypePtr & index_type = header.getByPosition(position).type;
                const auto & array_type = assert_cast<const DataTypeArray &>(*index_type);
                const auto & array_nested_type = array_type.getNestedType();
                const auto & converted_column = castColumn(ColumnWithTypeAndName{column, type, ""}, array_nested_type);
                out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithColumn(array_nested_type, converted_column, 0, row_size)));
            }
            else
            {
                return false;
            }

            if (function_name == "in"  || function_name == "globalIn")
                out.function = RPNElement::FUNCTION_IN;

            if (function_name == "notIn"  || function_name == "globalNotIn")
                out.function = RPNElement::FUNCTION_NOT_IN;

            return true;
        }
    }

    return false;
}

bool MergeTreeIndexConditionSuccinctRangeFilter::traverseTreeEquals(
    const String & function_name,
    const RPNBuilderTreeNode & key_node,
    const DataTypePtr & value_type,
    const Field & value_field,
    RPNElement & out/*,
    const RPNBuilderTreeNode * parent*/)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeEquals {} ", function_name);

    auto key_column_name = key_node.getColumnName();

    if (header.has(key_column_name))
    {
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeEquals header has column name {} ", key_column_name);
        size_t position = header.getPositionByName(key_column_name);
        const DataTypePtr & index_type = header.getByPosition(position).type;
        const auto * array_type = typeid_cast<const DataTypeArray *>(index_type.get());

        if (function_name == "has" || function_name == "indexOf")
        {
            if (!array_type)
                throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "First argument for function {} must be an array.", function_name);

            /// We can treat `indexOf` function similar to `has`.
            /// But it is little more cumbersome, compare: `has(arr, elem)` and `indexOf(arr, elem) != 0`.
            /// The `parent` in this context is expected to be function `!=` (`notEquals`).
            if (function_name == "has"/* || indexOfCanUseBloomFilter(parent)*/)
            {
                out.function = RPNElement::FUNCTION_HAS;
                const DataTypePtr actual_type = BloomFilter::getPrimitiveType(array_type->getNestedType());
                auto converted_field = convertFieldToType(value_field, *actual_type, value_type.get());
                if (converted_field.isNull())
                    return false;

                out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithField(actual_type.get(), converted_field)));
            }
        }
        else if (function_name == "hasAny" || function_name == "hasAll")
        {
            if (!array_type)
                throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "First argument for function {} must be an array.", function_name);

            if (value_field.getType() != Field::Types::Array)
                throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "Second argument for function {} must be an array.", function_name);

            const DataTypePtr actual_type = BloomFilter::getPrimitiveType(array_type->getNestedType());
            ColumnPtr column;
            {
                const bool is_nullable = actual_type->isNullable();
                auto mutable_column = actual_type->createColumn();

                for (const auto & f : value_field.safeGet<Array>())
                {
                    if ((f.isNull() && !is_nullable) || f.isDecimal(f.getType())) /// NOLINT(readability-static-accessed-through-instance)
                        return false;

                    auto converted = convertFieldToType(f, *actual_type);
                    if (converted.isNull())
                        return false;

                    mutable_column->insert(converted);
                }

                column = std::move(mutable_column);
            }

            out.function = function_name == "hasAny" ?
                RPNElement::FUNCTION_HAS_ANY :
                RPNElement::FUNCTION_HAS_ALL;
            out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithColumn(actual_type, column, 0, column->size())));
        }
        else if (function_name == "greater" || function_name == "less")
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "greater or less");
            out.function = function_name == "greater" ? RPNElement::FUNCTION_GREATER : RPNElement::FUNCTION_LESS;
            // out.predicate.emplace_back(std::make_pair(position, value_field));
            return true;
        }
        else if (function_name == "greaterOrEquals" || function_name == "lessOrEquals")
        {
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "between");
            out.function = function_name == "greaterOrEquals" ? RPNElement::FUNCTION_GREATER_OR_EQUALS : RPNElement::FUNCTION_LESS_OR_EQUALS;
            // out.predicate.emplace_back(std::make_pair(position, value_field));
            return true;
        }
        else
        {
            if (array_type)
                throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT,
                                "An array type of surf supports only has(), indexOf(), and hasAny() functions.");

            out.function = function_name == "equals" ? RPNElement::FUNCTION_EQUALS : RPNElement::FUNCTION_NOT_EQUALS;
            const DataTypePtr actual_type = BloomFilter::getPrimitiveType(index_type);
            auto converted_field = convertFieldToType(value_field, *actual_type, value_type.get());
            if (converted_field.isNull())
                return false;

            out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithField(actual_type.get(), converted_field)));
        }

        return true;
    }

    if (function_name == "mapContains" || function_name == "has")
    {
        auto map_keys_index_column_name = fmt::format("mapKeys({})", key_column_name);
        if (!header.has(map_keys_index_column_name))
            return false;

        size_t position = header.getPositionByName(map_keys_index_column_name);
        const DataTypePtr & index_type = header.getByPosition(position).type;
        const auto * array_type = typeid_cast<const DataTypeArray *>(index_type.get());

        if (!array_type)
            return false;

        out.function = RPNElement::FUNCTION_HAS;
        const DataTypePtr actual_type = BloomFilter::getPrimitiveType(array_type->getNestedType());
        auto converted_field = convertFieldToType(value_field, *actual_type, value_type.get());
        if (converted_field.isNull())
            return false;

        out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithField(actual_type.get(), converted_field)));
        return true;
    }

    if (key_node.isFunction())
    {
        WhichDataType which(value_type);

        auto key_node_function = key_node.toFunctionNode();
        auto key_node_function_name = key_node_function.getFunctionName();
        size_t key_node_function_arguments_size = key_node_function.getArgumentsSize();

        if (which.isTuple() && key_node_function_name == "tuple")
        {
            const Tuple & tuple = value_field.safeGet<const Tuple &>();
            const auto * value_tuple_data_type = typeid_cast<const DataTypeTuple *>(value_type.get());

            if (tuple.size() != key_node_function_arguments_size)
                throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "Illegal types of arguments of function {}", function_name);

            bool match_with_subtype = false;
            const DataTypes & subtypes = value_tuple_data_type->getElements();

            for (size_t index = 0; index < tuple.size(); ++index)
                match_with_subtype |= traverseTreeEquals(function_name, key_node_function.getArgumentAt(index), subtypes[index], tuple[index], out/*, &key_node*/);

            return match_with_subtype;
        }

        if (key_node_function_name == "arrayElement" && (function_name == "equals" || function_name == "notEquals"))
        {
            /** Try to parse arrayElement for mapKeys index.
              * It is important to ignore keys like column_map['Key'] = '' because if key does not exist in the map
              * we return default the value for arrayElement.
              *
              * We cannot skip keys that does not exist in map if comparison is with default type value because
              * that way we skip necessary granules where map key does not exist.
              */
            if (value_field == value_type->getDefault())
                return false;

            auto first_argument = key_node_function.getArgumentAt(0);
            const auto column_name = first_argument.getColumnName();

            auto map_keys_index_column_name = fmt::format("mapKeys({})", column_name);
            auto map_values_index_column_name = fmt::format("mapValues({})", column_name);

            size_t position = 0;
            Field const_value = value_field;
            DataTypePtr const_type;

            if (header.has(map_keys_index_column_name))
            {
                position = header.getPositionByName(map_keys_index_column_name);
                auto second_argument = key_node_function.getArgumentAt(1);

                if (!second_argument.tryGetConstant(const_value, const_type))
                    return false;
            }
            else if (header.has(map_values_index_column_name))
            {
                position = header.getPositionByName(map_values_index_column_name);
            }
            else
            {
                return false;
            }

            out.function = function_name == "equals" ? RPNElement::FUNCTION_EQUALS : RPNElement::FUNCTION_NOT_EQUALS;

            const auto & index_type = header.getByPosition(position).type;
            const auto actual_type = BloomFilter::getPrimitiveType(index_type);
            out.predicate.emplace_back(std::make_pair(position, BloomFilterHash::hashWithField(actual_type.get(), const_value)));

            return true;
        }
    }

    return false;
}

MergeTreeIndexSuccinctRangeFilter::MergeTreeIndexSuccinctRangeFilter(
    const IndexDescription & index_,
    size_t ds_ratio_,
    size_t num_columns_)
    : IMergeTreeIndex(index_)
    , num_columns(num_columns_)
    , ds_ratio(ds_ratio_)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexSuccinctRangeFilter {}", num_columns);
    assert(num_columns != 0);
    assert(ds_ratio != 0);
}

MergeTreeIndexGranulePtr MergeTreeIndexSuccinctRangeFilter::createIndexGranule() const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "createIndexGranule {}", num_columns);
    return std::make_shared<MergeTreeIndexGranuleSuccinctRangeFilter>(ds_ratio, index.column_names.size());
}

MergeTreeIndexAggregatorPtr MergeTreeIndexSuccinctRangeFilter::createIndexAggregator(const MergeTreeWriterSettings & /*settings*/) const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "createIndexAggregator {}", num_columns);
    return std::make_shared<MergeTreeIndexAggregatorSuccinctRangeFilter>(ds_ratio, index.column_names);
}

MergeTreeIndexConditionPtr MergeTreeIndexSuccinctRangeFilter::createIndexCondition(const ActionsDAG * filter_actions_dag, ContextPtr context) const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "createIndexCondition {}", num_columns);
    return std::make_shared<MergeTreeIndexConditionSuccinctRangeFilter>(filter_actions_dag, context, index.sample_block);
}

void MergeTreeIndexAggregatorSuccinctRangeFilter::update(const Block & block, size_t * pos, size_t limit) // This would work much better if the granules were in alphabetical order
{
    if (*pos >= block.rows())
        throw Exception(ErrorCodes::LOGICAL_ERROR, "The provided position is not less than the number of block rows. "
                        "Position: {}, Block rows: {}.", *pos, block.rows());

    Block granule_index_block;
    size_t max_read_rows = std::min(block.rows() - *pos, limit);
    
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "update {} {} {}", *pos, limit, block.rows());

    for (size_t column = 0; column < index_columns_name.size(); ++column)
    {
        const auto & column_and_type = block.getByName(index_columns_name[column]);

        // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "column name {}", index_columns_name[column]);
        // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "column type {}", column_and_type.type->getName());

        // for (size_t i = *pos; i < *pos + 2; ++i)
        // {
        //     Field res;
        //     column_and_type.column->get(i, res);
        //     LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "start {} size {}", toString(res), toString(res).size());
        // }
        // for (size_t i = *pos + max_read_rows - 1; i < *pos + max_read_rows + 1; ++i)
        // {
        //     Field res;
        //     column_and_type.column->get(i, res);
        //     LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "end {} size {}", toString(res), toString(res).size());
        // }

        // TrieNode node = root;
        for (size_t i = *pos; i < *pos + max_read_rows; ++i)
        {
            Field res;
            column_and_type.column->get(i, res); // Probably a better way of doing this

            const String & key = toString(res);

            TrieNode *current = &root;
            for (char c : key)
            {
                if (current->children.find(c) == current->children.end())
                {
                    current->children[c] = std::make_unique<TrieNode>();
                }
                // Descend into the child
                current = current->children[c].get();
            }

            current->is_terminal = true;
        }
    }

    *pos += max_read_rows;
    total_rows += max_read_rows;

}

// Create an empty trie (root node)
TrieNode createEmptyTrie()
{
    return TrieNode();
}


MergeTreeIndexAggregatorSuccinctRangeFilter::MergeTreeIndexAggregatorSuccinctRangeFilter(size_t ds_ratio_, const Names & columns_name_)
    : ds_ratio(ds_ratio_), index_columns_name(columns_name_)
{
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexAggregatorSuccinctRangeFilter");
    root = createEmptyTrie();
    assert(ds_ratio != 0);
}

bool MergeTreeIndexAggregatorSuccinctRangeFilter::empty() const
{
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "aggregator empty");
    return root.children.size() == 0;
}

MergeTreeIndexGranulePtr MergeTreeIndexAggregatorSuccinctRangeFilter::getGranuleAndReset()
{
    // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "getGranuleAndReset");
    const auto granule = std::make_shared<MergeTreeIndexGranuleSuccinctRangeFilter>(ds_ratio, root, total_rows);
    total_rows = 0;
    root.children.clear();
    root.is_terminal = false;
    return granule;
}

MergeTreeIndexPtr succinctRangeFilterIndexCreator(
    const IndexDescription & index)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "succinctRangeFilterIndexCreator assuming one column for now");
    size_t default_ds_ratio = 64;

    if (!index.arguments.empty())
    {
        const auto & argument = index.arguments[0];
        default_ds_ratio = std::min<size_t>(1024, std::max<size_t>(argument.safeGet<size_t>(), 0));
    }

    return std::make_shared<MergeTreeIndexSuccinctRangeFilter>(
        index, default_ds_ratio, 1);    /// We assume that the index is created for one column.
}

void succinctRangeFilterIndexValidator(const IndexDescription & index, bool attach)
{
    // assertIndexColumnsType(index.sample_block);

    // if (index.arguments.size() > 1)
    // {
    //     if (!attach) /// This is for backward compatibility.
    //         throw Exception(ErrorCodes::NUMBER_OF_ARGUMENTS_DOESNT_MATCH, "BloomFilter index cannot have more than one parameter.");
    // }

    // if (!index.arguments.empty())
    // {
    //     const auto & argument = index.arguments[0];

    //     if (!attach && (argument.getType() != Field::Types::Float64 || argument.safeGet<Float64>() < 0 || argument.safeGet<Float64>() > 1))
    //         throw Exception(ErrorCodes::BAD_ARGUMENTS, "The BloomFilter false positive must be a double number between 0 and 1.");
    // }
    if (attach)
        LOG_DEBUG(getLogger("succinctRangeFilterIndexValidator"), "for now we assume validity (attach) {}", index.arguments.size());
    else
        LOG_DEBUG(getLogger("succinctRangeFilterIndexValidator"), "for now we assume validity (not attach) {}", index.arguments.size());
}

}

