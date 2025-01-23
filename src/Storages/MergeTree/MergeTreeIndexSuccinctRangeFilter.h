#pragma once

#include <Columns/IColumn.h>
#include <Common/HashTable/HashSet.h>
#include <Interpreters/SuccinctRangeFilter.h>
#include <Storages/MergeTree/KeyCondition.h>
#include <Storages/MergeTree/MergeTreeIndices.h>

namespace DB
{

namespace ErrorCodes
{
    extern const int LOGICAL_ERROR;
}

struct LoudsDense
{
    std::vector<UInt64> d_labels;
    std::vector<UInt64> d_has_child;
    std::vector<UInt64> d_is_prefix_tree;
};

struct LoudsSparse
{
    std::vector<UInt64> s_labels;
    std::vector<UInt64> s_has_child;
    std::vector<UInt64> s_louds;
};

struct TrieNode
{
    std::unordered_map<char, std::unique_ptr<TrieNode>> children; // Map of child nodes
    bool is_terminal = false;                                    // Marks if the node is the end of a valid key
};

class MergeTreeIndexGranuleSuccinctRangeFilter final : public IMergeTreeIndexGranule
{
public:
    MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, const TrieNode & root);

    MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, size_t num_columns_);

    // MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, const std::vector<HashSet<UInt64>> & column_hashes);

    bool empty() const override;

    size_t findLargestDepth(const std::unique_ptr<TrieNode> & root, size_t ratio);
    std::pair<std::unique_ptr<TrieNode>, size_t> pruneSubtree(const TrieNode & old_node);

    void serializeBinary(WriteBuffer & ostr) const override;
    void deserializeBinary(ReadBuffer & istr, MergeTreeIndexVersion version) override;

    const std::vector<SuccinctRangeFilterPtr> & getFilters() const { return surfs; }

private:
    const size_t ds_ratio;
    size_t num_columns;
    size_t dense_depth;

    std::vector<SuccinctRangeFilterPtr> surfs;

    void fillingSuccinctRangeFilter(SuccinctRangeFilterPtr & surf) const;
};

class MergeTreeIndexConditionSuccinctRangeFilter final : public IMergeTreeIndexCondition, WithContext
{
public:
    struct RPNElement
    {
        enum Function
        {
            /// Atoms of a Boolean expression.
            // FUNCTION_GREATER,
            // FUNCTION_LESS,
            // FUNCTION_GREATER_OR_EQUALS,
            // FUNCTION_LESS_OR_EQUALS,
            FUNCTION_EQUALS, // Can I delete these?
            FUNCTION_NOT_EQUALS,
            FUNCTION_HAS,
            FUNCTION_HAS_ANY,
            FUNCTION_HAS_ALL,
            FUNCTION_IN,
            FUNCTION_NOT_IN,
            FUNCTION_UNKNOWN, /// Can take any value.
            /// Operators of the logical expression.
            FUNCTION_NOT,
            FUNCTION_AND,
            FUNCTION_OR,
            /// Constants
            ALWAYS_FALSE,
            ALWAYS_TRUE,
        };

        RPNElement(Function function_ = FUNCTION_UNKNOWN) : function(function_) {} /// NOLINT

        Function function = FUNCTION_UNKNOWN;
        std::vector<std::pair<size_t, ColumnPtr>> predicate;
    };

    MergeTreeIndexConditionSuccinctRangeFilter(const ActionsDAG * filter_actions_dag, ContextPtr context_, const Block & header_);

    bool alwaysUnknownOrTrue() const override;

    bool mayBeTrueOnGranule(MergeTreeIndexGranulePtr granule) const override
    {
        if (const auto & surf_granule = typeid_cast<const MergeTreeIndexGranuleSuccinctRangeFilter *>(granule.get()))
            return mayBeTrueOnGranule(surf_granule);

        throw Exception(ErrorCodes::LOGICAL_ERROR, "Requires succinct range filter index granule.");
    }

private:
    const Block & header;
    // const size_t hash_functions;
    std::vector<RPNElement> rpn;

    bool mayBeTrueOnGranule(const MergeTreeIndexGranuleSuccinctRangeFilter * granule) const;

    bool extractAtomFromTree(const RPNBuilderTreeNode & node, RPNElement & out);

    bool traverseFunction(const RPNBuilderTreeNode & node, RPNElement & out, const RPNBuilderTreeNode * parent);

    bool traverseTreeIn(
        const String & function_name,
        const RPNBuilderTreeNode & key_node,
        const ConstSetPtr & prepared_set,
        const DataTypePtr & type,
        const ColumnPtr & column,
        RPNElement & out);

    bool traverseTreeEquals(
        const String & function_name,
        const RPNBuilderTreeNode & key_node,
        const DataTypePtr & value_type,
        const Field & value_field,
        RPNElement & out,
        const RPNBuilderTreeNode * parent);
};

class MergeTreeIndexAggregatorSuccinctRangeFilter final : public IMergeTreeIndexAggregator
{
public:
    MergeTreeIndexAggregatorSuccinctRangeFilter(size_t ds_ratio, const Names & columns_name_);

    bool empty() const override;

    MergeTreeIndexGranulePtr getGranuleAndReset() override;

    void update(const Block & block, size_t * pos, size_t limit) override;

private:
    size_t ds_ratio;
    const Names index_columns_name;

    TrieNode root;
    size_t total_rows = 0;
};


class MergeTreeIndexSuccinctRangeFilter final : public IMergeTreeIndex
{
public:
    MergeTreeIndexSuccinctRangeFilter(
        const IndexDescription & index_,
        size_t ds_ratio,
        size_t num_columns);

    MergeTreeIndexGranulePtr createIndexGranule() const override;

    MergeTreeIndexAggregatorPtr createIndexAggregator(const MergeTreeWriterSettings & settings) const override;

    MergeTreeIndexConditionPtr createIndexCondition(const ActionsDAG * filter_actions_dag, ContextPtr context) const override;

private:
    size_t num_columns;
    size_t ds_ratio;
};

}

