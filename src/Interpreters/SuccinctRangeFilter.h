#pragma once

#include <vector>
#include <bitset>
#include <base/types.h>
#include <Core/Field.h>
#include <Common/PODArray.h>
#include <Common/Allocator.h>
#include <Columns/IColumn.h>
#include <Columns/ColumnVector.h>
#include <DataTypes/IDataType.h>


namespace DB
{
// struct SuccinctRangeFilterParameters
// {
//     SuccinctRangeFilterParameters(size_t ds_ratio);

//     /// LOUDS-DENSE to LOUDS-SPARSE size ratio, i.e. LOUDS_DENSE_size * ds_ratio = LOUDS_SPARSE_size. 64 by default
//     size_t ds_ratio;

// };

struct TrieNode
{
    std::unordered_map<char, std::unique_ptr<TrieNode>> children; // Map of child nodes
    bool is_terminal = false;                                    // Marks if the node is the end of a valid key
};

struct Iterator {
    // Trace from root to leaf: one position per level.
    std::vector<size_t> levelPositions;
    // Current level (leaf level)
    size_t currentLevel;
    // In a full implementation, we’d also store which part (dense vs sparse) we are in.
    bool valid = true;
};

struct BFSItem
{
    TrieNode * node;
    size_t depth;
};

struct KeyItem
{
    TrieNode * node;
    std::vector<char> key;
};

struct LOUDSDenseTrie
{
    std::vector<std::bitset<256>> d_labels;
    std::vector<std::bitset<256>> d_hasChild;
    std::vector<bool> d_isPrefixKey;
    std::vector<uint64_t> d_values;
};

struct LOUDSSparseTrie
{
    std::vector<uint8_t> s_labels;
    std::vector<bool> s_hasChild;
    std::vector<bool> s_LOUDS;
    std::vector<uint64_t> s_values;
};

struct LOUDSdsTrie
{
    LOUDSDenseTrie dense;
    LOUDSSparseTrie sparse;
};

class SuccinctRangeFilter
{

public:

    SuccinctRangeFilter(std::unique_ptr<TrieNode> root, size_t l_depth_);

    const LOUDSdsTrie & getFilter() const { return surf; }
    LOUDSdsTrie & getFilter() { return surf; }

    size_t getWriteSize();

private:

    size_t childPositionDense(size_t currentPos, uint8_t label, size_t level) const;
    size_t getSparseNodeStart(size_t pos) const;
    size_t getSparseNodeEnd(size_t pos) const;
    size_t rankTrue(const std::vector<bool>& vec, size_t pos) const;
    size_t selectFalse(const std::vector<bool>& vec, size_t count) const;
    size_t childPositionSparse(size_t currentPos) const;
    size_t findLowerBoundInDense(/*size_t currentPos, */size_t level, uint8_t target) const;
    bool canMoveNext(const Iterator & iter, size_t level) const;
    size_t nextPosition(size_t pos, size_t level) const;
    bool hasChild(size_t pos, size_t level) const;
    size_t leftMostChild(size_t pos, size_t level) const;

    LOUDSdsTrie surf;
    size_t l_depth;

public:

    std::optional<uint64_t> ExactKeySearch(const std::string & key) const;
    Iterator LowerBound(const std::string & key) const;
    bool MoveToNext(Iterator & iter) const;

};

using SuccinctRangeFilterPtr = std::shared_ptr<SuccinctRangeFilter>;

//bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b);


}
