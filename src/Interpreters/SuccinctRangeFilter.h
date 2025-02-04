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

struct BFSItem
{
    TrieNode * node;
    size_t depth;
};

struct LOUDSDenseTrie // Change variable names a some point
{
    std::vector<std::bitset<256>> d_labels;
    std::vector<std::bitset<256>> d_hasChild;
    std::vector<bool> d_isPrefixKey;
    std::vector<uint64_t> d_values;
};

struct LOUDSSparseTrie // Change variable names a some point
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


    // using UnderType = UInt64;
    // using Container = std::vector<UnderType>;

    SuccinctRangeFilter(std::unique_ptr<TrieNode> root, size_t l_depth);

    // bool find(const char * data, size_t len);
    // void add(const char * data, size_t len);
    // void clear();

    const LOUDSdsTrie & getFilter() const { return surf; }
    LOUDSdsTrie & getFilter() { return surf; }

    size_t getWriteSize();

    /// For debug.
//    UInt64 isEmpty() const;

    // friend bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b);
private:

    LOUDSdsTrie surf;

// public:
//     static ColumnPtr getPrimitiveColumn(const ColumnPtr & column);
//     static DataTypePtr getPrimitiveType(const DataTypePtr & data_type);
};

using SuccinctRangeFilterPtr = std::shared_ptr<SuccinctRangeFilter>;

//bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b);


}
