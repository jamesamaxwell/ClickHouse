// #pragma once

// #include <vector>
// #include <base/types.h>
// #include <Core/Field.h>
// #include <Common/PODArray.h>
// #include <Common/Allocator.h>
// #include <Columns/IColumn.h>
// #include <Columns/ColumnVector.h>
// #include <DataTypes/IDataType.h>


// namespace DB
// {
// // struct SuccinctRangeFilterParameters
// // {
// //     SuccinctRangeFilterParameters(size_t ds_ratio);

// //     /// LOUDS-DENSE to LOUDS-SPARSE size ratio, i.e. LOUDS_DENSE_size * ds_ratio = LOUDS_SPARSE_size. 64 by default
// //     size_t ds_ratio;

// // };

// struct LoudsDense
// {
//     std::vector<UnderType> d_labels;
//     std::vector<UnderType> d_has_child;
//     std::vector<UnderType> d_is_prefix_tree;
// };

// struct LoudsSparse
// {
//     std::vector<UnderType> s_labels;
//     std::vector<UnderType> s_has_child;
//     std::vector<UnderType> s_louds;
// };

// class SuccinctRangeFilter
// {

// public:
//     using UnderType = UInt64;
//     using Container = std::vector<UnderType>;

//     SuccinctRangeFilter(size_t size_);

//     void resize(size_t size_);
//     bool find(const char * data, size_t len);
//     void add(const char * data, size_t len);
//     void clear();

//     const Container & getFilter() const { return filter; }
//     Container & getFilter() { return filter; }

//     /// For debug.
//     UInt64 isEmpty() const;

//     friend bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b);
// private:

//     size_t ds_ratio;
//     Container filter;

// public:
//     static ColumnPtr getPrimitiveColumn(const ColumnPtr & column);
//     static DataTypePtr getPrimitiveType(const DataTypePtr & data_type);
// };

// using SuccinctRangeFilterPtr = std::shared_ptr<SuccinctRangeFilter>;

// bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b);


// }
