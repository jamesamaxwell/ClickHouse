#include <Interpreters/SuccinctRangeFilter.h>
#include <city.h>
#include <Columns/ColumnArray.h>
#include <Columns/ColumnNullable.h>
#include <Columns/ColumnLowCardinality.h>
#include <DataTypes/DataTypeArray.h>
#include <DataTypes/DataTypeNullable.h>
#include <DataTypes/DataTypeLowCardinality.h>

#include <Common/logger_useful.h>

namespace DB
{
namespace ErrorCodes
{
    extern const int BAD_ARGUMENTS;
}

SuccinctRangeFilter::SuccinctRangeFilter(size_t ds_ratio_)
{
    ds_ratio = ds_ratio_;
}

bool operator== (const SuccinctRangeFilter & a, const SuccinctRangeFilter & b)
{
    if (a.ds_ratio != b.ds_ratio)
        return false;
    else
        return true;
    
}

void SuccinctRangeFilter::add(const char * data, size_t len)
{
    LOG_DEBUG(getLogger("SuccinctRangeFilter"), "add {} {}", data, len);

}

bool SuccinctRangeFilter::find(const char * data, size_t len)
{
    LOG_DEBUG(getLogger("SuccinctRangeFilter"), "find {} {}", data, len);
    if (len > 3)
    {
        return true;
    }
    else
    {
        return false;
    }

}

// UInt64 isEmpty() const
// {
//     return 0;

// }

}
