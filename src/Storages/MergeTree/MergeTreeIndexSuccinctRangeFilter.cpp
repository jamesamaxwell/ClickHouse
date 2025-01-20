#include "MergeTreeIndexSuccinctRangeFilter.h"

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

MergeTreeIndexGranuleSuccinctRangeFilter::MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_)
    : ds_ratio(ds_ratio_)
{
    total_rows = 0;
    for (size_t column = 0; column < index_columns_; ++column)
        surfs[column] = std::make_shared<SuccinctRangeFilter>(ds_ratio);
}

bool MergeTreeIndexGranuleSuccinctRangeFilter::empty() const
{
    /// One possible definition:
    /// Return true if we haven't added any rows yet, or if no filters exist.
    return (total_rows == 0) || surfs.empty();
}

void MergeTreeIndexGranuleSuccinctRangeFilter::serializeBinary(WriteBuffer & ostr) const
{
    return;
}

void MergeTreeIndexGranuleSuccinctRangeFilter::deserializeBinary(ReadBuffer & istr, MergeTreeIndexVersion version)
{
    return;
}

void fillingSuccinctRangeFilter(SuccinctRangeFilterPtr & surf)
{
    return;
}

MergeTreeIndexConditionSuccinctRangeFilter::MergeTreeIndexConditionSuccinctRangeFilter(const ActionsDAG * filter_actions_dag, ContextPtr context_, const Block & header_)
{
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
    rpn_stack.emplace_back(true, true);

    LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ alwaysUnknownOrTrue ------------------------------------------------------ ");
    return rpn_stack[0];
}

bool MergeTreeIndexConditionSuccinctRangeFilter::mayBeTrueOnGranule(const MergeTreeIndexGranuleSuccinctRangeFilter * granule) const
{
    std::vector<BoolMask> rpn_stack;
    const auto & filters = granule->getFilters();

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
            || element.function == RPNElement::FUNCTION_HAS_ALL)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ first case ------------------------------------------------------ ");
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
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ FUNCTION_NOT ------------------------------------------------------ ");

//            rpn_stack.back() = !rpn_stack.back();
        }
        else if (element.function == RPNElement::FUNCTION_OR)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ FUNCTION_OR ------------------------------------------------------ ");

            // auto arg1 = rpn_stack.back();
            // rpn_stack.pop_back();
            // auto arg2 = rpn_stack.back();
            // rpn_stack.back() = arg1 | arg2;
        }
        else if (element.function == RPNElement::FUNCTION_AND)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ FUNCTION_AND ------------------------------------------------------ ");

            // auto arg1 = rpn_stack.back();
            // rpn_stack.pop_back();
            // auto arg2 = rpn_stack.back();
            // rpn_stack.back() = arg1 & arg2;
        }
        else if (element.function == RPNElement::ALWAYS_TRUE)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ ALWAYS_TRUE ------------------------------------------------------ ");

            // rpn_stack.emplace_back(true, false);
        }
        else if (element.function == RPNElement::ALWAYS_FALSE)
        {
            LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ ALWAYS_FALSE ------------------------------------------------------ ");

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
    LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ extractAtomFromTree ------------------------------------------------------ ");

    {
        Field const_value;
        DataTypePtr const_type;

        if (node.tryGetConstant(const_value, const_type))
        {
            if (const_value.getType() == Field::Types::UInt64)
            {
                out.function = const_value.safeGet<UInt64>() ? RPNElement::ALWAYS_TRUE : RPNElement::ALWAYS_FALSE;
                return true;
            }

            if (const_value.getType() == Field::Types::Int64)
            {
                out.function = const_value.safeGet<Int64>() ? RPNElement::ALWAYS_TRUE : RPNElement::ALWAYS_FALSE;
                return true;
            }

            if (const_value.getType() == Field::Types::Float64)
            {
                out.function = const_value.safeGet<Float64>() != 0.0 ? RPNElement::ALWAYS_TRUE : RPNElement::ALWAYS_FALSE;
                return true;
            }
        }
    }

    return traverseFunction(node, out, nullptr /*parent*/);
}

bool MergeTreeIndexConditionSuccinctRangeFilter::traverseFunction(const RPNBuilderTreeNode & node, RPNElement & out, const RPNBuilderTreeNode * parent)
{
    LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ traverseFunction ------------------------------------------------------ ");

    if (!node.isFunction())
        return false;

    const auto function = node.toFunctionNode();
    auto arguments_size = function.getArgumentsSize();
    auto function_name = function.getFunctionName();

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
        function_name == "hasAll")
    {
        Field const_value;
        DataTypePtr const_type;

        if (rhs_argument.tryGetConstant(const_value, const_type))
        {
            if (traverseTreeEquals(function_name, lhs_argument, const_type, const_value, out, parent))
                return true;
        }
        else if (lhs_argument.tryGetConstant(const_value, const_type) && (function_name == "equals" || function_name == "notEquals"))
        {
            if (traverseTreeEquals(function_name, rhs_argument, const_type, const_value, out, parent))
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
    LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ traverseTreeIn ------------------------------------------------------ ");

    return false;
}

bool MergeTreeIndexConditionSuccinctRangeFilter::traverseTreeEquals(
    const String & function_name,
    const RPNBuilderTreeNode & key_node,
    const DataTypePtr & value_type,
    const Field & value_field,
    RPNElement & out,
    const RPNBuilderTreeNode * parent)
{
    LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ traverseTreeEquals ------------------------------------------------------ ");
    return true;
}

MergeTreeIndexAggregatorSuccinctRangeFilter::MergeTreeIndexAggregatorSuccinctRangeFilter(size_t ds_ratio_, const Names & columns_name_)
    : ds_ratio(ds_ratio_), index_columns_name(columns_name_), column_hashes(columns_name_.size())
{
    assert(ds_ratio != 0);
}

bool MergeTreeIndexAggregatorSuccinctRangeFilter::empty() const
{
    return !total_rows;
}

MergeTreeIndexGranulePtr MergeTreeIndexAggregatorSuccinctRangeFilter::getGranuleAndReset()
{
    LOG_DEBUG(getLogger("MergeTreeIndexConditionSuccinctRangeFilter"), " ------------------------------------------------------ getGranuleAndReset ------------------------------------------------------ ");

    const auto granule = std::make_shared<MergeTreeIndexAggregatorSuccinctRangeFilter>(bits_per_row, hash_functions, column_hashes);
    total_rows = 0;
    column_hashes.clear();
    return granule;
}

}