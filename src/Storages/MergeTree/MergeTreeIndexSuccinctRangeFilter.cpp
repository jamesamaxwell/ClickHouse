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

MergeTreeIndexGranuleSuccinctRangeFilter::MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, const TrieNode & root_)
    : ds_ratio(ds_ratio_)
{
    num_columns = 0;
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter {}", root_.is_terminal);
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "ds_ratio {}", ds_ratio);
}

MergeTreeIndexGranuleSuccinctRangeFilter::MergeTreeIndexGranuleSuccinctRangeFilter(size_t ds_ratio_, size_t num_columns_)
    : ds_ratio(ds_ratio_), num_columns(num_columns_)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexGranuleSuccinctRangeFilter {}", num_columns);
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "ds_ratio {}", ds_ratio);
}

bool MergeTreeIndexGranuleSuccinctRangeFilter::empty() const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "empty granule");
    /// One possible definition:
    /// Return true if we haven't added any rows yet, or if no filters exist.
    return surfs.empty();
}

void MergeTreeIndexGranuleSuccinctRangeFilter::serializeBinary(WriteBuffer & ostr) const
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "serializeBinary");
    if (!surfs.empty())
        ostr.write(reinterpret_cast<const char *>(surfs[0]->getFilter().data()), 4);
}

void MergeTreeIndexGranuleSuccinctRangeFilter::deserializeBinary(ReadBuffer & istr, MergeTreeIndexVersion version)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "deserializeBinary");
    if (version != 1)
        throw Exception(ErrorCodes::LOGICAL_ERROR, "Unknown index version {}.", version);

    if (!surfs.empty())
        istr.readStrict(reinterpret_cast<char *>(surfs[0]->getFilter().data()), 4);
}

void fillingSuccinctRangeFilter(SuccinctRangeFilterPtr & surf)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "fillingSuccinctRangeFilter");
    surf->add("test", 5);
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
    rpn_stack.push_back(true);

    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "alwaysUnknownOrTrue");
    return rpn_stack[0];
}

bool MergeTreeIndexConditionSuccinctRangeFilter::mayBeTrueOnGranule(const MergeTreeIndexGranuleSuccinctRangeFilter * granule) const
{
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
            || element.function == RPNElement::FUNCTION_HAS_ALL)
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
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ extractAtomFromTree ------------------------------------------------------ ");

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
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " ------------------------------------------------------ traverseFunction ------------------------------------------------------ ");

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
                    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "traverseTreeIn");
                    // const auto prepared_info = getPreparedSetInfo(prepared_set);
                    // if (traverseTreeIn(function_name, lhs_argument, prepared_set, prepared_info.type, prepared_info.column, out))
                    //     return true;
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
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeIn {} ", function_name);
    auto key_node_column_name = key_node.getColumnName();
    if (header.has(key_node_column_name))
    {
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeIn ");
    }
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

    size_t row_size = column->size();
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeIn row size {} ", row_size);
    out.function = RPNElement::FUNCTION_IN;

    // if (value_field.getType() != Field::Types::Array)
    //     throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "Second argument for function {} must be an array.", function_name);

    return true;
}

bool MergeTreeIndexConditionSuccinctRangeFilter::traverseTreeEquals(
    const String & function_name,
    const RPNBuilderTreeNode & key_node,
    const DataTypePtr & value_type,
    const Field & value_field,
    RPNElement & out,
    const RPNBuilderTreeNode * parent)
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeEquals ");
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeEquals {} ", function_name);
    auto key_node_column_name = key_node.getColumnName();
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), " traverseTreeIn column name {} ", key_node_column_name);

    out.function = RPNElement::FUNCTION_HAS;
    const DataTypePtr actual_type = value_type;
    auto converted_field = convertFieldToType(value_field, *actual_type, value_type.get());
    if (converted_field.isNull())
        return false;

    if (value_field.getType() != Field::Types::Array)
        throw Exception(ErrorCodes::ILLEGAL_TYPE_OF_ARGUMENT, "Second argument for function {} must be an array.", function_name);

    out.function = RPNElement::FUNCTION_IN;

    if (!parent->isFunction())
        return false;

    return true;
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

void MergeTreeIndexAggregatorSuccinctRangeFilter::update(const Block & block, size_t * pos, size_t limit)
{
    if (*pos >= block.rows())
        throw Exception(ErrorCodes::LOGICAL_ERROR, "The provided position is not less than the number of block rows. "
                        "Position: {}, Block rows: {}.", *pos, block.rows());

    Block granule_index_block;
    size_t max_read_rows = std::min(block.rows() - *pos, limit);
    
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "update {} {} {}", pos[0], limit, block.rows());

    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "index_columns_name.size() {}", index_columns_name.size());
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "ds_ratio {}", ds_ratio);

    // size_t i = 0;
    // for (auto c : block)
    // {
    //     Field res;
    //     c.column->get(c.column->size() - 1, res);
    //     LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "string {} size {}", i, toString(res).size());
    //     ++i;
    //     if (i > 10)
    //         break;
    // }

    for (size_t column = 0; column < index_columns_name.size(); ++column)
    {
        const auto & column_and_type = block.getByName(index_columns_name[column]);

        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "column name {}", index_columns_name[column]);
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "column type {}", column_and_type.type->getName());

        for (size_t i = 0; i < 25; ++i)
        {
            Field res;
            column_and_type.column->get(i, res);
            LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "string {} size {}", toString(res), toString(res).size());
        }

        // TrieNode node = root;
        for (size_t i = 0; i < max_read_rows; ++i)
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
        // Field res;
        // column_and_type.column->get(column_and_type.column->size() - 1, res);
        // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "column {}", toString(res));
        // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "column size {}", column_and_type.column->size());
        // LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "string size {}", toString(res).size());

        // const auto & index_col = checkAndGetColumn<ColumnUInt64>(*index_column);
        // const auto & index_data = index_col.getData();
        // for (const auto & hash: index_data)
        //     column_hashes[column].insert(hash);
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
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "MergeTreeIndexAggregatorSuccinctRangeFilter");
    root = createEmptyTrie();
    assert(ds_ratio != 0);
}

bool MergeTreeIndexAggregatorSuccinctRangeFilter::empty() const
{
    if (!total_rows)
    {
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "empty");
    }
    else
    {
        LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "not empty");
    }
    return !total_rows;
}

MergeTreeIndexGranulePtr MergeTreeIndexAggregatorSuccinctRangeFilter::getGranuleAndReset()
{
    LOG_DEBUG(getLogger("MergeTreeIndexSuccinctRangeFilter"), "getGranuleAndReset");

    const auto granule = std::make_shared<MergeTreeIndexGranuleSuccinctRangeFilter>(ds_ratio, root);
    total_rows = 0;
    // column_hashes.clear();
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

