//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// index_test.cpp
//
// Identification: tests/index/index_test.cpp
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "gtest/gtest.h"
#include "harness.h"
#include "backend/common/logger.h"
#include "backend/index/index_factory.h"
#include "backend/index/index_key.h"
#include "backend/index/bwtree.h"
#include "backend/storage/tuple.h"

namespace peloton {
namespace test {

//===--------------------------------------------------------------------===//
// Index Tests
//===--------------------------------------------------------------------===//

#define TestKeyType index::GenericKey<256>
#define TestValueType ItemPointer
#define TestComparatorType index::GenericComparator<256>
#define TestEqualityChecker index::GenericEqualityChecker<256>

catalog::Schema *key_schema = nullptr;
catalog::Schema *tuple_schema = nullptr;

ItemPointer item0(120, 5);
ItemPointer item1(120, 7);
ItemPointer item2(123, 19);

IndexType index_type = INDEX_TYPE_BTREE;

enum INDEX_KEY_TYPE { UNIQUE_KEY = 0, NON_UNIQUE_KEY = 1 };

std::vector<INDEX_KEY_TYPE> index_types = {UNIQUE_KEY, NON_UNIQUE_KEY};

index::IndexMetadata *BuildIndexMetadata(INDEX_KEY_TYPE index_key_type) {
  // Build tuple and key schema
  std::vector<std::vector<std::string>> column_names;
  std::vector<catalog::Column> columns;
  std::vector<catalog::Schema *> schemas;
  index_type = INDEX_TYPE_BWTREE;

  catalog::Column column1(VALUE_TYPE_INTEGER, GetTypeSize(VALUE_TYPE_INTEGER),
                          "A", true);
  catalog::Column column2(VALUE_TYPE_VARCHAR, 1024, "B", true);
  catalog::Column column3(VALUE_TYPE_DOUBLE, GetTypeSize(VALUE_TYPE_DOUBLE),
                          "C", true);
  catalog::Column column4(VALUE_TYPE_INTEGER, GetTypeSize(VALUE_TYPE_INTEGER),
                          "D", true);

  columns.push_back(column1);
  columns.push_back(column2);

  // INDEX KEY SCHEMA -- {column1, column2}
  key_schema = new catalog::Schema(columns);
  key_schema->SetIndexedColumns({0, 1});

  columns.push_back(column3);
  columns.push_back(column4);

  // TABLE SCHEMA -- {column1, column2, column3, column4}
  tuple_schema = new catalog::Schema(columns);

  // Build index metadata
  bool unique_keys;
  if (index_key_type == NON_UNIQUE_KEY) {
    unique_keys = false;
  } else {
    unique_keys = true;
  }

  index::IndexMetadata *index_metadata = new index::IndexMetadata(
      "test_index", 125, index_type, INDEX_CONSTRAINT_TYPE_DEFAULT,
      tuple_schema, key_schema, unique_keys);

  return index_metadata;
}

index::BWTree<TestKeyType, TestValueType, TestComparatorType> *BuildBWTree(
    INDEX_KEY_TYPE index_key_type) {
  auto metadata = BuildIndexMetadata(index_key_type);
  auto *map = new index::BWTree<TestKeyType, TestValueType, TestComparatorType>(
      metadata);

  return map;
}

index::Index *BuildIndex(INDEX_KEY_TYPE index_key_type) {
  auto index_metadata = BuildIndexMetadata(index_key_type);
  // Build index
  index::Index *index = index::IndexFactory::GetInstance(index_metadata);
  EXPECT_TRUE(index != NULL);

  return index;
}

/*============================================================================
 * BLACK BOX TEST FOR BOTH BWTREE AND BTREE
 *===========================================================================*/
void BasicTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();

  std::vector<ItemPointer> locations;
  std::unique_ptr<index::Index> index(BuildIndex(index_key_type));

  // INDEX
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);

  // INSERT
  index->InsertEntry(key0.get(), item0);

  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 1);
  EXPECT_EQ(locations[0].block, item0.block);

  locations = index->ScanAllKeys();
  EXPECT_EQ(locations.size(), 1);

  // TEST SCAN
  std::vector<peloton::Value> values(2);
  std::vector<oid_t> key_column_ids(2);
  std::vector<ExpressionType> expr_types(2);
  ScanDirectionType direction = SCAN_DIRECTION_TYPE_FORWARD;

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("a");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 1);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("a");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_LESSTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 0);

  // DELETE
  index->DeleteEntry(key0.get(), item0);

  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 0);

  locations = index->ScanAllKeys();
  EXPECT_EQ(locations.size(), 0);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("a");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_LESSTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 0);

  delete tuple_schema;
}

TEST(IndexTests, BasicTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    BasicTestHelper(index_types[i]);
  }
}

// INSERT HELPER FUNCTION
// Loop based on scale factor
void InsertTest(index::Index *index, VarlenPool *pool, size_t scale_factor) {
  for (size_t scale_itr = 1; scale_itr <= scale_factor; scale_itr++) {
    // Insert a bunch of keys based on scale itr
    std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key4(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> keynonce(
        new storage::Tuple(key_schema, true));

    key0->SetValue(0, ValueFactory::GetIntegerValue(100 * scale_itr), pool);
    key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
    key1->SetValue(0, ValueFactory::GetIntegerValue(100 * scale_itr), pool);
    key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
    key2->SetValue(0, ValueFactory::GetIntegerValue(100 * scale_itr), pool);
    key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
    key3->SetValue(0, ValueFactory::GetIntegerValue(400 * scale_itr), pool);
    key3->SetValue(1, ValueFactory::GetStringValue("d"), pool);
    key4->SetValue(0, ValueFactory::GetIntegerValue(500 * scale_itr), pool);
    key4->SetValue(1, ValueFactory::GetStringValue(
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"),
                   pool);
    keynonce->SetValue(0, ValueFactory::GetIntegerValue(1000 * scale_itr),
                       pool);
    keynonce->SetValue(1, ValueFactory::GetStringValue("f"), pool);

    // INSERT
    index->InsertEntry(key0.get(), item0);
    index->InsertEntry(key1.get(), item1);
    index->InsertEntry(key1.get(), item2);
    index->InsertEntry(key1.get(), item1);
    index->InsertEntry(key1.get(), item1);
    index->InsertEntry(key1.get(), item0);

    index->InsertEntry(key2.get(), item1);
    index->InsertEntry(key3.get(), item1);
    index->InsertEntry(key4.get(), item1);
    //    index->InsertEntry(key4.get(), item0);
  }
}

// DELETE HELPER FUNCTION
void DeleteTest(index::Index *index, VarlenPool *pool, size_t scale_factor) {
  // Loop based on scale factor
  for (size_t scale_itr = 1; scale_itr <= scale_factor; scale_itr++) {
    // Delete a bunch of keys based on scale itr
    std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));
    std::unique_ptr<storage::Tuple> key4(new storage::Tuple(key_schema, true));

    key0->SetValue(0, ValueFactory::GetIntegerValue(100 * scale_itr), pool);
    key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
    key1->SetValue(0, ValueFactory::GetIntegerValue(100 * scale_itr), pool);
    key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
    key2->SetValue(0, ValueFactory::GetIntegerValue(100 * scale_itr), pool);
    key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
    key3->SetValue(0, ValueFactory::GetIntegerValue(400 * scale_itr), pool);
    key3->SetValue(1, ValueFactory::GetStringValue("d"), pool);
    key4->SetValue(0, ValueFactory::GetIntegerValue(500 * scale_itr), pool);
    key4->SetValue(1, ValueFactory::GetStringValue(
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                          "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"),
                   pool);

    // DELETE
    index->DeleteEntry(key0.get(), item0);
    index->DeleteEntry(key1.get(), item1);
    index->DeleteEntry(key2.get(), item2);
    index->DeleteEntry(key3.get(), item1);
    index->DeleteEntry(key4.get(), item1);
  }
}

void DeleteTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  std::vector<ItemPointer> locations;

  // INDEX
  std::unique_ptr<index::Index> index(BuildIndex(index_key_type));

  // Single threaded test
  size_t scale_factor = 1;
  LaunchParallelTest(1, InsertTest, index.get(), pool, scale_factor);

  // Checks
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);

  // Test insert
  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 1);

  locations = index->ScanKey(key1.get());
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 5);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  locations = index->ScanKey(key2.get());
  EXPECT_EQ(locations.size(), 1);

  locations = index->ScanAllKeys();
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 9 * scale_factor);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 5 * scale_factor);
  }

  // TEST SCAN
  std::vector<peloton::Value> values(2);
  std::vector<oid_t> key_column_ids(2);
  std::vector<ExpressionType> expr_types(2);
  ScanDirectionType direction = SCAN_DIRECTION_TYPE_FORWARD;

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("b");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 5);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  // setup values
  values[0] = ValueFactory::GetIntegerValue(99);
  values[1] = ValueFactory::GetStringValue("c");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_GREATERTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 1);

  LaunchParallelTest(1, DeleteTest, index.get(), pool, scale_factor);

  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 0);

  locations = index->ScanKey(key1.get());
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 2);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  locations = index->ScanKey(key2.get());
  EXPECT_EQ(locations.size(), 1);
  EXPECT_EQ(locations[0].block, item1.block);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("b");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 2);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  // THIS TEST IS FOR BWTREE ONLY. Not implemented for btree
  std::vector<ItemPointer> reverse_locations = index->Scan(
      values, key_column_ids, expr_types, SCAN_DIRECTION_TYPE_BACKWARD);
  if (index_key_type == NON_UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 2);
    EXPECT_EQ(locations[0].block, reverse_locations[1].block);
    EXPECT_EQ(locations[0].offset, reverse_locations[1].offset);
    EXPECT_EQ(locations[1].block, reverse_locations[0].block);
    EXPECT_EQ(locations[1].offset, reverse_locations[0].offset);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("c");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_NOTEQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 2);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("d");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_GREATERTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 0);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 0);
  }

  delete tuple_schema;
}

TEST(IndexTests, DeleteTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    DeleteTestHelper(index_types[i]);
  }
}

TEST(IndexTests, MultiThreadedTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  std::vector<ItemPointer> locations;

  // INDEX
  std::unique_ptr<index::Index> index(BuildIndex(NON_UNIQUE_KEY));

  // Parallel Test
  size_t num_threads = 24;
  size_t scale_factor = 1;
  LaunchParallelTest(num_threads, InsertTest, index.get(), pool, scale_factor);

  locations = index->ScanAllKeys();
  EXPECT_EQ(locations.size(), 9 * num_threads);

  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> keynonce(
      new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));

  keynonce->SetValue(0, ValueFactory::GetIntegerValue(1000), pool);
  keynonce->SetValue(1, ValueFactory::GetStringValue("f"), pool);

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);

  locations = index->ScanKey(keynonce.get());
  EXPECT_EQ(locations.size(), 0);

  locations = index->ScanKey(key1.get());
  EXPECT_EQ(locations.size(), 5 * num_threads);
  EXPECT_EQ(locations[0].block, item0.block);

  // TEST SCAN
  std::vector<peloton::Value> values(2);
  std::vector<oid_t> key_column_ids(2);
  std::vector<ExpressionType> expr_types(2);
  ScanDirectionType direction = SCAN_DIRECTION_TYPE_FORWARD;

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("b");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  // assume non_unique_key
  EXPECT_EQ(locations.size(), 5 * num_threads);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(99);
  values[1] = ValueFactory::GetStringValue("d");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_GREATERTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  // assume non_unique_key
  EXPECT_EQ(locations.size(), 1 * num_threads);

  // DELETE
  LaunchParallelTest(num_threads, DeleteTest, index.get(), pool, scale_factor);

  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 0);

  locations = index->ScanKey(key1.get());
  EXPECT_EQ(locations.size(), 2 * num_threads);

  locations = index->ScanKey(key2.get());
  EXPECT_EQ(locations.size(), 1 * num_threads);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("b");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 2 * num_threads);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("c");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 1 * num_threads);

  // setup values
  values[0] = ValueFactory::GetIntegerValue(99);
  values[1] = ValueFactory::GetStringValue("c");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = EXPRESSION_TYPE_COMPARE_GREATERTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_EQUAL;
  locations = index->Scan(values, key_column_ids, expr_types, direction);
  EXPECT_EQ(locations.size(), 1 * num_threads);

  delete tuple_schema;
}

/*============================================================================
 * WHITE BOX TEST SPECIFIC TO BWTREE
 *===========================================================================*/

TEST(IndexTests, BWTreeMappingTableTest) {
  // the values of the templates dont really matter;
  int size_to_test = 1025;
  auto initial_nodes =
      new index::BWTreeNode<TestKeyType, TestValueType,
                            TestComparatorType> *[size_to_test];
  for (int i = 0; i < size_to_test; i++) {
    initial_nodes[i] = reinterpret_cast<
        index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *>(
        new int(52));
  }
  auto LPIDs = new index::LPID[size_to_test];
  index::MappingTable<TestKeyType, TestValueType, TestComparatorType> map;
  for (int i = 0; i < size_to_test; i++) {
    LPIDs[i] = map.InstallPage(initial_nodes[i]);
  }
  for (int i = 0; i < size_to_test; i++) {
    EXPECT_EQ(initial_nodes[i], map.GetNode(LPIDs[i]));
    ;
  }
  auto swapped_nodes =
      new index::BWTreeNode<TestKeyType, TestValueType,
                            TestComparatorType> *[size_to_test];
  for (int i = 0; i < size_to_test; i++) {
    if (i % 2) {
      swapped_nodes[i] = reinterpret_cast<
          index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *>(
          new int(52));
      map.SwapNode(LPIDs[i], initial_nodes[i], swapped_nodes[i]);
    }
  }

  for (int i = 0; i < size_to_test; i++) {
    if (i % 2) {
      EXPECT_EQ(swapped_nodes[i], map.GetNode(LPIDs[i]));
    } else {
      EXPECT_EQ(initial_nodes[i], map.GetNode(LPIDs[i]));
    };
  }
  delete[] LPIDs;
  delete[] initial_nodes;
  delete[] swapped_nodes;
}

index::LNodeStateBuilder<TestKeyType, TestValueType, TestComparatorType> *
BuildLNodeStateBuilder(
    peloton::VarlenPool *pool, INDEX_KEY_TYPE index_key_type,
    index::BWTree<TestKeyType, TestValueType, TestComparatorType> *map) {
  std::pair<TestKeyType, TestValueType> item_locations[LPAGE_ARITY];
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
  key3->SetValue(0, ValueFactory::GetIntegerValue(400), pool);
  key3->SetValue(1, ValueFactory::GetStringValue("d"), pool);
  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  TestKeyType index_key3;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());
  index_key3.SetFromKey(key3.get());

  oid_t size;
  if (index_key_type == NON_UNIQUE_KEY) {
    item_locations[0] =
        std::pair<TestKeyType, TestValueType>(index_key0, item0);
    item_locations[1] =
        std::pair<TestKeyType, TestValueType>(index_key0, item1);
    item_locations[2] =
        std::pair<TestKeyType, TestValueType>(index_key0, item2);
    item_locations[3] =
        std::pair<TestKeyType, TestValueType>(index_key1, item0);
    item_locations[4] =
        std::pair<TestKeyType, TestValueType>(index_key1, item2);
    item_locations[5] =
        std::pair<TestKeyType, TestValueType>(index_key1, item1);
    item_locations[6] =
        std::pair<TestKeyType, TestValueType>(index_key2, item2);
    item_locations[7] =
        std::pair<TestKeyType, TestValueType>(index_key3, item2);

    size = 8;
  } else {
    item_locations[0] =
        std::pair<TestKeyType, TestValueType>(index_key0, item0);
    item_locations[1] =
        std::pair<TestKeyType, TestValueType>(index_key1, item1);
    item_locations[2] =
        std::pair<TestKeyType, TestValueType>(index_key2, item2);
    item_locations[3] =
        std::pair<TestKeyType, TestValueType>(index_key3, item2);
    size = 4;
  }

  index::LNodeStateBuilder<TestKeyType, TestValueType, TestComparatorType> *
      builder = new index::LNodeStateBuilder<TestKeyType, TestValueType,
                                             TestComparatorType>(
          INVALID_LPID, INVALID_LPID, item_locations, size, map);
  return builder;
}

void ScanKeyHelper(
    const TestKeyType &index_key0, INDEX_KEY_TYPE index_key_type,
    index::LPage<TestKeyType, TestValueType, TestComparatorType> *lpage) {
  // TEST SCAN KEY
  std::vector<TestValueType> locations;
  lpage->ScanKey(index_key0, locations);

  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 3);
  } else {
    EXPECT_EQ(locations.size(), 1);
  }
  return;
}

void ScanAllKeysHelper(
    INDEX_KEY_TYPE index_key_type,
    index::LPage<TestKeyType, TestValueType, TestComparatorType> *lpage) {
  // TEST SCAN ALL KEYS
  std::vector<TestValueType> locations;
  lpage->ScanAllKeys(locations);

  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 8);
  } else {
    EXPECT_EQ(locations.size(), 4);
  }
}
// TEST SCAN
void ScanHelper(
    INDEX_KEY_TYPE index_key_type,
    index::LPage<TestKeyType, TestValueType, TestComparatorType> *lpage,
    VarlenPool *pool) {
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());

  std::vector<TestValueType> locations;
  std::vector<peloton::Value> values(2);
  std::vector<oid_t> key_column_ids(2);
  std::vector<ExpressionType> expr_types(2);
  ScanDirectionType direction = SCAN_DIRECTION_TYPE_FORWARD;
  // TestValueType last_val_non_unique, last_val_unique;

  // setup values
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("a");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's
  expr_types[0] = (EXPRESSION_TYPE_COMPARE_EQUAL);
  expr_types[1] = (EXPRESSION_TYPE_COMPARE_EQUAL);
  locations.clear();
  lpage->Scan(values, key_column_ids, expr_types, direction, locations,
              &index_key0);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 3);
  } else {
    EXPECT_EQ(locations.size(), 1);
  }

  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("a");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  // setup expr's, col1 == 100 && col 2 != 'a'
  expr_types[0] = EXPRESSION_TYPE_COMPARE_EQUAL;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_NOTEQUAL;
  locations.clear();

  // Construct the lower bound key tuple
  std::unique_ptr<storage::Tuple> start_key;
  start_key.reset(new storage::Tuple(
      BuildIndexMetadata(index_key_type)->GetKeySchema(), true));
  start_key->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  start_key->SetValue(1, Value::GetMinValue(VALUE_TYPE_VARCHAR), pool);
  TestKeyType index_key_low;
  index_key_low.SetFromKey(start_key.get());

  lpage->Scan(values, key_column_ids, expr_types, direction, locations,
              &index_key_low);
  if (index_key_type == NON_UNIQUE_KEY) {
    // TODO fix this!!!
    EXPECT_EQ(locations.size(), 4);
    // last_val_non_unique = locations[3];
  } else {
    // EXPECT_EQ(locations.size(), 2);
    // last_val_unique =locations[1];
  }

  // col1 > 100 && col 2 != 'a'
  values[0] = ValueFactory::GetIntegerValue(100);
  values[1] = ValueFactory::GetStringValue("a");
  // setup column id's
  key_column_ids[0] = 0;
  key_column_ids[1] = 1;
  expr_types[0] = EXPRESSION_TYPE_COMPARE_GREATERTHAN;
  expr_types[1] = EXPRESSION_TYPE_COMPARE_NOTEQUAL;
  locations.clear();
  lpage->Scan(values, key_column_ids, expr_types, direction, locations,
              nullptr);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    EXPECT_EQ(locations.size(), 1);
  }
}

void LPageScanTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  index::BWTree<TestKeyType, TestValueType, TestComparatorType> *map =
      BuildBWTree(index_key_type);
  index::LNodeStateBuilder<TestKeyType, TestValueType, TestComparatorType> *
      builder = BuildLNodeStateBuilder(pool, index_key_type, map);
  index::LPage<TestKeyType, TestValueType, TestComparatorType> *lpage =
      reinterpret_cast<
          index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
          builder->GetPage());

  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  TestKeyType index_key0;
  index_key0.SetFromKey(key0.get());

  // TEST SCAN KEY
  ScanKeyHelper(index_key0, index_key_type, lpage);

  // TEST SCAN ALL KEYS
  ScanAllKeysHelper(index_key_type, lpage);

  // TEST SCAN
  ScanHelper(index_key_type, lpage, pool);

  delete builder;
  delete lpage;
}

void IPageScanTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  index::BWTree<TestKeyType, TestValueType, TestComparatorType> *map =
      BuildBWTree(index_key_type);
  index::LNodeStateBuilder<TestKeyType, TestValueType, TestComparatorType> *
      builder = BuildLNodeStateBuilder(pool, index_key_type, map);
  index::LPage<TestKeyType, TestValueType, TestComparatorType> *lpage =
      reinterpret_cast<
          index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
          builder->GetPage());

  // TODO TEST SCAN KEY

  // TODO TEST SCAN ALL KEYS

  // TODO TEST SCAN

  delete builder;
  delete lpage;
}

TEST(IndexTests, LPageScanTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    LPageScanTestHelper(index_types[i]);
  }
}

TEST(IndexTests, IPageScanTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    IPageScanTestHelper(index_types[i]);
  }
}

void BWTreeLPageDeltaConsilidationTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(index_key_type);
  auto baseNode =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(map);
  index::LPID lpid = map->GetMappingTable()->InstallPage(baseNode);
  std::vector<ItemPointer> locations;
  // Insert a bunch of keys based on scale itr

  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key4(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> keynonce(
      new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
  key3->SetValue(0, ValueFactory::GetIntegerValue(400), pool);
  key3->SetValue(1, ValueFactory::GetStringValue("d"), pool);
  key4->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key4->SetValue(1, ValueFactory::GetStringValue(
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"
                        "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"),
                 pool);
  keynonce->SetValue(0, ValueFactory::GetIntegerValue(1000), pool);
  keynonce->SetValue(1, ValueFactory::GetStringValue("f"), pool);

  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  TestKeyType index_key3;
  TestKeyType index_key4;
  TestKeyType index_nonce;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());
  index_key3.SetFromKey(key3.get());
  index_key4.SetFromKey(key4.get());
  index_nonce.SetFromKey(keynonce.get());

  locations.clear();
  baseNode->ScanKey(index_key0, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  baseNode->ScanKey(index_key1, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  baseNode->ScanKey(index_key2, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  baseNode->ScanKey(index_key3, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  baseNode->ScanKey(index_key4, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  baseNode->ScanKey(index_nonce, locations);
  EXPECT_EQ(locations.size(), 0);
  // perform many inserts
  index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *prev =
      baseNode;
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key0,
                                                         item0);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key1,
                                                         item1);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key1,
                                                         item2);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key1,
                                                         item1);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key1,
                                                         item1);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key1,
                                                         item0);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key2,
                                                         item1);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key3,
                                                         item1);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(map, prev, index_key4,
                                                         item1);
  locations.clear();
  prev->ScanKey(index_key0, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_key1, locations);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    EXPECT_EQ(locations.size(), 5);
  }

  locations.clear();
  prev->ScanKey(index_key2, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_key3, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_key4, locations);
  EXPECT_EQ(locations.size(), 1);

  EXPECT_TRUE(map->CompressDeltaChain(lpid, baseNode, prev));
  EXPECT_NE(map->GetMappingTable()->GetNode(lpid), prev);
  EXPECT_NE(map->GetMappingTable()->GetNode(lpid), baseNode);
  locations.clear();
  prev->ScanKey(index_key0, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_key1, locations);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    EXPECT_EQ(locations.size(), 5);
  }
  locations.clear();
  prev->ScanKey(index_key2, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_key3, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_key4, locations);
  EXPECT_EQ(locations.size(), 1);
  locations.clear();
  prev->ScanKey(index_nonce, locations);
  EXPECT_EQ(locations.size(), 0);

  // TODO destruct all prev
  index::LPID right_split_id = 1231921234;
  index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *
      new_base_node = map->GetMappingTable()->GetNode(lpid);

  index::LPageSplitDelta<TestKeyType, TestValueType, TestComparatorType> *
      split_delta = new index::LPageSplitDelta<TestKeyType, TestValueType,
                                               TestComparatorType>(
          map, new_base_node, index_key2, 3, right_split_id);

  EXPECT_TRUE(map->CompressDeltaChain(lpid, new_base_node, split_delta));
  auto compressed_node = map->GetMappingTable()->GetNode(lpid);

  auto compressed_lnode = reinterpret_cast<
      index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
      compressed_node);

  EXPECT_EQ(compressed_lnode->GetRightSiblingLPID(), right_split_id);

  EXPECT_NE(compressed_node, split_delta);
  locations.clear();
  prev->ScanKey(index_key0, locations);
  EXPECT_EQ(locations.size(), 1);
  // TODO test for duplicate key case (key1)
  locations.clear();
  compressed_node->ScanKey(index_key1, locations);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    // TODO this is not implemented yet, put back later
    EXPECT_EQ(locations.size(), 3);
  }
  locations.clear();
  compressed_node->ScanKey(index_key2, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  compressed_node->ScanKey(index_key3, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  compressed_node->ScanKey(index_key4, locations);
  EXPECT_EQ(locations.size(), 0);
  locations.clear();
  compressed_node->ScanKey(index_nonce, locations);
  EXPECT_EQ(locations.size(), 0);

  delete map;
}

TEST(IndexTests, BWTreeLPageDeltaConsilidationTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    BWTreeLPageDeltaConsilidationTestHelper(index_types[1]);
  }
}

TEST(IndexTests, BWTreeLPageSplitTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(UNIQUE_KEY);

  index::IPage<TestKeyType, TestValueType, TestComparatorType> *originalRoot =
      reinterpret_cast<
          index::IPage<TestKeyType, TestValueType, TestComparatorType> *>(
          map->GetMappingTable()->GetNode(0));
  index::LPage<TestKeyType, TestValueType, TestComparatorType> *baseNode =
      reinterpret_cast<
          index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
          map->GetMappingTable()->GetNode(1));

  // Insert a bunch of keys based on scale itr
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);

  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());

  for (int i = 0; i < LPAGE_ARITY; i++) {
    switch (i % 3) {
      case 0:
        baseNode->GetLocationsArray()[i].first = index_key0;
        baseNode->GetLocationsArray()[i].second = item0;
        break;
      case 1:
        baseNode->GetLocationsArray()[i].first = index_key1;
        baseNode->GetLocationsArray()[i].second = item1;
        break;
      case 2:
        baseNode->GetLocationsArray()[i].first = index_key2;
        baseNode->GetLocationsArray()[i].second = item2;
        break;
    }
  }

  baseNode->SetSize(LPAGE_ARITY);

  baseNode->SplitNodes(1, 0);

  // If this passes, it means that the IPageUpdateDelta has been successfully
  // posted
  EXPECT_NE(originalRoot, map->GetMappingTable()->GetNode(0));

  // If this passes, it means that the LPageSplitDelta has been successfully
  // posted
  EXPECT_NE(baseNode, map->GetMappingTable()->GetNode(1));

  index::IPageUpdateDelta<TestKeyType, TestValueType, TestComparatorType> *
      ipage_update_delta =
          reinterpret_cast<index::IPageUpdateDelta<TestKeyType, TestValueType,
                                                   TestComparatorType> *>(
              map->GetMappingTable()->GetNode(0));

  EXPECT_EQ(originalRoot, ipage_update_delta->GetModifiedNode());

  index::LPageSplitDelta<TestKeyType, TestValueType, TestComparatorType> *
      split_delta =
          reinterpret_cast<index::LPageSplitDelta<TestKeyType, TestValueType,
                                                  TestComparatorType> *>(
              map->GetMappingTable()->GetNode(1));

  EXPECT_EQ(baseNode, split_delta->GetModifiedNode());

  EXPECT_EQ(split_delta->GetRightSplitPageLPID(), 2);

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *
      right_split_node = reinterpret_cast<
          index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
          map->GetMappingTable()->GetNode(2));

  int expected_size_right_child;

  expected_size_right_child = LPAGE_ARITY - (LPAGE_ARITY / 2) - 1;

  EXPECT_EQ(right_split_node->GetSize(), expected_size_right_child);

  for (int i = 0; i < expected_size_right_child; i++) {
    switch ((i + (LPAGE_ARITY - expected_size_right_child)) % 3) {
      case 0:
        // EXPECT_EQ(right_split_node->GetLocationsArray()[i].first,
        // index_key0);
        EXPECT_EQ(right_split_node->GetLocationsArray()[i].second.block,
                  item0.block);
        EXPECT_EQ(right_split_node->GetLocationsArray()[i].second.offset,
                  item0.offset);
        break;
      case 1:
        // EXPECT_EQ(right_split_node->GetLocationsArray()[i].first,
        // index_key1);
        EXPECT_EQ(right_split_node->GetLocationsArray()[i].second.block,
                  item1.block);
        EXPECT_EQ(right_split_node->GetLocationsArray()[i].second.offset,
                  item1.offset);
        break;
      case 2:
        // EXPECT_EQ(right_split_node->GetLocationsArray()[i].first,
        // index_key2);
        EXPECT_EQ(right_split_node->GetLocationsArray()[i].second.block,
                  item2.block);
        EXPECT_EQ(right_split_node->GetLocationsArray()[i].second.offset,
                  item2.offset);
        break;
    }
  }

  delete map;
}

}  // End test namespace
}  // End peloton namespace
