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

#include "backend/storage/tuple.h"
#define private public
#define protected public
#include "backend/index/bwtree.h"
#include "backend/index/bwtree.cpp"

namespace peloton {
namespace test {

//===--------------------------------------------------------------------===//
// Index Tests
//===--------------------------------------------------------------------===//

#define TestKeyType index::GenericKey<12>
#define TestValueType ItemPointer
#define TestComparatorType index::GenericComparator<12>
#define TestEqualityChecker index::GenericEqualityChecker<12>

catalog::Schema *key_schema = nullptr;
catalog::Schema *tuple_schema = nullptr;
index::IndexMetadata *metadata_ptr;

ItemPointer item0(120, 5);
ItemPointer item1(120, 7);
ItemPointer item2(123, 19);
ItemPointer item3(124, 20);
ItemPointer item4(125, 25);

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

  metadata_ptr = index_metadata;
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

  size_t footprint = 0;
  footprint +=
      sizeof(index::IPage<TestKeyType, TestValueType, TestComparatorType>);
  footprint += sizeof(
      index::LPageUpdateDelta<TestKeyType, TestValueType, TestComparatorType>);
  footprint +=
      sizeof(index::LPage<TestKeyType, TestValueType, TestComparatorType>);
  footprint += sizeof(
      index::MappingTable<TestKeyType, TestValueType, TestComparatorType>);
  footprint +=
      sizeof(
          index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *) *
      MAPPING_TABLE_INITIAL_CAP;
  EXPECT_EQ(index->GetMemoryFootprint(), footprint);

  index->BWTreeCheck();

  // DELETE
  index->DeleteEntry(key0.get(), item0);

  // index->Debug();

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
  assert(locations.size() == 0);

  // Test the mem foot print
  footprint += sizeof(
      index::LPageUpdateDelta<TestKeyType, TestValueType, TestComparatorType>);
  EXPECT_EQ(index->GetMemoryFootprint(), footprint);
  index->BWTreeCheck();

  // Clean up
  index->Cleanup();
  footprint -= 2 * sizeof(index::LPageUpdateDelta<TestKeyType, TestValueType,
                                                  TestComparatorType>);
  EXPECT_EQ(index->GetMemoryFootprint(), footprint);
  index->BWTreeCheck();

  delete tuple_schema;
  //  delete values;
}

// TEST(IndexTests, BasicUniqueTest) {
//    BasicTestHelper(UNIQUE_KEY);
//}
//
// TEST(IndexTests, BasicNonUniqueTest) { BasicTestHelper(NON_UNIQUE_KEY); }

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
    key4->SetValue(1, ValueFactory::GetStringValue("ee"), pool);
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
    key4->SetValue(1, ValueFactory::GetStringValue("ee"), pool);

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
  size_t footprint_before_insert = index->GetMemoryFootprint();

  // Single threaded test
  size_t scale_factor = 200;
  LaunchParallelTest(1, InsertTest, index.get(), pool, scale_factor);

  index->BWTreeCheck();
  size_t footprint_after_insert = index->GetMemoryFootprint();
  EXPECT_EQ(footprint_after_insert > footprint_before_insert, true);

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

  int key1_count = 5;

  locations = index->ScanKey(key1.get());
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), key1_count);
  } else if (index_key_type == UNIQUE_KEY && index_type == INDEX_TYPE_BWTREE) {
    EXPECT_EQ(locations.size(), 1);
  }

  locations = index->ScanKey(key2.get());
  EXPECT_EQ(locations.size(), 1);

  locations = index->ScanAllKeys();
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), (key1_count + 4) * scale_factor);
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
    EXPECT_EQ(locations.size(), key1_count);
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
  EXPECT_EQ(locations.size(), 1 * scale_factor);

  LaunchParallelTest(1, DeleteTest, index.get(), pool, scale_factor);

  index->BWTreeCheck();
  size_t footprint_after_delete = index->GetMemoryFootprint();
  EXPECT_EQ(footprint_after_delete > footprint_after_insert, true);

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

  index->Cleanup();
  index->BWTreeCheck();
  size_t footprint_after_compress = index->GetMemoryFootprint();
  EXPECT_EQ(footprint_after_delete > footprint_after_compress, true);

  delete tuple_schema;
}

TEST(IndexTests, EpochManagerTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(NON_UNIQUE_KEY);

  index::EpochManager<TestKeyType, TestValueType,
                      TestComparatorType> *epoch_manager_ =
      new index::EpochManager<TestKeyType, TestValueType, TestComparatorType>;

  int curr_epoch = epoch_manager_->GetCurrentEpoch();
  LOG_INFO("Current epoch is %d", (int)curr_epoch);
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key4(new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
  key3->SetValue(0, ValueFactory::GetIntegerValue(400), pool);
  key3->SetValue(1, ValueFactory::GetStringValue("d"), pool);
  key4->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key4->SetValue(1, ValueFactory::GetStringValue("e"), pool);

  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  TestKeyType index_key3;
  TestKeyType index_key4;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());
  index_key3.SetFromKey(key3.get());
  index_key4.SetFromKey(key4.get());

  for (int i = 0; i < 100; i++) {
    auto baseNode =
        new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
            map, index_key4, true);

    index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *prev =
        baseNode;
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key0, item0, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item1, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item2, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item1, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item1, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item0, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key2, item1, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key3, item1, index_key4, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key4, item1, index_key4, true, INVALID_LPID);

    epoch_manager_->AddNodeToEpoch(prev);
  }

  epoch_manager_->ReleaseEpoch(curr_epoch);

  delete epoch_manager_;

  LOG_INFO("Sleeping..");
  std::this_thread::sleep_for(std::chrono::milliseconds(EPOCH_LENGTH_MILLIS));

  LOG_INFO("Woke up.. Nodes should hopefully have been cleaned up");
  //

  //  delete map->GetMappingTable()->GetNode(0);
  //  delete map->GetMappingTable()->GetNode(1);
  delete metadata_ptr;

  // delete key_schema;

  delete tuple_schema;
  delete map;
}

TEST(IndexTests, DeleteTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    DeleteTestHelper(index_types[1]);
  }
}

TEST(IndexTests, IPageMergeTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(NON_UNIQUE_KEY);

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
  key4->SetValue(1, ValueFactory::GetStringValue("e"), pool);

  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  TestKeyType index_key3;
  TestKeyType index_key4;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());
  index_key3.SetFromKey(key3.get());
  index_key4.SetFromKey(key4.get());

  index::IPage<TestKeyType, TestValueType, TestComparatorType> *base_ipage =
      new index::IPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key4, true);

  index::IPage<TestKeyType, TestValueType, TestComparatorType> *left_ipage =
      new index::IPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key1, true);

  index::IPage<TestKeyType, TestValueType, TestComparatorType> *center_ipage =
      new index::IPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key2, true);

  index::IPage<TestKeyType, TestValueType, TestComparatorType> *right_ipage =
      new index::IPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key4, true);

  map->GetMappingTable()->InstallPage(base_ipage);
  int left_lpid = map->GetMappingTable()->InstallPage(left_ipage);
  int center_lpid = map->GetMappingTable()->InstallPage(center_ipage);
  int right_lpid = map->GetMappingTable()->InstallPage(right_ipage);

  //  LOG_INFO("base: %d, left: %d, center: %d, right: %d", base_lpid,
  //  left_lpid,
  //           center_lpid, right_lpid);

  base_ipage->children_[0].first = index_key1;
  base_ipage->children_[0].second = left_lpid;
  base_ipage->children_[1].first = index_key2;
  base_ipage->children_[1].second = center_lpid;
  base_ipage->children_[2].first = index_key4;
  base_ipage->children_[2].second = right_lpid;

  left_ipage->children_[0].first = index_key0;
  left_ipage->children_[0].second = 100;
  left_ipage->children_[1].first = index_key1;
  left_ipage->children_[1].second = 200;

  center_ipage->children_[0].first = index_key2;
  center_ipage->children_[0].second = 300;

  right_ipage->children_[0].first = index_key3;
  right_ipage->children_[0].second = 400;
  right_ipage->children_[1].first = index_key4;
  right_ipage->children_[1].second = 500;

  base_ipage->size_ = 3;
  left_ipage->size_ = 2;
  center_ipage->size_ = 1;
  right_ipage->size_ = 2;

  //  base_ipage->Cleanup();

  //  std::cout << base_ipage->Debug(0, base_lpid).c_str();
  //  std::cout << map->GetMappingTable()
  //                   ->GetNode(left_lpid)
  //                   ->Debug(1, left_lpid)
  //                   .c_str();
  //  std::cout << map->GetMappingTable()
  //                   ->GetNode(center_lpid)
  //                   ->Debug(1, center_lpid)
  //                   .c_str();
  //  std::cout << map->GetMappingTable()
  //                   ->GetNode(right_lpid)
  //                   ->Debug(1, right_lpid)
  //                   .c_str();
}

TEST(IndexTests, LPageMergeTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(NON_UNIQUE_KEY);

  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key4(new storage::Tuple(key_schema, true));

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);
  key2->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key2->SetValue(1, ValueFactory::GetStringValue("c"), pool);
  key3->SetValue(0, ValueFactory::GetIntegerValue(400), pool);
  key3->SetValue(1, ValueFactory::GetStringValue("d"), pool);
  key4->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key4->SetValue(1, ValueFactory::GetStringValue("e"), pool);

  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  TestKeyType index_key3;
  TestKeyType index_key4;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());
  index_key3.SetFromKey(key3.get());
  index_key4.SetFromKey(key4.get());

  index::IPage<TestKeyType, TestValueType, TestComparatorType> *base_ipage =
      new index::IPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key4, true);

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *left_lpage =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key1, false);

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *center_lpage1 =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key2, false);

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *center_lpage2 =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key2, false);

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *center_lpage3 =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key2, false);

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *right_lpage =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_key4, true);

  map->GetMappingTable()->InstallPage(base_ipage);
  int left_lpid = map->GetMappingTable()->InstallPage(left_lpage);
  int center_lpid1 = map->GetMappingTable()->InstallPage(center_lpage1);
  int center_lpid2 = map->GetMappingTable()->InstallPage(center_lpage2);
  int center_lpid3 = map->GetMappingTable()->InstallPage(center_lpage3);
  int right_lpid = map->GetMappingTable()->InstallPage(right_lpage);

  left_lpage->right_sib_ = center_lpid1;
  center_lpage1->right_sib_ = center_lpid2;
  center_lpage2->right_sib_ = center_lpid3;
  center_lpage3->right_sib_ = right_lpid;
  right_lpage->right_sib_ = INVALID_LPID;

  //  LOG_INFO(
  //      "base: %d, left: %d, center1: %d, center2: %d, center3: %d, right:
  //      %d",
  //      base_lpid, left_lpid, center_lpid1, center_lpid2, center_lpid3,
  //      right_lpid);

  base_ipage->children_[0].first = index_key1;
  base_ipage->children_[0].second = left_lpid;
  base_ipage->children_[1].first = index_key2;
  base_ipage->children_[1].second = center_lpid1;
  base_ipage->children_[2].first = index_key4;
  base_ipage->children_[2].second = right_lpid;

  left_lpage->locations_[0].first = index_key0;
  left_lpage->locations_[0].second = item0;
  left_lpage->locations_[1].first = index_key1;
  left_lpage->locations_[1].second = item1;
  left_lpage->locations_[2].first = index_key1;
  left_lpage->locations_[2].second = item1;

  center_lpage1->locations_[0].first = index_key2;
  center_lpage1->locations_[0].second = item2;
  center_lpage2->locations_[0].first = index_key2;
  center_lpage2->locations_[0].second = item2;
  center_lpage3->locations_[0].first = index_key2;
  center_lpage3->locations_[0].second = item2;

  right_lpage->locations_[0].first = index_key3;
  right_lpage->locations_[0].second = item3;
  right_lpage->locations_[1].first = index_key4;
  right_lpage->locations_[1].second = item4;

  base_ipage->size_ = 3;

  left_lpage->size_ = 3;
  center_lpage1->size_ = 1;
  center_lpage2->size_ = 1;
  center_lpage3->size_ = 1;
  right_lpage->size_ = 2;

  //  base_ipage->Cleanup();

  //  std::cout << base_ipage->Debug(0, base_lpid).c_str();
  //
  //  std::cout << map->GetMappingTable()
  //                   ->GetNode(left_lpid)
  //                   ->Debug(1, left_lpid)
  //                   .c_str();
  //  std::cout << map->GetMappingTable()
  //                   ->GetNode(center_lpid)
  //                   ->Debug(1, center_lpid)
  //                   .c_str();
  //  std::cout << map->GetMappingTable()
  //                   ->GetNode(right_lpid)
  //                   ->Debug(1, right_lpid)
  //                   .c_str();
}

// TEST(IndexTests, DeleteTest) { DeleteTestHelper(index_types[1]); }

TEST(IndexTests, MultiThreadedTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  std::vector<ItemPointer> locations;

  // INDEX
  std::unique_ptr<index::Index> index(BuildIndex(NON_UNIQUE_KEY));

  // Parallel Test
  size_t num_threads = 1;
  size_t scale_factor = 2048;
  LaunchParallelTest(num_threads, InsertTest, index.get(), pool, scale_factor);

  locations = index->ScanAllKeys();
  EXPECT_EQ(locations.size(), 9 * num_threads * scale_factor);

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
  EXPECT_EQ(locations.size(), 1 * num_threads * scale_factor);

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
  EXPECT_EQ(locations.size(), 1 * num_threads * scale_factor);

  delete tuple_schema;
}

TEST(IndexTests, SingleThreadedSplitTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  std::vector<ItemPointer> locations;

  // INDEX
  std::unique_ptr<index::Index> index(BuildIndex(NON_UNIQUE_KEY));

  // Parallel Test
  size_t num_threads = 24;
  size_t scale_factor = 1;
  // LaunchParallelTest(num_threads, InsertTest, index.get(), pool,
  // scale_factor);

  for (int i = 0; i < (int)num_threads; i++) {
    for (size_t scale_itr = 1; scale_itr <= scale_factor; scale_itr++) {
      // Insert a bunch of keys based on scale itr
      std::unique_ptr<storage::Tuple> key0(
          new storage::Tuple(key_schema, true));
      std::unique_ptr<storage::Tuple> key1(
          new storage::Tuple(key_schema, true));
      std::unique_ptr<storage::Tuple> key2(
          new storage::Tuple(key_schema, true));
      std::unique_ptr<storage::Tuple> key3(
          new storage::Tuple(key_schema, true));
      std::unique_ptr<storage::Tuple> key4(
          new storage::Tuple(key_schema, true));
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
      key4->SetValue(1,
                     ValueFactory::GetStringValue(
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

  LOG_INFO("Before the first ScanAllKeys");
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
  // clean up test
  for (int i = 0; i < size_to_test; i++) {
    map.mapping_table_[i] = nullptr;
  }
  printf("here");
  delete[] LPIDs;
  delete[] initial_nodes;
  delete[] swapped_nodes;
}

/*============================================================================
 * WHITE BOX TEST SPECIFIC TO BWTREE
 *===========================================================================*/

// TEST(IndexTests, BWTreeMappingTableTest) {
//  // the values of the templates dont really matter;
//  int size_to_test = 1025;
//  auto initial_nodes =
//      new index::BWTreeNode<TestKeyType, TestValueType,
//                            TestComparatorType> *[size_to_test];
//  for (int i = 0; i < size_to_test; i++) {
//    initial_nodes[i] = reinterpret_cast<
//        index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *>(
//        new int(52));
//  }
//  auto LPIDs = new index::LPID[size_to_test];
//  index::MappingTable<TestKeyType, TestValueType, TestComparatorType> map;
//  for (int i = 0; i < size_to_test; i++) {
//    LPIDs[i] = map.InstallPage(initial_nodes[i]);
//  }
//  for (int i = 0; i < size_to_test; i++) {
//    EXPECT_EQ(initial_nodes[i], map.GetNode(LPIDs[i]));
//    ;
//  }
//  auto swapped_nodes =
//      new index::BWTreeNode<TestKeyType, TestValueType,
//                            TestComparatorType> *[size_to_test];
//  for (int i = 0; i < size_to_test; i++) {
//    if (i % 2) {
//      swapped_nodes[i] = reinterpret_cast<
//          index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType>
//          *>(
//          new int(52));
//      map.SwapNode(LPIDs[i], initial_nodes[i], swapped_nodes[i]);
//    }
//  }
//
//  for (int i = 0; i < size_to_test; i++) {
//    if (i % 2) {
//      EXPECT_EQ(swapped_nodes[i], map.GetNode(LPIDs[i]));
//    } else {
//      EXPECT_EQ(initial_nodes[i], map.GetNode(LPIDs[i]));
//    };
//  }
//  delete[] LPIDs;
//  delete[] initial_nodes;
//  delete[] swapped_nodes;
//}

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
          INVALID_LPID, INVALID_LPID, item_locations, size, map, index_key0,
          true);
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
  delete map;
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
  delete map;
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

void BWTreeIPageDeltaConsilidationTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(index_key_type);

  map->epoch_manager_.GetCurrentEpoch();
  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key2(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key3(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key4(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key5(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key6(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key7(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key8(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key9(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key10(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key11(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key12(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key13(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key14(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key15(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key16(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key17(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key18(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key19(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> key20(new storage::Tuple(key_schema, true));
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
  key4->SetValue(1, ValueFactory::GetStringValue("e"), pool);
  key5->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key5->SetValue(1, ValueFactory::GetStringValue("f"), pool);
  key6->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key6->SetValue(1, ValueFactory::GetStringValue("g"), pool);
  key7->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key7->SetValue(1, ValueFactory::GetStringValue("h"), pool);
  key8->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key8->SetValue(1, ValueFactory::GetStringValue("i"), pool);
  key9->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key9->SetValue(1, ValueFactory::GetStringValue("j"), pool);
  key10->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key10->SetValue(1, ValueFactory::GetStringValue("k"), pool);
  key11->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key11->SetValue(1, ValueFactory::GetStringValue("l"), pool);
  key12->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key12->SetValue(1, ValueFactory::GetStringValue("m"), pool);
  key13->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key13->SetValue(1, ValueFactory::GetStringValue("n"), pool);
  key14->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key14->SetValue(1, ValueFactory::GetStringValue("o"), pool);
  key15->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key15->SetValue(1, ValueFactory::GetStringValue("p"), pool);
  key16->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key16->SetValue(1, ValueFactory::GetStringValue("q"), pool);
  key17->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key17->SetValue(1, ValueFactory::GetStringValue("r"), pool);
  key18->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key18->SetValue(1, ValueFactory::GetStringValue("s"), pool);
  key19->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key19->SetValue(1, ValueFactory::GetStringValue("t"), pool);
  key20->SetValue(0, ValueFactory::GetIntegerValue(500), pool);
  key20->SetValue(1, ValueFactory::GetStringValue("u"), pool);
  keynonce->SetValue(0, ValueFactory::GetIntegerValue(10000000), pool);
  keynonce->SetValue(1, ValueFactory::GetStringValue("z"), pool);

  TestKeyType index_key0;
  TestKeyType index_key1;
  TestKeyType index_key2;
  TestKeyType index_key3;
  TestKeyType index_key4;
  TestKeyType index_key5;
  TestKeyType index_key6;
  TestKeyType index_key7;
  TestKeyType index_key8;
  TestKeyType index_key9;
  TestKeyType index_key10;
  TestKeyType index_key11;
  TestKeyType index_key12;
  TestKeyType index_key13;
  TestKeyType index_key14;
  TestKeyType index_key15;
  TestKeyType index_key16;
  TestKeyType index_key17;
  TestKeyType index_key18;
  TestKeyType index_key19;
  TestKeyType index_key20;
  TestKeyType index_nonce;

  index_key0.SetFromKey(key0.get());
  index_key1.SetFromKey(key1.get());
  index_key2.SetFromKey(key2.get());
  index_key3.SetFromKey(key3.get());
  index_key4.SetFromKey(key4.get());
  index_key5.SetFromKey(key5.get());
  index_key6.SetFromKey(key6.get());
  index_key7.SetFromKey(key7.get());
  index_key8.SetFromKey(key8.get());
  index_key9.SetFromKey(key9.get());
  index_key10.SetFromKey(key10.get());
  index_key11.SetFromKey(key11.get());
  index_key12.SetFromKey(key12.get());
  index_key13.SetFromKey(key13.get());
  index_key14.SetFromKey(key14.get());
  index_key15.SetFromKey(key15.get());
  index_key16.SetFromKey(key16.get());
  index_key17.SetFromKey(key17.get());
  index_key18.SetFromKey(key18.get());
  index_key19.SetFromKey(key19.get());
  index_key20.SetFromKey(key20.get());
  index_nonce.SetFromKey(keynonce.get());

  auto baseNode =
      new index::IPage<TestKeyType, TestValueType, TestComparatorType>(
          map, index_nonce, true);
  baseNode->children_[0].second = 2000;
  index::LPID lpid = map->GetMappingTable()->InstallPage(baseNode);
  std::vector<ItemPointer> locations;
  // Insert a bunch of keys based on scale itr

  EXPECT_EQ(baseNode->children_[map->GetChild(index_key0, baseNode->children_,
                                              baseNode->size_)].second,
            2000);
  EXPECT_EQ(baseNode->children_[map->GetChild(index_key1, baseNode->children_,
                                              baseNode->size_)].second,
            2000);
  EXPECT_EQ(baseNode->children_[map->GetChild(index_key2, baseNode->children_,
                                              baseNode->size_)].second,
            2000);
  EXPECT_EQ(baseNode->children_[map->GetChild(index_key3, baseNode->children_,
                                              baseNode->size_)].second,
            2000);
  EXPECT_EQ(baseNode->children_[map->GetChild(index_key4, baseNode->children_,
                                              baseNode->size_)].second,
            2000);
  EXPECT_EQ(baseNode->children_[map->GetChild(index_nonce, baseNode->children_,
                                              baseNode->size_)].second,
            2000);
  // perform many inserts
  index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *prev =
      baseNode;
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key0, index_nonce, true, 2000, 2004, false,
      baseNode->GetRightMostKey(), baseNode->IsInifinity());
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key4, index_nonce, true, 2004, 2008, false,
      baseNode->GetRightMostKey(), baseNode->IsInifinity());
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key8, index_nonce, true, 2008, 2012, false,
      baseNode->GetRightMostKey(), baseNode->IsInifinity());
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key12, index_nonce, true, 2012, 2016, false,
      baseNode->GetRightMostKey(), baseNode->IsInifinity());
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key16, index_nonce, true, 2016, 2020, false,
      baseNode->GetRightMostKey(), baseNode->IsInifinity());
  EXPECT_TRUE(map->CompressDeltaChain(lpid, baseNode, prev));
  auto newBaseNode = reinterpret_cast<
      index::IPage<TestKeyType, TestValueType, TestComparatorType> *>(
      map->GetMappingTable()->GetNode(lpid));
  EXPECT_EQ(newBaseNode->size_, 6);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key0, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key1, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2004);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key4, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2004);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key8, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2008);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key12, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2012);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key13, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2016);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key14, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2016);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key15, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2016);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key16, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2016);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_nonce, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2020);
  prev = newBaseNode;
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key18, index_nonce, true, 2020, 2021, false,
      newBaseNode->GetRightMostKey(), newBaseNode->IsInifinity());

  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key2, index_key4, false, 2004, 3000, false,
      newBaseNode->GetRightMostKey(), newBaseNode->IsInifinity());

  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key14, index_key16, false, 2016, 2014, false,
      newBaseNode->GetRightMostKey(), newBaseNode->IsInifinity());

  EXPECT_TRUE(map->CompressDeltaChain(lpid, newBaseNode, prev));
  newBaseNode = reinterpret_cast<
      index::IPage<TestKeyType, TestValueType, TestComparatorType> *>(
      map->GetMappingTable()->GetNode(lpid));
  EXPECT_EQ(newBaseNode->size_, 9);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key0, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key1, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2004);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key2, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2004);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key3, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      3000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key4, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      3000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key5, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2008);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key8, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2008);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key12, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2012);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key13, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2016);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key14, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2016);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key15, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2014);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key16, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2014);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key17, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2020);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key18, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2020);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key19, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2021);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_nonce, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2021);
  for (unsigned int i = 0; i < newBaseNode->size_; i++) {
    LOG_INFO("%s , %lu", map->ToString(newBaseNode->children_[i].first).c_str(),
             newBaseNode->children_[i].second);
  }
  prev = newBaseNode;
  // this should not happen because it should be redireceted to the split lpid
  //    prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
  //            TestComparatorType>(
  // map, prev, index_key15, index_key16, false, 2016, 4000, false,
  // newBaseNode->GetRightMostKey(), newBaseNode->IsInifinity());

  EXPECT_EQ(map->CompareKey(newBaseNode->children_[3].first, index_key8), 0);
  prev = new index::IPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key6, index_key8, false, 2008, 4006, false,
      newBaseNode->GetRightMostKey(), newBaseNode->IsInifinity());

  prev = new index::IPageSplitDelta<TestKeyType, TestValueType,
                                    TestComparatorType>(
      map, prev, index_key8, 4001, 3, newBaseNode->GetRightMostKey(),
      newBaseNode->IsInifinity());

  EXPECT_TRUE(map->CompressDeltaChain(lpid, newBaseNode, prev));

  newBaseNode = reinterpret_cast<
      index::IPage<TestKeyType, TestValueType, TestComparatorType> *>(
      map->GetMappingTable()->GetNode(lpid));
  for (unsigned int i = 0; i < newBaseNode->size_; i++) {
    LOG_INFO("%s , %lu", map->ToString(newBaseNode->children_[i].first).c_str(),
             newBaseNode->children_[i].second);
  }
  EXPECT_EQ(newBaseNode->size_, 5);
  EXPECT_EQ(map->CompareKey(newBaseNode->GetRightMostKey(), index_key8), 0);
  EXPECT_FALSE(newBaseNode->IsInifinity());

  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key0, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key1, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2004);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key2, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2004);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key3, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      3000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key4, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      3000);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key5, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2008);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key6, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      2008);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key7, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      4006);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key8, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      4006);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_key9, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      4006);
  EXPECT_EQ(
      newBaseNode->children_[map->GetChild(index_nonce, newBaseNode->children_,
                                           newBaseNode->size_)].second,
      4006);

  size_t prev_mem_footprint = map->GetMemoryFootprint();
  map->GetMappingTable()->RemovePage(lpid);
  size_t current_mem_footprint = map->GetMemoryFootprint();
  EXPECT_EQ(
      prev_mem_footprint - current_mem_footprint,
      sizeof(index::IPage<TestKeyType, TestValueType, TestComparatorType>));

  delete map;
}

TEST(IndexTests, BWTreeIPageDeltaConsilidationTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    BWTreeIPageDeltaConsilidationTestHelper(index_types[1]);
  }
}

void BWTreeLPageDeltaConsilidationTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  auto map = BuildBWTree(index_key_type);
  map->epoch_manager_.GetCurrentEpoch();
  TestKeyType dummy_key_type;
  auto baseNode =
      new index::LPage<TestKeyType, TestValueType, TestComparatorType>(
          map, dummy_key_type, true);
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

  // First test many LPageUpdateDeltas on top of each other
  index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *prev =
      baseNode;
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key0, item0, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key1, item1, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key1, item2, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key1, item1, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key1, item1, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key1, item0, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key2, item1, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key3, item1, dummy_key_type, true, INVALID_LPID);
  prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                     TestComparatorType>(
      map, prev, index_key4, item1, dummy_key_type, true, INVALID_LPID);
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

  index::LPID right_split_id = 1231921234;
  index::BWTreeNode<TestKeyType, TestValueType, TestComparatorType> *
      new_base_node = map->GetMappingTable()->GetNode(lpid);

  new_base_node->BWTreeCheck();

  index::LPage<TestKeyType, TestValueType, TestComparatorType> *new_base_lpage =
      reinterpret_cast<
          index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
          new_base_node);

  LOG_INFO("After compressing, the size of the new LPage is %d",
           (int)new_base_lpage->size_);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(new_base_lpage->size_, 5);
  } else {
    EXPECT_EQ(new_base_lpage->size_, 9);
  }

  // Now test a single SplitDelta on top of LPage
  index::LPageSplitDelta<TestKeyType, TestValueType, TestComparatorType> *
      split_delta;
  if (index_key_type == UNIQUE_KEY) {
    split_delta = new index::LPageSplitDelta<TestKeyType, TestValueType,
                                             TestComparatorType>(
        map, new_base_node, index_key1, 2, right_split_id, dummy_key_type, true,
        INVALID_LPID);
  } else {
    split_delta = new index::LPageSplitDelta<TestKeyType, TestValueType,
                                             TestComparatorType>(
        map, new_base_node, index_key1, 4, right_split_id, dummy_key_type, true,
        INVALID_LPID);
  }

  EXPECT_TRUE(map->CompressDeltaChain(lpid, new_base_node, split_delta));
  auto compressed_node = map->GetMappingTable()->GetNode(lpid);

  compressed_node->BWTreeCheck();

  auto compressed_lnode = reinterpret_cast<
      index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
      compressed_node);

  // Make sure the consolidation update the right sibling pointer
  EXPECT_EQ(compressed_lnode->GetRightSiblingLPID(), right_split_id);

  EXPECT_NE(compressed_node, split_delta);
  EXPECT_NE(compressed_node, new_base_node);

  // Set the right sibling pointer to nullptr in order to do scan
  compressed_lnode->right_sib_ = INVALID_LPID;
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
    EXPECT_EQ(locations.size(), 4);
  }
  locations.clear();
  compressed_node->ScanKey(index_key2, locations);
  if (index_key_type == UNIQUE_KEY)
    EXPECT_EQ(locations.size(), 1);
  else
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

  // Now test several update deltas below a single split delta
  if (index_key_type == NON_UNIQUE_KEY) {
    prev = new_base_lpage;

    new_base_lpage->right_sib_ = INVALID_LPID;
    new_base_lpage->right_most_key = index_key4;
    right_split_id = new_base_lpage->GetRightSiblingLPID();

    lpid = map->GetMappingTable()->InstallPage(new_base_lpage);
    LOG_INFO("Again, reiterating, size of the new base page is %d",
             (int)new_base_lpage->size_);

    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key0, item0, dummy_key_type, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item1, dummy_key_type, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item2, dummy_key_type, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key2, item1, dummy_key_type, true, INVALID_LPID);

    split_delta = new index::LPageSplitDelta<TestKeyType, TestValueType,
                                             TestComparatorType>(
        map, prev, index_key2, 6, INVALID_LPID, dummy_key_type, true,
        INVALID_LPID);

    EXPECT_TRUE(map->CompressDeltaChain(lpid, new_base_lpage, split_delta));
    auto compressed_node = map->GetMappingTable()->GetNode(lpid);
    compressed_node->BWTreeCheck();

    auto compressed_lnode = reinterpret_cast<
        index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
        compressed_node);

    // Make sure the consolidation update the right sibling pointer
    EXPECT_EQ(compressed_lnode->GetRightSiblingLPID(), right_split_id);

    EXPECT_NE(compressed_node, split_delta);
    EXPECT_NE(compressed_node, new_base_node);

    // std::cout << compressed_lnode->Debug(1, lpid);
    EXPECT_EQ(compressed_lnode->size_, 11);

    locations.clear();
    compressed_lnode->ScanKey(index_key0, locations);
    EXPECT_EQ(locations.size(), 2);
    locations.clear();
    compressed_lnode->ScanKey(index_key1, locations);
    EXPECT_EQ(locations.size(), 7);
    locations.clear();
    compressed_lnode->ScanKey(index_key2, locations);
    EXPECT_EQ(locations.size(), 2);
    locations.clear();
    compressed_lnode->ScanKey(index_key3, locations);
    EXPECT_EQ(locations.size(), 0);
    locations.clear();
    compressed_lnode->ScanKey(index_key4, locations);
    EXPECT_EQ(locations.size(), 0);

    std::cout << compressed_lnode->Debug(1, lpid);

    // Now must check stuff about responsibility
    EXPECT_EQ(map->CompareKey(compressed_lnode->right_most_key, index_key2), 0);
  } else {
    prev = new_base_lpage;

    new_base_lpage->right_sib_ = INVALID_LPID;
    new_base_lpage->right_most_key = index_key4;
    right_split_id = new_base_lpage->GetRightSiblingLPID();

    lpid = map->GetMappingTable()->InstallPage(new_base_lpage);
    LOG_INFO("Again, reiterating, size of the new base page is %d",
             (int)new_base_lpage->size_);

    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key0, item0, dummy_key_type, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item1, dummy_key_type, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key1, item2, dummy_key_type, true, INVALID_LPID);
    prev = new index::LPageUpdateDelta<TestKeyType, TestValueType,
                                       TestComparatorType>(
        map, prev, index_key2, item1, dummy_key_type, true, INVALID_LPID);

    split_delta = new index::LPageSplitDelta<TestKeyType, TestValueType,
                                             TestComparatorType>(
        map, prev, index_key2, 2, INVALID_LPID, dummy_key_type, true,
        INVALID_LPID);

    EXPECT_TRUE(map->CompressDeltaChain(lpid, new_base_lpage, split_delta));
    auto compressed_node = map->GetMappingTable()->GetNode(lpid);
    compressed_node->BWTreeCheck();

    auto compressed_lnode = reinterpret_cast<
        index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
        compressed_node);

    // Make sure the consolidation update the right sibling pointer
    EXPECT_EQ(compressed_lnode->GetRightSiblingLPID(), right_split_id);

    EXPECT_NE(compressed_node, split_delta);
    EXPECT_NE(compressed_node, new_base_node);

    // std::cout << compressed_lnode->Debug(1, lpid);
    EXPECT_EQ(compressed_lnode->size_, 3);

    locations.clear();
    compressed_lnode->ScanKey(index_key0, locations);
    EXPECT_EQ(locations.size(), 1);
    locations.clear();
    compressed_lnode->ScanKey(index_key1, locations);
    EXPECT_EQ(locations.size(), 1);
    locations.clear();
    compressed_lnode->ScanKey(index_key2, locations);
    EXPECT_EQ(locations.size(), 1);
    locations.clear();
    compressed_lnode->ScanKey(index_key3, locations);
    EXPECT_EQ(locations.size(), 0);
    locations.clear();
    compressed_lnode->ScanKey(index_key4, locations);
    EXPECT_EQ(locations.size(), 0);

    std::cout << compressed_lnode->Debug(1, lpid);

    // Now must check stuff about responsibility
    EXPECT_EQ(map->CompareKey(compressed_lnode->right_most_key, index_key2), 0);
  }

  //  delete map;
}

TEST(IndexTests, BWTreeLPageDeltaConsilidationTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    BWTreeLPageDeltaConsilidationTestHelper(index_types[i]);
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

  baseNode->SplitNodes(1);

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
