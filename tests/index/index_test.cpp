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
#include "backend/index/bwtree.h"
#include "backend/index/index_key.h"
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

enum INDEX_KEY_TYPE { UNIQUE_KEY = 0, NON_UNIQUE_KEY = 1 };

std::vector<INDEX_KEY_TYPE> index_types = {UNIQUE_KEY, NON_UNIQUE_KEY};

index::IndexMetadata *BuildIndexMetadata(INDEX_KEY_TYPE index_key_type) {
  // Build tuple and key schema
  std::vector<std::vector<std::string>> column_names;
  std::vector<catalog::Column> columns;
  std::vector<catalog::Schema *> schemas;
  IndexType index_type = INDEX_TYPE_BTREE;
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
  bool unique_keys;
  if (index_key_type == NON_UNIQUE_KEY) {
    unique_keys = false;
  } else {
    unique_keys = true;
  }

  auto metadata = BuildIndexMetadata(index_key_type);
  TestComparatorType comparator(metadata);
  auto *map = new index::BWTree<TestKeyType, TestValueType, TestComparatorType>(
      unique_keys, comparator);

  return map;
}

index::Index *BuildIndex(INDEX_KEY_TYPE index_key_type) {
  auto index_metadata = BuildIndexMetadata(index_key_type);
  // Build index
  index::Index *index = index::IndexFactory::GetInstance(index_metadata);
  EXPECT_TRUE(index != NULL);

  return index;
}

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

  // DELETE
  index->DeleteEntry(key0.get(), item0);

  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 0);

  delete tuple_schema;
}

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
  } else if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  }

  locations = index->ScanKey(key2.get());
  EXPECT_EQ(locations.size(), 1);

  LaunchParallelTest(1, DeleteTest, index.get(), pool, scale_factor);

  locations = index->ScanKey(key0.get());
  EXPECT_EQ(locations.size(), 0);

  locations = index->ScanKey(key1.get());
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 2);
  } else if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  }

  locations = index->ScanKey(key2.get());
  EXPECT_EQ(locations.size(), 1);
  EXPECT_EQ(locations[0].block, item1.block);

  delete tuple_schema;
}

TEST(IndexTests, DeleteTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    DeleteTestHelper(index_types[i]);
  }
}

TEST(IndexTests, MultiThreadedInsertTest) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  std::vector<ItemPointer> locations;

  // INDEX
  std::unique_ptr<index::Index> index(BuildIndex(NON_UNIQUE_KEY));

  // Parallel Test
  size_t num_threads = 4;
  size_t scale_factor = 6;
  LaunchParallelTest(num_threads, InsertTest, index.get(), pool, scale_factor);

  //  locations = index->ScanAllKeys();
  //  EXPECT_EQ(locations.size(), 9 * num_threads);

  std::unique_ptr<storage::Tuple> key0(new storage::Tuple(key_schema, true));
  std::unique_ptr<storage::Tuple> keynonce(
      new storage::Tuple(key_schema, true));

  std::unique_ptr<storage::Tuple> key1(new storage::Tuple(key_schema, true));

  keynonce->SetValue(0, ValueFactory::GetIntegerValue(1000), pool);
  keynonce->SetValue(1, ValueFactory::GetStringValue("f"), pool);

  key0->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key0->SetValue(1, ValueFactory::GetStringValue("a"), pool);
  key1->SetValue(0, ValueFactory::GetIntegerValue(100), pool);
  key1->SetValue(1, ValueFactory::GetStringValue("b"), pool);

  locations = index->ScanKey(keynonce.get());
  EXPECT_EQ(locations.size(), 0);

  locations = index->ScanKey(key1.get());
  EXPECT_EQ(locations.size(), 5 * num_threads);
  EXPECT_EQ(locations[0].block, item0.block);

  delete tuple_schema;
}

void LPageScanTestHelper(INDEX_KEY_TYPE index_key_type) {
  auto pool = TestingHarness::GetInstance().GetTestingPool();
  index::BWTree<TestKeyType, TestValueType, TestComparatorType> *map =
      BuildBWTree(index_key_type);
  index::LPage<TestKeyType, TestValueType, TestComparatorType> *lpage = nullptr;

  std::pair<TestKeyType, TestValueType> item_locations[LPAGE_ARITY];
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

    size = 6;
  } else {
    item_locations[0] =
        std::pair<TestKeyType, TestValueType>(index_key0, item0);
    item_locations[1] =
        std::pair<TestKeyType, TestValueType>(index_key1, item1);
    item_locations[2] =
        std::pair<TestKeyType, TestValueType>(index_key2, item2);
    size = 3;
  }

  index::LNodeStateBuilder<TestKeyType, TestValueType, TestComparatorType>
      builder(INVALID_LPID, INVALID_LPID, item_locations, size, map);
  lpage = reinterpret_cast<
      index::LPage<TestKeyType, TestValueType, TestComparatorType> *>(
      builder.GetPage());

  std::vector<TestValueType> locations;
  locations = lpage->ScanKey(index_key0);
  if (index_key_type == NON_UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 3);
  } else {
    EXPECT_EQ(locations.size(), 1);
  }
}

TEST(IndexTests, LPageScanTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    LPageScanTestHelper(index_types[1]);
  }
}

}  // End test namespace
}  // End peloton namespace
