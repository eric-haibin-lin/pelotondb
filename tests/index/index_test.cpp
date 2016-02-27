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

// namespace peloton{
// namespace index{
//
////===--------------------------------------------------------------------===//
//// LPageSplitDelta Methods Begin
////===--------------------------------------------------------------------===//
// template <typename KeyType, typename ValueType, class KeyComparator>
// std::vector<ValueType>
// LPageSplitDelta < KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
//  std::vector<ValueType> result;
//  assert(this->modified_node != nullptr);
//  LOG_INFO("LPageSplitDelta::ScanKey");
//
//  bool greater_than_left_key = this->map->CompareKey(key, modified_key_) > 0;
//
//  if (greater_than_left_key) {
//    LOG_INFO(
//        "LPageSplitDelta::ScanKey Found a matching key for right split page");
//    result = this->map->GetMappingTable()
//                 ->GetNode(right_split_page_lpid_)->ScanKey(key);
//  } else {
//    // Scan the modified node
//    result = this->modified_node->ScanKey(key);
//  }
//  return result;
//};
//
// template <typename KeyType, typename ValueType, class KeyComparator>
// NodeStateBuilder < KeyType, ValueType, KeyComparator> *
// LPageSplitDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
//  // Children of IPageDelta always return a INodeStateBuilder
//  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
//      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator>
//      *>(
//          this->modified_node->BuildNodeState());
//  assert(builder != nullptr);
//  builder->SeparateFromKey(modified_key_, modified_key_location_,
//                           right_split_page_lpid_);
//
//  return builder;
//}
//
// template <typename KeyType, typename ValueType, class KeyComparator>
// bool LPageSplitDelta<KeyType, ValueType, KeyComparator>::InsertEntry(
//    KeyType key, ValueType location, LPID self) {
//  if (this->map->CompareKey(key, modified_key_) ==
//      1)  // this key is greater than modified_key_
//  {
//    return this->map->GetMappingTable()->GetNode(right_split_page_lpid_)
//        ->InsertEntry(key, location, right_split_page_lpid_);
//  }
//
//  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
//      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
//                                                              key, location);
//  bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
//  if (!status) {
//    delete new_delta;
//  }
//  return status;
//};
//
// template <typename KeyType, typename ValueType, class KeyComparator>
// bool LPageSplitDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(
//    KeyType key, ValueType location, LPID self) {
//  if (this->map->CompareKey(key, modified_key_) ==
//      1)  // this key is greater than modified_key_
//  {
//    return this->map->GetMappingTable()->GetNode(right_split_page_lpid_)
//        ->InsertEntry(key, location, right_split_page_lpid_);
//  }
//
//  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
//      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
//                                                              key, location);
//  new_delta->SetDeleteFlag();
//  bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
//  if (!status) {
//    delete new_delta;
//  }
//  return status;
//};
////===--------------------------------------------------------------------===//
//// LPageSplitDelta Methods End
////===--------------------------------------------------------------------===//
//
// template <typename KeyType, typename ValueType, class KeyComparator>
// void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::SeparateFromKey(
//    KeyType separator_key, ValueType location, LPID split_new_page_id) {
//  assert(locations_ != nullptr);
//
//  int index = this->map->BinarySearch(separator_key, locations_,
//      this->size);
//  assert(index < this->size && index >= 0);
//  // assume we include the key at the split page
//  // decrement size
//  this->size = index + 1;
//  // update separator info
//  this->is_separated = true;
//  this->separator_key = separator_key;
//  separator_location_ = location;
//  this->split_new_page_id = split_new_page_id;
//}

//}}
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
  size_t num_threads = 24;
  size_t scale_factor = 1;
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

  locations = baseNode->ScanKey(index_key0);
  EXPECT_EQ(locations.size(), 0);
  locations = baseNode->ScanKey(index_key1);
  EXPECT_EQ(locations.size(), 0);
  locations = baseNode->ScanKey(index_key2);
  EXPECT_EQ(locations.size(), 0);
  locations = baseNode->ScanKey(index_key3);
  EXPECT_EQ(locations.size(), 0);
  locations = baseNode->ScanKey(index_key4);
  EXPECT_EQ(locations.size(), 0);
  locations = baseNode->ScanKey(index_nonce);
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

  locations = prev->ScanKey(index_key0);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_key1);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    EXPECT_EQ(locations.size(), 5);
  }

  locations = prev->ScanKey(index_key2);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_key3);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_key4);
  EXPECT_EQ(locations.size(), 1);

  EXPECT_TRUE(map->CompressDeltaChain(lpid, baseNode, prev));
  EXPECT_NE(map->GetMappingTable()->GetNode(lpid), prev);
  EXPECT_NE(map->GetMappingTable()->GetNode(lpid), baseNode);
  locations = prev->ScanKey(index_key0);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_key1);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    EXPECT_EQ(locations.size(), 5);
  }

  locations = prev->ScanKey(index_key2);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_key3);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_key4);
  EXPECT_EQ(locations.size(), 1);
  locations = prev->ScanKey(index_nonce);
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
  EXPECT_NE(compressed_node, split_delta);
  locations = prev->ScanKey(index_key0);
  EXPECT_EQ(locations.size(), 1);
  // TODO test for duplicate key case (key1)
  locations = compressed_node->ScanKey(index_key1);
  if (index_key_type == UNIQUE_KEY) {
    EXPECT_EQ(locations.size(), 1);
  } else {
    // TODO this is not implemented yet, put back later
    EXPECT_EQ(locations.size(), 3);
  }

  locations = compressed_node->ScanKey(index_key2);
  EXPECT_EQ(locations.size(), 0);
  locations = compressed_node->ScanKey(index_key3);
  EXPECT_EQ(locations.size(), 0);
  locations = compressed_node->ScanKey(index_key4);
  EXPECT_EQ(locations.size(), 0);
  locations = compressed_node->ScanKey(index_nonce);
  EXPECT_EQ(locations.size(), 0);

  delete map;
}

TEST(IndexTests, BWTreeLPageDeltaConsilidationTest) {
  for (unsigned int i = 0; i < index_types.size(); i++) {
    BWTreeLPageDeltaConsilidationTestHelper(index_types[1]);
  }
}

}  // End test namespace
}  // End peloton namespace
