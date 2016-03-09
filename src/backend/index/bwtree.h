//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// BWTree.h
//
// Identification: src/backend/index/BWTree.h
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once
#include "backend/common/types.h"
#include "backend/common/logger.h"
#include "backend/index/index.h"
#include "backend/storage/tuple.h"
#include <map>
#include <vector>
#include <climits>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <thread>
#include <iostream>

namespace peloton {
namespace index {

//===--------------------------------------------------------------------===//
// Types and Enums
//===--------------------------------------------------------------------===//
typedef uint64_t LPID;

enum BWTreeNodeType {
  TYPE_LPAGE = 0,
  TYPE_IPAGE = 1,
  TYPE_OTHER = 2,
};

// The default capacity of each page
#define IPAGE_ARITY 1024
#define LPAGE_ARITY 1024

// The threshold above which pages get split
#define LPAGE_SPLIT_THRESHOLD 1000
#define IPAGE_SPLIT_THRESHOLD 1000

// The threshold under which pages get merged
#define IPAGE_MERGE_THRESHOLD 800
#define LPAGE_MERGE_THRESHOLD 800

// The initial capacity of the mapping table
#define MAPPING_TABLE_INITIAL_CAP 128
#define INVALID_LPID ULLONG_MAX

// The length limit of the delta chain to be compressed
#define LPAGE_DELTA_CHAIN_LIMIT 6
#define IPAGE_DELTA_CHAIN_LIMIT 4

// The size of each page which epoch manager manages
#define EPOCH_PAGE_SIZE 1024
#define MAX_ACTIVE_EPOCHS 2
#define EPOCH_LENGTH_MILLIS 40
#define TEMPLATE_TYPE KeyType, ValueType, KeyComparator

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode;

template <typename KeyType, typename ValueType, class KeyComparator>
class Delta;

template <typename KeyType, typename ValueType, class KeyComparator>
class IPage;

template <typename KeyType, typename ValueType, class KeyComparator>
class IPageUpdateDelta;

template <typename KeyType, typename ValueType, class KeyComparator>
class IPageSplitDelta;

template <typename KeyType, typename ValueType, class KeyComparator>
class LPage;

template <typename KeyType, typename ValueType, class KeyComparator>
class EpochManager;

//===--------------------------------------------------------------------===//
// NodeStateBuilder
//===--------------------------------------------------------------------===//
/**
 * Builder for a node state, to be used with Delta Compression and
 * scans on indexes
 */
template <typename KeyType, typename ValueType, class KeyComparator>
class NodeStateBuilder {
 protected:
  // number of k-v pairs. everything after index (size - 1) is ignored
  oid_t size;

  // pointer to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  // whether this page is separated from a certain key
  bool is_separated = false;

  // the key which splits the original page
  KeyType separator_key;

  // the LPID of the new split page
  LPID split_new_page_id = 0;

  // the right most key of the page to build
  KeyType right_most_key;

  // whether the page to build has infinity as right most key
  bool infinity = false;

 public:
  // construct an empty node state
  NodeStateBuilder(oid_t size, BWTree<KeyType, ValueType, KeyComparator> *map,
                   KeyType right_most_key, bool infinity)
      : size(size),
        map(map),
        right_most_key(right_most_key),
        infinity(infinity){};

  // materialize the page to build
  virtual BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage() = 0;

  // perform scan operation on a page builder
  virtual LPID Scan(const std::vector<Value> &values,
                    const std::vector<oid_t> &key_column_ids,
                    const std::vector<ExpressionType> &expr_types,
                    const ScanDirectionType &scan_direction,
                    std::vector<ValueType> &result,
                    const KeyType *index_key) = 0;

  // perform scan all keys operation on a page builder
  virtual LPID ScanAllKeys(std::vector<ValueType> &result) = 0;

  // set whether this page to build has right most key as infinity
  inline void SetInfinity(bool infinity) { this->infinity = infinity; }

  // set the right most key
  inline void SetRightMostKey(KeyType right_most_key) {
    this->right_most_key = right_most_key;
  }

  // whether this page was split
  inline bool IsSeparated() { return is_separated; }

  // get the LPID of the new split right page
  inline LPID GetSplitNewPageId() { return split_new_page_id; }

  virtual ~NodeStateBuilder() {}

  friend class LPage<KeyType, ValueType, KeyComparator>;
  friend class IPage<KeyType, ValueType, KeyComparator>;
};

/**
 * Builder for the node state of an IPage
 */
template <typename KeyType, typename ValueType, class KeyComparator>
class INodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // IPage children nodes
  std::pair<KeyType, LPID> children_[IPAGE_ARITY + IPAGE_DELTA_CHAIN_LIMIT];

 public:
  // Construct an IPage builder from an IPage
  INodeStateBuilder(std::pair<KeyType, LPID> *children, int children_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map,
                    KeyType right_most_key, bool infinity)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(
            children_len, map, right_most_key, infinity) {
    memcpy(children_, children,
           sizeof(std::pair<KeyType, LPID>) * children_len);
  };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage();

  ~INodeStateBuilder() {}

  //***************************************************
  // IPage Methods
  //***************************************************

  // add a new <key, LPID> pair to its children if key is not
  // found. Otherwise overwrite the pair with the same key
  void AddChild(std::pair<KeyType, LPID> &new_pair);

  // replace the last <key, LPID> pair from the IPage
  void ReplaceLastChild(LPID &inf_LPID);

  // remove an <key, LPID> pair from the IPage
  void RemoveChild(KeyType &key_to_remove);

  // remove the last <key, LPID> pair from the IPage
  void RemoveLastChild();

  LPID Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  LPID ScanAllKeys(std::vector<ValueType> &result);

  friend class LPage<KeyType, ValueType, KeyComparator>;
  friend class IPage<KeyType, ValueType, KeyComparator>;
  friend class IPageUpdateDelta<KeyType, ValueType, KeyComparator>;
  friend class IPageSplitDelta<KeyType, ValueType, KeyComparator>;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class LNodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // the LPID of the left sibling
  LPID left_sibling_ = INVALID_LPID;

  // the LPID of the right sibling
  LPID right_sibling_ = INVALID_LPID;

  // the index of the separator key
  int separator_index_ = -1;

  // <key, value> pairs
  std::pair<KeyType, ValueType>
      locations_[IPAGE_ARITY + LPAGE_DELTA_CHAIN_LIMIT];

 public:
  // construct a builder from an LPage
  LNodeStateBuilder(LPID left_sibling, LPID right_sibling,
                    std::pair<KeyType, ValueType> *locations, int location_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map,
                    KeyType right_most_key, bool infinity)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(
            location_len, map, right_most_key, infinity),
        left_sibling_(left_sibling),
        right_sibling_(right_sibling) {
    memcpy(locations_, locations,
           sizeof(std::pair<KeyType, ValueType>) * location_len);
  };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage();

  ~LNodeStateBuilder(){};

  //***************************************************
  // LPage Methods
  //***************************************************

  // set the value of left sibling LPID
  void UpdateLeftSib(LPID new_left_sib) { left_sibling_ = new_left_sib; }

  // set the value of right sibling LPID
  void UpdateRightSib(LPID new_right_sib) { right_sibling_ = new_right_sib; }

  // add a <key, value> pair
  void AddLeafData(std::pair<KeyType, ValueType> &new_pair);

  // remove a <key, value> pair
  void RemoveLeafData(std::pair<KeyType, ValueType> &entry_to_remove);

  LPID Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,

            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  LPID ScanAllKeys(std::vector<ValueType> &result);

  // perform scan key on the page to build
  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

 private:
  // returns true if the two item pointers are equal
  bool ItemPointerEquals(ValueType v1, ValueType v2);

  friend class LPage<KeyType, ValueType, KeyComparator>;
};

//===--------------------------------------------------------------------===//
// ReadWriteLatch is an implementation of read write lock using compare and
// swap intrinsics
//===--------------------------------------------------------------------===//
class ReadWriteLatch {
 public:
  // acquire the read lock
  inline void AquireRead() {
    LOG_INFO("Aquiring read Lock");
    while (true) {
      while (current_writers_ == 1)
        ;
      __sync_add_and_fetch(&current_readers_, 1);
      if (current_writers_ == 0)
        break;
      else
        __sync_add_and_fetch(&current_readers_, -1);
    }
  }

  // release the read lock
  inline void ReleaseRead() {
    LOG_INFO("Releasing read Lock");
    __sync_add_and_fetch(&current_readers_, -1);
  }

  // acquire the write lock
  inline void AquireWrite() {
    LOG_INFO("Aquiring write Lock");
    while (!__sync_bool_compare_and_swap(&current_writers_, 0, 1))
      ;
    while (current_readers_ > 0)
      ;
  }

  // release the write lock
  inline void ReleaseWrite() {
    LOG_INFO("Releasing write Lock");
    __sync_bool_compare_and_swap(&current_writers_, 1, 0);
  }

 private:
  // number of readers
  int current_readers_ = 0;

  // number of writers
  int current_writers_ = 0;
};

/*
 * The mapping table from LPID to BWTreeNode pointers. Updating the table
 * requires the acquisition of the read write latch. The table is initialized
 * with MAPPING_TABLE_INITIAL_CAP and double its capacity when its full
 */
template <typename KeyType, typename ValueType, class KeyComparator>
class MappingTable {
 private:
  // the capacity of mapping table
  unsigned int mapping_table_cap_ = MAPPING_TABLE_INITIAL_CAP;

  // the size of the mapping table
  unsigned int mapping_table_size_ = 0;

  // the content of the table
  BWTreeNode<KeyType, ValueType, KeyComparator> **mapping_table_;

  // next available LPID
  LPID next_LPID_ = 0;

  // the read write latch
  ReadWriteLatch latch_;

 public:
  // initialize an empty mapping table
  MappingTable() {
    LOG_INFO("Constructing Mapping Table with initial capacity %d",
             mapping_table_cap_);
    mapping_table_ = new BWTreeNode<KeyType, ValueType,
                                    KeyComparator> *[mapping_table_cap_]();
  };

  // clean up all the bwtree nodes in the table
  ~MappingTable() {
    for (int i = 0; i < (long)this->mapping_table_cap_; i++) {
      if (mapping_table_[i] != nullptr) {
        mapping_table_[i]->SetCleanUpChildren();
        delete mapping_table_[i];
      }
    }
    delete[] mapping_table_;
  }

  // install a new page with the next available LPID to the mapping table
  LPID InstallPage(BWTreeNode<KeyType, ValueType, KeyComparator> *node) {
    LOG_INFO("Installing page in mapping table");
    LPID newLPID = __sync_fetch_and_add(&next_LPID_, 1);

    // table grew too large, expand it
    while (newLPID >= mapping_table_cap_) {
      LOG_INFO("mapping table has grown too large");

      // only one thread should expand the table
      latch_.AquireWrite();
      if (newLPID < mapping_table_cap_) {
        latch_.ReleaseWrite();
        break;
      }
      // doubling the capacity
      int new_mapping_table_cap = mapping_table_cap_ * 2;
      LOG_INFO("doubling size of mapping table capacity from %d to %d",
               mapping_table_cap_, new_mapping_table_cap);
      auto new_mapping_table =
          new BWTreeNode<KeyType, ValueType,
                         KeyComparator> *[new_mapping_table_cap]();
      memcpy(new_mapping_table, mapping_table_,
             mapping_table_cap_ *
                 sizeof(new BWTreeNode<KeyType, ValueType, KeyComparator> *));
      delete[] mapping_table_;
      mapping_table_ = new_mapping_table;
      mapping_table_cap_ = new_mapping_table_cap;
      latch_.ReleaseWrite();
    }
    latch_.AquireRead();
    LOG_INFO("adding LPID: %lu to mapping table", newLPID);
    mapping_table_[newLPID] = node;
    latch_.ReleaseRead();
    return newLPID;
  }

  // remove a bwtree node associated with the given LPID
  void RemovePage(LPID id) {
    LOG_INFO("Inside RemovePage for LPID %d", (int)id);
    delete mapping_table_[id];
    mapping_table_[id] = nullptr;
  }

  // swap an entry in the mapping table using CAS
  bool SwapNode(LPID id, BWTreeNode<KeyType, ValueType, KeyComparator> *oldNode,
                BWTreeNode<KeyType, ValueType, KeyComparator> *newNode) {
    LOG_INFO("swapping node for LPID: %lu into mapping table", id);
    latch_.AquireRead();
    bool ret =
        __sync_bool_compare_and_swap(mapping_table_ + id, oldNode, newNode);
    latch_.ReleaseRead();
    return ret;
  }

  // return the bwtree node pointer of the given LPID
  BWTreeNode<KeyType, ValueType, KeyComparator> *GetNode(LPID id) {
    LOG_INFO("getting node for LPID: %lu from mapping table", id);
    latch_.AquireRead();
    auto ret = mapping_table_[id];
    latch_.ReleaseRead();
    return ret;
  }

  // calculating the memory footprint by aggregating all bwtree nodes in
  // the mapping table
  size_t GetMemoryFootprint() {
    LOG_INFO("MappingTable::GetMemoryFootprint");
    size_t total = 0;
    for (int i = 0; i < (long)this->mapping_table_cap_; i++) {
      if (mapping_table_[i] != nullptr) {
        total += mapping_table_[i]->GetMemoryFootprint();
      }
    }
    // also take into account of the size of the mapping table
    total += sizeof(MappingTable<KeyType, ValueType, KeyComparator>);
    total += mapping_table_cap_ *
             sizeof(BWTreeNode<KeyType, ValueType, KeyComparator> *);
    return total;
  }

  friend class BWTree<KeyType, ValueType, KeyComparator>;
};

// an internal structure of the epoch linked list, where each list node
// is an array of the bwtree node pointer to be free'd
template <typename KeyType, typename ValueType, class KeyComparator>
struct EpochPage {
  BWTreeNode<KeyType, ValueType, KeyComparator> *pointers[EPOCH_PAGE_SIZE];
  EpochPage *next_page = nullptr;
};

//===--------------------------------------------------------------------===//
// EpochManager keeps track of all the pages to be collected and the number
// of active thread per epoch. The pages are garbage collected once there're
// no active thread in the current epoch. Epoch page is the internal structure
// to keep track of all the pages.
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class EpochManager {
 public:
  // initialize resources for epoch manager
  EpochManager() {
    // should start GC pthread
    current_epoch_ = 0;
    destructor_called_ = false;
    memset(&epoch_size_, 0, MAX_ACTIVE_EPOCHS * sizeof(int));
    memset(&epoch_users_, 0, MAX_ACTIVE_EPOCHS * sizeof(int));
    for (int i = 0; i < MAX_ACTIVE_EPOCHS; i++) {
      first_pages_[i] = new EpochPage<KeyType, ValueType, KeyComparator>();
    }
    pthread_create(&management_pthread_, nullptr, &epoch_management, this);
  }

  // join the garbage collection thread
  ~EpochManager() {
    destructor_called_ = true;
    LOG_INFO("set destructor called flag");
    pthread_join(management_pthread_, nullptr);
  }

  // get the current epoch number
  unsigned long GetCurrentEpoch() {
    auto current_epoch = current_epoch_;
    __sync_add_and_fetch(&epoch_users_[current_epoch], 1);
    LOG_INFO("Got current epoch %lu", current_epoch);
    return current_epoch;
  }

  // release the current epoch
  void ReleaseEpoch(unsigned long i) {
    __sync_add_and_fetch(&epoch_users_[i], -1);
  }

  // add a page to current epoch table
  void AddNodeToEpoch(
      BWTreeNode<KeyType, ValueType, KeyComparator> *node_to_destroy) {
    LOG_INFO("Adding node to epoch");
    auto curr_epoch = current_epoch_;
    auto write_pos = __sync_fetch_and_add(&epoch_size_[curr_epoch], 1);
    // find correct page
    LOG_INFO("write_pos is %d", (int)write_pos);

    auto curr_page = first_pages_[curr_epoch];
    for (int i = 0; i < (write_pos / EPOCH_PAGE_SIZE); i++) {
      // LOG_INFO("Epoch Page Full, getting child page");
      LOG_INFO("At page number: %d", (int)i);
      if (curr_page->next_page == nullptr) {
        // allocate a new epoch page
        EpochPage<KeyType, ValueType, KeyComparator> *new_page =
            new EpochPage<KeyType, ValueType, KeyComparator>();
        if (!__sync_bool_compare_and_swap(&(curr_page->next_page), nullptr,
                                          new_page)) {
          delete new_page;
        }
      }
      curr_page = curr_page->next_page;
    }
    // add the node to the epoch page
    curr_page->pointers[write_pos % EPOCH_PAGE_SIZE] = node_to_destroy;
  }

 private:
  EpochPage<KeyType, ValueType, KeyComparator> *first_pages_[MAX_ACTIVE_EPOCHS];
  int epoch_size_[MAX_ACTIVE_EPOCHS];

  // the number of active users in each epoch
  int epoch_users_[MAX_ACTIVE_EPOCHS];

  // the current epoch number
  unsigned long current_epoch_;

  bool destructor_called_;

  // the thread responsible for garbage collection
  pthread_t management_pthread_;

  // private internal methods for GC pthread
  static inline void CollectGarbageForPage(
      EpochPage<KeyType, ValueType, KeyComparator> *top_page) {
    auto curr_epoch_page = top_page;
    while (curr_epoch_page != nullptr) {
      LOG_INFO("Garbage Collector cleaning page");
      // free all pages
      for (int i = 0; i < EPOCH_PAGE_SIZE; i++) {
        auto currptr = curr_epoch_page->pointers[i];
        if (currptr != nullptr) {
          currptr->SetCleanUpChildren();
          delete currptr;
        }
      }
      auto next_epoch_page = curr_epoch_page->next_page;
      delete curr_epoch_page;
      curr_epoch_page = next_epoch_page;
    }
  }

  static void *epoch_management(__attribute__((unused)) void *args) {
    auto manager =
        reinterpret_cast<EpochManager<KeyType, ValueType, KeyComparator> *>(
            args);
    LOG_INFO("Epoch Management Thread initialized");
    while (!manager->destructor_called_) {
      // give it time to start up
      std::this_thread::sleep_for(
          std::chrono::milliseconds(EPOCH_LENGTH_MILLIS));
      LOG_INFO("Collecting Garbage start");
      auto old_epoch = manager->current_epoch_;
      auto new_epoch = (old_epoch + 1) % MAX_ACTIVE_EPOCHS;
      manager->current_epoch_ = new_epoch;

      // wait for other threads to see the new epoch
      while (manager->epoch_size_[new_epoch] == 0 &&
             !manager->destructor_called_)
        ;

      // wait for users of the old epoch to complete
      while (manager->epoch_users_[old_epoch] > 0 &&
             !manager->destructor_called_)
        ;
      LOG_INFO("Collecting Garbage for page %lu", old_epoch);
      CollectGarbageForPage(manager->first_pages_[old_epoch]);
      manager->first_pages_[old_epoch] =
          new EpochPage<KeyType, ValueType, KeyComparator>();
      manager->epoch_size_[old_epoch] = 0;
    }
    // we are exiting collect all garbage
    // we assume no threads are currently adding_stuff
    LOG_INFO("Epoch Management Thread cleaning up remaining garbage");
    for (int i = 0; i < MAX_ACTIVE_EPOCHS; i++) {
      CollectGarbageForPage(manager->first_pages_[i]);
    }
    LOG_INFO("Epoch Management Thread exiting");
    return nullptr;
  }
};

//===--------------------------------------------------------------------===//
// BWTree holds the mapping table and schema information of the index. For
// each scan/insert/delete operation, it pass the operation to the root.
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree {
 private:
  // the logical page id for the root node
  LPID root_;

  // the metadata of the index
  IndexMetadata *metadata_;

  // the key comparator
  KeyComparator comparator_;

  // the mapping table
  MappingTable<KeyType, ValueType, KeyComparator> *mapping_table_;

  // the epoch manager
  EpochManager<KeyType, ValueType, KeyComparator> epoch_manager_;

 public:
  /*
   * On construction, BWTree will create a IPage and an empty LPage.
   * The IPage has only one pointer, which points to the empty LPage.
   * The IPage serves as the root of all other nodes.
   */
  BWTree(IndexMetadata *metadata)
      : metadata_(metadata),
        comparator_(KeyComparator(metadata)),
        unique_keys(metadata->unique_keys) {
    LOG_INFO("Inside BWTree Constructor");
    mapping_table_ = new MappingTable<KeyType, ValueType, KeyComparator>();
    // create the root IPage
    KeyType dummy_key;
    IPage<KeyType, ValueType, KeyComparator> *root =
        new IPage<KeyType, ValueType, KeyComparator>(this, dummy_key, true);
    // install root IPage
    root_ = GetMappingTable()->InstallPage(root);
    LOG_INFO("Root got LPID: %lu", root_);

    // create the first LPage
    LPage<KeyType, ValueType, KeyComparator> *first_lpage =
        new LPage<KeyType, ValueType, KeyComparator>(this, dummy_key, true);

    std::pair<KeyType, LPID> first_lpage_pair;
    LPID first_lpage_lpid;

    // install the first LPage
    first_lpage_lpid = GetMappingTable()->InstallPage(first_lpage);
    LOG_INFO("The first LPage got LPID: %d", (int)first_lpage_lpid);

    first_lpage_pair.second = first_lpage_lpid;
    root->GetChildren()[0] = first_lpage_pair;

    LOG_INFO("Leaving BWTree constructor");
  };

  ~BWTree() { delete mapping_table_; }

  // get the index of the child at next level, which contains the given key
  int GetChild(KeyType key, std::pair<KeyType, LPID> *children, oid_t len);

  template <typename PairSecond>
  // get the index of the first occurrence of the given key
  // positive index indicates found, negative indicates not found.
  // 0 could be either case
  int BinarySearch(KeyType key, std::pair<KeyType, PairSecond> *locations,
                   oid_t len);

  // insert an entry to bwtree
  bool InsertEntry(KeyType key, ValueType location);

  // delete an entry from bwtree
  bool DeleteEntry(KeyType key, ValueType location) {
    LPID child_lpid;
    bool complete = false;
    // retry delete until it succeeds
    while (!complete) {
      LOG_INFO("BWtree::DeleteEntry - %s", ToString(key).c_str());
      auto epochNum = epoch_manager_.GetCurrentEpoch();
      child_lpid = root_;
      complete = GetMappingTable()
                     ->GetNode(child_lpid)
                     ->DeleteEntry(key, location, child_lpid);
      epoch_manager_.ReleaseEpoch(epochNum);
    }
    return true;
  };

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction,
                              const KeyType *index_key);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  // install a delta on a given node
  bool InstallParentDelta(
      IPageUpdateDelta<KeyType, ValueType, KeyComparator> *delta,
      KeyType right_most_key, bool right_most_key_is_infinity,
      LPID search_lpid);

  // print debug function
  void Debug();

  // check the integrity of the bwtree
  void BWTreeCheck();

  // perform cleanup to reduce memory footprint
  bool Cleanup() {
    auto epochnum = epoch_manager_.GetCurrentEpoch();
    std::cout << "Cleanup attempt 1:" << std::endl;
    GetMappingTable()->GetNode(root_)->Cleanup();

    std::cout << "Cleanup attempt 2:" << std::endl;
    GetMappingTable()->GetNode(root_)->Cleanup();

    // for now, do twice to handle if multiple adjacent nodes can be merged
    epoch_manager_.ReleaseEpoch(epochnum);
    return true;
  }

  // get the memory footprint of the bwtree
  size_t GetMemoryFootprint();

  // compress all the delta chains in the bwtree
  void CompressAllPages();

  // get the LPID of the root page
  inline LPID GetRootLPID() { return root_; }

 public:
  // whether unique key is required
  bool unique_keys;

  inline MappingTable<KeyType, ValueType, KeyComparator> *GetMappingTable() {
    return mapping_table_;
  }

  inline const catalog::Schema *GetKeySchema() {
    return metadata_->GetKeySchema();
  }

  inline std::string GetStringPrefix(std::string str, int len) {
    std::string prefix;
    for (int i = 0; i < len; i++) {
      prefix += str;
    }
    return prefix;
  }

  // get a string representation of the item pointer
  inline std::string ToString(ValueType val) {
    return "(" + std::to_string(val.block) + "," + std::to_string(val.offset) +
           ")";
  }

  // get a string representation of the key
  inline std::string ToString(KeyType key) {
    // assume it's a generic key, with an integer type as first column,
    // and varchar type as second column
    auto tuple = key.GetTupleForComparison(GetKeySchema());
    const Value &k1 = tuple.GetValue(0);
    const Value &k2 = tuple.GetValue(1);
    auto k1_str = k1.Debug();
    k1_str.erase(0, 9);

    auto k2_str = k2.Debug();
    k2_str.erase(0, 9);
    // remove type information
    k2_str = k2_str.substr(0, k2_str.find_last_of('['));
    k2_str = k2_str.substr(k2_str.find_first_of(']') + 1);
    if (k2_str.size() == 0) {
      k2_str = "?";
    }
    return "(" + k1_str + "," + k2_str + ")";
  }

  inline bool CompareValue(ValueType left, ValueType right) {
    return left.block == right.block && left.offset == right.offset;
  }

  // return 0 if equal, -1 if left < right, 1 otherwise
  inline int CompareKey(KeyType left, KeyType right) {
    bool less_than_right = comparator_(left, right);
    bool greater_than_right = comparator_(right, left);
    if (!less_than_right && !greater_than_right) {
      return 0;
    } else if (less_than_right) {
      return -1;
    }
    return 1;
  }

  // returns true if the key is not greater than the right most key
  inline bool KeyNotGreaterThan(KeyType key, KeyType right_most_key, bool inf) {
    return inf || CompareKey(key, right_most_key) <= 0;
  }

  // compress delta chain
  bool CompressDeltaChain(
      LPID page_to_compress,
      BWTreeNode<KeyType, ValueType, KeyComparator> *old_node_ptr,
      BWTreeNode<KeyType, ValueType, KeyComparator> *new_node) {
    LOG_INFO("Compressing delta chain for LPID: %lu", page_to_compress);
    auto new_node_state = new_node->BuildNodeState();
    auto new_node_ptr = new_node_state->GetPage();

    // delete the temporary state
    delete new_node_state;
    bool completed = GetMappingTable()->SwapNode(page_to_compress, old_node_ptr,
                                                 new_node_ptr);
    // if we didn't failed to install we should clean up the page we created
    if (completed) {
      this->epoch_manager_.AddNodeToEpoch(old_node_ptr);
    } else {
      // if we failed we should clean up the new page
      delete new_node_ptr;
    }

    return completed;
  }

  LPID *GetRootLPIDAddress() { return &root_; }

  std::vector<oid_t> ScanKeyInternal(KeyType key,
                                     std::pair<KeyType, ValueType> *locations,
                                     oid_t size);

  // helper function to perform scan on LPage
  LPID ScanHelper(const std::vector<Value> &values,
                  const std::vector<oid_t> &key_column_ids,
                  const std::vector<ExpressionType> &expr_types,
                  __attribute__((unused))
                  const ScanDirectionType &scan_direction,
                  const KeyType *index_key, std::vector<ValueType> &result,
                  std::pair<KeyType, ValueType> *locations, oid_t size,
                  LPID right_sibling);

  // helper function to perform scan on IPage
  LPID ScanHelper(const std::vector<Value> &values,
                  const std::vector<oid_t> &key_column_ids,
                  const std::vector<ExpressionType> &expr_types,
                  __attribute__((unused))
                  const ScanDirectionType &scan_direction,
                  const KeyType *index_key, std::vector<ValueType> &result,
                  std::pair<KeyType, LPID> *children, oid_t size);

  // helper function to perform scan all keys
  LPID ScanAllKeysHelper(oid_t size, std::pair<KeyType, ValueType> *locations,
                         oid_t right_sibling, std::vector<ValueType> &result);

  // helper function to perform scan key
  LPID ScanKeyHelper(KeyType key, oid_t size,
                     std::pair<KeyType, ValueType> *locations,
                     oid_t right_sibling, std::vector<ValueType> &result,
                     bool page_is_infinity, KeyType page_right_most_key);

 private:
  // check if the leading column is equal to the given key
  bool MatchLeadingColumn(const AbstractTuple &index_key,
                          const std::vector<oid_t> &key_column_ids,
                          const std::vector<Value> &values);
};

//===--------------------------------------------------------------------===//
// BWTreeNode is the abstract class for all tree nodes in the bwtree
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode {
 public:
  BWTreeNode(BWTree<KeyType, ValueType, KeyComparator> *map,
             int delta_chain_len, KeyType right_most_key, bool infinity)
      : map(map),
        delta_chain_len_(delta_chain_len),
        right_most_key(right_most_key),
        infinity(infinity){};

  virtual bool InsertEntry(KeyType key, ValueType location, LPID self) = 0;

  virtual bool DeleteEntry(KeyType key, ValueType location, LPID self,
                           bool update_parent = true) = 0;

  virtual bool AddINodeEntry(
      LPID, IPageUpdateDelta<KeyType, ValueType, KeyComparator> *) {
    return false;
  };

  virtual LPID Scan(const std::vector<Value> &values,
                    const std::vector<oid_t> &key_column_ids,
                    const std::vector<ExpressionType> &expr_types,
                    const ScanDirectionType &scan_direction,
                    std::vector<ValueType> &result,
                    const KeyType *index_key) = 0;

  virtual LPID ScanAllKeys(std::vector<ValueType> &result) = 0;

  virtual LPID ScanKey(KeyType key, std::vector<ValueType> &result) = 0;

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState() {
    return this->BuildNodeState(-1);
  }

  // build the new state of the node given the LPage/IPage/Delta
  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index) = 0;

  void SetCleanUpChildren() { clean_up_children_ = true; }

  virtual std::string Debug(int depth, LPID self) {
    return "Undefined type" + std::string(depth, '=') + std::to_string(self) +
           "\n";
  }

  virtual void BWTreeCheck() = 0;

  virtual size_t GetMemoryFootprint() = 0;

  // Each sub-class will have to implement this function to return their type
  virtual BWTreeNodeType GetTreeNodeType() const = 0;

  virtual ~BWTreeNode(){};

  inline int GetDeltaChainLen() {
    LOG_INFO("Got delta chain len %d", delta_chain_len_);
    return delta_chain_len_;
  }

  inline bool IsInifinity() { return this->infinity; }

  inline KeyType GetRightMostKey() { return this->right_most_key; }

  virtual bool Cleanup() { return true; }

  oid_t GetSize() { return 0; }

  virtual bool InstallParentDelta(
      LPID self, IPageUpdateDelta<KeyType, ValueType, KeyComparator> *delta,
      KeyType right_most_key, bool right_most_key_is_infinity,
      LPID search_lpid) = 0;

 protected:
  // perform delta insertion
  bool PerformDeltaInsert(
      LPID my_lpid, Delta<KeyType, ValueType, KeyComparator> *new_delta,
      BWTreeNode<KeyType, ValueType, KeyComparator> *old_delta) {
    bool status;
    LOG_INFO("Inside PerformDeltaChainInsert with new_delta len = %d",
             new_delta->GetDeltaChainLen());

    // check if it's over the delta chain length limit
    if (new_delta->GetDeltaChainLen() > new_delta->GetDeltaChainLimit()) {
      status = this->map->CompressDeltaChain(my_lpid, old_delta, new_delta);
      if (status) {
        delete new_delta;
      }
    } else {
      status =
          this->map->GetMappingTable()->SwapNode(my_lpid, old_delta, new_delta);
    }
    return status;
  }

  bool PerformDeltaInsert(LPID my_lpid,
                          Delta<KeyType, ValueType, KeyComparator> *new_delta) {
    return this->PerformDeltaInsert(my_lpid, new_delta, this);
  }

  // the handler to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  // length of delta chain
  int delta_chain_len_;

  bool clean_up_children_ = false;

  KeyType right_most_key;

  bool infinity = false;
};

//===--------------------------------------------------------------------===//
// IPage
// The IPage hold pointers to all its children. An IPage with n + 1 keys
// k1, k2, .... kn has (n + 1) pointers to its children, where
// p1 for (-infinity, k1), p2 for [k1, k2), p3 for [k2, k3) ...
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  IPage(BWTree<KeyType, ValueType, KeyComparator> *map, KeyType right_most_key,
        bool infinity)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0, right_most_key,
                                                      infinity) {
    size_ = 1;
    should_split_ = false;
  };

  // construct an IPage from a IPage Builder
  IPage(BWTree<KeyType, ValueType, KeyComparator> *map,
        NodeStateBuilder<KeyType, ValueType, KeyComparator> *state)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(
            map, 0, state->right_most_key, state->infinity) {
    size_ = state->size;
    INodeStateBuilder<KeyType, ValueType, KeyComparator> *ibuilder =
        reinterpret_cast<
            INodeStateBuilder<KeyType, ValueType, KeyComparator> *>(state);
    assert(size_ <= IPAGE_ARITY);
    // copy the children over
    memcpy(children_, ibuilder->children_,
           size_ * sizeof(std::pair<KeyType, LPID>));
    // set the split flag if threshold is exceeded
    should_split_ = size_ > IPAGE_SPLIT_THRESHOLD;
  };

  ~IPage() {}

  bool AddINodeEntry(
      LPID self,
      IPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta);

  bool InsertEntry(KeyType key, ValueType location, LPID self);

  bool DeleteEntry(KeyType key, ValueType location, LPID self, bool);

  LPID Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  LPID ScanAllKeys(std::vector<ValueType> &result);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  bool InstallParentDelta(
      LPID self, IPageUpdateDelta<KeyType, ValueType, KeyComparator> *delta,
      KeyType right_most_key, bool right_most_key_is_infinity,
      LPID search_lpid);

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  std::string Debug(int depth, LPID self);

  void BWTreeCheck();

  size_t GetMemoryFootprint() {
    LOG_INFO("IPage::GetMemoryFootprint");
    auto size = sizeof(IPage<KeyType, ValueType, KeyComparator>);
    return size;
  }

  // get the index of the child at next level, which contains the given key
  int GetChild(KeyType key, std::pair<KeyType, LPID> *children, oid_t len);

  inline std::pair<KeyType, LPID> *GetChildren() { return children_; }

  // split the IPage by creating a new IPage at right, installing a delta on
  // parent
  // and a split delta on top of itself
  void SplitNodes(LPID self);

  bool Cleanup() {
    int left_size, right_size;
    BWTreeNode<KeyType, ValueType, KeyComparator> *left_node, *right_node;
    LPage<KeyType, ValueType, KeyComparator> *left_lnode, *right_lnode;
    IPage<KeyType, ValueType, KeyComparator> *left_inode, *right_inode;
    LPID left_lpid, right_lpid;
    int i = 0;
    LOG_INFO("Cleanup invoked, size is %d", (int)this->size_);
    std::pair<KeyType, LPID> new_children[IPAGE_ARITY];
    int new_children_size = 0;

    for (int i = 0; i < (long)this->size_; i++) {
      this->map->GetMappingTable()
          ->GetNode(this->children_[i].second)
          ->Cleanup();
    }

    if (size_ >= 2)  // merges possible
    {
      LOG_INFO("Merge possible, size_ = %d", (int)size_);
      switch (this->map->GetMappingTable()
                  ->GetNode(this->children_[0].second)
                  ->GetTreeNodeType()) {
        case TYPE_LPAGE:
          LOG_INFO("Children are LPages");
          for (i = 0; i < (long)this->size_ - 1; i++) {
            left_lpid = this->children_[i].second;
            right_lpid = this->children_[i + 1].second;

            LOG_INFO("Left lpid is %d, right lpid is %d", (int)left_lpid,
                     (int)right_lpid);

            left_node = this->map->GetMappingTable()->GetNode(left_lpid);
            right_node = this->map->GetMappingTable()->GetNode(right_lpid);

            left_lnode =
                reinterpret_cast<LPage<KeyType, ValueType, KeyComparator> *>(
                    left_node);
            right_lnode =
                reinterpret_cast<LPage<KeyType, ValueType, KeyComparator> *>(
                    right_node);

            left_size = left_lnode->GetSizeForCleanup();
            right_size = right_lnode->GetSizeForCleanup();

            LOG_INFO("Left Lnode size is %d", (int)left_size);
            LOG_INFO("Right Lnode size is %d", (int)right_size);

            if (left_size + right_size < LPAGE_MERGE_THRESHOLD) {
              LOG_INFO("Can merge!");
              LPage<KeyType, ValueType, KeyComparator> *new_merged_lpage =
                  new LPage<KeyType, ValueType, KeyComparator>(
                      this->map, right_lnode->GetRightMostKey(),
                      right_lnode->IsInifinity());

              std::vector<std::pair<KeyType, ValueType>> left_node_list;
              std::vector<std::pair<KeyType, ValueType>> right_node_list;

              left_lnode->GetLocationsList(left_node_list);
              right_lnode->GetLocationsList(right_node_list);

              int new_lpage_size = 0;
              for (int index = 0; index < (long)left_node_list.size(); index++)
                new_merged_lpage->GetLocationsArray()[new_lpage_size++] =
                    left_node_list[index];

              for (int index = 0; index < (long)right_node_list.size(); index++)
                new_merged_lpage->GetLocationsArray()[new_lpage_size++] =
                    right_node_list[index];

              new_merged_lpage->SetSize(new_lpage_size);

              // now set the right sibling of this newly created LPage
              new_merged_lpage->SetRightSiblingLPID(
                  right_lnode->GetRightMostSibling());

              LOG_INFO("Size of new merged lpage is %d", (int)new_lpage_size);
              LOG_INFO("It's right most sibling is %d",
                       (int)new_merged_lpage->GetRightSiblingLPID());

              new_children[new_children_size].first =
                  right_lnode->GetRightMostKey();
              new_children[new_children_size].second = left_lpid;
              new_children_size++;

              this->map->GetMappingTable()->SwapNode(left_lpid, left_lnode,
                                                     new_merged_lpage);
              // TODO delete in chain, the two old lpages
              // TODO reclaim right lpid!
              i++;  // skip the merged node!

              left_lnode->DestroyRightSiblings();
              right_lnode->DestroyRightSiblings();
              this->map->GetMappingTable()->RemovePage(right_lpid);
              delete left_lnode;

            } else {
              new_children[new_children_size].first =
                  left_lnode->GetRightMostKey();
              new_children[new_children_size].second = left_lpid;
              new_children_size++;
            }
          }
          break;

        case TYPE_IPAGE:
          LOG_INFO("Children are IPages");
          for (i = 0; i < (long)this->size_ - 1; i++) {
            left_lpid = this->children_[i].second;
            right_lpid = this->children_[i + 1].second;

            LOG_INFO("Left lpid is %d", (int)left_lpid);
            LOG_INFO("Right lpid is %d", (int)right_lpid);

            left_node = this->map->GetMappingTable()->GetNode(left_lpid);

            right_node = this->map->GetMappingTable()->GetNode(right_lpid);

            left_inode =
                reinterpret_cast<IPage<KeyType, ValueType, KeyComparator> *>(
                    left_node);
            right_inode =
                reinterpret_cast<IPage<KeyType, ValueType, KeyComparator> *>(
                    right_node);

            left_size = left_inode->GetSize();
            right_size = right_inode->GetSize();

            LOG_INFO("Left size is %d", (int)left_size);
            LOG_INFO("Right size is %d", (int)right_size);

            if (left_size + right_size < IPAGE_MERGE_THRESHOLD) {
              LOG_INFO("Can Merge!");
              IPage<KeyType, ValueType, KeyComparator> *new_merged_ipage =
                  new IPage<KeyType, ValueType, KeyComparator>(
                      this->map, right_inode->GetRightMostKey(),
                      right_inode->IsInifinity());

              int new_ipage_size = 0;
              for (int index = 0; index < (long)left_inode->size_; index++)
                new_merged_ipage->children_[new_ipage_size++] =
                    left_inode->children_[index];

              for (int index = 0; index < (long)right_inode->size_; index++)
                new_merged_ipage->children_[new_ipage_size++] =
                    right_inode->children_[index];

              new_merged_ipage->size_ = new_ipage_size;
              LOG_INFO("New IPage size is %d", (int)new_merged_ipage->size_);

              new_children[new_children_size].first =
                  right_inode->GetRightMostKey();
              new_children[new_children_size].second = left_lpid;
              new_children_size++;
              LOG_INFO("New children size is %d", (int)new_children_size);

              this->map->GetMappingTable()->SwapNode(left_lpid, left_inode,
                                                     new_merged_ipage);

              delete left_inode;
              // delete right_inode;
              i++;  // skip the merged node!
              this->map->GetMappingTable()->RemovePage(right_lpid);

            } else {
              new_children[new_children_size].first =
                  left_inode->GetRightMostKey();
              new_children[new_children_size].second = left_lpid;
              new_children_size++;
            }
          }
          break;

        default:
          LOG_INFO(
              "Compression at some place has not been done! Can't merge. "
              "Self-Destruct...");
          assert(1 == 0);
          break;
      }
      if (i == (long)size_ - 1)
        new_children[new_children_size++] = children_[size_ - 1];

      LOG_INFO("All merges done! Size of new children array is %d",
               (int)new_children_size);
      for (int i = 0; i < new_children_size; i++) {
        this->children_[i] = new_children[i];
      }
      this->size_ = new_children_size;
    }

    return true;
  }

  oid_t GetSize() {
    LOG_INFO("Inside GetSize");
    return this->size_;
  }

 private:
  // the <key, LPID> pair of IPage's children
  std::pair<KeyType, LPID> children_[IPAGE_ARITY];

  oid_t size_;

  // a flag indicating whether split should be triggered on this IPage
  bool should_split_;
};

//===--------------------------------------------------------------------===//
// Delta is the abstract class for all deltas.
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class Delta : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  Delta(BWTree<KeyType, ValueType, KeyComparator> *map,
        BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
        KeyType right_most_key, bool infinity)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0, right_most_key,
                                                      infinity),
        modified_node(modified_node) {
    if (this->modified_node != nullptr) {
      this->delta_chain_len_ = modified_node->GetDeltaChainLen() + 1;
    }
  };

  LPID Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  LPID ScanAllKeys(std::vector<ValueType> &result);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  // set the linked modified node and increment the chain length
  void SetModifiedNode(
      BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node) {
    this->modified_node = modified_node;
    this->delta_chain_len_ = modified_node->GetDeltaChainLen() + 1;
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetModifiedNode() {
    return this->modified_node;
  }
  virtual ~Delta() {
    LOG_TRACE("Inside Delta Destructor");
    if (this->clean_up_children_) {
      LOG_TRACE("Have to clean up children, passing along the delete below");
      this->modified_node->SetCleanUpChildren();
      delete modified_node;
    }
  };

  virtual int GetDeltaChainLimit() = 0;

 protected:
  // the modified node could either be a LPage or IPage or Delta
  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node;

  void SetCleanUpChildren() {
    BWTreeNode<KeyType, ValueType, KeyComparator>::SetCleanUpChildren();
    modified_node->SetCleanUpChildren();
  }
};

//===--------------------------------------------------------------------===//
// IPageDelta is the abstract class of all delta's on IPages
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  IPageDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
             BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
             KeyType right_most_key, bool infinity)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node,
                                                 right_most_key, infinity){};

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  inline int GetDeltaChainLimit() { return IPAGE_DELTA_CHAIN_LIMIT; }
};

//===--------------------------------------------------------------------===//
// IPageSplitDelta stores the information of split on an IPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageSplitDelta : public IPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  IPageSplitDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType key, LPID value, int modified_index,
                  KeyType right_most_key, bool infinity)
      : IPageDelta<KeyType, ValueType, KeyComparator>(map, modified_node,
                                                      right_most_key, infinity),
        modified_key_(key),
        modified_val_(value),
        modified_index_(modified_index),
        split_completed_(false){};

  bool InsertEntry(KeyType key, ValueType location, LPID self);

  bool AddINodeEntry(
      LPID self,
      IPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta);

  bool DeleteEntry(KeyType key, ValueType location, LPID self, bool);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  bool InstallParentDelta(
      LPID self, IPageUpdateDelta<KeyType, ValueType, KeyComparator> *delta,
      KeyType right_most_key, bool right_most_key_is_infinity,
      LPID search_lpid);

  void SetSplitCompleted() { split_completed_ = true; }

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  std::string Debug(int depth, LPID self);

  void BWTreeCheck();

  size_t GetMemoryFootprint() {
    LOG_INFO("IPageSplitDelta::GetMemoryFootprint");
    size_t size = sizeof(IPageSplitDelta<KeyType, ValueType, KeyComparator>);
    size += this->modified_node->GetMemoryFootprint();
    return size;
  }

 private:
  // This key excluded in left child, included in the right child
  KeyType modified_key_;

  // The LPID of the new LPage
  LPID modified_val_;

  // index of last key in the left child;
  int modified_index_;

  bool split_completed_;
};

//===--------------------------------------------------------------------===//
// LPageDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  LPageDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
             BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
             KeyType right_most_key, bool infinity, LPID right_sibling)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node,
                                                 right_most_key, infinity),
        right_sibling(right_sibling){};

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  inline int GetDeltaChainLimit() { return LPAGE_DELTA_CHAIN_LIMIT; }

  inline bool InstallParentDelta(
      LPID, IPageUpdateDelta<KeyType, ValueType, KeyComparator> *, KeyType,
      bool, LPID) {
    return false;
  }

 protected:
  LPID right_sibling = INVALID_LPID;
};
//===--------------------------------------------------------------------===//
// LPageSplitDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageSplitDelta : public LPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  LPageSplitDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType splitterKey, int modified_index, LPID rightSplitPage,
                  KeyType right_most_key, bool infinity, LPID right_sibling)
      : LPageDelta<KeyType, ValueType, KeyComparator>(
            map, modified_node, right_most_key, infinity, right_sibling),
        modified_key_(splitterKey),
        modified_key_index_(modified_index),
        right_split_page_lpid_(rightSplitPage),
        split_completed_(false){};

  // This method will either try to create a delta on top of itself, if the
  // current key is less
  // than or equal to the modified_key_, or it will simply call InsertEntry on
  // the LPID of the
  // newly created right_split_page_lpid
  bool InsertEntry(KeyType key, ValueType location, LPID self);

  bool DeleteEntry(KeyType key, ValueType location, LPID self, bool);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  std::string Debug(int depth, LPID self);

  void BWTreeCheck();

  size_t GetMemoryFootprint() {
    LOG_INFO("LPageSplitDelta::GetMemoryFootprint");
    size_t size = sizeof(LPageSplitDelta<KeyType, ValueType, KeyComparator>);
    // calculate left page only
    size += this->modified_node->GetMemoryFootprint();
    return size;
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetModifiedNode() {
    return this->modified_node;
  }

  LPID GetRightSplitPageLPID() { return right_split_page_lpid_; }

  void SetSplitCompleted() { split_completed_ = true; }

 private:
  // This key included in left child, excluded in the right child
  KeyType modified_key_;
  int modified_key_index_;

  // The LPID of the new LPage
  LPID right_split_page_lpid_;

  bool split_completed_;
};

//===--------------------------------------------------------------------===//
// LPageUpdateDelta
// LPageUpdateDelta either represents a insert delta or delete delta on LPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageUpdateDelta : public LPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  LPageUpdateDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                   BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                   KeyType key, ValueType value, KeyType right_most_key,
                   bool infinity, LPID right_sibling)
      : LPageDelta<KeyType, ValueType, KeyComparator>(
            map, modified_node, right_most_key, infinity, right_sibling),
        modified_key_(key),
        modified_val_(value) {
    LOG_INFO("Inside LPageUpdateDelta Constructor");
  };

  bool InsertEntry(KeyType key, ValueType location, LPID my_lpid);

  bool DeleteEntry(KeyType key, ValueType location, LPID my_lpid, bool);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  void SetKey(KeyType modified_key) { modified_key_ = modified_key; }

  void SetValue(ValueType modified_val) { modified_val_ = modified_val; }

  void SetDeleteFlag() { is_delete_ = true; };

  std::string Debug(int depth, LPID self);

  void BWTreeCheck();

  size_t GetMemoryFootprint() {
    LOG_INFO("LPageUpdateDelta::GetMemoryFootprint");
    size_t size = sizeof(LPageUpdateDelta<KeyType, ValueType, KeyComparator>);
    size += this->modified_node->GetMemoryFootprint();
    return size;
  }

  inline void SetShouldSplit(bool should_split) {
    this->should_split_ = should_split;
  }

 private:
  // The key which is modified
  KeyType modified_key_;

  // The item pointer of the updated key
  ValueType modified_val_;

  // Whether it's a delete delta
  bool is_delete_ = false;

  bool should_split_ = false;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class LPageRemoveDelta : public LPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  LPageRemoveDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                   KeyType right_most_key, bool infinity, LPID right_sibling)
      : LPageDelta<KeyType, ValueType, KeyComparator>(map, right_most_key,
                                                      infinity, right_sibling) {
    LOG_INFO("Inside LPageRemoveDelta Constructor");
  };

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID my_lpid);

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location, LPID my_lpid,
                   bool);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  size_t GetMemoryFootprint() {
    // assume we never call GetFootprint on Remove Delta
    return 0;
  }

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

 private:
  // Whether it's a delete delta
  bool is_delete_ = false;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class LPageMergeDelta : public LPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  LPageMergeDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType right_most_key, bool infinity, LPID right_sibling)
      : LPageDelta<KeyType, ValueType, KeyComparator>(
            map, modified_node, right_most_key, infinity, right_sibling) {
    LOG_INFO("Inside LPageMergeDelta Constructor");
  };

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID my_lpid);

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location, LPID my_lpid,
                   bool);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  size_t GetMemoryFootprint() {
    // assume we never call GetFootprint on Merge Delta
    return 0;
  }

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

 private:
  // Whether it's a delete delta
  bool is_delete_ = false;
};

//===--------------------------------------------------------------------===//
// IPageUpdateDelta
//
// IPageUpdateDelta either represents a insert delta or delete delta on IPage
// The insert delta is essentially a separator in the case of a split SMO,
// where the new LPID has to be inserted. The delete delta is created in the
// case of a merge SMO, where an old LPID has to be removed.
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageUpdateDelta : public IPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  IPageUpdateDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                   BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                   KeyType max_key_left_split_node,
                   KeyType max_key_right_split_node,
                   bool right_node_is_infinity, LPID left_split_node_lpid,
                   LPID right_split_node_lpid, bool is_delete,
                   KeyType right_most_key, bool infinity)
      : IPageDelta<KeyType, ValueType, KeyComparator>(map, modified_node,
                                                      right_most_key, infinity),
        max_key_left_split_node_(max_key_left_split_node),
        max_key_right_split_node_(max_key_right_split_node),
        left_split_node_lpid_(left_split_node_lpid),
        right_split_node_lpid_(right_split_node_lpid),
        right_node_is_infinity_(right_node_is_infinity),
        is_delete_(is_delete) {
    LOG_INFO("Inside IPageUpdateDelta Constructor");
  };

  bool AddINodeEntry(
      LPID self,
      IPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta);

  bool InsertEntry(KeyType key, ValueType location, LPID self);

  bool DeleteEntry(KeyType key, ValueType location, LPID self, bool);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  std::string Debug(int depth, LPID self);

  void BWTreeCheck();

  size_t GetMemoryFootprint() {
    LOG_INFO("IPageUpdateDelta::GetMemoryFootprint");
    size_t size = sizeof(IPageUpdateDelta<KeyType, ValueType, KeyComparator>);
    size += this->modified_node->GetMemoryFootprint();
    return size;
  }

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetModifiedNode() {
    return this->modified_node;
  }

  void SetRightMostKey(KeyType right_most_key) {
    this->right_most_key = right_most_key;
  }

  void SetInfinity(bool infinity) { this->infinity = infinity; }

  bool InstallParentDelta(
      LPID self, IPageUpdateDelta<KeyType, ValueType, KeyComparator> *delta,
      KeyType right_most_key, bool right_most_key_is_infinity,
      LPID search_lpid);

 private:
  // The key which is modified
  KeyType max_key_left_split_node_, max_key_right_split_node_;

  // two members
  LPID left_split_node_lpid_, right_split_node_lpid_;

  bool right_node_is_infinity_;
  // Whether it's a delete delta
  bool is_delete_ = false;
};

// TODO More delta classes such as
// merge, remove_page

//===--------------------------------------------------------------------===//
// LPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  // get the index of the first occurrence of the given <k, v> pair
  // used when deleting a <k-v> entry from non-unique keys
  static int BinarySearch(std::pair<KeyType, ValueType> pair,
                          std::pair<KeyType, ValueType> *locations, oid_t len);

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map, KeyType right_most_key,
        bool infinity)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0, right_most_key,
                                                      infinity) {
    left_sib_ = INVALID_LPID;
    right_sib_ = INVALID_LPID;
    size_ = 0;
  };

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map,
        NodeStateBuilder<KeyType, ValueType, KeyComparator> *state)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(
            map, 0, state->right_most_key, state->infinity) {
    LOG_INFO(" ");
    size_ = state->size;
    LNodeStateBuilder<KeyType, ValueType, KeyComparator> *lstate =
        reinterpret_cast<
            LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(state);
    assert(size_ < LPAGE_ARITY);
    memcpy(locations_, lstate->locations_,
           sizeof(std::pair<KeyType, ValueType>) * size_);
    left_sib_ = lstate->left_sibling_;
    if (state->IsSeparated()) {
      right_sib_ = state->split_new_page_id;
    } else {
      right_sib_ = lstate->right_sibling_;
    }
  };

  ~LPage(){
      // LOG_INFO("Destroying an LPage");
  };

  void GetLocationsList(
      std::vector<std::pair<KeyType, ValueType>> &children_list) {
    for (int i = 0; i < (long)size_; i++) {
      children_list.push_back(locations_[i]);
    }

    LPID right_sibling_lpid = right_sib_;
    BWTreeNode<KeyType, ValueType, KeyComparator> *next_node;
    LPage<KeyType, ValueType, KeyComparator> *next_lpage;

    while (right_sibling_lpid != INVALID_LPID) {
      next_node = this->map->GetMappingTable()->GetNode(right_sibling_lpid);
      next_lpage = reinterpret_cast<LPage<KeyType, ValueType, KeyComparator> *>(
          next_node);

      if (this->map->CompareKey(next_lpage->GetRightMostKey(),
                                this->right_most_key) != 0)
        break;
      for (int i = 0; i < next_lpage->GetSize(); i++) {
        children_list.push_back(next_lpage->locations_[i]);
      }
      right_sibling_lpid = next_lpage->right_sib_;
    }
  }

  bool InsertEntry(KeyType key, ValueType location, LPID self);

  bool DeleteEntry(KeyType key, ValueType location, LPID self, bool);

  LPID Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  LPID ScanAllKeys(std::vector<ValueType> &result);

  LPID ScanKey(KeyType key, std::vector<ValueType> &result);

  std::string Debug(int depth, LPID self);

  void BWTreeCheck();

  size_t GetMemoryFootprint() {
    LOG_INFO("LPage::GetMemoryFootprint");
    auto size = sizeof(LPage<KeyType, ValueType, KeyComparator>);
    return size;
  }

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  bool SplitNodes(LPID self);

  void MergeNodes(LPID self, LPID right_sibling_lpid);

  bool InstallParentDelta(LPID,
                          IPageUpdateDelta<KeyType, ValueType, KeyComparator> *,
                          KeyType, bool, LPID) {
    return false;
  }

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  inline LPID GetLeftSiblingLPID() { return left_sib_; }

  inline LPID GetRightSiblingLPID() { return right_sib_; }

  void SetRightSiblingLPID(LPID new_right_sib) { right_sib_ = new_right_sib; }

  std::pair<KeyType, ValueType> *GetLocationsArray() { return locations_; }

  void SetSize(int size) { size_ = size; }

  int GetSize() { return size_; }

  int GetSizeForCleanup() {
    int total_size;
    total_size = size_;
    LPID right_sibling_lpid = right_sib_;

    BWTreeNode<KeyType, ValueType, KeyComparator> *next_node;
    LPage<KeyType, ValueType, KeyComparator> *next_lpage;

    while (right_sibling_lpid != INVALID_LPID) {
      next_node = this->map->GetMappingTable()->GetNode(right_sibling_lpid);
      next_lpage = reinterpret_cast<LPage<KeyType, ValueType, KeyComparator> *>(
          next_node);
      if (this->map->CompareKey(next_lpage->GetRightMostKey(),
                                this->right_most_key) != 0)
        break;
      total_size += next_lpage->GetSize();
      right_sibling_lpid = next_lpage->right_sib_;
    }
    return total_size;
  }

  LPID GetRightMostSibling() {
    LPID curr_right_sibling_lpid;
    curr_right_sibling_lpid = right_sib_;
    LPage<KeyType, ValueType, KeyComparator> *sibling;

    while (curr_right_sibling_lpid != INVALID_LPID) {
      sibling = reinterpret_cast<LPage<KeyType, ValueType, KeyComparator> *>(
          this->map->GetMappingTable()->GetNode(curr_right_sibling_lpid));

      if (this->map->CompareKey(this->right_most_key,
                                sibling->right_most_key) == 0)
        curr_right_sibling_lpid = sibling->right_sib_;
      else
        break;
    }
    return curr_right_sibling_lpid;
  }

  void DestroyRightSiblings() {
    LPID curr_right_sibling_lpid, old_right_sibling;
    curr_right_sibling_lpid = right_sib_;
    LPage<KeyType, ValueType, KeyComparator> *sibling;

    while (curr_right_sibling_lpid != INVALID_LPID) {
      sibling = reinterpret_cast<LPage<KeyType, ValueType, KeyComparator> *>(
          this->map->GetMappingTable()->GetNode(curr_right_sibling_lpid));

      if (this->map->CompareKey(this->right_most_key,
                                sibling->right_most_key) == 0) {
        old_right_sibling = curr_right_sibling_lpid;
        curr_right_sibling_lpid = sibling->right_sib_;
        this->map->GetMappingTable()->RemovePage(old_right_sibling);
      } else
        break;
    }
  }

 private:
  // return a vector of indices of the matched slots
  std::vector<oid_t> ScanKeyInternal(KeyType key);

  // left_sibling pointer is used for reverse scan
  LPID left_sib_;
  LPID right_sib_;

  // the key-value pairs
  std::pair<KeyType, ValueType> locations_[LPAGE_ARITY];

  // the size of stored locations
  oid_t size_;
};

}  // End index namespace
}  // End peloton namespace
