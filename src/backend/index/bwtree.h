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
#include <map>
#include <vector>
#include <climits>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <thread>

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

#define IPAGE_ARITY 256
#define LPAGE_ARITY 256

#define LPAGE_SPLIT_THRESHOLD 64

#define INVALID_LPID ULLONG_MAX

#define LPAGE_DELTA_CHAIN_LIMIT 5
#define IPAGE_DELTA_CHAIN_LIMIT 5

#define EPOCH_PAGE_SIZE 256
#define MAX_ACTIVE_EPOCHS 2
#define EPOCH_LENGTH_MILLIS 40

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode;

template <typename KeyType, typename ValueType, class KeyComparator>
class IPage;

template <typename KeyType, typename ValueType, class KeyComparator>
class LPage;

template <typename KeyType, typename ValueType, class KeyComparator>
class EpochManager;

//===--------------------------------------------------------------------===//
// NodeStateBuilder
//===--------------------------------------------------------------------===//
// TODO add access methods for LNode scan
// TODO add methods for use by split deltas, etc
// TODO add constructors to LNode and INode to build node based on state
// TODO add methods on each node to build the NodeState
/**
 * Builder for a node state, to be used with Delta Compression and
 * scans on indexes with multiple keys
 */
template <typename KeyType, typename ValueType, class KeyComparator>
class NodeStateBuilder {
 protected:
  // number of k-v pairs. everything after size is ignored
  oid_t size;

  // pointer to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  bool is_separated = false;

  KeyType separator_key;
  LPID split_new_page_id = 0;

 public:
  NodeStateBuilder(oid_t size, BWTree<KeyType, ValueType, KeyComparator> *map)
      : size(size), map(map){};

  virtual BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage() = 0;

  virtual void Scan(__attribute__((unused)) const std::vector<Value> &values,
                    __attribute__((unused))
                    const std::vector<oid_t> &key_column_ids,

                    __attribute__((unused))
                    const std::vector<ExpressionType> &expr_types,
                    __attribute__((unused))
                    const ScanDirectionType &scan_direction,
                    __attribute__((unused)) std::vector<ValueType> &result,
                    __attribute__((unused)) const KeyType *index_key){};

  virtual void ScanAllKeys(std::vector<ValueType> &result) = 0;

  virtual void ScanKey(__attribute__((unused)) KeyType key,
                       __attribute__((unused))
                       std::vector<ValueType> &result){};

  inline bool IsSeparated() { return is_separated; }

  inline KeyType GetSeparatorKey() { return separator_key; }

  inline LPID GetSplitNewPageId() { return split_new_page_id; }

  virtual ~NodeStateBuilder() { ; };

  friend class LPage<KeyType, ValueType, KeyComparator>;
  friend class IPage<KeyType, ValueType, KeyComparator>;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class INodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // IPage children nodes
  std::pair<KeyType, LPID> children_[IPAGE_ARITY + IPAGE_DELTA_CHAIN_LIMIT];

 public:
  // IPage constructor
  INodeStateBuilder(std::pair<KeyType, LPID> *children, int children_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(children_len, map) {
    memcpy(children_, children,
           sizeof(std::pair<KeyType, LPID>) * children_len);
  };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage();

  ~INodeStateBuilder(){};

  //***************************************************
  // IPage Methods
  //***************************************************

  /* AddChild adds a new <key, LPID> pair to its children if key is not
   * found. Otherwise overwrite the pair with the same key */
  void AddChild(std::pair<KeyType, LPID> &new_pair);

  void RemoveChild(KeyType &key_to_remove);

  void SeparateFromKey(KeyType separator_key, LPID split_new_page_id);

  void ScanAllKeys(std::vector<ValueType> &result);

  friend class LPage<KeyType, ValueType, KeyComparator>;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class LNodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // LPage members
  LPID left_sibling_ = INVALID_LPID;
  LPID right_sibling_ = INVALID_LPID;
  int separator_index_ = -1;

  // LPage members
  std::pair<KeyType, ValueType>
      locations_[IPAGE_ARITY + LPAGE_DELTA_CHAIN_LIMIT];

 public:
  // LPage constructor
  LNodeStateBuilder(LPID left_sibling, LPID right_sibling,
                    std::pair<KeyType, ValueType> *locations, int location_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(location_len, map),
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

  void UpdateLeftSib(LPID new_left_sib) { left_sibling_ = new_left_sib; }

  void UpdateRightSib(LPID new_right_sib) { right_sibling_ = new_right_sib; }

  void AddLeafData(std::pair<KeyType, ValueType> &new_pair);

  void RemoveLeafData(std::pair<KeyType, ValueType> &entry_to_remove);

  void SeparateFromKey(KeyType separator_key, int index,
                       LPID split_new_page_id);

  /*
   * Methods for Scan
   */
  void Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,

            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  void ScanAllKeys(std::vector<ValueType> &result);

  void ScanKey(KeyType key, std::vector<ValueType> &result);

 private:
  bool ItemPointerEquals(ValueType v1, ValueType v2);

  friend class LPage<KeyType, ValueType, KeyComparator>;
};

//===--------------------------------------------------------------------===//
// ReadWriteLatch
//===--------------------------------------------------------------------===//
class ReadWriteLatch {
 public:
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
  inline void ReleaseRead() {
    LOG_INFO("Releasing read Lock");
    __sync_add_and_fetch(&current_readers_, -1);
  }
  inline void AquireWrite() {
    LOG_INFO("Aquiring write Lock");
    while (!__sync_bool_compare_and_swap(&current_writers_, 0, 1))
      ;
    while (current_readers_ > 0)
      ;
  }
  inline void ReleaseWrite() {
    LOG_INFO("Releasing write Lock");
    assert(__sync_bool_compare_and_swap(&current_writers_, 1, 0));
  }

 private:
  int current_readers_ = 0;
  int current_writers_ = 0;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class MappingTable {
 private:
  // the mapping table
  unsigned int mapping_table_cap_ = 128;
  unsigned int mapping_table_size_ = 0;

  BWTreeNode<KeyType, ValueType, KeyComparator> **mapping_table_;
  LPID next_LPID_ = 0;
  ReadWriteLatch latch_;

 public:
  MappingTable() {
    LOG_INFO("Constructing Mapping Table with initial capacity %d",
             mapping_table_cap_);
    mapping_table_ =
        new BWTreeNode<KeyType, ValueType, KeyComparator> *[mapping_table_cap_];
  };
  ~MappingTable() { delete[] mapping_table_; }
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

      int new_mapping_table_cap = mapping_table_cap_ * 2;
      LOG_INFO("doubling size of mapping table capacity from %d to %d",
               mapping_table_cap_, new_mapping_table_cap);
      auto new_mapping_table =
          new BWTreeNode<KeyType, ValueType,
                         KeyComparator> *[new_mapping_table_cap];
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

  bool SwapNode(LPID id, BWTreeNode<KeyType, ValueType, KeyComparator> *oldNode,
                BWTreeNode<KeyType, ValueType, KeyComparator> *newNode) {
    LOG_INFO("swapping node for LPID: %lu into mapping table", id);
    latch_.AquireRead();
    bool ret =
        __sync_bool_compare_and_swap(mapping_table_ + id, oldNode, newNode);
    latch_.ReleaseRead();
    return ret;
  }

  // assumes that LPID is valid
  BWTreeNode<KeyType, ValueType, KeyComparator> *GetNode(LPID id) {
    LOG_INFO("getting node for LPID: %lu frommapping table", id);
    latch_.AquireRead();
    auto ret = mapping_table_[id];
    latch_.ReleaseRead();
    return ret;
  }
};

template <typename KeyType, typename ValueType, class KeyComparator>
struct EpochPage {
  BWTreeNode<KeyType, ValueType, KeyComparator> *pointers[EPOCH_PAGE_SIZE];
  EpochPage *next_page = nullptr;
};

//===--------------------------------------------------------------------===//
// EpochManager
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class EpochManager {
 public:
  EpochManager() {
    // should start GC pthread
    current_epoch_ = 0;
    destructor_called_ = false;
    memset(&epoch_size_, 0, MAX_ACTIVE_EPOCHS * sizeof(int));
    memset(&epoch_users_, 0, MAX_ACTIVE_EPOCHS * sizeof(int));
    pthread_create(&management_pthread_, nullptr, &epoch_management, this);
    for (int i = 0; i < MAX_ACTIVE_EPOCHS; i++) {
      first_pages_[i] = new EpochPage<KeyType, ValueType, KeyComparator>();
    }
  }

  ~EpochManager() {
    destructor_called_ = true;
    LOG_INFO("set destructor called flag");
    pthread_join(management_pthread_, nullptr);
  }

  unsigned long GetCurrentEpoch() {
    auto current_epoch = current_epoch_;
    __sync_add_and_fetch(&epoch_users_[current_epoch], 1);
    LOG_INFO("Got current epoch %lu", current_epoch);
    return current_epoch;
  }
  void ReleaseEpoch(unsigned long i) {
    __sync_add_and_fetch(&epoch_users_[i], -1);
  }

  void AddNodeToEpoch(
      BWTreeNode<KeyType, ValueType, KeyComparator> *node_to_destroy) {
    LOG_INFO("Adding node to epoch");
    auto curr_epoch = current_epoch_;
    auto write_pos = __sync_fetch_and_add(&epoch_size_[curr_epoch], 1);
    // find correct page
    auto curr_page = first_pages_[curr_epoch];
    for (int i = 0; i < (write_pos / EPOCH_PAGE_SIZE) + 1; i++) {
      LOG_INFO("Epoch Page Full, getting child page");
      if (curr_page->next_page == nullptr) {
        EpochPage<KeyType, ValueType, KeyComparator> *new_page =
            new EpochPage<KeyType, ValueType, KeyComparator>();
        if (!__sync_bool_compare_and_swap(&(curr_page->next_page), nullptr,
                                          new_page)) {
          delete new_page;
        }
      }
      curr_page = curr_page->next_page;
    }
    curr_page->pointers[write_pos] = node_to_destroy;
  }

 private:
  EpochPage<KeyType, ValueType, KeyComparator> *first_pages_[MAX_ACTIVE_EPOCHS];
  int epoch_size_[MAX_ACTIVE_EPOCHS];
  int epoch_users_[MAX_ACTIVE_EPOCHS];
  unsigned long current_epoch_;
  bool destructor_called_;
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
    LOG_INFO("Epoch Management Thread cleaing up remaining garbage");
    for (int i = 0; i < MAX_ACTIVE_EPOCHS; i++) {
      CollectGarbageForPage(manager->first_pages_[i]);
    }
    LOG_INFO("Epoch Management Thread exiting");
    return nullptr;
  }
};

//===--------------------------------------------------------------------===//
// BWTree
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree {
 private:
  // the logical page id for the root node
  LPID root_;
  IndexMetadata *metadata_;
  KeyComparator comparator_;
  MappingTable<KeyType, ValueType, KeyComparator> *mapping_table_;
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

    IPage<KeyType, ValueType, KeyComparator> *root =
        new IPage<KeyType, ValueType, KeyComparator>(this);

    root_ = GetMappingTable()->InstallPage(root);

    LOG_INFO("Root got LPID: %lu", root_);

    LPage<KeyType, ValueType, KeyComparator> *first_lpage =
        new LPage<KeyType, ValueType, KeyComparator>(this);

    std::pair<KeyType, LPID> first_lpage_pair;

    LPID first_lpage_lpid;

    first_lpage_lpid = GetMappingTable()->InstallPage(first_lpage);

    LOG_INFO("The first LPage got LPID: %d", (int)first_lpage_lpid);
    first_lpage_pair.second = first_lpage_lpid;

    root->GetChildren()[0] = first_lpage_pair;

    LOG_INFO("Leaving BWTree constructor");

    // with the given comparator
  };

  ~BWTree() { delete mapping_table_; }

  // get the index of the first occurrence of the given key
  template <typename PairSecond>

  // positive index indicates found, negative indicates not found. 0 could be
  // either case
  int BinarySearch(KeyType key,

                   std::pair<KeyType, PairSecond> *locations, oid_t len);

  bool InsertEntry(KeyType key, ValueType location);

  bool DeleteEntry(KeyType key, ValueType location) {
    LPID child_lpid;
    bool complete = false;
    while (!complete) {
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

 public:
  // whether unique key is required
  bool unique_keys;

  // return 0 if the page install is not successful

  inline MappingTable<KeyType, ValueType, KeyComparator> *GetMappingTable() {
    return mapping_table_;
  }

  inline const catalog::Schema *GetKeySchema() {
    return metadata_->GetKeySchema();
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
      // TODO add this back when other memory problems fixed

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

  void ScanHelper(const std::vector<Value> &values,
                  const std::vector<oid_t> &key_column_ids,
                  const std::vector<ExpressionType> &expr_types,
                  const ScanDirectionType &scan_direction,
                  const KeyType *index_key, std::vector<ValueType> &result,
                  std::pair<KeyType, ValueType> *locations, oid_t size,
                  LPID right_sibling);

  void ScanAllKeysHelper(oid_t size, std::pair<KeyType, ValueType> *locations,
                         oid_t right_sibling, std::vector<ValueType> &result);

  void ScanKeyHelper(KeyType key, oid_t size,
                     std::pair<KeyType, ValueType> *locations,
                     oid_t right_sibling, std::vector<ValueType> &result);
};

//===--------------------------------------------------------------------===//
// BWTreeNode
// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode {
 public:
  BWTreeNode(BWTree<KeyType, ValueType, KeyComparator> *map,
             int delta_chain_len)
      : map(map), delta_chain_len_(delta_chain_len){};

  // These are virtual methods which child classes have to implement.
  // They also have to be redeclared in the child classes
  virtual bool InsertEntry(KeyType key, ValueType location, LPID self,
                           LPID parent) = 0;

  virtual bool DeleteEntry(KeyType key, ValueType location, LPID self) = 0;

  virtual void Scan(const std::vector<Value> &values,
                    const std::vector<oid_t> &key_column_ids,
                    const std::vector<ExpressionType> &expr_types,
                    const ScanDirectionType &scan_direction,
                    std::vector<ValueType> &result,
                    const KeyType *index_key) = 0;

  virtual void ScanAllKeys(std::vector<ValueType> &result) = 0;

  virtual void ScanKey(KeyType key, std::vector<ValueType> &result) = 0;

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState() {
    return this->BuildNodeState(-1);
  }

  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index) = 0;

  // for scanKey, we have to build the node state as well. But we only care
  // about the keys we want to scan
  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      __attribute__((unused)) KeyType key) {
    return nullptr;
  };

  void SetCleanUpChildren() { clean_up_children_ = true; }

  // TODO to be implemented in the future
  //    virtual NodeStateBuilder<KeyType, ValueType, KeyComparator>
  //    *BuildScanState(
  //        __attribute__((unused)) const std::vector<Value> &values,
  //        __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
  //        __attribute__((unused)) const std::vector<ExpressionType>
  //        &expr_types,
  //        __attribute__((unused)) const ScanDirectionType &scan_direction) {
  //      return nullptr;
  //    };

  // Each sub-class will have to implement this function to return their type
  virtual BWTreeNodeType GetTreeNodeType() const = 0;

  virtual ~BWTreeNode(){

  };

  inline int GetDeltaChainLen() {
    LOG_INFO("got delta chain len %d\n", delta_chain_len_);
    return delta_chain_len_;
  }

 protected:
  // the handler to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  // length of delta chain
  int delta_chain_len_;

  bool clean_up_children_ = false;
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
  IPage(BWTree<KeyType, ValueType, KeyComparator> *map)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    size_ = 0;
    // children_ = new std::pair<KeyType, LPID>();
  };

  ~IPage(){};

  bool InsertEntry(KeyType key, ValueType location, LPID self, LPID parent);

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    int child_index = GetChild(key, children_, size_);
    LPID child_lpid = this->children_[child_index].second;
    return this->map->GetMappingTable()
        ->GetNode(child_lpid)
        ->DeleteEntry(key, location, child_lpid);
  };

  void Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  void ScanAllKeys(std::vector<ValueType> &result);

  void ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  // get the index of the child at next level, which contains the given key
  int GetChild(__attribute__((unused)) KeyType key,
               __attribute__((unused)) std::pair<KeyType, LPID> *children,
               __attribute__((unused)) oid_t len);

  std::pair<KeyType, LPID> *GetChildren() { return children_; }

  void SplitNodes(LPID self, LPID parent);

 private:
  std::pair<KeyType, LPID> children_[IPAGE_ARITY];

  oid_t size_;
};

//===--------------------------------------------------------------------===//
// Delta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class Delta : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  Delta(BWTree<KeyType, ValueType, KeyComparator> *map,
        BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(
            map, modified_node->GetDeltaChainLen() + 1),
        modified_node(modified_node){};

  void Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  void ScanAllKeys(std::vector<ValueType> &result);

  void ScanKey(KeyType key, std::vector<ValueType> &result);

  ~Delta() {
    if (this->clean_up_children_) {
      delete modified_node;
    }
  };

 protected:
  // the modified node could either be a LPage or IPage or Delta
  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node;

  virtual int GetDeltaChainLimit() = 0;

  bool PerformDeltaInsert(LPID my_lpid,
                          Delta<KeyType, ValueType, KeyComparator> *new_delta) {
    bool status;
    if (new_delta->GetDeltaChainLen() > this->GetDeltaChainLimit()) {
      status = this->map->CompressDeltaChain(my_lpid, this, new_delta);
    } else {
      status = this->map->GetMappingTable()->SwapNode(my_lpid, this, new_delta);
    }
    return status;
  }

  void SetCleanUpChildren() {
    BWTreeNode<KeyType, ValueType, KeyComparator>::SetCleanUpChildren();
    modified_node->SetCleanUpChildren();
  }
};

//===--------------------------------------------------------------------===//
// IPageDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  IPageDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
             BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node){};

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  inline int GetDeltaChainLimit() { return IPAGE_DELTA_CHAIN_LIMIT; }
};

//===--------------------------------------------------------------------===//
// IPageSplitDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageSplitDelta : public IPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  IPageSplitDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType key, LPID value, int modified_index)
      : IPageDelta<KeyType, ValueType, KeyComparator>(map, modified_node),
        modified_key_(key),
        modified_val_(value),
        modified_index_(modified_index){};

  bool InsertEntry(KeyType key, ValueType location, LPID self, LPID parent);

  bool DeleteEntry(KeyType key, ValueType location, LPID self);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

 private:
  // This key excluded in left child, included in the right child
  KeyType modified_key_;

  // The LPID of the new LPage
  LPID modified_val_;
  // index of last key in the left child;
  int modified_index_;
};

//===--------------------------------------------------------------------===//
// IPageDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  LPageDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
             BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node){};

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  inline int GetDeltaChainLimit() { return LPAGE_DELTA_CHAIN_LIMIT; }
};
//===--------------------------------------------------------------------===//
// LPageSplitDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageSplitDelta : public LPageDelta<KeyType, ValueType, KeyComparator> {
 public:
  LPageSplitDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType splitterKey, int modified_index, LPID rightSplitPage)
      : LPageDelta<KeyType, ValueType, KeyComparator>(map, modified_node),
        modified_key_(splitterKey),
        modified_key_index_(modified_index),
        right_split_page_lpid_(rightSplitPage){};

  // This method will either try to create a delta on top of itself, if the
  // current key is less
  // than or equal to the modified_key_, or it will simply call InsertEntry on
  // the LPID of the
  // newly created right_split_page_lpid
  bool InsertEntry(KeyType key, ValueType location, LPID self, LPID parent);

  bool DeleteEntry(KeyType key, ValueType location, LPID self);

  void ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetModifiedNode() {
    return this->modified_node;
  }

  LPID GetRightSplitPageLPID() { return right_split_page_lpid_; }

 private:
  // This key included in left child, excluded in the right child
  KeyType modified_key_;
  int modified_key_index_;

  // The LPID of the new LPage
  LPID right_split_page_lpid_;
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
                   KeyType key, ValueType value)
      : LPageDelta<KeyType, ValueType, KeyComparator>(map, modified_node),
        modified_key_(key),
        modified_val_(value) {
    LOG_INFO("Inside LPageUpdateDelta Constructor");
  };

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID my_lpid,
                   __attribute__((unused)) LPID parent) {
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
        new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
                                                                key, location);
    bool status = this->PerformDeltaInsert(my_lpid, new_delta);
    if (!status) {
      delete new_delta;
    }
    return status;
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location, LPID my_lpid) {
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
        new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
                                                                key, location);

    new_delta->SetDeleteFlag();
    bool status = this->PerformDeltaInsert(my_lpid, new_delta);
    if (!status) {
      delete new_delta;
    }
    return status;
  };

  void ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      __attribute__((unused)) KeyType key) {
    return nullptr;
  };

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      __attribute__((unused)) const std::vector<Value> &values,
      __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
      __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
      __attribute__((unused)) const ScanDirectionType &scan_direction) {
    return nullptr;
  };

  void SetKey(KeyType modified_key) { modified_key_ = modified_key; }

  void SetValue(ValueType modified_val) { modified_val_ = modified_val; }

  void SetDeleteFlag() { is_delete_ = true; };

 private:
  // The key which is modified
  KeyType modified_key_;

  // The item pointer of the updated key
  ValueType modified_val_;

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
                   KeyType max_key_right_split_node, LPID left_split_node_lpid,
                   LPID right_split_node_lpid, bool is_delete)
      : IPageDelta<KeyType, ValueType, KeyComparator>(map, modified_node),
        max_key_left_split_node_(max_key_left_split_node),
        max_key_right_split_node_(max_key_right_split_node),
        left_split_node_lpid_(left_split_node_lpid),
        right_split_node_lpid_(right_split_node_lpid),
        is_delete_(is_delete) {
    LOG_INFO("Inside IPageUpdateDelta Constructor");
  };

  bool InsertEntry(KeyType key, ValueType location, LPID self, LPID parent);

  bool DeleteEntry(KeyType key, ValueType location, LPID self);

  void ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetModifiedNode() {
    return this->modified_node;
  }

 private:
  // The key which is modified
  KeyType max_key_left_split_node_, max_key_right_split_node_;

  // two members
  LPID left_split_node_lpid_, right_split_node_lpid_;

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

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    left_sib_ = INVALID_LPID;
    right_sib_ = INVALID_LPID;
    size_ = 0;
  };

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map,
        NodeStateBuilder<KeyType, ValueType, KeyComparator> *state)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    LOG_INFO(" ");
    size_ = state->size;
    LNodeStateBuilder<KeyType, ValueType, KeyComparator> *lstate =
        reinterpret_cast<
            LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(state);

    memcpy(locations_, lstate->locations_,
           sizeof(std::pair<KeyType, ValueType>) * size_);
    left_sib_ = lstate->left_sibling_;
    if (state->IsSeparated()) {
      right_sib_ = state->split_new_page_id;
    } else {
      right_sib_ = lstate->right_sibling_;
    }
  };

  ~LPage(){};

  bool InsertEntry(KeyType key, ValueType location, LPID self, LPID parent);

  bool DeleteEntry(KeyType key, ValueType location, LPID self);

  void Scan(const std::vector<Value> &values,
            const std::vector<oid_t> &key_column_ids,
            const std::vector<ExpressionType> &expr_types,
            const ScanDirectionType &scan_direction,
            std::vector<ValueType> &result, const KeyType *index_key);

  void ScanAllKeys(std::vector<ValueType> &result);

  void ScanKey(KeyType key, std::vector<ValueType> &result);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState(
      int max_index);

  bool SplitNodes(LPID self, LPID parent);

  void MergeNodes(LPID self, LPID right_sibling_lpid);

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  inline LPID GetLeftSiblingLPID() { return left_sib_; }

  inline LPID GetRightSiblingLPID() { return right_sib_; }

  std::pair<KeyType, ValueType> *GetLocationsArray() { return locations_; }

  void SetSize(int size) { size_ = size; }

  int GetSize() { return size_; }

 private:
  // return a vector of indices of the matched slots
  std::vector<oid_t> ScanKeyInternal(KeyType key);

  // return true if a given ItemPointer is invalid
  //  bool IsInvalidItemPointer(ValueType val);

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
