/*-------------------------------------------------------------------------
*
* backend.h
* file description
*
* Copyright(c) 2015, CMU
*
* /n-store/src/storage/backend.h
*
*-------------------------------------------------------------------------
*/

#pragma once

#include <string>

#include "common/types.h"

namespace nstore {
namespace storage {

//===--------------------------------------------------------------------===//
// Backend (for physical storage)
//===--------------------------------------------------------------------===//

/// Represents a storage backend. Can reside on MM or NVM.
class Backend {

public:
	virtual ~Backend(){};

	//===--------------------------------------------------------------------===//
	// Interface
	//===--------------------------------------------------------------------===//

	virtual void* Allocate(size_t size) = 0;

	virtual void Free(void* ptr) = 0;

	virtual void Sync(void* ptr) = 0;

	virtual std::string GetBackendType() const = 0;
};

} // End storage namespace
} // End nstore namespace