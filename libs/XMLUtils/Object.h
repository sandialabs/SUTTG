/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

//                               -*- C++ -*-

#ifndef _base_class_Object_
#define _base_class_Object_

#include "RefPtr.h"

/**
 * \class Object
 * Base class for reference counted objects.
 *
 * Part of the Particle Simulation Toolkit (pst)
 *
 * The "friend" classes "RefPtr" and "ConstRefPtr" take care of the
 * reference counting and garbage collection.  This means that it 
 * should be safe to create an array of reference counted objects,
 * as long as you *do not* assign a reference counted pointer to 
 * any at the entries in the array at any time.
 */
class Object {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Construct a new reference counted object with a zero reference count
  Object() : refs_(0) {
  }
  
  /// Destroy this object
  virtual ~Object() {
  }
  
  /// Returns the number of references that are held to this object
  long int reference_count() const {
    return refs_;
  }
  
protected:
  /// Enables the friends of the class to increment and decrement the
  /// reference count.
  long int reference_grab() const {
    return ++refs_;
  }
  
  long int reference_release() const {
    return --refs_;
  }    
  
private:
  mutable long int refs_;
};

#endif //_utility_ref_Object_
