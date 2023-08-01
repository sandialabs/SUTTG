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

#ifndef _utility_ref_RefPtr_
#define _utility_ref_RefPtr_

#include "MyException.h"
#include <typeinfo>  // for dynamic_cast
#include <unistd.h>
#include <stddef.h>

/**
 * \class RefPtr
 * Reference counted pointer that gives the holder 
 * modification privileges to the pointee.
 *
 * Part of the Particle Simulation Toolkit (pst)
 */
template <class T>
class RefPtr {
public:
  /// Make the typename that this pointer holds accessible to other objects.
  typedef T Type;
  
  /// Construct a new RefPtr and initialize the pointee to NULL.
  RefPtr() : ptr_(NULL) {}
  
  /// Construct a new RefPtr and initialize the pointee as given.
  RefPtr(T* p) : ptr_(p) {
    grab();
  }
  
  /// Construct a new RefPtr and initialize to the given RefPtr pointee.
  RefPtr(const RefPtr<T>& p) : ptr_(p.ptr_) {
    grab();
  }
  
  /// Perform a static cast to initialize this pointee.
  /// This cast is only valid if T is a parent class of Other
  template <class Other>
  RefPtr(RefPtr<Other> p) : ptr_(static_cast<T*>(p.pointee())) {
    grab();
  }
  
  /// Destroy this RefPtr.
  ~RefPtr() {
    release();
  }
  
  /// Assign the value of this RefPtr to the given pointee.
  RefPtr<T>& operator=(T* p) {
    if(p != ptr_) {
      release();  // release the old pointer
      ptr_ = p;   // assign our value to this one
      grab();     // grab this pointer
    }
    return *this;
  }
  
  /// Assign the value of this RefPtr to the pointee of the given RefPtr.
  RefPtr<T>& operator=(const RefPtr<T>& p) {
    if(p.ptr_ != ptr_) {
      release();
      ptr_ = p.ptr_;
      grab();
    }
    return *this;
  }
  
  /// Use dynamic_cast to set the pointee to the
  /// pointer that was passed in, and return *this.
  /// The returned value is NULL if the cast fails.
  template <class Other>
  RefPtr<T>& cast(Other* p) {
    if(p.ptr_ != ptr_) {
      release();
      
      //std::cout << "DEBUG: Dynamic cast from type " << typeid(p).name()
      //	<< " to " << typeid(T).name() << std::endl; 
      
      ptr_ = dynamic_cast<T*>(p);
      if(p != NULL && ptr_ == NULL) {
	throw MyException
	  (std::string("RefPtr::cast(Other):  Failed dynamic cast from ")
	   + std::string(typeid(Other).name()) + std::string(" to ") +
	   std::string(typeid(Type).name()));
	
      }
      grab();
    }
    return *this;
  }
  
  /// Use dynamic_cast to set the pointee to the
  /// pointee of the RefPtr given, and return *this.
  /// The returned value is NULL if the cast fails.
  template <class Other>
  RefPtr<T>& cast(RefPtr<Other> p) {
    if(p.ptr_ != ptr_) {
      release();
      
      //std::cout << "DEBUG: Dynamic cast from type " 
      //	<< typeid(p.pointee()).name()
      //	<< " to " << typeid(T).name() << std::endl; 
      
      ptr_ = dynamic_cast<T*>(p.pointee());
      if(p != NULL && ptr_ == NULL) {
	throw MyException
	  (std::string("RefPtr::cast(Other):  Failed dynamic cast from ")
	   + std::string(typeid(Other).name()) + std::string(" to ") +
	   std::string(typeid(Type).name()));
      }
      grab();
    }
    return *this;
  }
  
  /// Return the pointee of this RefPtr.
  /// This will throw an exception if the pointee is NULL.
  T* operator->() const {
    if(ptr_ == NULL) {
      std::cerr << "RefPtr<" << typeid(T).name() 
		<< ">::operator->() const invoked on a null pointer\n";
      throw MyException("RefPtr::operator->() const");
    }
    return ptr_;
  }
  
  /// Return a reference to the pointee of this RefPtr.
  /// This will not work right if the pointee is NULL.
  T& operator*() const {
    if(ptr_ == NULL) {
      std::cerr << "RefPtr<" << typeid(T).name()
		<< ">::operator*() const invoked on a null pointer\n";
      throw MyException("RefPtr::operator*() const");
    } 
    return *ptr_;
  }
  
  /// Return the pointee of this RefPtr.
  T* pointee() {
    return ptr_;
  }
  
  /// Return the pointee of this RefPtr in a const context.
  const T* pointee() const {
    return ptr_;
  }
  
  /// Compare the pointee of this RefPtr with the given pointer.
  bool operator==(const T* p) const {
    return ptr_ == p;
  }
  
  /// Compare the value of this pointee with the pointee of the given RefPtr.
  bool operator==(const RefPtr<T>& p) const {
    return ptr_ == p.pointee();
  }
  
  /// Test inequality.
  bool operator!=(const T* p) const {
    return ptr_ != p;
  }
  
  /// Test inequality.
  bool operator!=(const RefPtr<T>& p) const {
    return ptr_ != p.ptr_;
  }
  
  /// Convenience routine to sort pointer values in standard containers.
  inline bool operator<(const RefPtr<T>& p) const {
    return ptr_ < p.ptr_;
  }
  
  /// Convenience routine to sort pointer values in standard containers.
  template <class Other>
  bool operator<(const RefPtr<Other>& p) const {
    return ptr_ < p.pointee();
  }
  
private:
  T* ptr_;
  
  /// Grab a reference to the current pointee if it is not NULL.
  inline void grab() {
    if(ptr_ != NULL) 
      ptr_->reference_grab();
  }
  
  /// Release the reference to the current pointee if it is not NULL.
  /// If this results in the reference count of the pointee dropping to zero,
  /// delete the object pointed to.
  inline void release() {
    if(ptr_ != NULL) {
      if(ptr_->reference_release() == 0)
	delete ptr_;
    }
  }
};

#endif //_utility_ref_RefPtr_
