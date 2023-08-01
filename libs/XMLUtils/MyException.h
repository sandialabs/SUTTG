/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

//                          -*- C++ -*-

#ifndef _MyException_
#define _MyException_

#include <iostream>
#include <exception>
#include <string.h>

/**
 * Just an example exception - feel free to override this.
 */
class MyException : public std::exception {
public:
  /// Construct an exception using a C-style character string.
  MyException(const char* errormessage) {
    std::cerr << "ERROR:  " << errormessage << "\n";
    error_ = std::string("MyException: ") + errormessage;
  }

  /// Construct an exception using a C++-style string
  MyException(const std::string& errormessage) {
    std::cerr << "ERROR:  " << errormessage << "\n";
    error_ = std::string("MyException: ") + errormessage;
  }

  /// Destroy.
  virtual ~MyException() throw() {
  }

  /// What's going on?
  const char* what() const throw() {
    try {
      return error_.c_str();
    } catch(...) {
      ;/// This function is not permitted to throw exceptions.
    }
    return error_.c_str();
  }
  
private:
  std::string error_;
};

#endif // _MyException_
