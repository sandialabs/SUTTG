/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

//                                 -*- C++ -*-

#ifndef _util_xml_class_XMLAttributeList_
#define _util_xml_class_XMLAttributeList_

#include "Object.h"
#include <map>

/**
 * \class XMLAttributeList
 *
 * The implementation of a container for XML attributes.
 *
 * Attributes are stored as key-value pairs, both of which are stored
 * as text strings.  The keys are stored and handled in a case-sensitive
 * manner, whereas the values are stored exactly as given.
 *
 * The order in which attributes were added is not preserved.
 * I am not sure whether this is a problem or not.
 *
 * Boolean values are a special case.  They are treated in a
 * case-insensitive manner.  True boolean values are returned for
 * the strings 'true', 'yes', and non-zero numerical values.  False boolean
 * values are returned for the strings 'false', 'no', and zero values.
 * All other strings are considered unacceptable boolean values.
 * 
 * Written by Helgi Adalsteinsson
 * Modified by Bert Debusschere, 3/13/08 to make keys case-sensitive
 */
class XMLAttributeList : public Object {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// The container type used to hold the attributes.
  typedef std::map< std::string, std::string > Map_t;
  
  /// The iterator type returned by this implementation.
  typedef Map_t::iterator iterator;
  typedef Map_t::const_iterator const_iterator;
  
  /// Construct a blank attribute list.
  XMLAttributeList();
  
private:
  /// Blocked copy constructor.  Throws an exception.
  XMLAttributeList(const XMLAttributeList&);

public:  
  /// Destroy this list.
  virtual ~XMLAttributeList();
  
private:
  /// Blocked assignment operator.  
  /// \throw MyException
  XMLAttributeList& operator=(const XMLAttributeList&);

public:
  /// Get the number of attributes in the list.
  int size() const;
  
  /// Return true if the given key is defined.
  bool has(const std::string&) const;
  
  /// Get the attribute associated with the given key.
  /// \throw MyException if the key is not defined.
  const std::string& get(const std::string&) const;
  
  /// Get the attribute associated with the given key
  /// _or_ return the given default value if the key is not defined.
  std::string get(const std::string&, const std::string&) const;
  
  /// Get the given attribute as an integer value.
  /// \throw MyException if the key is not defined.
  /// \throw MyException if the value is not a valid integer.
  int get_int(const std::string&) const;
  
  /// Get the given attribute as an integer value
  /// _or_ return the given default value if the key is not defined.
  /// \throw MyException if the value is set and is not a valid integer.
  int get_int(const std::string&, int) const;
  
  /// Get the given attribute as a real value.
  /// \throw MyException if the key is not defined.
  /// \throw MyException if the value is not a valid number.
  double get_double(const std::string&) const;
  
  /// Get the given attribute as a real value
  /// _or_ return the given default value if the key is not defined.
  /// \throw MyException if the value is set and is not a valid number.
  double get_double(const std::string&, double) const;
  
  /// Get the given attribute as a boolean value.
  /// \throw MyException if the key is not defined.
  /// \throw MyException if the value is not a valid boolean value.
  /// True boolean values are "yes", "true", and 'non-zero' numerical values.
  /// False boolean values are "no" "false", and 'zero' numerical values.
  bool get_bool(const std::string&) const;
  
  /// Get the given attribute as a boolean value
  /// _or_ return the given default value if the key is not defined.
  /// \throw MyException if the value is not a valid boolean value.
  bool get_bool(const std::string&, bool) const;
  
  /// Assign a text attribute to the given key.
  void set(const std::string&, const std::string&);
  
  /// Assign an integer value to the given key.
  void set_int(const std::string&, int);
  
  /// Assign a numerical value to the given key.
  void set_double(const std::string&, double);
  
  /// Assign a boolean attribute to the given key.
  /// True boolean values are added as "true", false as "false"
  void set_bool(const std::string&, bool);
  
  /// Get an iterator to the first element
  iterator begin();
  
  /// Get an iterator past the last element.
  iterator end();
  
  /// Get an iterator to the first element in a const context.
  const_iterator begin() const;
  
  /// Get an iterator past the last element in a const context.
  const_iterator end() const;
  
private:
  /// The attributes.
  Map_t attribute_;
  
  /// Convert a string to lower case (conversion in place).
  /// (needed to handle boolean values)
  void make_lower_case(std::string&) const;
  
  /// Get an iterator pointing to the location of the given string.
  iterator get_location(const std::string&);
  
  /// Get an iterator to a location in a const context.
  const_iterator get_location(const std::string&) const;
  
  /// Return the boolean value of the given string.
  /// \throw MyException if the string is not a valid boolean value.
  bool boolean_value(const std::string&, const char* where) const;
};

#endif // _util_xml_class_XMLAttributeList_
