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

#ifndef _util_xml_class_XMLParser_
#define _util_xml_class_XMLParser_

#include "XMLElement.h"
#include <iostream>

/**
 * \class XMLParser
 * A pure abstract base class for parsers that read data from
 * an XML file and return the top node of a parse tree.
 * The parse tree node is a RefPtr< XMLElement >.
 */
class XMLParser : virtual public Object {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Default constructor.  Intended for derived classes.
  XMLParser();
  
  /// Destructor.
  virtual ~XMLParser();
  
  /// Parse the given input buffer and return the parse tree.
  virtual RefPtr<XMLElement> parse(std::istream&) = 0;

  /// Parse the given input string and return the parse tree.
  virtual RefPtr<XMLElement> parse(const std::string&) = 0;
};

#endif // _util_xml_class_XMLParser_
