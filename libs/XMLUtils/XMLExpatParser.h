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

#ifndef _util_xml_class_XMLExpatParser_
#define _util_xml_class_XMLExpatParser_

#include "XMLParser.h"
#include <iostream>
#include <vector>
#include <expat.h>

/**
 * An XML parser that uses the Expat library to handle the gruntwork.
 * This class requires that the Expat be installed on your system.
 *
 * Expat is available <a href="http://www.jclark.com/xml/expat.html">
 * at the Expat site</a>
 *
 * This class may not be fully exception safe, since there is no
 * good way of enforcing that the Expat parser is destroyed cleanly.
 */
class XMLExpatParser : public XMLParser {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Construct a new parser.
  XMLExpatParser();
  
private:
  /// Blocked copy constructor.  It is not safe to copy this object
  /// since the Expat parser may save state which cannot be duplicated.
  XMLExpatParser(const XMLExpatParser&);
  
public:
  /// Destructor.
  virtual ~XMLExpatParser() throw();
  
private:
  /// Blocked assignment operator.  Not for public consumption.
  XMLExpatParser& operator=(const XMLExpatParser&);
  
public:
  /// Parse the given input buffer and return a parse tree.
  RefPtr<XMLElement> parse(std::istream&);

public:
  /// Parse the given input string and return a parse tree.
  RefPtr<XMLElement> parse(const std::string&);
  
private:
  /// The Expat parser.
  XML_Parser parser_;
  
  /// The path that we have traversed so far in building the tree.
  std::vector< RefPtr<XMLElement> > path_;
  
  /// The current leaf of the parse tree.
  RefPtr<XMLElement> leaf_;
  
  /// The method used to parse the start tag.
  void do_start(const XML_Char*, const XML_Char**);
  
  /// The method used to parse the end tag.
  void do_end(const XML_Char*);
  
  /// The method used to parse character (content) data.
  void do_character_data(const XML_Char*, int);
  
  /// Initialize the state of the parser.
  void init();
  
public:
  /// Static wrapper method used as a callback to get the 'start' tag.
  /// This method is for internal use only.  Calling this method
  /// directly will most likely result in a segmentation fault.
  static void start_(void*, const XML_Char*, const XML_Char**);
  
  /// Static wrapper method used as a callback to get the 'end' tag.
  /// This method is for internal use only.  Calling this method
  /// directly will most likely result in a segmentation fault.
  static void end_(void*, const XML_Char*);
  
  /// Static wrapper method used as a callback to get character data.
  /// This method is for internal use only.  Calling this method
  /// directly will most likely result in a segmentation fault.
  static void character_data_(void*, const XML_Char*, int);
};

#endif // _util_xml_class_XMLExpatParser_
