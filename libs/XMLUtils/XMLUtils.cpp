/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include "XMLUtils.h"

//
// Dump out XML tree to a file or output stream
//
void dump_xml_tree(RefPtr<XMLElement> tree, const std::string& indentation,std::ostream& outfile) {
  // Print the 'header' for this node.
  outfile << indentation << "<" << tree->label();

  // Dump all attributes associated with this node.
  RefPtr<XMLAttributeList> att = tree->attributes();
  XMLAttributeList::iterator it, end = att->end();
  for(it = att->begin(); it != end; ++it)
    outfile << " " << it->first << "=\"" << it->second << "\"";

  if(tree->count_children() == 0) {
    // If this node has no children, we're done.
    outfile << " />\n";
  }
  else {
    // ... otherwise, we recursively dump all the child nodes.
    outfile << ">\n";
    int children = tree->count_children();
    std::string new_indent = indentation + "  ";
    for(int kid = 0; kid < children; ++kid)
      dump_xml_tree(tree->get_child(kid), new_indent, outfile);
    // .. and print the closing bracket for this element.
    outfile << indentation << "</" << tree->label() << ">\n";
  }
}

