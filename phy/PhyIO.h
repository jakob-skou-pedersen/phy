/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __PhyIO_h
#define __PhyIO_h

#include "phy/NewickParser.h"
#include "phy/Tree.h"
#include "phy/utils.h"
#include <boost/graph/graphviz.hpp>
#include <fstream>

namespace phy {

  /** Parse a newick tree string and return the corresponding TreeInfo struct.*/
  TreeInfo newickToTreeInfo(string const & str);

  /** Read and return a newick tree string from instream str (which is not closed).*/
  string readNewickStr(istream & inFile);

  /** Read newick tree string from fileName. */
  string readNewickStr(string const & fileName);

  /** Dump tree in dot format to string. (Can be visualized using the dot tool. */
  string writeDot(TreeInfo const & treeInfo);

  /** Dump graph in dot format to fileName. */
  void writeDot(TreeInfo const & treeInfo, string const &fileName);

  /** Write phylo tree in Newick format. */
  string writeNewick(TreeInfo const & treeInfo);

  /** Read treeInfo and return a boost graph (of type defined in Tree.h). */
  tree_t mkBoostGraph(TreeInfo const & treeInfo);


  // deprecated functions based on graph structure defined in Tree.h.

  /** Parse a newick tree string and return a reference to the corresponding graph. */
  tree_t readNewick(string const & str);

  /** Parse a tree_t graph and write it to a string in Newick format. */
  string writeNewick(tree_t const & tree);

  /** Dump graph in dot format to string. */
  string writeDot(tree_t const & tree);

  /** Dump graph in dot format to fileName. */
  void writeDot(tree_t const & tree, string const & fileName);

} //namespace phy

#endif  //__PhyIO_h
