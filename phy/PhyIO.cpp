/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/PhyIO.h"

namespace phy {

  /** Copy information from Newick node to tree_t vertex and define edge from parent (p) to current vertex (v) */
  vertex_t copyNode1(NewickNode *np, TreeInfo & treeInfo, unsigned p)
  {
    treeInfo.nodeMap.push_back( np->getName() );
    unsigned ndIndex = treeInfo.nodeMap.size() - 1;
    if (ndIndex != 0) {  // not root
      treeInfo.branchMap.push_back( np->getBranchLength() );
      vector<unsigned> neighbors;
      neighbors.push_back(p);
      neighbors.push_back(ndIndex);
      treeInfo.branchNeighbors.push_back(neighbors);
    }
    return ndIndex;
  }


  void copyTreeRec1(NewickNode *np, TreeInfo & treeInfo, unsigned p)
  {
    p = copyNode1(np, treeInfo, p);
    unsigned n = np->getChildrenNumber();
    if (n) 
      for (unsigned i = 0; i < n; i++)
	copyTreeRec1(np->getChild(i), treeInfo, p);
  }


  TreeInfo newickToTreeInfo(string const &str)
  {
    NewickParser parser;
    NewickNode * root = parser.parseTree(str);
    TreeInfo treeInfo;
    copyTreeRec1(root, treeInfo, 0);
    treeInfo.root = 0;
    return treeInfo;
  }

  string readNewickStr(istream & str)
  {
  string s;
  getline(str, s, ';');
  if ( str.eof() )
    errorAbort("End of file reached when reading newick string from file.");
  s += ";";
  return s;
  }

  string readNewickStr(string const & fileName)
  {
    ifstream f;
    openInFile(f, fileName);
    string treeStr = readNewickStr(f);
    f.close();
    return treeStr;
  }


  tree_t mkBoostGraph(TreeInfo const & treeInfo)
  {
    tree_t tree;
    // add nodes
    for (unsigned i = 0; i < treeInfo.nodeMap.size(); i++) {
      vertex_t v = add_vertex(tree);
      tree[v].name = treeInfo.nodeMap[i];
    }
    // add branches
    for (unsigned i = 0; i < treeInfo.branchMap.size(); i++) {
      vector<unsigned> const & nei = treeInfo.branchNeighbors[i];
      edge_t e = add_edge(nei[0], nei[1], tree).first;
      tree[e].length = treeInfo.branchMap[i];
    }
    return tree;
  }


  /** Dump tree in dot format to string. (Can be visualized using the dot tool. */
  string writeDot(TreeInfo const & treeInfo)
  {
    tree_t tree = mkBoostGraph(treeInfo);
    return writeDot(tree);
  }


  /** Dump graph in dot format to fileName. */
  void writeDot(TreeInfo const & treeInfo, string const &fileName)
  {
    tree_t tree = mkBoostGraph(treeInfo);
    writeDot(tree, fileName);
  }


  /** Write phylo tree in Newick format. */
  string writeNewick(TreeInfo const & treeInfo)
  {
    tree_t tree = mkBoostGraph(treeInfo);
    return writeNewick(tree);
  }


  // deprecated functions (although the above are defined in terms of these, they would not normally be used directly -- i.e, they are now implementation details)

  /** Copy information from Newick node to tree_t vertex and define edge from parent (p) to current vertex (v) */
  vertex_t copyNode(NewickNode *np, tree_t & tree, vertex_t p)
  {
    vertex_t v = add_vertex(tree);
    tree[v].name = np->getName();
    if (treeTraits_t::null_vertex() != p) {
      edge_t e = add_edge(p, v, tree).first;
      tree[e].length = np->getBranchLength();
    }
    return v;
  }

  void copyTreeRec(NewickNode *np, tree_t & tree, vertex_t p)
  {
    p = copyNode(np, tree, p);
    int n = np->getChildrenNumber();
    if (n) 
      for (int i = 0; i < n; i++)
	copyTreeRec(np->getChild(i), tree, p);
  }

  /** Functor that assigns ids to edges of tree  */
  struct PreVisitorEdgeIds {
  
    PreVisitorEdgeIds(tree_t & tree, int initEdgeId=0) : tree_(tree), id_(initEdgeId) {};
  
    template <class Vertex>
    void operator()(const Vertex &v) 
    {
      if (in_degree(v, tree_)) 
	{
	  if (in_degree(v, tree_) > 1)
	    errorAbort("Graph is not a tree");
	  treeTraits_t::in_edge_iterator ii, ii_end;
	  tie(ii, ii_end) = in_edges(v, tree_);
	  tree_[*ii].id = id_;    
	  id_++;
	}
    }
    tree_t & tree_;
    /** edge id */
    int id_;
  };

  
  tree_t readNewick(string const &str)
  {
    NewickParser parser;
    NewickNode * root = parser.parseTree(str);
    tree_t tree = tree_t();
    copyTreeRec(root, tree, treeTraits_t::null_vertex() );

    // assign edge IDs
    PreVisitorEdgeIds preVisit(tree, 0);
    treeTraversal(tree, 0, preVisit);

    return tree;
  }


  void writeNewickRec(tree_t const & tree, vertex_t const &v, stringstream & ss)
  {
    // recursively list children
    treeTraits_t::adjacency_iterator ai, ai_end;
    tie(ai, ai_end) = adjacent_vertices(v, tree);
    if (ai != ai_end) {  // has children
      ss << "(";    
      for (; ai != ai_end; ++ai) {
	writeNewickRec(tree, *ai, ss);      
	if (ai + 1 != ai_end)
	  ss << ",";
      }
      ss << ")";    
    }    

    // write node name and (parent)branchlength
    if (tree[v].name.size() != 0)
      ss << "'" << tree[v].name << "'";

    if (in_degree(v, tree)) {
      if (in_degree(v, tree) > 1)
	errorAbort("Graph is not a tree");
      treeTraits_t::in_edge_iterator ii, ii_end;
      tie(ii, ii_end) = in_edges(v, tree);
      ss << ':' << fixed << tree[*ii].length;
    }
  }

  string writeNewick(tree_t const &tree)
  {
    stringstream ss;
    writeNewickRec(tree, 0, ss);
    ss << ";" << endl;
    return ss.str();
  }


  /** functor for use with write_graphviz. Writes out internal properties (name) */
  struct phyNodeWriter 
  {
    phyNodeWriter (tree_t const & tree) : tree_ (tree) {};

    template <class Vertex>
    void operator()(std::ostream& out, const Vertex &v) 
    {
      out << " [label=\"" << tree_[v].name << "\"]" << std::endl;
    }

    tree_t const & tree_;
  };

  /** functor for use with write_graphviz. Writes out internal edge (name) */
  struct phyEdgeWriter 
  {
    phyEdgeWriter (tree_t const & tree) : tree_ (tree) {};

    template <class Edge>
    void operator()(std::ostream& out, const Edge &e) 
    {
      out << " [label=\"" << tree_[e].length << "\"]" << std::endl;
    }

    tree_t const & tree_;
  };

  string writeDot(tree_t const &tree)
  {
    stringstream ss;
    write_graphviz(ss, tree, phyNodeWriter(tree), phyEdgeWriter(tree));
    return ss.str();
  }

  void writeDot(tree_t const &tree, string const &fileName)
  {
    ofstream f(fileName.c_str(), ios::out);
    if (!f)
      errorAbort("Cannot open file: " + fileName + "\n");
  
    f << writeDot(tree);
    f.close();
  }

} // namespace phy
