/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "phy/PhyIO.h"
#include <algorithm>
#include <functional>

using boost::unit_test::test_suite;

string const IN_TREE_STRING = "('Bovine':0.693950,('Hylobates':0.360790,('Pongo':0.336360,('G._Gorilla':0.171470,('P._paniscus':0.192680,'H._sapiens':0.119270):0.083860):0.061240):0.150570):0.549390,'Rodent':1.214600);\n";

BOOST_AUTO_TEST_CASE(PhyIO_TreeInfo_readNewick_1) 
{
  phy::TreeInfo treeInfo = phy::newickToTreeInfo(IN_TREE_STRING);
  string treeStr = writeNewick(treeInfo);
  BOOST_CHECK(treeStr == IN_TREE_STRING);

  //  cout << treeStr << endl;
  //  writeDot(treeInfo, "out.dot");
}

// test of deprecated tree structure follows below.

BOOST_AUTO_TEST_CASE(PhyIO_boosGraph_readNewick_1) 
{
  phy::tree_t tree = phy::readNewick(IN_TREE_STRING);

  // test ids of root edges
  phy::treeTraits_t::out_edge_iterator ii, iiEnd;
  int const rootVertex = 0;
  tie(ii, iiEnd) = out_edges(rootVertex, tree);
  for (; ii < iiEnd; ++ii)
    {
    BOOST_CHECK( tree[*ii].id >= 0 );
    //    BOOST_TEST_MESSAGE( tree[*ii].id );
    }
}

// compare IDs of edges of graph
struct EdgeIdComp : public std::binary_function<phy::edge_t, phy::edge_t, bool> {
  EdgeIdComp(phy::tree_t const & tree) : tree_(tree) {};
  bool operator()( phy::edge_t const & a, phy::edge_t const & b) { return tree_[a].id < tree_[b].id; }
  phy::tree_t const & tree_;
};

// check that largest edge id equals number of edges -1 (zero based);
BOOST_AUTO_TEST_CASE(PhyIO_boosGraph_readNewick_2) 
{
  phy::tree_t tree = phy::readNewick(IN_TREE_STRING);
  phy::treeTraits_t::edge_iterator ei, eiBegin, eiEnd;
  tie(eiBegin, eiEnd) = edges(tree);
  EdgeIdComp cmp(tree);
  int maxId = tree[ *std::max_element(eiBegin, eiEnd, cmp) ].id;
  //  BOOST_TEST_MESSAGE( maxId );
  //  BOOST_TEST_MESSAGE( eiEnd - eiBegin - 1 );
  BOOST_CHECK( maxId == (eiEnd - eiBegin - 1) );
}

BOOST_AUTO_TEST_CASE(PhyIO_boosGraph_readWrite_1) 
{
  phy::tree_t tree = phy::readNewick(IN_TREE_STRING);
  string outTree = writeNewick(tree);
  BOOST_CHECK( IN_TREE_STRING == outTree );
}
