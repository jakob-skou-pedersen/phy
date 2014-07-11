/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __PhyloTreeIO_h
#define __PhyloTreeIO_h

#include "phy/PhyloTree.h"

using namespace std;

namespace phy {

  // read / write of time
  // reversible phylogenetic 

  PhyloTree readPhyloTree(string const file, string const & treeString = "");
  PhyloTree readPhyloTree(istream & str, string const & treeString = "");

  void writePhyloTree(string const file, PhyloTree const & phyloTree);
  void writePhyloTree(ostream & str, PhyloTree const & phyloTree);

} // end namespace phy

#endif  //__PhyloTreeIO_h
