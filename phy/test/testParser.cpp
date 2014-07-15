/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include <string>
#include <iostream>
#include "phy/NewickParser.h"
#include "phy/PhyIO.h"

void printTreeNewickFormatRec(NewickNode *np, stringstream &sStream)
{
  if (np->getChildrenNumber() != 0)
    sStream << "(";    

  for (int i = 0; i < np->getChildrenNumber(); i++){
    printTreeNewickFormatRec(np->getChild(i), sStream);
    if (i + 1 < np->getChildrenNumber() )
      sStream << ",";
  }

  if (np->getChildrenNumber() != 0)
    sStream << ")";


  if (np->getName().size() != 0) {
    sStream << "'";
    sStream << np->getName();
    sStream << "'";
  }

  if (np->getBranchLength() != 0)
    sStream << ':' << fixed << np->getBranchLength();
}

string const printTreeNewickFormat(NewickNode * root)
{
  stringstream sStream;

  sStream << "(";
  for (int i = 0; i < root->getChildrenNumber(); i++){
    printTreeNewickFormatRec(root->getChild(i), sStream);
    if (i + 1 < root->getChildrenNumber() )
	sStream << ",";
  }
  sStream << ")";      

  if (root->getName().size() != 0) {
    sStream << "'";
    sStream << root->getName();
    sStream << "'";
  }

  if (root->getBranchLength() != 0)
    sStream << ':' << root->getBranchLength();

  sStream << ";";

  return sStream.str();
}

int main() {
  string treeString = "(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);";

  phy::tree_t tree = phy::readNewick(treeString);
  //cout << writeDot(tree);

  cout << writeNewick(tree);
  writeDot(tree, "out.dot");

}
