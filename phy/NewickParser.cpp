/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "NewickParser.h"

NewickNode *NewickParser::makeNode(string nm, 
				   double bl, 
				   NewickNode *p, 
				   vector<NewickNode *> *dl)
{
  NewickNode *np;
  np = new NewickNode(p, nm, bl);
  if (dl)
    for (unsigned int i = 0; i < dl->size(); i++)
      np->addChild((*dl)[i]);
  return np;
}


// undefined functions from NewickParserBase
NewickNode *NewickParserBase::parseTree(string const &tree)
{
  istringstream treeStream(tree);
  lexer.switch_streams(&treeStream, &cerr);
  yyparse();	  
  setParentPointer(theTree, 0);

  return theTree;
}				  


void NewickParserBase::setParentPointer(NewickNode *np, NewickNode *pp)
{
  np->setParent(pp);
  for(int i = 0; i < np->getChildrenNumber(); i++) {
    setParentPointer(np->getChild(i), np);
  }
}
