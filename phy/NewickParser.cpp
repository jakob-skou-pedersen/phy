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


# ifndef NO_FLEX
// undefined functions from NewickParserBase
NewickNode *NewickParserBase::parseTree(string const &tree)
{
  istringstream treeStream(tree);
  lexer.switch_streams(&treeStream, &cerr);
  yyparse();	  
  setParentPointer(theTree, 0);

  return theTree;
}				  

# else
NewickNode *NewickParserBase::parseTree(string const &tree)
{
  cerr << "NewickParserBase::parseTree cannot be called when compiled with NO_FLEX preprocessor flag. Install needed libraries (Flex and Bison) and recompile witout the NO_FLEX flag. Program terminating." << endl;
  exit(0);
  
  NewickNode * np = NULL;
  return np;
}				  
# endif /* NO_FLEX */

void NewickParserBase::setParentPointer(NewickNode *np, NewickNode *pp)
{
  np->setParent(pp);
  for(int i = 0; i < np->getChildrenNumber(); i++) {
    setParentPointer(np->getChild(i), np);
  }
}
