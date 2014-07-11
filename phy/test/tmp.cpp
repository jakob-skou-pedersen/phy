/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/Grammar.h"

using namespace phy;

int main(void)
{
  Grammar g = readGrammar("data/grammarOneState.txt");

  // emit models
  xvector_t missingData(4, 1);
  vector<xvector_t> sngEmit(1, missingData);

  g.resetEmissions(sngEmit);

  cout << g.inside() << endl;
}
