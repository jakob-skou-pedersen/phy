/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __EmitWrapIO_h
#define __EmitWrapIO_h

#include "phy/EmitWrapPhylo.h"
#include "phy/EmitWrapDfg.h"
#include <iostream>

namespace phy {

  /** Overloaded input operator for EmitWrapper. */
  void readEmitWrappers(istream & str, EmitWrappers & ew, string const & filePrefix = "");

  /** Read EmitWrapper from file. */
  void readEmitWrappers(string const & file, EmitWrappers & ew);
  EmitWrappers readEmitWrappers(string const & file);

  /** Version of readEmitWrappers that allow a global tree (a parameter to phloModels) to be specified */
  void readEmitWrappers(string const & file, EmitWrappers & ew, string const & treeString);
  EmitWrappers readEmitWrappers(string const & file, string const & treeString);


} // end namespace phy

#endif  // __EmitWrapIO_hh
