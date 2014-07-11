/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __GrammarWrap_h
#define __GrammarWrap_h

#include "phy/Grammar.h"
#include "phy/Factors.h"
#include "phy/SeqIO.h"
#include "phy/EmitWrap.h"

////////////////////////////////////////////////////////////////
// The header defines convenience functions for dealing with the
// Grammar class. 
////////////////////////////////////////////////////////////////

namespace phy {

  namespace grammar {

    /** Calculate emissions and reset them in the grammar. */ 
    void calcAndResetEmissions(SeqData const & seqData, Grammar & g, EmitWrappers & ew, EmitAnnoMap & eam, string const & annoName);

  } // end namespace grammar

} // end namespace phy

#endif  // __GrammarWrap_h
