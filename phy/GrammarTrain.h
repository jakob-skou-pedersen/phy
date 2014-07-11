/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __GrammarTrain_h
#define __GrammarTrain_h

#include "phy/Grammar.h"
#include "phy/GrammarWrap.h"
#include "phy/GrammarIO.h"
#include "phy/Factors.h"
#include "phy/SeqIO.h"
#include "phy/EmitWrap.h"

namespace phy {

  namespace grammar {

    /** Expectation maximization training of transition probs of grammar. */
    void grammarEm(vector<SeqData> const & seqDataVec, Grammar & g, EmitWrappers & ew, EmitAnnoMap & eam, string const & annoName, number_t totalPseudoCounts, number_t minDeltaLogLik, unsigned maxIter, string const & logFile, string const & tmpOutputFile);

    /** Creates a compositeFactorSet with a factor for each
	transitionDistribution of the grammar. pseudoCounts defines
	the total number of pseudoCounts for each transition
	distribution. */
    CompositeFactorSet mkTransFactorSet(Grammar & g, number_t totalPseudoCounts);


  } // end namespace grammar

} // end namespace phy

#endif  // __GrammarTrain_h
