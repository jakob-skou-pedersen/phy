/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/GrammarWrap.h"

namespace phy {

  namespace grammar {

    void calcAndResetEmissions(SeqData const & seqData, Grammar & g, EmitWrappers & ew, EmitAnnoMap & eam, string const & annoName)
    {
      // mk emissions
      unsigned n = seqData.seqSize();
      xvector_t tmpV(n);
      reset(tmpV);
      vector<xvector_t> vv(ew.vecWrapNames().size(), tmpV);
      ew.vectorEmissions(vv, seqData);
      
      xmatrix_t tmpM(n, n);
      reset(tmpM);
      vector<xmatrix_t> vm(ew.matWrapNames().size(), tmpM);
      ew.matrixEmissions(vm, seqData);
      
      // mask according to anno
      if ( eam.size() and annoName.size() ) {
	string const & anno = seqData.getAnno(annoName);
	maskEmissions(vv, ew.vecWrapNames(), anno, eam);
	maskEmissions(vm, ew.matWrapNames(), anno, eam);
      }
      
      g.resetEmissions(vv, vm);
    }

  } // end namespace grammar

} // end namespace phy
