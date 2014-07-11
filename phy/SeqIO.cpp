/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "phy/SeqIO.h"

namespace phy {

  void SeqData::toUpper()
  {
    for (unsigned i = 0; i < sequences.size(); i++) {
      string & s = sequences[i];
      transform(s.begin(), s.end(), s.begin(), (int(*)(int)) toupper);
    }
  }
 

  SeqData amaToSeqData(ama const & a, bool toUpper, string const & idFeatureKey)
  {
    SeqData sd;

    // id
    sd.entryId = a.getFeature(idFeatureKey);

    // sequences
    for (unsigned i = 0; i < a.sequences.size(); i++)
      sd.addSeq(a.sequences[i].src, a.sequences[i].text);
    
//    // lSequences are taken from ama lAnno
//    for (unsigned i = 0; i < a.lAnno.size(); i++)
//      sd.addLSeq(a.lAnno[i].first, a.lAnno[i].second);

    // annotations (only ama anno are included, not ama lAnno)
    for (unsigned i = 0; i < a.anno.size(); i++)
      sd.addAnno(a.anno[i].first, a.anno[i].second);

    // to upper case
    if (toUpper)
      sd.toUpper();

    return sd;
  }


  vector<SeqData> amaToSeqData(vector<ama> const & v, bool toUpper, string const & idFeatureKey)
  {
    vector<SeqData> u;
    for (unsigned i = 0; i < v.size(); i++) 
      u.push_back(amaToSeqData(v[i], toUpper, idFeatureKey) );
    return u;
  }



}
