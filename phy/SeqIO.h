/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __SeqIO_h
#define __SeqIO_h

#include "phy/utils.h"
#include "phy/ama.h"
#include <fstream>


namespace phy {

  using namespace std;

  // struct for sequential data
  struct SeqData {
    
    string entryId;             // a unique id 
    vector<string> sequences;   // sequences
    vector<string> seqNames;    // names of sequences (same order as sequences)
    vector<vector<string> > lSequences;   // sequences of long, multi-char symbols
    vector<string> lSeqNames;    // names of sequences (same order as sequences)
    vector<string> annotations; // annotation sequences
    vector<string> annoNames;   // names of annotations (same order as annotation sequences)

    /** Convert all sequences to upper case */
    void toUpper();

    /** Return length of sequences[0] if it exists and zero otherwise. */
    unsigned seqSize() const {return (sequences.size() > 0) ? sequences[0].size() : 0;}

    // convenience getters and setters
    inline void addSeq(string const & name, string const & seq) {seqNames.push_back(name); sequences.push_back(seq);}
    //    inline void addLSeq(string const & name, vector<string> const & lSeq) {lSeqNames.push_back(name); lSequences.push_back(lSeq);}
    inline void addAnno(string const & name, string const & anno) {annoNames.push_back(name); annotations.push_back(anno);}

    inline bool hasSeq(string const & name) const {return hasElement(seqNames, name);}
    //    inline bool hasLSeq(string const & name) const {return hasElement(lSeqNames, name);}
    inline bool hasAnno(string const & name) const {return hasElement(annoNames, name);}

    inline string const & getSeq(string const & name) const {return sequences[ getIndex(seqNames, name) ];}
    //    inline vector<string> const & getLSeq(string const & name) const {return lSequences[ getIndex(lSeqNames, name) ];}
    inline string const & getAnno(string const & name) const {return annotations[ getIndex(annoNames, name) ];}
  };

  /** copy id, sequences, and 'anno' (not 'lAnno') annotation from ama
      to SeqData. If toUpper is true, then sequences will be converted to
      upper case. */
  SeqData amaToSeqData(ama const & a, bool toUpper = true, string const & idFeatureKey = "ENTRY");
  vector<SeqData> amaToSeqData(vector<ama> const & v, bool toUpper = true, string const & idFeatureKey = "ENTRY");


} // namespace phy

#endif //  __SeqIO_h
