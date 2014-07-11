
(1) SeqCppFiles.tar : contains the modified SeqIO.cpp and SeqIO.h. These files are modified to read probe data from ama files. Please see the .ama files included

(2) dfgSpec.tar     : contains 3 emit models in the directory ./dfgSpec/...., for example dfgSpec/pairEmitModel/, dfgSpec/singleEmitModel/  and dfgSpec/stackingEmitModel/ 
                      Each of these EmitModels contain 4 different files for dfg inputs.

(3) gramarFiles.tar : contains grammarStacking.txt, grammarStackingAnnoMap.txt and grammarStackingEmitModels.txt. In grammarStackingEmitModels.txt - the above emitModels are specified through the directory /dfgSpec/..EmitModels/

(4) SeqDataToDfg.txt: file for mapping input sequence.
