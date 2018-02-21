#ifndef CLSREADSMAPPING_H
#define CLSREADSMAPPING_H
#include "clsparsestdcircrnainfo.h"
#include "../../ShareLibrary/clsgtfparse.h"

enum En_ReadsType{rtPairEnd, rtSingle, rtMax};

class ClsReadsMapping
{
public:
    ClsReadsMapping();

public:
    void MapReadsToCircRNADb(string strRefPath, string strReads1Path, string strReads2Path, En_ReadsType enReadsType);
};

#endif // CLSREADSMAPPING_H
