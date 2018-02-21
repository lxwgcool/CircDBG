#ifndef CLSPREEVALUATION_H
#define CLSPREEVALUATION_H

#include "clsparsestdcircrnainfo.h"
#include "clsreadsmapping.h"
#include "clsconfig.h"

using namespace std;

class ClsPreEvaluation
{
public:
    ClsPreEvaluation(St_Config& stConfig);
    ~ClsPreEvaluation();

private:
    St_Config m_stConfig;
    En_CellType m_enCellType;

public:
    void GeneralEvaluate();
    void GetValidStdCandiByChormoson(int iChromIndex);
};

#endif // CLSPREEVALUATION_H
