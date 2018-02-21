#ifndef CLSCOMPAREBASEFUNCTION_H
#define CLSCOMPAREBASEFUNCTION_H

#include "clscircrnadetection.h"
#include "clsparsestdcircrnainfo.h"

struct St_CandiSoftwareSupport
{
    St_Candidate stCandi;
    int iSoftwareCount;

    St_CandiSoftwareSupport()
    {
        iSoftwareCount = 0;
    }

    St_CandiSoftwareSupport(St_Candidate stV1, int iV2)
    {
        stCandi = stV1;
        iSoftwareCount = iV2;
    }
};

class ClsCompareBaseFunction
{
public:
    ClsCompareBaseFunction();
    ~ClsCompareBaseFunction();

private:
    int m_iReadLen;    
    En_CellType m_enCellType;

public:
    void Init(int iReadLen, string strTissue);

    // Brief.txt (same directory as executable file)
    // 所有的结果都将会转换成这个格式
    void ParseBrief(string strCurPath, vector<St_Candidate>& vCandi); // Brief is my current result
                                                                      // It is also the common format of comparison

    void ParseSimulationResult(string strSimulatePath, vector<St_Candidate>& vCandi);

    //**********parse the result of public databse, including: CircBase**********
    //Parse the result of other tools,including: CIRI, Find-circ, CIRCExplorer,
    //***************************************************************************
    /// Databse
    void ParseCircBaseResult(string strCircBasePath, vector<St_Candidate>& vCandi);
    /// CIRI
    void ParseCIRIResult(string strCIRIPath, vector<St_Candidate>& vCandi);
    /// Find-circ
    void ParseFindCircResult(string strFindCircPath, vector<St_Candidate>& vCandi);
    /// CIRCExplorer
    void ParseCircExplorerResult(string strCIRCexplorerPath, vector<St_Candidate>& vCandi);
    /// circRNA_finder
    void ParseCircRNAFinderResult(string strCircRNAFinderResult, vector<St_Candidate>& vCandi);

    /// Parse Brief circRNAdb
    void ParseCircRNADbBrief(string strFilePath, vector<St_Candidate>& vCandi);

    //*********************The basic comparision function************************
    //***************************************************************************
    int CheckBothHitExonBoundary(vector<St_Row_Chrom>& vChrom, St_Candidate* pCandi);

    bool CompareTwoCandi(St_Candidate* pCandi1, St_Candidate* pCandi2);

    bool CheckIfBoundaryFromSameExon(vector<St_Row_Chrom>& vChrom, St_Candidate* pCandi, bool& bShortExon);

    void GetCandiForSpeciChrom(vector<St_Row_Chrom>& vChrom, vector<St_Candidate>& vCandiAll,
                               vector<St_Candidate>& vCandiForSpeciChrom, string strChromName);

    void SaveHitCandi(vector<St_Candidate>& vCandiCIRI, string strCIIRIHitPath);
    void SaveCandi(vector<St_Candidate>& vCandi, string strCandiPath);

    string CreateConsensusBanchmark(string strMyProgram, string strCIRI, string strFind_Circ, string strCIRCExplorer);
    void CandiMerge(vector<St_Candidate>& vCandi, vector<St_CandiSoftwareSupport>& vSum);    

    ///set the circular rna type for each candidate (the result from other tools)
    void SetCircTag(vector<St_Row_Chrom>& vChrom, vector<St_Candidate>& vCandi);

};

#endif // CLSCOMPAREBASEFUNCTION_H
