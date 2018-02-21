#ifndef CLSPARSESTDCIRCRNAINFO_H
#define CLSPARSESTDCIRCRNAINFO_H
#include "../../ShareLibrary/clsgtfparse.h"
#include "clscircrnadetection.h"
//#include <string>
//#include <vector>
//using namespace std;

struct St_StdCircRNAInfo
{
    string strName;
    string strChrom; //unsigned char ucChromIndex;
    int iStartPos;
    int iEndPos;
    bool bRC;
    En_CircType enType; // Self or regular
    vector<string> vOrgCell;

    //for record the total value of cell types
    //string strCellTypes;

    St_StdCircRNAInfo():strName(""), strChrom(""), iStartPos(-1), iEndPos(-1),
                        bRC(false), enType(ctMax)
    {}
};

enum En_CellType{ctH9, ctGlioblastoma, ctBrain, ctMaxx};

class ClsParseStdCircRNAInfo
{
public:
    ClsParseStdCircRNAInfo();
    ~ClsParseStdCircRNAInfo();

public:
    void Init(string strTissue);
    //For circRNADb:
    ///1: Parse circRNA record
    void ParseCircRNADb(string strFilePath, vector<St_StdCircRNAInfo>& vInfo);

    ///2: Parse circRNADb style gtf file
    void ParseGTFWithCircRNADbStyle(string strGTFPath, vector<St_Row_Chrom>& vChrom);

    ///3: Check the irregular circular RNA (do not located in the edge of the real exon)
    void CheckIrregularCircRNA(vector<St_StdCircRNAInfo>& vInfo,
                               vector<St_Row_Chrom>& vChrom); // 将所有不规则的CircRNA都输出


    ///3: Create Reference File
    string GenerateCircularRelatedRef(vector<St_StdCircRNAInfo>& vInfo,
                                      vector<St_Row_Chrom>& vChrom,
                                      string strRefPath);
    bool CreatCircRef(vector<St_Raw_Exon>& vExon, string& strRef,
                      bool bRC, string strCircRNAName, En_CircType& enType,
                      ofstream& ofs, ofstream& ofsExon, int& iSelfCircCase, int& iRegCircCase);

    void FilterByCell(vector<St_StdCircRNAInfo>& vInfo, En_CellType enCellType);

    //Get The most wanted valid standard CircRNA
    void GetValidStdCircRNAByChrom(vector<St_StdCircRNAInfo>& vTargetValidStdInfo,
                                   vector<St_Row_Chrom>& vChrom,
                                   string strDbPath, int iChromIndex);
    void PrintStdHitInfoByBriefCandiStyle(vector<St_StdCircRNAInfo>& vTargetValidStdInfo,
                                          vector<St_Row_Chrom>& vChrom, string strDbPath, int iChromIndex);

private:
    En_CellType m_enCellType;
};



#endif // CLSPARSESTDCIRCRNAINFO_H
