#ifndef CLSBLAST_H
#define CLSBLAST_H

#include "clsbasealgorithm.h"
#include <vector>

struct St_SegUnit
{
    int iStart;
    int iEnd;
    St_SegUnit():iStart(-1), iEnd(-1)
    {}
    St_SegUnit(int iStartValue, int iEndValue)
    {
        iStart = iStartValue;
        iEnd = iEndValue;
    }
    void Init()
    {
        iStart = -1;
        iEnd = -1;
    }

};

struct St_BlastAlignUnit
{
    string strName;
    int iAlignNum;
    float fAlignPercent;
    int iStartPos; //This is for subject (means reference)
    int iEndPos;   //This is for subject (means reference)
    float fIdentifyRatio;
    float fGapRatio;
    string strAlignInfo;
    int iQueryStartPos;
    int iQueryEndPos;

    St_BlastAlignUnit()
    {
        Init();
    }

    void Init()
    {
        strName = "";
        fAlignPercent = -1;
        iStartPos = -1;
        iEndPos = -1;
        iAlignNum = -1;
        fIdentifyRatio = 0;
        fGapRatio = 0;
        string strAlignInfo = "";
    }
    bool IsRC() // If it is reverse complementary
    {
        if(iStartPos < 0 || iEndPos < 0)
            return false;
        if(iStartPos > iEndPos)
            return true;
        else
            return false;
    }
};

struct St_BlastResult
{
    vector<St_BlastAlignUnit> vBlastAlign;

    void Init()
    {
        vBlastAlign.clear();
    }
};

enum En_FastaType{ftQuery=0, ftRef, ftMax};

class ClsBlast
{
public:
    ClsBlast();
   ~ClsBlast();

public:
    void Init(string strBlastRootFolder);

    void TwoSeqBlast(string& strQueryPath, string& strSubjectPath, bool bShowFlank);

    void GetTwoSeqAlignResult(string& strQueryPath, string& strSubjectPath,
                              vector<St_BlastResult>& vBlastResult,
                              bool bRecord, bool bSaveAlignDetail);

    //Create Fa File
    string CreatFaFile(string strName, string& strSeq, En_FastaType enFastaType);

private:
    //Make comparison by blastn
    string TwoFastaFileAlign(string strQuerySeqPath, string strTargetSeqPath,
                             int iSeqSizeType=0); //0: normal size type, 1: small size type
                                                //2: large size type

    St_BlastResult ParseResult(string strRstPath, vector<St_BlastResult>& vBlastResult,
                               bool bRecord, bool bSaveAlignDetail=false);

private:
    string m_strBlastRootFolder;
};

#endif // CLSBLAST_H
