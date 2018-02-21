#ifndef CLSCLASSIFICAITON_H
#define CLSCLASSIFICAITON_H

#include "clsbasealgorithm.h"
#include "clsgtfparse.h"
#include "clsblast.h"

enum En_CircType{ctSelf, ctRegular, ctMax};

enum En_MapStatus{msImbalance, msLowQuality,
                  msAdditionalPart, msAdditionalPart_Ref,
                  msMultiSupportR,
                  msMultiSupportS, msMultiSupportSRDiff,
                  msMultiSupportSRChimeric,
                  msGood, msBad,
                  msMax}; // 8 Cases

struct St_Candidate
{
    string strChromName;
    int iStartPos;
    int iEndPos;
    //-->
    int iSupportNum; // Especially For My result
    //<--
    bool bRC;
    En_CircType enType;
    string strTag;
    bool bHit;

    St_Raw_Exon stExonDonor;  //This is Donor
    St_Raw_Exon stExonAcceptor; // Only For Reg Case

    //---->For Blast
    vector<string> vReadsSeq;
    vector<En_MapStatus> vMapStatus;

    string strBlastRef;
    string strSelfDonorRef;
    string strSelfAcceptorRef;

    int iBreakPoint;
    int iSelfDonorBP;
    int iSelfAcceptorBP;
    En_MapStatus enMajorMapStatus;  //来自于vMapStatus
    //<-----

    St_Candidate()
    {
        Clear();
    }

    void Clear()
    {
        strChromName = "";
        iStartPos = 0;
        iEndPos = 0;
        iSupportNum = 0;
        bRC = false;
        enType = ctMax;
        strTag = "xx";
        bHit = false;
        vReadsSeq.clear();
        strBlastRef = "";
        strSelfDonorRef = "";
        strSelfAcceptorRef = "";
        iBreakPoint = -1;
        enMajorMapStatus = msMax;
        iSelfDonorBP = -1;
        iSelfAcceptorBP = -1;
    }

    St_Candidate(string strV1, int iV2, int iV3, int iV4=1, bool bV5=false,
                 En_CircType enV6=ctMax, string strV7="xx")
    {
        strChromName = strV1;
        iStartPos = iV2;
        iEndPos = iV3;
        iSupportNum = iV4;
        bRC = bV5;
        enType = enV6;
        strTag = strV7;
    }


    void SetChromStartEnd(string strV1, int iV2, int iV3, int iV4=1)
    {
        this->strChromName = strV1;
        this->iStartPos = iV2;
        this->iEndPos = iV3;
        this->iSupportNum = iV4;
    }

    void SetCircType(En_CircType enType)
    {
        this->enType = enType;
    }

    int GetLength()
    {
        return abs(iEndPos - iStartPos);
    }

    string GetTypeString()
    {
        if(enType == ctMax)
            return "M"; //Max or Mixture
        else if(enType == ctSelf)
            return "S"; //Single
        else if(enType == ctRegular)
            return "R"; //Regular
        else
            return "N"; //Nil
    }

    bool operator == (const St_Candidate& rhs) const // 我们在这里不对direction进行比较
    {
        if( this->strChromName == rhs.strChromName &&
            this->iStartPos == rhs.iStartPos &&
            this->iEndPos == rhs.iEndPos )
            return true;
        else
            return false;
    }

    bool CheckOneSideSimilarity(St_Candidate& rhs)
    {
        int iMaxDiff = 8;
        if( this->strChromName == rhs.strChromName &&
            (abs(this->iStartPos - rhs.iStartPos) < iMaxDiff ||
             abs(this->iEndPos - rhs.iEndPos) < iMaxDiff))
            return true;
        else
            return false;
    }

    bool CheckSimilarity(St_Candidate& rhs)
    {
        int iMaxDiff = 8;
        if( this->strChromName == rhs.strChromName &&
            abs(this->iStartPos - rhs.iStartPos) < iMaxDiff &&
            abs(this->iEndPos - rhs.iEndPos) < iMaxDiff)
            return true;
        else
            return false;
    }

    bool operator < (const St_Candidate& rhs) const
    {
        return this->iStartPos < rhs.iStartPos;
    }
};


class ClsClassificaiton
{
public:
    ClsClassificaiton();
    ~ClsClassificaiton();

public:
    void ClassifyCirc();
    void Init(string strV1, string strV2,
              string strV3, string strV4, string strV5,
              int iV5, int iV6, bool bV7);

private:
    void ParseDBGResult();
    void ParseBrief(string strCurPath, vector<St_Candidate>& vCandi);

    void FilterCircCandi();

    void SetBlastRefForCandi();//(vector<St_Row_Chrom>& vChrom);
    void CreateNormalRef(string& strExonDonorSeq, string& strExonAcceptorSeq,  St_Candidate& stCandi);
    void CreateDonorSelfRef(string& strExonDonorSeq, St_Candidate& stCandi);
    void CreateAcceptorSelfRef(string& strExonAcceptorSeq, St_Candidate& stCandi);

    void CheckCircCategory();

    void CheckBlastResult(St_BlastResult& stBlastResult, St_Candidate& stCandi,
                          string& strCurReads, ofstream& ofs);
    void GetSpecialMappingCase();

    void CheckChimeric();
    bool CheckDonorSelfAlign(string strQueryFa, St_Candidate& stCandi, ClsBlast* pBlast);
    bool CheckAcceptorSelfAlign(string strQueryFa, St_Candidate& stCandi, ClsBlast* pBlast);

private:
    string m_strDBGResult;
    string m_strGTF;
    string m_strRef;
    string m_strCircCandi;
    string m_strBlastRootFolder;
    int m_iReadLen;
    int m_iKmerLen;
    bool m_bShowFlank;

    vector<St_Candidate> m_vCandi;
};

#endif // CLSCLASSIFICAITON_H
