#ifndef CLSCIRCRNADETECTION_H
#define CLSCIRCRNADETECTION_H

#include "clsconfig.h"
#include <vector>
#include "../../ShareLibrary/clsgtfparse.h"
#include "../../ShareLibrary/clsfastqreader.h"
#include "clsdebruijngraph.h"

using namespace std;

enum En_CircType{ctSelf=0, ctRegular, ctMax};
enum En_HitSeed{hsPre=0, hsNext, hsMax};

struct St_Candidate
{
    //unsigned char ucChromIndex;
    string strChromName; // it could be the integer, it also could be the letters
    int iStartPos;
    int iEndPos;
    //-->
    int iSupportNum; // Especially For My result
    //<--
    bool bRC;
    En_CircType enType;
    string strTag;
    bool bHit;

    St_Raw_Exon* pExon;  //This is Donor
    St_Raw_Exon* pExonAcceptor; // Only For Reg Case

    St_Candidate()
    {
        //ucChromIndex = 0;
        strChromName = "";
        iStartPos = 0;
        iEndPos = 0;
        iSupportNum = 0;
        bRC = false;
        enType = ctMax;
        strTag = "xx";
        bHit = false;
        pExon = NULL;
        pExonAcceptor = NULL;
    }

    St_Candidate(string strV1, int iV2, int iV3, int iV4=1, bool bV5=false,
                 En_CircType enV6=ctMax, string strV7="xx", St_Raw_Exon* pV8=NULL, St_Raw_Exon* pV9=NULL)
    {
        strChromName = strV1;
        iStartPos = iV2;
        iEndPos = iV3;
        iSupportNum = iV4;
        bRC = bV5;
        enType = enV6;
        strTag = strV7;
        pExon = pV8;
        pExonAcceptor = pV9;
    }

    void SetExon(St_Raw_Exon* pExon)
    {
        this->pExon = pExon;
    }

    void SetExonAcceptor(St_Raw_Exon* pExon)
    {
        this->pExonAcceptor = pExon;
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

    int GetDonorLength()
    {
        return pExon == NULL ? -1 : pExon->GetLength();
    }

    int GetAcceptorLength()
    {
        return pExonAcceptor == NULL ? -1 : pExon->GetLength();
    }

    int GetLength()
    {
        if(this->pExon != NULL)
        {
            if(this->pExonAcceptor != NULL) //Regular case: both valid
            {
                return this->pExon->GetLength() + this->pExonAcceptor->GetLength();
            }
            else //self case
                return this->pExon->GetLength();
        }
        else
        {
            if(this->pExonAcceptor != NULL) //self case
            {
                return this->pExonAcceptor->GetLength();
            }
            else //error case
                return -1;
        }
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

struct St_HitNode
{
    unsigned int uiValue;
    St_Node* pNode;
    int iPreGapLen;

    St_HitNode()
    {
        uiValue = 0;
        pNode = NULL;
        iPreGapLen = 0;
    }

    St_HitNode(unsigned int uiV1, St_Node* pV2, int iV3)
    {
        uiValue = uiV1;
        pNode = pV2;
        iPreGapLen = iV3;
    }
};


struct St_ExonPosCount
{
    St_ExonPos stExonPos;
    unsigned char ucCount; // maximum value = 256

    St_ExonPosCount()
    {
        ucCount = 0;
    }

    St_ExonPosCount(St_ExonPos stV1, unsigned char ucV2)
    {
        stExonPos = stV1;
        ucCount = ucV2;
    }

    bool operator == (St_ExonPosCount stV1)
    {
        if(this->stExonPos == stV1.stExonPos)
            return true;
        else
            return false;

    }
};

struct St_CandiExAssistSelf
{
    St_Candidate stCandi;
    St_ExonPos stExonPos;
    int iCount;
    int iAddNodHitNum;

    St_CandiExAssistSelf():iCount(0), iAddNodHitNum(0)
    {}

    St_CandiExAssistSelf(St_Candidate& stCandi, St_ExonPos& stExonPos, int iCount)
    {
        this->stCandi = stCandi;
        this->stExonPos = stExonPos;
        this->iCount = iCount;
    }

    int GetLength()
    {
        return stCandi.GetLength();
    }
};

struct St_CandiExAssistReg
{
    St_Candidate stCandi;
    St_ExonPos stExonPosDonor; // 就是后面的一个exon
    St_ExonPos stExonPosAcceptor; //就是前面的一个exon
    int iCount;
    int iAddNodHitNum;

    St_CandiExAssistReg():iCount(0), iAddNodHitNum(0)
    {}

    St_CandiExAssistReg(St_Candidate& stCandi, St_ExonPos& stExonPs1, St_ExonPos& stExonPos2, int iCount)
    {
        this->stCandi = stCandi;
        this->stExonPosDonor = stExonPs1;
        this->stExonPosAcceptor = stExonPos2;
        this->iCount = iCount;
        this->iAddNodHitNum = 0;
    }

    int GetLength()
    {
        return stCandi.GetLength();
    }

    bool IsNeighbor()
    {
        int iDiff = abs((int)this->stExonPosDonor.ucExonIndex - (int)this->stExonPosAcceptor.ucExonIndex);
        if(iDiff == 1)
            return true;
        else
            return false;
    }
};

class ClsCircRNADetection
{
public:
    ClsCircRNADetection();
    ~ClsCircRNADetection();

private:
    string m_strReads1Path;
    string m_strReads2Path;
    int m_iKmerLen;
    float m_fKmerRatio;
    int m_iReadLen;
    int m_iMinHitNum;
    int m_iMaxSupportNum;
    vector<St_Row_Chrom>* m_pvChrom;

private:
    //vector<St_Candidate> m_vSelfCircCandi;
    //vector<St_Candidate> m_vRegCircCandi;

public:
    void Init(St_Config& stConfig, vector<St_Row_Chrom>& vChrom);
    void FindCirc(unordered_map<unsigned int, St_Node>& mpSDBG,
                  St_Fasta* pCurRefFa,
                  vector<St_Fastq>& vReads, string strChromName);

private:
    int CheckSampling(unsigned int* arrySamplingKmer,
                      unordered_map<unsigned int, St_Node>& mpSDBG);
    void CheckHittingCase(string strSeq,
                          unordered_map<unsigned int, St_Node>& mpSDBG, bool bSeqRC);
    bool CheckSamePos(vector< vector<St_ExonPosCount> >& vvExonPosCount, vector<St_ExonInfo>& vSENode);
    bool FindMaxValidHit(vector<St_ExonPosCount>& vExonPosCount, St_ExonPos& stExonPos);

    bool CheckSelfCircularCase(vector<St_HitNode>& vHitNode, St_ExonPos& stExonPos, St_Candidate& stValidCandi);

    void PrintCircCandidate(vector<St_Candidate>& vSelfCircCandi,
                            vector<St_Candidate>& vRegCircCandi, string strChromName);

    //New: for the method based on DBG
    void CheckHittingCaseByDBG(string strSeq,
                               unordered_map<unsigned int, St_Node>& mpSDBG,
                               St_Fasta* pCurRefFa,
                               bool bSeqRC, string strChromName,
                               vector<St_Candidate>& vSelfCircCandi, vector<St_Candidate>& vRegCircCandi,
                               ofstream& ofsDebug, ofstream& ofsDBGResult);

    bool CheckDBGSeedHit(vector<St_ExonInfo>& vEI, vector<int>& vCount,
                         string& strSeq, int iSeedPos,
                         unordered_map<unsigned int, St_Node>& mpSDBG, En_HitSeed enHitSeed);

    int KeepMaxCount(vector<St_ExonInfo>& vEI, vector<int>& vCount);
    bool CheckDuplicateInEI(vector<St_ExonInfo>& vEI);

    bool CheckCurSelfCaseByDBG(St_Candidate& stProbCandi, string& strSeq, St_ExonPos& stProbExonPos,
                               unordered_map<unsigned int, St_Node>& mpSDBG,
                               vector<St_Candidate>& vSelfCircCandi);

    bool CheckCurRegCaseByDBG(St_Candidate& stProbCandi, string& strSeq,
                              St_ExonPos& stProbExonPosDonor, St_ExonPos& stProbExonPosAcceptor,
                              unordered_map<unsigned int, St_Node>& mpSDBG,
                              St_Candidate& stCurCandi, vector<St_Candidate>& vRegCircCandi);

    int CheckDBGTail(St_Link& stCurLink, int& iPathFullLen,
                     vector<unsigned int>& vKmerNode, vector<char>& vTagTail,
                     St_ExonPos& stProbExonPos, unsigned int uiPreSeed,
                     unordered_map<unsigned int, St_Node>& mpSDBG);
    int CheckDBGHead(St_Link& stCurLink, int& iPathFullLen,
                     vector<unsigned int>& vKmerNode, vector<char>& vTagHead,
                     St_ExonPos& stProbExonPos, unsigned int uiNextSeed,
                     unordered_map<unsigned int, St_Node>& mpSDBG);
    int CheckIfExonPosInEI(St_ExonPos& stExonPos, vector<St_ExonInfo>& vEI);
    int LooseMatch(string strSeq1, string strSeq2, int iLen, int iPos);

    void CollectCandiAssistSelf(vector<St_ExonInfo>::iterator itrEI,
                                vector<St_CandiExAssistSelf>& vCandiAssistSelf,
                                string strChromName);

    void CollectCandiAssistReg(vector<St_ExonInfo>::iterator itrEIHead,
                               vector<St_ExonInfo>::iterator itrEITail,
                               vector<St_CandiExAssistReg>& vCandiAssistReg,
                               string strChromName);

    void PurifySelfCase(vector<St_CandiExAssistSelf>& vCandiAssistSelf,
                        unordered_map<unsigned int, St_Node>& mpSDBG,
                        St_Fasta* pCurRefFa,
                        string& strSeq, bool bSeqRC,
                        unordered_map<unsigned int, St_Node>::iterator& itrAdd1,
                        unordered_map<unsigned int, St_Node>::iterator& itrAdd2,
                        vector<St_Candidate>& vSelfCircCandi,
                        ofstream& ofsDBGResult);

    void PurifyRegCase(vector<St_CandiExAssistReg>& vCandiAssistReg,
                       unordered_map<unsigned int, St_Node>& mpSDBG,
                       St_Fasta* pCurRefFa,
                       string& strSeq, bool bSeqRC,
                       unordered_map<unsigned int, St_Node>::iterator& itrAdd1,
                       unordered_map<unsigned int, St_Node>::iterator& itrAdd2,
                       vector<St_Candidate>& vRegCircCandi,
                       ofstream& ofsDBGResult);
};

#endif // CLSCIRCRNADETECTION_H
