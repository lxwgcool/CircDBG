#ifndef CLSDEBRUIJNGRAPH_H
#define CLSDEBRUIJNGRAPH_H

#include "../../ShareLibrary/clsgtfparse.h"
#include "../../ShareLibrary/clskmeralgorithm.h"
#include "../../ShareLibrary/clsfastareader.h"
#include "clsconfig.h"
#include <stdint.h>
#include <unordered_map>

using namespace std;


struct St_ExonPos  // For each
{
    unsigned char ucChromIndex;      // For Human: maximum 256
    uint16_t uiGeneIndex;              // For Human: maximum less than 32768
    unsigned char ucTranscriptIndex; // For Human: maximum should less than 256
    unsigned char ucExonIndex;       // For Human: maximum less than 256    
    bool bRC;

    St_ExonPos()
    {
        ucChromIndex = 0;
        uiGeneIndex = 0;
        ucTranscriptIndex = 0;
        ucExonIndex = 0;
        bRC = false;
    }

    void Clear()
    {
        ucChromIndex = 0;
        uiGeneIndex = 0;
        ucTranscriptIndex = 0;
        ucExonIndex = 0;
        bRC = false;
    }

    bool operator == (const St_ExonPos& stRightExonPos)
    {
        if(this->ucChromIndex == stRightExonPos.ucChromIndex &&
           this->uiGeneIndex == stRightExonPos.uiGeneIndex &&
           this->ucTranscriptIndex == stRightExonPos.ucTranscriptIndex &&
           this->ucExonIndex == stRightExonPos.ucExonIndex)
        {
            return true;
        }
        else
            return false;
    }

    bool operator < (const St_ExonPos& stRightExonPos) const
    {
        if(this->ucChromIndex <= stRightExonPos.ucChromIndex)
            return true;
        else if(this->ucChromIndex <= stRightExonPos.uiGeneIndex)
            return true;
        else if(this->ucChromIndex <= stRightExonPos.ucTranscriptIndex)
            return true;
        else if(this->ucChromIndex <= stRightExonPos.ucExonIndex)
            return true;
        else
            return false;
    }
};

struct St_Link  // How to know if it is head node or none head node? --> Go!!
{
    unsigned int uiPreKmer;
    unsigned int uiNextKmer;

    /* This tag means if it is:
     * "S": In the beginning part of current exon: Start to KmerRatio * ReadLen
     * "M": In the middle part of current exon: Otherwise
     * "E": In the ending part of current exon: (ExonLen  - KmerRatio * ReadLen) to End
     * 我们需要新增一些tag,来标示real head 和 read tail 的部分：
     *   "H": the first 10 Nodes;
     *   "T": The last 10 Nodes;
     *   有效的hit 必须包含至少三个 --> Cool --> I think 30% is a premium percentage for threshould identification
     */
    char cCurTag;  // Current kmer Tag
    char cPreTag;  // Pre-kmer Tag
    char cNextTag; // Next-Kmer Tag

    St_Link()
    {
        uiPreKmer = 0;
        uiNextKmer = 0;
        cCurTag = 'U';
        cPreTag = 'U';
        cNextTag = 'U';
    }
};


//*********************************************************
//************* Simple DBG (De Bruijn Graph) **************
//*********************************************************
struct St_ExonInfo // SE: means single Exon Node
{
    //unsigned int uiValue; // Current Kmer Sequence Value
    vector<St_Link> vLink; //This is link
    St_ExonPos stExonPos;    
    St_ExonInfo()
    {}

    bool CheckIsSameTranscript(St_ExonInfo& stEI)
    {
        if(this->stExonPos.ucChromIndex == stEI.stExonPos.ucChromIndex &&
           this->stExonPos.uiGeneIndex == stEI.stExonPos.uiGeneIndex &&
           this->stExonPos.ucTranscriptIndex == stEI.stExonPos.ucTranscriptIndex)
        {
            return true;
        }
        else
            return false;
    }

    //--> 我们默认，对于输入而言，当前的为Back
    // 我们需要 当前的：
    // (1) 在正常的方向上在后面
    // (2) 在RC的方向上在前面
    bool CheckIsBackToBegin(St_ExonInfo& stEI) // RC: Reverse Complementary
    {
        if(CheckIsSameTranscript(stEI))
        {
            if(!this->stExonPos.bRC) // Normal Direction
            {
                if(this->stExonPos.ucExonIndex > stEI.stExonPos.ucExonIndex)
                    return true;
                else
                    return false;
            }
            else // Reverse Complementory
            {
                if(this->stExonPos.ucExonIndex < stEI.stExonPos.ucExonIndex)
                    return true;
                else
                    return false;
            }
        }
        else
            return false;
    }
};

struct St_Node //Normal Node: 这个是在于合并之后的从全局(也就是所有的exon merge之后的角度出发)的相应的node
{
    vector<St_ExonInfo> vEI; // EI means Exon Info
};

//struct St_DBG
//{
//    int iChromIndex;
//    string strChromName;
//    unordered_map<unsigned int, St_Node> mpDBG;

//    St_DBG():iChromIndex(-1), strChromName("")
//    {}

//    void Reset()
//    {
//        iChromIndex = -1;
//        strChromName = "";
//        mpDBG.clear();
//    }
//};

class ClsDeBruijnGraph
{
public:
    ClsDeBruijnGraph();
    ~ClsDeBruijnGraph();

private:
    //To record the configuration value
    string m_strDNARefPath;
    int m_iKmerLen;
    int m_iReadLen;
    float m_fKmerRatio;
    //For DBG
    //vector<St_DBG> m_vDBG;

public:
    void Init(St_Config& stConfig);

    //-->Build DBG Graph
    void BuildGraph(unordered_map<unsigned int, St_Node>& mpDBG,
                    St_Row_Chrom* pChrom, St_Fasta* pCurRefFa, int iCurChromIndex);
//    vector<St_DBG>& GetDBG();

private:
    void BuildDBGForSingleExon(unordered_map<unsigned int, St_Node>& mpDBG,
                               St_Raw_Exon& stExon, St_ExonPos& stExonPos,
                               string& strRef, ofstream& ofsExon);

    void BuildByShortExon(unordered_map<unsigned int, St_ExonInfo>& mpSEDBG,
                          string& strExonSeq, St_ExonPos& stExonPos);

    void BuildByRegularExon(unordered_map<unsigned int, St_ExonInfo>& mpSEDBG,
                            string& strHeadPartSeq, string& strTailPartSeq, St_ExonPos& stExonPos);

    void UpdateSingleExonDBG(unordered_map<unsigned int, St_ExonInfo>& mpSEDBG,
                             St_ExonPos& stExonPos, char cTag, unsigned int uiKmer,
                             unsigned int uiPreKmer, int iItrIndex);    
};

#endif // CLSDEBRUIJNGRAPH_H
