#include "clscircrnadetection.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <algorithm>
#include <map>

//Some Macro-Variates
const int SAMPLINGNUM = 8;
const float ARRYPOSRATIO[SAMPLINGNUM] = {.1, .2, .3, .4, .5, .6, .7, .8};

const int CSAMPLINGHITMIN = 2;
const int CSDBGFULLLENMIN = 5;
const int MAXHITDIFF = 5; // the hitting length difference (exon compare with exon)
const int PREMIUMLEN = 10;

#define CIRC_CANDI_PRINT

//#define _DEBUG

ClsCircRNADetection::ClsCircRNADetection()
{
}

ClsCircRNADetection::~ClsCircRNADetection()
{
}

void ClsCircRNADetection::Init(St_Config& stConfig, vector<St_Row_Chrom>& vChrom)
{
    m_strReads1Path = stConfig.strReads1Path;
    m_strReads2Path = stConfig.strReads2Path;
    m_iKmerLen = stConfig.iKmerLen;
    m_iReadLen = stConfig.iReadLen;
    m_fKmerRatio = stConfig.fKmerRatio;
    m_iMaxSupportNum = stConfig.iMaxSupportNum;

    m_iMinHitNum = (m_iReadLen-m_iKmerLen + 1) * m_fKmerRatio; //* 2 * .6;  --> we need to simple the logic here

    m_pvChrom = &vChrom; // this is a const pointer
}

void ClsCircRNADetection::FindCirc(unordered_map<unsigned int, St_Node>& mpSDBG,
                                   St_Fasta* pCurRefFa,
                                   vector<St_Fastq>& vReads, string strChromName)
{    
    if(vReads.empty())
        return;

    //Clear candidate container
    vector<St_Candidate> vSelfCircCandi;
    vector<St_Candidate> vRegCircCandi;

    //cout << "Do real detection" << endl;
    ofstream ofsDebug;

#ifdef _DEBUG
    ofsDebug.open(("./Detection_Result/_debug_" + strChromName + ".txt").c_str());
    ofsDebug << "Do it" << endl;
#endif

    //2: Search CircRNA through DBG
    unsigned int arrySamplingKmer[SAMPLINGNUM];
    ofstream ofsDBGResult;
    ofsDBGResult.open(("./Detection_Result/DBG_Result_" + strChromName + ".txt").c_str());

    for(vector<St_Fastq>::iterator itr = vReads.begin(); itr != vReads.end(); itr++)
    {
        if((int)itr->strSeq.length() < m_iReadLen - 5) // 5 diffs torlerant
            continue;

        if(m_iReadLen < m_iKmerLen*2)
            continue;

        //(1) delete all of reads contain letter 'N'
        if(itr->strSeq.find('N') != string::npos) // 这个表示找到了
            continue;

        m_iReadLen = itr->strSeq.length();

        //如果没有N
        //(2) evenly sampling reads
        string strCurSeq = ""; // If this is valid reads could be used for further detection

        // a) normal direction
        //   i) get kmer ID
        bool bSeqRC = false;
        for(int i = 0; i < SAMPLINGNUM; i++)
        {
            string strTmp = itr->strSeq.substr(itr->strSeq.length()*ARRYPOSRATIO[i], m_iKmerLen);
            arrySamplingKmer[i] = ConvertKmerToNum32(strTmp);
        }

        int iOrgHitNumber = CheckSampling(arrySamplingKmer, mpSDBG);
        if(iOrgHitNumber >= SAMPLINGNUM * .3) //足够多  --> 也不要多到50%啦e
        {
            strCurSeq = itr->strSeq;
        }
        else
        {
            for(int i = 0; i < SAMPLINGNUM; i++)
            {
                string strTmp = GetReverseCompelement(itr->strSeq.substr(itr->strSeq.length()*ARRYPOSRATIO[i],
                                                                         m_iKmerLen));
                arrySamplingKmer[i] = ConvertKmerToNum32(strTmp);
            }
            int iRCHitNumber = CheckSampling(arrySamplingKmer, mpSDBG);

            if(iOrgHitNumber >= iRCHitNumber) //防止错误的，但是侥幸对上了一些
            {
                if(iOrgHitNumber >= CSAMPLINGHITMIN)
                    strCurSeq = itr->strSeq;
            }
            else
            {
                if(iRCHitNumber >= CSAMPLINGHITMIN)
                {
                    strCurSeq = GetReverseCompelement(itr->strSeq);
                    bSeqRC = true;
                }
            }
        }

//        //   ii) Check how may could be found back in vChrom

//        if(CheckSampling(arrySamplingKmer, mpSDBG))
//        {
//            strCurSeq = itr->strSeq;
//        }
//        else // b) RC direction
//        {
//            for(int i = 0; i < SAMPLINGNUM; i++)
//            {
//                string strTmp = GetReverseCompelement(itr->strSeq.substr(itr->strSeq.length()*ARRYPOSRATIO[i],
//                                                                         m_iKmerLen));
//                arrySamplingKmer[i] = ConvertKmerToNum32(strTmp);
//            }
//            if(CheckSampling(arrySamplingKmer, mpSDBG))
//            {
//                strCurSeq = GetReverseCompelement(itr->strSeq);
//                bSeqRC = true;
//            }
//        }

        //(3)Go Further Detection
        if(strCurSeq != "")
        {
//#ifdef CIRC_CANDI_PRINT
//         cout << strCurSeq << endl << endl;
//#endif
            /// 2:Find Hitting Case
            // Go Further Deteciton ->Go
            //cout << strCurSeq << endl;
            //CheckHittingCase(strCurSeq, mpSDBG, bSeqRC);
            CheckHittingCaseByDBG(strCurSeq, mpSDBG, pCurRefFa,
                                  bSeqRC, strChromName,
                                  vSelfCircCandi, vRegCircCandi,
                                  ofsDebug, ofsDBGResult); //纯靠DBG的方法来干这个事儿
        }
    }

    //3: know the order of those exons. (based on the hitting position from reads)
    //Output all of m_vSelfCircCandi and Output all of m_vRefCircCandi
    //cout << "Print Circ Candidate" << endl;
    PrintCircCandidate(vSelfCircCandi, vRegCircCandi, strChromName);

#ifdef _DEBUG
    ofsDebug.close();
#endif
}

int ClsCircRNADetection::CheckSampling(unsigned int* arrySamplingKmer,
                                       unordered_map<unsigned int, St_Node>& mpSDBG)
{
    int iHitNumber = 0;
    for(int i = 0; i< SAMPLINGNUM; i++)
    {
        if(mpSDBG.find(arrySamplingKmer[i]) == mpSDBG.end()) // do not find
            continue;
        iHitNumber++;
    }

    return iHitNumber;
//    if(iHitNumber < CSAMPLINGHITMIN)
//        return false;
//    return true;
}

//Check Hitting Status
bool Sort_Large_To_Small_ExonHitCount(vector<St_ExonPosCount> vValue1, vector<St_ExonPosCount> vValue2)
{
    if(vValue1.size() > vValue2.size())
        return true;
    else
        return false;
}

void ClsCircRNADetection::CheckHittingCase(string strSeq,
                                           unordered_map<unsigned int, St_Node>& mpSDBG,
                                           bool bSeqRC)
{
    /*
    //Now, we need to do the core part of checking hitting
    //1: Find All of Hitting Node:
    vector<St_HitNode> vHitNode; // 这里我们并不合并重复的点
    int iPreHitPos = 0;
    for(unsigned int i = 0; i < strSeq.length() - m_iKmerLen + 1; i++)
    {
        unsigned int uiKmer = ConvertKmerToNum32(strSeq.substr(i, m_iKmerLen));
        if(mpSDBG.find(uiKmer) != mpSDBG.end()) // 重复的没有关系，直接放进去不用担心
        {
            ////Debug Print --> <--

            int iPreGapLen = i - iPreHitPos;
            vHitNode.push_back(St_HitNode(uiKmer, &(mpSDBG[uiKmer]), iPreGapLen));

            iPreHitPos = i + 1;
        }
    }

    ////Debug Print --> <--

    //2: Merge and pick up the most hit exons.
    /// (1) 将前后一致的给合并起来 (这里不见的一定是要相邻的)  --> 这里主要着眼于找相应的exon的值 --> comes from the same exon
    vector< vector<St_ExonPosCount> > vvExonHitCount;  // use vector here is correct!!

    for(vector<St_HitNode>::iterator itr = vHitNode.begin(); itr != vHitNode.end(); itr++)
    {
        if(CheckSamePos(vvExonHitCount, itr->pNode->vEI)) // 这里的same pos不是指的一个地方，而是很多个地方的hit,
                                                           // 然后hit的情况都是一样的
        {
            continue;
        }
        else // 这是不一样的hit的位置(们)
        {
            //for each hit node
            vector<St_ExonPosCount> vEPCount; //Exon Pos Count
            for(vector<St_ExonInfo>::iterator subItr = itr->pNode->vEI.begin();
                subItr != itr->pNode->vEI.end(); subItr++)
            {
                vEPCount.push_back(St_ExonPosCount(subItr->stExonPos, 1));
            }
            vvExonHitCount.push_back(vEPCount);
        }
    }

    if(vvExonHitCount.empty())
        return;

    sort(vvExonHitCount.begin(), vvExonHitCount.end(), Sort_Large_To_Small_ExonHitCount);

    ////Debug Print --> <--

    ///2: 我们再找到最好的点，然后再把他们merge到一起  --> later move on   --> 还有个逻辑忘记了，当hit数目一样的时候要找length最短的那个exon
    /// 将子集逐一合并进去  --> 从后往前进行遍历 --> 我觉得还是有多少融多少
    for(vector< vector<St_ExonPosCount> >::iterator itr = vvExonHitCount.end()-1;
        itr >= vvExonHitCount.begin(); itr--)
    {
        for(vector< vector<St_ExonPosCount> >::iterator subItr = itr-1;
            subItr >= vvExonHitCount.begin(); subItr--)
        {
            //如果是完全的子集，那么就融进去  --> 有一部分我们也要融进去 --> Go!!
            for(vector<St_ExonPosCount>::iterator itrMain = subItr->begin();
                itrMain != subItr->end(); itrMain++)
            {
                if(itr->empty())
                    break;
                vector<St_ExonPosCount>::iterator itrFind = find(itr->begin(), itr->end(), *itrMain);
                if(itrFind != itr->end()) // 证明找到了，那么这一下我们就要融进去了
                {
                    //Update itrMain
                    itrMain->ucCount += itrFind->ucCount;
                    //erase the old one
                    itr->erase(itrFind);
                }
            }
            if(itr->empty())
                break;
        }
        if(itr->empty())
        {
            vvExonHitCount.erase(itr);
        }
    }

    ////Debug Print --> <--

    ///3: 将一些hit比较少的case 给滤掉
    int iCountMinThreshould = (m_iReadLen - m_iKmerLen + 1) * .2;
    for(vector< vector<St_ExonPosCount> >::iterator itr = vvExonHitCount.end()-1; itr >= vvExonHitCount.begin(); itr--)
    {
        for(vector<St_ExonPosCount>::iterator subItr = itr->end() - 1; subItr >= itr->begin(); subItr--)
        {
            if(subItr->ucCount < iCountMinThreshould)
            {
                itr->erase(subItr);
            }
        }
        if(itr->empty()) // 全部不符合要求
        {
            vvExonHitCount.erase(itr);
        }
    }

    ///4: For single exon hitting case: pick out the hitting case with number of kmer > min threshould
    vector<St_ExonPosCount> vPotentialSelfCirc;
    if(vvExonHitCount.size() == 1)
    {
        vPotentialSelfCirc = vvExonHitCount[0];
    }
    else
    {
        for(vector< vector<St_ExonPosCount> >::iterator itr = vvExonHitCount.begin();
            itr != vvExonHitCount.end(); itr++)
        {
            for(vector<St_ExonPosCount>::iterator subItr = itr->begin(); subItr != itr->end(); subItr++)
            {
                if(subItr->ucCount > m_iMinHitNum)
                {
                    vPotentialSelfCirc.push_back(*subItr);
                }
            }
        }
    }

    /// Find the top 5: those top 3 should satisfy some conditions  --> this is just try to consider all of cases without any filter.
    if(!vPotentialSelfCirc.empty()) // 找到这个里面 最大的，然后这个case肯定就是self circle  --> Self Circ Case
    {
        St_ExonPos stExonPos;

        if(FindMaxValidHit(vPotentialSelfCirc, stExonPos))
        {
            //--> 基于这个去查找相应的back splicing
            St_Candidate stValidCandi;
            if(CheckSelfCircularCase(vHitNode, stExonPos, stValidCandi))
            {
#ifdef CIRC_CANDI_PRINT
                cout << "Start & End: " << IntToStr(stValidCandi.iStartPos) << " : "
                     << IntToStr(stValidCandi.iEndPos) << endl;

                cout << "RC: " << (stValidCandi.bRC ? "Yes" : "No") << endl;

                cout << "Current Reads: " << (bSeqRC ? "RC" : "Reg") << endl
                     << strSeq << endl << endl;
#endif
            }
        }

//        vector<St_ExonPosCount> vTmp;
//        for(vector<St_ExonPosCount>::iterator itr = vPotentialSelfCirc.begin();
//            itr != vPotentialSelfCirc.end(); itr++)
//        {
//            vTmp.push_back(*itr);
//            if(FindMaxValidHit(vTmp, stExonPos))
//            {
//                //--> 基于这个去查找相应的back splicing
//                St_Candidate stValidCandi;
//                if(CheckSelfCircularCase(vHitNode, stExonPos, stValidCandi))
//                {
////                    cout << "Start & End: " << IntToStr(stValidCandi.iStartPos) << " : " << IntToStr(stValidCandi.iEndPos) << endl;
////                    cout << "RC: " << (stValidCandi.bRC ? "Yes" : "No") << endl;
////                    cout << "Current Reads: " << (bSeqRC ? "RC" : "Reg") << endl
////                         << strSeq << endl << endl;
//                }
//            }
//            vTmp.clear();
//        }
    }
    else if(vvExonHitCount.size() == 2)  // Move on
    {
        //cout << "Size 2" << endl;

    }
    else if(vvExonHitCount.size() > 2)
    {
        //cout << "Size > 2" << endl;
    }
    else
    {}*/
}



bool ClsCircRNADetection::CheckSamePos(vector< vector<St_ExonPosCount> >& vvExonPosCount,
                                       vector<St_ExonInfo>& vEI)  //target: combine the same types of EI (everything is the same)
{
    bool bSame = false;
    for(vector< vector<St_ExonPosCount> >::iterator itr = vvExonPosCount.begin();
        itr != vvExonPosCount.end(); itr++)
    {
        if(vEI.size() != itr->size())
            continue;
        else
        {
            int iCount = 0;
            for(vector<St_ExonPosCount>::iterator subItr = itr->begin(); subItr != itr->end(); subItr++)
            {
                for(vector<St_ExonInfo>::iterator itrEI = vEI.begin();
                    itrEI != vEI.end(); itrEI++)
                {
                    if(subItr->stExonPos == itrEI->stExonPos)
                    {
                        iCount++;
                        break;
                    }
                }
            }
            if(iCount == (int)itr->size())
            {
                //Update Count Number
                for(vector<St_ExonPosCount>::iterator subItr = itr->begin(); subItr != itr->end(); subItr++)
                {
                    subItr->ucCount++;
                }
                bSame = true;
                break;
            }
        }
    }
    return bSame;
}

bool ClsCircRNADetection::FindMaxValidHit(vector<St_ExonPosCount>& vExonPosCount, St_ExonPos& stExonPos)
{
    // 如果hit一样那么我们选择4，短的exon
    int iHit = 0;
    St_ExonPos stTmpExonPos;
    for(vector<St_ExonPosCount>::iterator itr = vExonPosCount.begin(); itr != vExonPosCount.end(); itr++)
    {
        if(itr->ucCount > iHit) // 如果hit的多，那么我们就选择hit最多的
        {
            stTmpExonPos = itr->stExonPos;
            iHit = itr->ucCount;
        }
        else if( abs(itr->ucCount - iHit) <= MAXHITDIFF ) // 选择最短的
        {
            if((*m_pvChrom)[stTmpExonPos.ucChromIndex].vRG[stTmpExonPos.uiGeneIndex].vRT[stTmpExonPos.ucTranscriptIndex].vRExon[stTmpExonPos.ucExonIndex].GetLength() >
               (*m_pvChrom)[itr->stExonPos.ucChromIndex].vRG[itr->stExonPos.uiGeneIndex].vRT[itr->stExonPos.ucTranscriptIndex].vRExon[itr->stExonPos.ucExonIndex].GetLength() )
            {
                stTmpExonPos = itr->stExonPos; // 我们选择更短的exon
                iHit = itr->ucCount;
            }
        }
    }

    if(iHit >= (m_iReadLen - m_iKmerLen + 1) * .2) // 不能太小
    {
        stExonPos = stTmpExonPos;
        return true;
    }
    else
        return false;
}

//Check Slef Circular Case --> 今天把self circular case 给做好
struct St_HitEI // Single Exon Node
{
    unsigned int uiValue;
    St_ExonInfo stEI;
    int iPreGapLen;

    St_HitEI()
    {}

    St_HitEI(unsigned int uiV1, St_ExonInfo& stV2, int iV3)
    {
        uiValue = uiV1;
        stEI = stV2;
        iPreGapLen = iV3;
    }
};

enum En_SeqPart{spHead=0, spStart, spEnd, spTail, spMax};

bool ClsCircRNADetection::CheckSelfCircularCase(vector<St_HitNode>& vHitNode, St_ExonPos& stExonPos,
                                                St_Candidate& stValidCandi)
{
    /*
    St_Raw_Exon& stExon = (*m_pvChrom)[stExonPos.ucChromIndex].vRG[stExonPos.uiGeneIndex].\
                          vRT[stExonPos.ucTranscriptIndex].vRExon[stExonPos.ucExonIndex];

    ///here: we find some of exons which contain circular rna do not have the typical 2bps signal around boundary
    if(!stExon.GetIsSupportCircRNA()) //GetIsBothSupportCircRNA())
    {
        return false;
    }

    // 我们主要在这里看，这个exon是不是相应的circular rna case
    //1: 把这些single exon information (EI) 拿出来
    vector<St_HitEI> vHitEI;
    bool bValid = false;
    int iJumpNode = 0;
    //2： 将不属于该exon的相应的点依照次序过滤掉
    for(vector<St_HitNode>::iterator itr = vHitNode.begin(); itr !=  vHitNode.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_ExonInfo>::iterator subItr = itr->pNode->vEI.begin();
            subItr != itr->pNode->vEI.end(); subItr++)
        {
            if(subItr->stExonPos == stExonPos)
            {
                bFind = true;
                vHitEI.push_back(St_HitEI(itr->uiValue, *subItr, itr->iPreGapLen + iJumpNode));

                //string strCurKmer = ConvertNum32ToKmer(itr->uiValue, m_iKmerLen);

                iJumpNode = 0;

//                if(vHitSENode.size() == 30)
//                {
//                    int i = 0;
//                }

                break;
            }
        }

        if(!bFind)
          iJumpNode += itr->iPreGapLen + 1;
    }

    //---> Temporary
//    for(vector<St_HitSENode>::iterator itr = vHitSENode.begin(); itr != vHitSENode.end(); itr++)
//    {
//        if(itr->uiValue == ConvertKmerToNum32("TGCGTTTCCTGGCTG"))
//            cout << "1 " << itr->stSENode.vTag[0] << endl;

//        if(itr->uiValue == ConvertKmerToNum32("GCGTTTCCTGGCTGG"))
//            cout << "2 " << itr->stSENode.vTag[0] << endl;

//        if(itr->uiValue == ConvertKmerToNum32("CGTTTCCTGGCTGGC"))
//            cout << "3 " << itr->stSENode.vTag[0] << endl;

//        if(itr->uiValue == ConvertKmerToNum32("GTTTCCTGGCTGGCT"))
//            cout << "4 " << itr->stSENode.vTag[0] << endl;

//        if(itr->uiValue == ConvertKmerToNum32("TTTCCTGGCTGGCTT"))
//            cout << "5 " << itr->stSENode.vTag[0] << endl;
//    }

    //<---

    //3: 顺序查找两个最长路径：因为vector中存储的是顺序 从前到后的匹配结果，因此在这里，我们从前往后寻找相应的最长路径，
    //   然后看看他们hit的相应的tag的变化趋势
    //4: 根据以往的经验，tag都应该是E to S
    //   Notice: 因为现在增加了新的case： H 和 T

    vector<char> vLongestPath;
    vector< vector<char> > vvValidPath;
    //-->
    vector<int> vLPGap; // 两个最长路径之间的距离  --> longest parth gap: notice: all of longest path will has their own gap length,
                                 // even if such longest path start at the beginning.
                        // 这个gap包含两个方面的东西 -->
                        // (1) 跳过的hit的点
                        // (2) 本身跳过的 un-hit 的 node.
    int iGapLen = 0;
    int iInnerJumpLen = 0;  // Try to record the distance betwee two valid long path
    //<--

    for(vector<St_HitEI>::iterator itr = vHitEI.begin(); itr != vHitEI.end(); itr++)
    {
        vLongestPath.clear();
        iGapLen += itr->iPreGapLen;
        //try to find the path
        iInnerJumpLen = 0;

        for(vector<St_HitEI>::iterator subItr = itr+1; subItr < vHitEI.end(); subItr++)
        {
            iInnerJumpLen += subItr->iPreGapLen;

            //先判断是不是最后一个
            if(subItr == vHitEI.end() - 1)
            {
                //check if current case could be found in the "next node list" of previous node
                for(unsigned int i = 0; i < itr->stEI.vLink.size(); i++)
                {
                    if(subItr->uiValue == itr->stEI.vLink[i].uiNextKmer)
                    {
                        char cTag = itr->stEI.vLink[i].cCurTag;
                        vLongestPath.push_back(cTag); //这个是itr related kmer当前的tag

                        //add the last as well, since this is the last node
                        cTag = (subItr->stEI.vLink.end()-1)->cCurTag;
                        vLongestPath.push_back(cTag); //这个是当前的tag

                        itr = subItr;
                        break;
                    }
                }

                if(vLongestPath.size() >= 3)
                {
                    vvValidPath.push_back(vLongestPath);
                    vLPGap.push_back(iGapLen);
                    iGapLen = 0;
                }
                else // 不满足条件
                {
                    iGapLen += iInnerJumpLen;
                }

                itr = subItr;
                break;
            }

            /// Go later--> For the majority case: does not the last one
            bool bFindNext = false;
            for(unsigned int i = 0; i < itr->stEI.vLink.size(); i++)
            {
                if(subItr->uiValue == itr->stEI.vLink[i].uiNextKmer)
                {
                    char cTag = itr->stEI.vLink[i].cCurTag;
                    vLongestPath.push_back(cTag); //这个是itr's 当前的tag

                    itr = subItr;
                    bFindNext = true;
                    break;
                }
            }

            if(!bFindNext) // 如果没找到next
            {
                if(vLongestPath.size() > 3) // 那么我们认为是靠谱的  --> 将所有valid 的longest path （path length > 3）
                {
                    vvValidPath.push_back(vLongestPath);
                    vLPGap.push_back(iGapLen);
                    iGapLen = 0;
                }
                else
                {
                    iGapLen += iInnerJumpLen - subItr->iPreGapLen;
                }
                break;
            }
        }
    }

    //将chain 转成相应的position tag (E, S or M)
    // New Important tag           (H and T)
    vector<char> vOrderedTag;

    // 记录在原始的longest path的路径中，tail path开始的index 以及 head path 开始的index
    // Notice: this is for current parth, rather than the original sequence
    int iHeadPathIndex = -1;  // head 应该是越靠前越好
    bool bHeadPathFind = false;
    int iTailPathIndex = -1;  // tail 应该是越靠后越好

    int iNumHead = 0;
    int iNumTail = 0;
    //在这里我们要记录相应的路径所对应的前后碰到的head 和 tail的情况--这里只要去计数就可以了  --> later go!!
    //我们首先要看多少head node 和多少 tail node 是我们所期待的->
    //Get current Exon
    int iMinPremiumNodeNum = PREMIUMLEN * .3;
    if(stExon.GetLength() - m_iKmerLen + 1 < PREMIUMLEN * 2) // 这么写也不存在错误，因为如果是长的length，那么这个肯定符合条件，
                                                             // 如果是小的length的exon那么我们直接用它的length去计算也是正确的
    {
        iMinPremiumNodeNum = (stExon.GetLength() - m_iKmerLen + 1) / 2;
    }
    int iValidHitSum = 0;
    for(vector< vector<char> >::iterator itr = vvValidPath.begin(); itr != vvValidPath.end(); itr++)
    {
        //---> 统计有效路径的长度总和
        iValidHitSum += itr->size();
        //<---

        map<char, int> mpTagCount;
        for(vector<char>::iterator subItr = itr->begin(); subItr != itr->end(); subItr++)
        {
            mpTagCount[*subItr]++;
        }
        //Find the maximum count of tag
        int iHSETCount[spMax] = {0, 0, 0, 0};

        for(map<char, int>::iterator subItr = mpTagCount.begin(); subItr != mpTagCount.end(); subItr++)
        {
            if(subItr->first == 'H')
                iHSETCount[spHead] = subItr->second;
            else if(subItr->first == 'S')
                iHSETCount[spStart] = subItr->second;
            else if(subItr->first == 'E')
                iHSETCount[spEnd] = subItr->second;
            else if(subItr->first == 'T')
                iHSETCount[spTail] = subItr->second;
            else
            {}
        }

        if(iHSETCount[spHead] != 0 && !bHeadPathFind)
        {
            iHeadPathIndex = itr - vvValidPath.begin();
            bHeadPathFind = true;
        }

        if(iHSETCount[spTail] != 0)
        {
            iTailPathIndex = itr - vvValidPath.begin();
        }

        iNumHead += iHSETCount[spHead];
        iNumTail += iHSETCount[spTail];

        char cTag = 'S';
        if( iHSETCount[spHead] + iHSETCount[spStart] < iHSETCount[spEnd] + iHSETCount[spTail] )
            cTag = 'E';

        vOrderedTag.push_back(cTag);
    }

    //Check Head Number and Tail Number
    //if(iNumHead < iMinPremiumNodeNum || iNumTail < iMinPremiumNodeNum)
    if(iNumHead < 1 || iNumTail < 1)
    {
        vHitEI.clear();
        return bValid;
    }

    // 这里先看看head 和 tail 的 index 关系能不能满足要求  ---->
    if(iHeadPathIndex == -1 || iTailPathIndex == -1)
    {
        vHitEI.clear();
        return bValid;
    }
    else
    {
        if(vvValidPath.size() == vLPGap.size())  // this is only for debug
        {
            //这里我们的距离是想从一个path的尾端 计算 到 另一个 parth的起始端，这么一来也就是说我们并不是考虑所有跨度的kmer
            //因此现在这里我们要减去 (kmerLen - 1). For example:

//                 =========             ==========
//                          -------------  -> This is what we want
//                    --------------------- --> 这种kmer包含了起始和终止位置的kmer，并不会四我们想要的！！！


            int iExpectGapLen = vLPGap[iHeadPathIndex] - (m_iKmerLen - 1);
            if(iExpectGapLen <= m_iKmerLen)   //注意，这里的这个逻辑，我们在实际上并不需要进行清除，我们可以作为一个分类
                                      //通过这个我们可以发现，这种case在find circ的结果中完全没有体现 101 - 54 ==> 95 - 54
                                      //pre value:12, Now, i set the value equal with kmer length
            {}  // 这个是符合标准的
            else
            {
                vHitEI.clear();
                return bValid;
            }
        }
        else
        {
            cout << "Logic Error !" << endl;
            vHitEI.clear();
            return bValid;
        }
    }
    //<------

//    //我们要在这里解决覆盖多次的问题，感觉还是有点复杂的  --> 我觉得直接用kmer hit 的数量作为threshold的判定是最好的 -->
//    int iBoundaryLen = m_iReadLen * m_fKmerRatio;
//    int iTotalExtractLen = iBoundaryLen * 2 + 2 * (m_iKmerLen - 1); // 这里我们是需要分开考虑的
//    int iKmerNum = (iBoundaryLen - 1) * 2;
//    if(stExon.GetLength() < iTotalExtractLen)
//        iKmerNum = stExon.GetLength() - m_iKmerLen;

//    if(iValidHitSum > iKmerNum * 1.2)
//    {
//        bValid = true;

//        //Put this case into the single self circular case
//        St_Candidate stCurCandi(stExonPos.ucChromIndex, stExon.iEnd, stExon.iStart, 1, stExon.bRC, ctSelf);

//        //Find if it is exsited
//        bool bFind = false;
//        for(vector<St_Candidate>::iterator itr = m_vSelfCircCandi.begin(); itr != m_vSelfCircCandi.end(); itr++)
//        {
//            if(*itr == stCurCandi)
//            {
//                itr->iSupportNum++;
//                bFind = true;
//                break;
//            }
//        }
//        if(!bFind) //找不到
//        {
//            m_vSelfCircCandi.push_back(stCurCandi);  //这里的理解方式应该为，从后往前，所以start是大的index，end是小的index
//        }

//        vHitSENode.clear();
//        return bValid;
//    }
//    //<----


    //Merge the same part: here we already know the all of the part are the valid path (length > 3)
    vector<char> vCompressedTag;

    for(vector<char>::iterator itr = vOrderedTag.begin(); itr != vOrderedTag.end(); itr++)
    {
        char cTag = *itr;
        for(vector<char>::iterator subItr = itr+1; subItr != vOrderedTag.end(); subItr++)
        {
            if(cTag == *subItr)
            {
                itr = subItr;
            }
            else
            {
                break;
            }
        }
        vCompressedTag.push_back(cTag);
    }

    //Now 现在可以来进行相应的判断了
    if(vCompressedTag.size() == 2 &&
       vCompressedTag[0] == 'E' && vCompressedTag[1] == 'S')
    {
        bValid = true;

        //Put this case into the single self circular case
        St_Candidate stCurCandi(stExonPos.ucChromIndex, stExon.iEnd, stExon.iStart, 1, stExon.bRC, ctSelf);
        stValidCandi = stCurCandi;

        //Find if it is exsited
        bool bFind = false;
        for(vector<St_Candidate>::iterator itr = m_vSelfCircCandi.begin(); itr != m_vSelfCircCandi.end(); itr++)
        {
            if(*itr == stCurCandi)
            {
                itr->iSupportNum++;
                bFind = true;
                break;
            }
        }
        if(!bFind) //找不到
        {
            m_vSelfCircCandi.push_back(stCurCandi);  //这里的理解方式应该为，从后往前，所以start是大的index，end是小的index
        }

#ifdef _DEBUG
#endif
    }
    else
    {}

    vHitEI.clear();  // 这里要debug一下 结果不正确
    return bValid;*/
}

void ClsCircRNADetection::PrintCircCandidate(vector<St_Candidate>& vSelfCircCandi,
                                             vector<St_Candidate>& vRegCircCandi,
                                             string strChromName)
{
    ofstream ofs;
    ofs.open(("./Detection_Result/Candidate_" + strChromName + ".txt").c_str());

    ofstream ofsBrief; // 这个是为了解析起来比较方便的文件
    ofsBrief.open("./Detection_Result/Brief_" + strChromName + ".txt");

    //1: Statistic Result
    ofs << "=============Statistic Number=============" <<endl;
    ofs << "Self Circular RNA   : " << IntToStr(vSelfCircCandi.size()) << endl;
    ofs << "Regular Circular RNA: " << IntToStr(vRegCircCandi.size()) << endl;
    ofs << "Total RNA Number    : " << IntToStr(vSelfCircCandi.size() + vRegCircCandi.size()) << endl << endl;

    map<unsigned int, int> mpDuplicateTest;

    //2: For the full set of Self Circular RNA
    ofs << "=============Self Circular RNA=============" <<endl;
    for(vector<St_Candidate>::iterator itr = vSelfCircCandi.begin(); itr != vSelfCircCandi.end(); itr++)
    {        
        ofs << strChromName << "\t" //IntToStr(itr->ucChromIndex) << "\t"
             << "<" << IntToStr(itr->iStartPos) << ", " << IntToStr(itr->iEndPos) << ">"
             << "\t" << IntToStr(itr->iSupportNum) << endl;
        ofsBrief << strChromName << " " // << IntToStr(itr->ucChromIndex) << " "
                 << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                 << IntToStr(itr->iSupportNum) << " "
                 //<< itr->strTag << " "
                 << itr->GetStrand() << " "
                 << "S" << endl;

        mpDuplicateTest[itr->iStartPos + itr->iEndPos] = itr->iSupportNum;
    }

    //3: For the full set of Regular Circular RNA
    int iDuplicatNum = 0;
    ofs << endl << "=============Regular Circular RNA=============" <<endl;
    for(vector<St_Candidate>::iterator itr = vRegCircCandi.begin(); itr != vRegCircCandi.end(); itr++)
    {
        if(mpDuplicateTest.find(itr->iStartPos + itr->iEndPos) != mpDuplicateTest.end())
        {
            iDuplicatNum++;
            continue;
        }

        ofs << strChromName << "\t" //IntToStr(itr->ucChromIndex) << "\t"
             << "<" << IntToStr(itr->iStartPos) << ", " << IntToStr(itr->iEndPos) << ">"
             << "\t" << IntToStr(itr->iSupportNum) << endl;

        ofsBrief << strChromName << " " // << IntToStr(itr->ucChromIndex) << " "
                 << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                 << IntToStr(itr->iSupportNum) << " "
                 //<< itr->strTag << " "
                 << itr->GetStrand() << " "
                 << "R" << endl;
    }

    ofs << endl << "Duplicated Num: " << iDuplicatNum << endl;

    ofs.close();
    ofsBrief.close();
}

//看看这个是在干嘛 --> 这里其实可以增加count, 从而进一步选取最优的  --> 因为都是premium的点，因此这里也需要用相应的similarity的处理方式
//这个逻辑有点奇怪 --> 我觉得应该需要更改 --> 应该是一个一个的加起来
bool ClsCircRNADetection::CheckDBGSeedHit(vector<St_ExonInfo>& vEI, vector<int>& vCount,
                                          string& strSeq, int iSeedPos,
                                          unordered_map<unsigned int, St_Node>& mpSDBG,
                                          En_HitSeed enHitSeed)
{
    string strSeed = strSeq.substr(iSeedPos, m_iKmerLen);
    //cout << strSeed << endl;
    unsigned int uiSeed = ConvertKmerToNum32(strSeed);
    unordered_map<unsigned int, St_Node>::iterator itrDBG = mpSDBG.find(uiSeed);

    if(itrDBG != mpSDBG.end())
    {
        vector<St_Link>::iterator itrLink = itrDBG->second.vEI.begin()->vLink.begin();
        //再用一个来进行定位
        string strDBGSeed = "";
        bool bFind = false;

        if(itrLink != itrDBG->second.vEI.begin()->vLink.end())
        {
            switch((int)enHitSeed)
            {
                case hsPre:
                    strDBGSeed = ConvertNum32ToKmer(itrLink->uiPreKmer, m_iKmerLen); //Check PreNode
                    if(strSeq.rfind(strDBGSeed, iSeedPos) != string::npos) //从后往前找 (rfind)
                        bFind = true;
                    break;
                case hsNext:
                    strDBGSeed = ConvertNum32ToKmer(itrLink->uiNextKmer, m_iKmerLen); //Check NextNode
                    if(strSeq.find(strDBGSeed, iSeedPos) != string::npos) //从前往后找 (find)
                        bFind = true;
                    break;
            }
        }

        if(bFind)
        {
            if(vEI.empty())
            {
                vEI.insert(vEI.end(), itrDBG->second.vEI.begin(), itrDBG->second.vEI.end());
                vCount.resize(itrDBG->second.vEI.size(), 1);
                //vEI.push_back(*itrEI);
            }
            else // 将新产生的EI加进去 --> 不应该是都加进去
            {
                //int iEI = vEI.size();
                for(vector<St_ExonInfo>::iterator itr = itrDBG->second.vEI.begin();
                    itr != itrDBG->second.vEI.end(); itr++)
                {
                    bool bFind = false;                    
                    for(int i=0; i<(int)vEI.size(); i++)
                    {
                        if(vEI[i].stExonPos == itr->stExonPos)
                        {
                            vCount[i]++;
                            bFind = true;
                            break;
                        }
                    }
                    if(!bFind)
                    {
                        vEI.push_back(*itr);
                        vCount.push_back(1);
                    }
                }
            }
        }
        return true;
    }
    else
        return false;

//    if(itrDBG != mpSDBG.end())
//    {
//        //Go
//        for(vector<St_ExonInfo>::iterator itrEI = itrDBG->second.vEI.begin();
//            itrEI != itrDBG->second.vEI.end(); itrEI++)
//        {
//            for(vector<St_Link>::iterator itrLink = itrEI->vLink.begin();
//                itrLink != itrEI->vLink.end(); itrLink++)
//            {
//                string strDBGSeed = "";
//                bool bFind = false;
//                switch(enHitSeed)
//                {
//                    case hsPre:
//                        strDBGSeed = ConvertNum32ToKmer(itrLink->uiPreKmer, m_iKmerLen); //Check PreNode
//                        if(strSeq.rfind(strDBGSeed, iSeedPos) != string::npos) //从后往前找 (rfind)
//                            bFind = true;
//                        break;
//                    case hsNext:
//                        strDBGSeed = ConvertNum32ToKmer(itrLink->uiNextKmer, m_iKmerLen); //Check NextNode
//                        if(strSeq.find(strDBGSeed, iSeedPos) != string::npos) //从前往后找 (find)
//                            bFind = true;
//                        break;
//                    default:
//                        break;
//                }
//                if(bFind)
//                {
////                    //------> For Debug
////                    St_ExonPos& stExonPos = itrEI->stExonPos;
////                    St_Raw_Exon& stExon = (*m_pvChrom)[stExonPos.ucChromIndex].vRG[stExonPos.uiGeneIndex].\
////                                               vRT[stExonPos.ucTranscriptIndex].vRExon[stExonPos.ucExonIndex];

////                    if(stExon.iEnd == 121116815)
////                    {
////                        int i =0;
////                        i++;
////                    }

////                    if(stExon.iEnd == 206567030)
////                    {
////                        int j =0;
////                        j++;
////                    }
////                    //<------

//                    if(vEI.empty())
//                    {
//                        vEI.insert(vEI.end(), itrDBG->second.vEI.begin(), itrDBG->second.vEI.end());
//                        //vEI.push_back(*itrEI);
//                    }
//                    else // 将新产生的EI加进去 --> 不应该是都加进去
//                    {
//                        int iEI = vEI.size();
//                        for(vector<St_ExonInfo>::iterator itr = itrDBG->second.vEI.begin();
//                            itr != itrDBG->second.vEI.end(); itr++)
//                        {
//                            bool bFind = false;
//                            for(int i=0; i<iEI; i++)
//                            {
//                                if(vEI[i].stExonPos == itr->stExonPos)
//                                {
//                                    bFind = true;
//                                    break;
//                                }
//                            }
//                            if(!bFind)
//                            {
//                                vEI.push_back(*itr);
//                            }
//                        }
////                        int iEI = vEI.size();
////                        bool bFind = false;
////                        for(int i=0; i<iEI; i++)
////                        {
////                            if(vEI[i].stExonPos == itrEI->stExonPos)
////                            {
////                                bFind = true;
////                                break;
////                            }
////                        }
////                        if(!bFind)
////                        {
////                            vEI.push_back(*itrEI);
////                        }
//                    }
//                    return true;
//                    //break;
//                }
//                else
//                {}
//                return false;
//            }
//        }
//        return true;
//    }
//    else
//        return false;
}

bool ClsCircRNADetection::CheckDuplicateInEI(vector<St_ExonInfo>& vEI)
{
    if(vEI.size() == 1)
        return false;

    vector<St_ExonInfo> vTmpExonInfo;

    for(vector<St_ExonInfo>::iterator itr = vEI.begin(); itr != vEI.end(); itr++)
    {
        St_ExonPos& stItrPos = itr->stExonPos;
        St_Raw_Exon& stItrExon = (*m_pvChrom)[stItrPos.ucChromIndex].vRG[stItrPos.uiGeneIndex].\
                                   vRT[stItrPos.ucTranscriptIndex].vRExon[stItrPos.ucExonIndex];

        for(vector<St_ExonInfo>::iterator subItr = vTmpExonInfo.begin();
            subItr != vTmpExonInfo.end(); subItr++)
        {
            if(itr->CheckIsSameTranscript(*subItr))
            {
                St_ExonPos& stSubItrPos = itr->stExonPos;
                St_Raw_Exon& stSubItrExon = (*m_pvChrom)[stSubItrPos.ucChromIndex].vRG[stSubItrPos.uiGeneIndex].\
                                           vRT[stSubItrPos.ucTranscriptIndex].vRExon[stSubItrPos.ucExonIndex];

                if(stItrExon.GetLength() == stSubItrExon.GetLength())
                    return true;
            }
        }
        vTmpExonInfo.push_back(*itr);
    }
    return false;
}

int ClsCircRNADetection::KeepMaxCount(vector<St_ExonInfo>& vEI, vector<int>& vCount)
{
    if(vEI.empty() || vCount.empty())
        return -1;

    int iMaxCount = 0;
    for(vector<int>::iterator itr = vCount.begin(); itr != vCount.end(); itr++)
    {
        if(iMaxCount < *itr)
            iMaxCount = *itr;
    }

    if(iMaxCount <= 1) // 需要有两个
    {
       vEI.clear();
       return -1;
    }

    int i = vEI.size() - 1;
    for(vector<St_ExonInfo>::iterator itr = vEI.end() - 1; itr >= vEI.begin(); itr--, i--)
    {
        if(vCount[i] != iMaxCount)
            vEI.erase(itr);
    }

    return iMaxCount;
}

void ClsCircRNADetection::CheckHittingCaseByDBG(string strSeq,
                                                unordered_map<unsigned int, St_Node>& mpSDBG,
                                                St_Fasta* pCurRefFa,
                                                bool bSeqRC, string strChromName,
                                                vector<St_Candidate>& vSelfCircCandi,
                                                vector<St_Candidate>& vRegCircCandi,
                                                ofstream& ofsDebug, ofstream& ofsDBGResult)
{
#ifdef _DEBUG
    bool bFindDebug = false;
    if(strSeq == "CTTCAAGTACTACATCCATGACCTATCTGACCTTATTGATGTCATGAAGACATATCACATGTACAATGCCGACAGCATCAGTGCTCAGAGCAAACTAAAG")
    {
        ofsDebug << "Find it --> RC" << endl;
        bFindDebug = true;
    }
    //else
    //    ofsDebug << "--- Sorry ---" << endl;
#endif    

    //1:首先完成的是DBG Filter
    int iLastKmerPos = m_iReadLen - m_iKmerLen;
    const int CHECKNUM = 4;
    ///前面CHECKNUM个 (从前往后)
    int ArryDBGPosCheckHead[CHECKNUM] = {1, 5,
                                         m_iKmerLen, 1.5*m_iKmerLen}; // {2, 5, 8}
//    const int CHECKNUM = 3;
//    ///前面三个 (从前往后)
//    int ArryDBGPosCheckHead[CHECKNUM] = {2, 5, 8};


    vector<St_ExonInfo> vEIHead;
    vector<int> vCount;
    //cout << ">>>> vEIHead" << endl;
    for(int i=0; i<CHECKNUM; i++)
    {
        CheckDBGSeedHit(vEIHead, vCount, strSeq, ArryDBGPosCheckHead[i], mpSDBG, hsNext);
    }
    int iMaxCount = KeepMaxCount(vEIHead, vCount);
    if(iMaxCount >= 2)
    {
        if(CheckDuplicateInEI(vEIHead))
            return;
    }

    ///后面CHECKNUM个 (从后往前)
    int ArryDBGPosCheckTail[CHECKNUM] = {iLastKmerPos - 1, iLastKmerPos - 5,
                                         iLastKmerPos - m_iKmerLen, iLastKmerPos - 1.5*m_iKmerLen};
//    int ArryDBGPosCheckTail[CHECKNUM] = {iLastKmerPos - 2, iLastKmerPos - 5,
//                                         iLastKmerPos - 8};

    vector<St_ExonInfo> vEITail;
    vCount.clear();
    //cout << ">>>> vEITail" << endl;
    for(int i=0; i<CHECKNUM; i++)
    {
        CheckDBGSeedHit(vEITail, vCount, strSeq, ArryDBGPosCheckTail[i], mpSDBG, hsPre);
    }
    iMaxCount = KeepMaxCount(vEITail, vCount);
    if(iMaxCount >= 2)
    {
        if(CheckDuplicateInEI(vEITail))
            return;
    }

    if(vEIHead.empty() || vEITail.empty())
        return;

    //--> Check more two --> Just Need to know if they existed  or not -->    
    unordered_map<unsigned int, St_Node>::iterator itrAdd1;
    unordered_map<unsigned int, St_Node>::iterator itrAdd2;

    ///Get the suitable DBG Sampling Position
    float ArryPosRatio[2] = {.15, .85}; //分别取15% 位置 和 85%位置的点 // 管他tail 或者 head ---> 只是用于之后的筛选
                                       //这里取值的原因是因为我们需要尽可能的是这样的点不掉在breakpoint的地方


    ///For itrAdd1
    string strSeed = strSeq.substr(iLastKmerPos*ArryPosRatio[0], m_iKmerLen);
    unsigned int uiSeed = ConvertKmerToNum32(strSeed);
    itrAdd1 = mpSDBG.find(uiSeed);
    ///For itrAdd2
    strSeed = strSeq.substr(iLastKmerPos*ArryPosRatio[1], m_iKmerLen);
    uiSeed = ConvertKmerToNum32(strSeed);
    itrAdd2 = mpSDBG.find(uiSeed);
    //<--

    //有料，那么我们就看看那结果如何  --> Map 这里不靠谱 ---> 不要用复杂结构体作key
    ///用于 collect self curcular case
    vector<St_CandiExAssistSelf> vCandiAssistSelf;
    ///用于 collect regular circular Case
    vector<St_CandiExAssistReg> vCandiAssistReg;

    for(vector<St_ExonInfo>::iterator itrEIHead = vEIHead.begin(); itrEIHead != vEIHead.end(); itrEIHead++)
    {
        if(itrEIHead->vLink.empty())
            continue;

//        //---------->
//        St_ExonPos& stExonPosDonor = itrEIHead->stExonPos;
//        St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
//                                   vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];
//        if(stExonDonor.iEnd == 121116815)
//        {
//            int i =0;
//            i++;
//        }

//        if(stExonDonor.iEnd == 206567030)
//        {
//            int j =0;
//            j++;
//        }
//        //<----

        for(vector<St_Link>::iterator itrEIHeadLink = itrEIHead->vLink.begin();
            itrEIHeadLink != itrEIHead->vLink.end(); itrEIHeadLink++)
        {
            bool bHeadSeqHitEnd = false;
            if(itrEIHeadLink->cCurTag == 'T' || itrEIHeadLink->cCurTag == 'E') //head link需要碰到的是尾巴--> 这个是针对self case来说的
            {  //do the following logic
               bHeadSeqHitEnd = true;
            }
            else
            {}

            //Compare with Tail                        
            for(vector<St_ExonInfo>::iterator itrEITail = vEITail.begin(); itrEITail != vEITail.end(); itrEITail++)
            {
                if(itrEITail->vLink.empty())
                    continue;

                bool bProbSelfCase = false;
                bool bProbRefCase = false;
                //这里还真要做修改--> 我们现在的代码不能handle一个reads同时支持self circular以及regular circular
                for(vector<St_Link>::iterator itrEITailLink = itrEITail->vLink.begin();
                    itrEITailLink != itrEITail->vLink.end(); itrEITailLink++)
                {                    
                    if(itrEITailLink->cCurTag == 'H' || itrEITailLink->cCurTag == 'S')  // 我们这里需要的是premium node！！！
                    {
                        //This is the case for self circular roughly candidate collection
                        if(itrEIHead->stExonPos == itrEITail->stExonPos)
                        {
                            bProbSelfCase = true;
                            //Do by the packaged function
                            CollectCandiAssistSelf(itrEIHead, vCandiAssistSelf, strChromName);  // 一次就够了!!!
                            break;
                        }
                        else // the regular circ case belongs to the exon with long length
                        {
                            bProbRefCase = true;
                            CollectCandiAssistReg(itrEIHead, itrEITail, vCandiAssistReg, strChromName);
                            break;
                        }
                    }                   
                }

                //我们在这里真的应该考虑所有的link 组合，有点 push 啊!
                if(!bProbSelfCase && !bProbRefCase) //其他的case都应该是潜在的regular circular的 case,因此这里我们需要好好的进行check
                {
                    //相当于都没搞的时候，我们看看是怎么回事儿！
                    //这里我们要check长度
                    St_ExonPos& stExonPosDonor = itrEIHead->stExonPos;
                    St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
                                               vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];

                    //Step 2: Set Acceptor
                    St_ExonPos& stExonPosAcceptor = itrEITail->stExonPos;
                    St_Raw_Exon& stExonAcceptor = (*m_pvChrom)[stExonPosAcceptor.ucChromIndex].vRG[stExonPosAcceptor.uiGeneIndex].\
                                               vRT[stExonPosAcceptor.ucTranscriptIndex].vRExon[stExonPosAcceptor.ucExonIndex];

                    // 当exon过短，导致reads的cover，cover到了整个exon的时候，这种情况也是需要考虑的
                    // 这个情况就是下面的两个逻辑的组合
                    if( ((itrEIHeadLink->cCurTag == 'H' || itrEIHeadLink->cCurTag == 'S') &&
                         stExonDonor.GetLength() <= (m_iReadLen - m_iKmerLen + 1)*2) ||
                        ((itrEITail->vLink.begin()->cCurTag == 'T' || itrEITail->vLink.begin()->cCurTag == 'E') &&
                         stExonAcceptor.GetLength() <= (m_iReadLen - m_iKmerLen + 1)*2) )
                    {
                        if(stExonDonor.GetIsCRNATailExon() && stExonAcceptor.GetIsCRNAHeadExon())
                        {
                            CollectCandiAssistReg(itrEIHead, itrEITail, vCandiAssistReg, strChromName);
                        }
                    }
                }
            }
        }
    }

#ifdef _DEBUG
    if(bFindDebug)
    {
        ofsDebug << "**********************vCandiAssistReg size: " << IntToStr(vCandiAssistReg.size()) << endl;
        for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.begin();
            itr != vCandiAssistReg.end(); itr++)
        {
            St_ExonPos& stExonPosDonor = itr->stExonPosDonor;
            St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
                                       vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];

            St_ExonPos& stExonPosAcceptor = itr->stExonPosAcceptor;
            St_Raw_Exon& stExonAcceptor = (*m_pvChrom)[stExonPosAcceptor.ucChromIndex].vRG[stExonPosAcceptor.uiGeneIndex].\
                                       vRT[stExonPosAcceptor.ucTranscriptIndex].vRExon[stExonPosAcceptor.ucExonIndex];

            ofsDebug << "ExonDonor &  stExonPosAcceptor: "
                 << "(" << IntToStr(stExonDonor.iStart) << ", " << IntToStr(stExonDonor.iEnd) << ")"
                 << ", "
                 << "(" << IntToStr(stExonAcceptor.iStart) << ", " << IntToStr(stExonAcceptor.iEnd) << ")"
                 << endl;
            ofsDebug << "iCount: " << IntToStr(itr->iCount) << endl;
            ofsDebug << "iAddCount: " << IntToStr(itr->iAddNodHitNum) << endl << endl;
        }
    }
#endif

//    //--> For debug
//    // Check if those case include the specific raw candidate
//    for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.begin(); itr != vCandiAssistReg.end(); itr++)
//    {
//        St_ExonPos& stExonPosDonor = itr->stExonPosDonor;
//        St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
//                                   vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];

//        St_ExonPos& stExonPosAcceptor = itr->stExonPosAcceptor;
//        St_Raw_Exon& stExonAcceptor = (*m_pvChrom)[stExonPosAcceptor.ucChromIndex].vRG[stExonPosAcceptor.uiGeneIndex].\
//                                   vRT[stExonPosAcceptor.ucTranscriptIndex].vRExon[stExonPosAcceptor.ucExonIndex];
//        if(stExonDonor.iEnd == 33099711 && stExonAcceptor.iStart == 33097428)
//        {
//            cout << strSeq << endl;
//            break;
//        }
//    }
//    //<--------

    //这里看看 --> 是不是连续的重复
    if(vCandiAssistReg.size() == 1 && vCandiAssistSelf.size() == 1)
    {
        if(vCandiAssistReg[0].IsNeighbor())
        {
            if(vCandiAssistReg[0].stCandi.GetDonorLength() == vCandiAssistReg[0].stCandi.GetAcceptorLength())
                return;
        }
    }

    if(vCandiAssistSelf.size() == 1)
    {
        //Check if it is as same as its neighbor
        //Check its previous:
        St_ExonPos& stCurExonPos = vCandiAssistSelf[0].stExonPos;
        St_Raw_Exon& stCurExon = (*m_pvChrom)[stCurExonPos.ucChromIndex].vRG[stCurExonPos.uiGeneIndex].\
                                             vRT[stCurExonPos.ucTranscriptIndex].vRExon[stCurExonPos.ucExonIndex];
        int iMaxExonIndex = (*m_pvChrom)[stCurExonPos.ucChromIndex].vRG[stCurExonPos.uiGeneIndex].\
                vRT[stCurExonPos.ucTranscriptIndex].vRExon.size() - 1;

        if(stCurExonPos.ucExonIndex > 0)
        {
            St_Raw_Exon& stExonPre = (*m_pvChrom)[stCurExonPos.ucChromIndex].vRG[stCurExonPos.uiGeneIndex].\
                                              vRT[stCurExonPos.ucTranscriptIndex].vRExon[stCurExonPos.ucExonIndex - 1];
            if(stCurExon.GetLength() == stExonPre.GetLength())
                return;
        }
        //Check it later:
        if(stCurExonPos.ucExonIndex < iMaxExonIndex)
        {
            St_Raw_Exon& stExonNext = (*m_pvChrom)[stCurExonPos.ucChromIndex].vRG[stCurExonPos.uiGeneIndex].\
                                              vRT[stCurExonPos.ucTranscriptIndex].vRExon[stCurExonPos.ucExonIndex + 1];
            if(stCurExon.GetLength() == stExonNext.GetLength())
                return;
        }
    }

    PurifySelfCase(vCandiAssistSelf, mpSDBG, pCurRefFa,
                   strSeq, bSeqRC, itrAdd1, itrAdd2, vSelfCircCandi, ofsDBGResult);
    PurifyRegCase(vCandiAssistReg, mpSDBG, pCurRefFa,
                  strSeq, bSeqRC, itrAdd1, itrAdd2, vRegCircCandi, ofsDBGResult);

//    /* 这里我们首先着眼于解决regular circle 的 case
//     * 也就是所谓的三次握手，搞起！！！因为之前数据结构铺垫的好，这个逻辑结构实现起来不是难事儿
//     * 因为这个不是
//     */
//    vector<St_HitNode> vHitNode; // 这里我们并不合并重复的点
//    for(unsigned int i = 0; i < strSeq.length() - m_iKmerLen + 1; i++)
//    {
//        unsigned int uiKmer = ConvertKmerToNum32(strSeq.substr(i, m_iKmerLen));
//        if(mpSDBG.find(uiKmer) != mpSDBG.end()) // 重复的没有关系，直接放进去不用担心
//        {
//            //第一次通信
//            St_Node& stNode = mpSDBG[uiKmer];
//            for(vector<St_ExonInfo>::iterator itrEI = stNode.vEI.begin(); itrEI != stNode.vEI.end(); itrEI++)
//            {
//                bool bFind = false;
//                unsigned int uiCur = string::npos;
//                unsigned int uiPre = string::npos;
//                for(vector<St_Link>::iterator itrLink = itrEI->vLink.begin();
//                    itrLink != itrEI->vLink.end(); itrLink++)
//                {
//                    if(strSeq.find(itrLink->uiNextKmer, i+1) != string::npos)
//                    {
//                        bFind = true;
//                        uiCur = itrLink->uiNextKmer;
//                        uiPre = uiKmer;
//                        break;
//                    }
//                }
//                if(bFind) // 看这个exon
//                {
//                    //第二次通信完成
//                    //开始第三次 寻找end  --> 通过看tag来知道是不是第一个节点还是最后一个！！ 明天搞起
//                    int iMaxLen = m_iReadLen;
//                    int iOffSet = 0;
//                    while(uiCur != && iOffSet < iMaxLen)
//                    {

//                    }

//                }
//            }
//        }
//    }
}

void ClsCircRNADetection::PurifySelfCase(vector<St_CandiExAssistSelf>& vCandiAssistSelf,
                                         unordered_map<unsigned int, St_Node>& mpSDBG,
                                         St_Fasta* pCurRefFa,
                                         string& strSeq, bool bSeqRC,
                                         unordered_map<unsigned int, St_Node>::iterator& itrAdd1,
                                         unordered_map<unsigned int, St_Node>::iterator& itrAdd2,
                                         vector<St_Candidate>& vSelfCircCandi,
                                         ofstream& ofsDBGResult)
{
    //get max hit
    if(vCandiAssistSelf.empty())
    {}
    else
    {
        if(vCandiAssistSelf.size() == 1)
        {}
        else
        {
            //Update the additional Node count
            for(vector<St_CandiExAssistSelf>::iterator itr = vCandiAssistSelf.begin(); itr != vCandiAssistSelf.end(); itr++)
            {
                if(itrAdd1 != mpSDBG.end())  // 有意义才做这个事儿
                {
                    itr->iAddNodHitNum = itr->iAddNodHitNum + CheckIfExonPosInEI(itr->stExonPos, itrAdd1->second.vEI);
                }
                if(itrAdd2 != mpSDBG.end()) // 有意义才做这个事儿
                {
                    itr->iAddNodHitNum = itr->iAddNodHitNum + CheckIfExonPosInEI(itr->stExonPos, itrAdd2->second.vEI);
                }
            }

            //找到最大的count
            int iMaxCount = 0;
            for(vector<St_CandiExAssistSelf>::iterator itr = vCandiAssistSelf.begin(); itr != vCandiAssistSelf.end(); itr++)
            {
                if(iMaxCount < itr->iCount)  //选择 majority
                {
                    iMaxCount = itr->iCount;
                }
            }

            //再找 Add Number Hit Number 最大的
            int iMaxAddNum = 0;
            for(vector<St_CandiExAssistSelf>::iterator itr = vCandiAssistSelf.begin(); itr != vCandiAssistSelf.end(); itr++)
            {
                if(itr->iCount == iMaxCount)  //选择 majority
                {
                     if(iMaxAddNum < itr->iAddNodHitNum)
                     {
                        iMaxAddNum =  itr->iAddNodHitNum;
                     }
                }
            }

            //Collect The Final Case
            for(vector<St_CandiExAssistSelf>::iterator itr = vCandiAssistSelf.end() - 1; itr >= vCandiAssistSelf.begin(); itr--)
            {
                if(itr->iCount == iMaxCount && itr->iAddNodHitNum == iMaxAddNum)
                {}
                else
                    vCandiAssistSelf.erase(itr);
            }
        }      

        //<---------
        if(vCandiAssistSelf.size() > 2)
            return;

        if( vCandiAssistSelf.size() == 2 &&
            vCandiAssistSelf[0].iCount == vCandiAssistSelf[1].iCount &&
            vCandiAssistSelf[0].iAddNodHitNum == vCandiAssistSelf[1].iAddNodHitNum &&
            vCandiAssistSelf[0].GetLength() == vCandiAssistSelf[1].GetLength())
        {
            return;
        }

        for(vector<St_CandiExAssistSelf>::iterator itr = vCandiAssistSelf.begin(); itr != vCandiAssistSelf.end(); itr++)
        {
            //Check this Exon!!-->
            if(CheckCurSelfCaseByDBG(itr->stCandi, strSeq, itr->stExonPos, mpSDBG, vSelfCircCandi))
            {
#ifdef CIRC_CANDI_PRINT
                St_Candidate& stCurCandi = itr->stCandi;
                //Get two exon: donor Exon and accepter Exon, then print the start and ending position out
                St_ExonPos& stExonPos = itr->stExonPos;
                St_Raw_Exon& stExon = (*m_pvChrom)[stExonPos.ucChromIndex].vRG[stExonPos.uiGeneIndex].\
                                           vRT[stExonPos.ucTranscriptIndex].vRExon[stExonPos.ucExonIndex];

                ofsDBGResult << "Donor(Start, End) & Acceptor(Start, End): "
                     << "(" << IntToStr(stExon.iStart) << ", " << IntToStr(stExon.iEnd) << ")"
                     << " -- "
                     << "(" << IntToStr(stExon.iStart) << ", " << IntToStr(stExon.iEnd) << ")" << endl;

                ofsDBGResult << "Candi Start & End: "<< IntToStr(stCurCandi.iStartPos) << " : "
                     << IntToStr(stCurCandi.iEndPos) << endl;

                ofsDBGResult << "RC: " << (stCurCandi.bRC ? "Yes" : "No") << endl;

                ofsDBGResult << "Chrom: " << stCurCandi.strChromName << endl;

                ofsDBGResult << "Current Reads: " << (bSeqRC ? "RC" : "Reg") << endl;

                ofsDBGResult << "CircType: " << IntToStr((int)stCurCandi.enType) << endl;

                ofsDBGResult << strSeq << endl;

                ofsDBGResult <<  "***" << endl;
//                cout << "Start & End: " << IntToStr(itr->stCandi.iStartPos) << " : "
//                     << IntToStr(itr->stCandi.iEndPos) << endl;

//                cout << "RC: " << (itr->stCandi.bRC ? "Yes" : "No") << endl;

//                cout << "Current Reads: " << (bSeqRC ? "RC" : "Reg") << endl
//                     << strSeq << endl << endl;
#endif
            }
        }
    }
}

void ClsCircRNADetection::PurifyRegCase(vector<St_CandiExAssistReg>& vCandiAssistReg,
                                        unordered_map<unsigned int, St_Node>& mpSDBG,
                                        St_Fasta* pCurRefFa,
                                        string& strSeq, bool bSeqRC,
                                        unordered_map<unsigned int, St_Node>::iterator& itrAdd1,
                                        unordered_map<unsigned int, St_Node>::iterator& itrAdd2,
                                        vector<St_Candidate>& vRegCircCandi,
                                        ofstream& ofsDBGResult)
{
#ifdef DEBUG
    bool bDebugTarget = false;
    if(strSeq == "CTTCAAGTACTACATCCATGACCTATCTGACCTTATTGATGTCATGAAGACATATCACATGTACAATGCCGACAGCATCAGTGCTCAGAGCAAACTAAAG")
    {
        cout << "PurifyRegCase vCandiAssistReg size: " << IntToStr(vCandiAssistReg.size()) << endl;
        bDebugTarget = true;
    }
#endif

    if(vCandiAssistReg.empty())
        return;
    else if(vCandiAssistReg.size() == 1)
    {
//        St_ExonPos& stExonPosDonor = vCandiAssistReg[0].stExonPosDonor;
//        St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
//                                   vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];
//        int i =0;

    }
    else // 对于大于一个candidate --> 我们需要过滤，实际上我们最理想的情况是一个reads只对应一个潜在的candidate
    {
        //Update the additional Node count
        for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.begin(); itr != vCandiAssistReg.end(); itr++)
        {
//            St_ExonPos& stExonPosDonor = itr->stExonPosDonor;
//            St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
//                                       vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];

            if(itrAdd1 != mpSDBG.end())  // 有意义才做这个事儿
            {
                int iTopUp = CheckIfExonPosInEI(itr->stExonPosDonor, itrAdd1->second.vEI);
                if(iTopUp == 0)
                    iTopUp = CheckIfExonPosInEI(itr->stExonPosAcceptor, itrAdd1->second.vEI);

                itr->iAddNodHitNum += iTopUp;
            }
            if(itrAdd2 != mpSDBG.end()) // 有意义才做这个事儿
            {
                int iTopUp = CheckIfExonPosInEI(itr->stExonPosDonor, itrAdd2->second.vEI);
                if(iTopUp == 0)
                    iTopUp = CheckIfExonPosInEI(itr->stExonPosAcceptor, itrAdd2->second.vEI);

                itr->iAddNodHitNum += iTopUp;
            }
        }

#ifdef DEBUG
        if(strSeq == "CTTCAAGTACTACATCCATGACCTATCTGACCTTATTGATGTCATGAAGACATATCACATGTACAATGCCGACAGCATCAGTGCTCAGAGCAAACTAAAG")
        {
            cout << "PurifyRegCase vCandiAssistReg size: " << IntToStr(vCandiAssistReg.size())
                 << endl
                 << IntToStr(vCandiAssistReg[0].stCandi.iStartPos) << ", " << IntToStr(vCandiAssistReg[1].stCandi.iStartPos)
                 << endl
                 << IntToStr(vCandiAssistReg[0].stCandi.iEndPos) << ", " << IntToStr(vCandiAssistReg[1].stCandi.iEndPos)
                 << endl
                 << IntToStr(vCandiAssistReg[0].iAddNodHitNum) << endl
                 << IntToStr(vCandiAssistReg[1].iAddNodHitNum) << endl
                 << IntToStr(abs(vCandiAssistReg[0].iCount - vCandiAssistReg[1].iCount)) << endl;
        }
#endif

        if( vCandiAssistReg.size() == 2 &&
            (vCandiAssistReg[0].stCandi.iStartPos != vCandiAssistReg[1].stCandi.iStartPos ||
             vCandiAssistReg[0].stCandi.iEndPos != vCandiAssistReg[1].stCandi.iEndPos) &&
            vCandiAssistReg[0].iAddNodHitNum == 2 &&
            vCandiAssistReg[1].iAddNodHitNum == 2 &&
            abs(vCandiAssistReg[0].iCount - vCandiAssistReg[1].iCount) <= 2 )
        {} // 不做过滤
        else
        {
            //找到最大的count
            int iMaxCount = 0;
            for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.begin(); itr != vCandiAssistReg.end(); itr++)
            {
                if(iMaxCount < itr->iCount)  //选择 majority
                {
                    iMaxCount = itr->iCount;
                }
            }

            //再找 Add Number Hit Number 最大的
            int iMaxAddNum = 0;
            for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.begin(); itr != vCandiAssistReg.end(); itr++)
            {
                if(itr->iCount == iMaxCount)  //选择 majority
                {
                     if(iMaxAddNum < itr->iAddNodHitNum)
                     {
                        iMaxAddNum =  itr->iAddNodHitNum;
                     }
                }
            }

            //Collect The Final Case
            for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.end() - 1; itr >= vCandiAssistReg.begin(); itr--)
            {
                if(itr->iCount == iMaxCount && itr->iAddNodHitNum == iMaxAddNum)
                {}
                else
                    vCandiAssistReg.erase(itr);
            }
        }
    }    

    if(vCandiAssistReg.size() > 2) //支持太多了，反而不真实了，所以应该去除掉
        return;

    for(vector<St_CandiExAssistReg>::iterator itr = vCandiAssistReg.begin(); itr != vCandiAssistReg.end(); itr++)
    {
        //Check if they are the same sequence --> go
        if(itr->stCandi.GetDonorLength() == itr->stCandi.GetAcceptorLength())
        {
            St_Raw_Exon* pDonor = itr->stCandi.pExon;
            St_Raw_Exon* pAcceptor = itr->stCandi.pExonAcceptor;
            if(pCurRefFa->strSeq.substr(pDonor->iStart - 1, abs(pDonor->iEnd - pDonor->iStart) + 1) ==
               pCurRefFa->strSeq.substr(pAcceptor->iStart - 1, abs(pAcceptor->iEnd - pAcceptor->iStart) + 1))
            {
                return;
            }
        }
        //<--

        //Check this Exon!!-->
        St_Candidate stCurCandi;
        if(CheckCurRegCaseByDBG(itr->stCandi, strSeq, itr->stExonPosDonor, itr->stExonPosAcceptor, mpSDBG,
                                stCurCandi, vRegCircCandi))
        {
#ifdef CIRC_CANDI_PRINT
            //Get two exon: donor Exon and accepter Exon, then print the start and ending position out
            St_ExonPos& stExonPosDonor = itr->stExonPosDonor;
            St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
                                       vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];

            //Step 2: Set Acceptor
            St_ExonPos& stExonPosAcceptor = itr->stExonPosAcceptor;
            St_Raw_Exon& stExonAcceptor = (*m_pvChrom)[stExonPosAcceptor.ucChromIndex].vRG[stExonPosAcceptor.uiGeneIndex].\
                                       vRT[stExonPosAcceptor.ucTranscriptIndex].vRExon[stExonPosAcceptor.ucExonIndex];

            ofsDBGResult << "Donor(Start, End) & Acceptor(Start, End): "
                 << "(" << IntToStr(stExonDonor.iStart) << ", " << IntToStr(stExonDonor.iEnd) << ")"
                 << " -- "
                 << "(" << IntToStr(stExonAcceptor.iStart) << ", " << IntToStr(stExonAcceptor.iEnd) << ")" << endl;

            ofsDBGResult << "Candi Start & End: "<< IntToStr(stCurCandi.iStartPos) << " : "
                 << IntToStr(stCurCandi.iEndPos) << endl;

            ofsDBGResult << "RC: " << (stCurCandi.bRC ? "Yes" : "No") << endl;

            ofsDBGResult << "Chrom: " << stCurCandi.strChromName << endl;

            ofsDBGResult << "Current Reads: " << (bSeqRC ? "RC" : "Reg") << endl;

            ofsDBGResult << "CircType: " << IntToStr((int)stCurCandi.enType) << endl;

            ofsDBGResult << strSeq << endl;

            ofsDBGResult <<  "***" << endl;
#endif
        }
    }
}

void ClsCircRNADetection::CollectCandiAssistSelf(vector<St_ExonInfo>::iterator itrEI,
                                                 vector<St_CandiExAssistSelf>& vCandiAssistSelf,
                                                 string strChromName)
{
    //Get it
    St_ExonPos& stExonPos = itrEI->stExonPos;
    St_Raw_Exon& stExon = (*m_pvChrom)[stExonPos.ucChromIndex].vRG[stExonPos.uiGeneIndex].\
                          vRT[stExonPos.ucTranscriptIndex].vRExon[stExonPos.ucExonIndex];
    St_Candidate stCurCandi(strChromName, stExon.iEnd, stExon.iStart, 1,
                            stExon.bRC, ctSelf, "xx", &stExon);
    bool bFind = false;
    for(vector<St_CandiExAssistSelf>::iterator itrCA = vCandiAssistSelf.begin();
        itrCA != vCandiAssistSelf.end(); itrCA++) // CA: Candi Assistant
    {
        if(itrCA->stCandi == stCurCandi)
        {
            itrCA->iCount++;
            bFind= true;
            break;
        }
    }
    if(!bFind)
    {
        St_CandiExAssistSelf tmpCandiAssist(stCurCandi, stExonPos, 1);
        vCandiAssistSelf.push_back(tmpCandiAssist);
    }
//    if(stTempCandi.enType == ctMax)
//        stTempCandi = stCurCandi;
//    else
//    {
//        //取最短的
//        if(stTempCandi.GetLength() < stCurCandi.GetLength())
//            stTempCandi = stCurCandi;
//    }
}

void ClsCircRNADetection::CollectCandiAssistReg(vector<St_ExonInfo>::iterator itrEIHead,
                                                vector<St_ExonInfo>::iterator itrEITail,
                                                vector<St_CandiExAssistReg>& vCandiAssistReg,
                                                string strChromName)
{
    /*
     * here we should know:
     * EIHead --> 对应的应该是 Donor
     * EITail --> 对应的应该是 Acceptor
     *
     * 这里的Head是说它对应的reads heead的片段
     * 这里的Tail是说它对应的是的 reads tail的片段
     */
    //这里主要是看看他们是不是属于同一个transcript,并且是不是一后一前 -->
    if(itrEIHead->CheckIsBackToBegin(*itrEITail))
    {
        //如果这个是真的，那么我们就干一个事儿
        //Step 1: Set Donor
        St_ExonPos& stExonPosDonor = itrEIHead->stExonPos;
        St_Raw_Exon& stExonDonor = (*m_pvChrom)[stExonPosDonor.ucChromIndex].vRG[stExonPosDonor.uiGeneIndex].\
                                   vRT[stExonPosDonor.ucTranscriptIndex].vRExon[stExonPosDonor.ucExonIndex];
        int iPosDonor = stExonDonor.iEnd;

        //Step 2: Set Acceptor
        St_ExonPos& stExonPosAcceptor = itrEITail->stExonPos;
        St_Raw_Exon& stExonAcceptor = (*m_pvChrom)[stExonPosAcceptor.ucChromIndex].vRG[stExonPosAcceptor.uiGeneIndex].\
                                   vRT[stExonPosAcceptor.ucTranscriptIndex].vRExon[stExonPosAcceptor.ucExonIndex];
        int iPosAcceptor = stExonAcceptor.iStart;

        //Step 3: Initial A new Candidate
        St_Candidate stCurCandi(strChromName, iPosDonor, iPosAcceptor, 1,
                                stExonPosDonor.bRC, ctRegular, "xx", &stExonDonor, &stExonAcceptor);

        //Merge this Candidate with the container
        bool bFind = false;
        for(vector<St_CandiExAssistReg>::iterator itrCA = vCandiAssistReg.begin();
            itrCA != vCandiAssistReg.end(); itrCA++) // CA: Candi Assistant
        {
            if(itrCA->stCandi == stCurCandi)
            {
                itrCA->iCount++;
                bFind= true;
                break;
            }
        }
        if(!bFind)
        {
            St_CandiExAssistReg tmpCandiAssist(stCurCandi, stExonPosDonor, stExonPosAcceptor, 1);
            vCandiAssistReg.push_back(tmpCandiAssist);
        }
    }
}

int ClsCircRNADetection::CheckIfExonPosInEI(St_ExonPos& stExonPos, vector<St_ExonInfo>& vEI)
{
    bool bFind = false;
    for(vector<St_ExonInfo>::iterator itr = vEI.begin(); itr != vEI.end(); itr++)
    {
        if(itr->stExonPos == stExonPos)
        {
            bFind = true;
            break;
        }
    }
    if(bFind)
        return 1;
    else
        return 0;
}

bool ClsCircRNADetection::CheckCurSelfCaseByDBG(St_Candidate& stProbCandi,
                                                string& strSeq, St_ExonPos& stProbExonPos,
                                                unordered_map<unsigned int, St_Node>& mpSDBG,
                                                vector<St_Candidate>& vSelfCircCandi)
{
    /* 我们的策略如下:
     * (1) 从sequence的两端入手
     *     a) 左端: 从左到右，使用next进行: 边走边走快 --> 到最后Premiun的点的时候停下来  --> 记录两个数组，然后进行整合
     *        *有个要注意，如果多个没有搞定，那么我们就要后递增的走知道找到相应的好的点
     *     b) 右端: 从右到左，使用pre进行:
     */

    //from start to end --> 我们在这里期待的是: 能够有一半的到到达相应的尾部
    //--> 这三个变量之后会用到
    unordered_map<unsigned int, St_Node>::iterator itrDBGSeed;
    St_Link stCurLink;
    unsigned int uiSeed = 0;
    int iOffSetTail = 0;
    //<--

    for(int i = 0; i < (int)strSeq.length() - m_iKmerLen + 1; i++)
    {
        string strKmer = strSeq.substr(i, m_iKmerLen);
        uiSeed = ConvertKmerToNum32(strKmer);
        itrDBGSeed = mpSDBG.find(uiSeed);
        if(itrDBGSeed != mpSDBG.end()) //证明找到了
        {
            bool bFindProbExon = false;
            for(vector<St_ExonInfo>::iterator itr = itrDBGSeed->second.vEI.begin();
                itr != itrDBGSeed->second.vEI.end(); itr++)
            {
                if(itr->stExonPos == stProbExonPos)
                {
                    bFindProbExon = true;
                    stCurLink = (*itr->vLink.begin());
                    break;
                }
            }
            if(bFindProbExon)
                break;
            else
                iOffSetTail++;
        }
        else
            iOffSetTail++;
    }

    if(itrDBGSeed == mpSDBG.end())
        return false;

    //--> 然后我们来遍历
    ///先从左边开始遍历 ---> 目的是tail
    /// 干脆将节点都拿出来 -->

    // 2个条件
    // (1) 到了tail处
    // (2) walk的距离大于一个reads的长度
    ///To Tail Part

    vector<unsigned int> vKmerNodeToTail;
    vector<char> vTagTail;
    int iPathFullLen = 0;
    int iResultCodeTail = CheckDBGTail(stCurLink, iPathFullLen,
                                       vKmerNodeToTail, vTagTail, stProbExonPos, uiSeed, mpSDBG);

    if(iPathFullLen < CSDBGFULLLENMIN)
        return false;

    ///From Begin Part --> Go! -->  Going to head part ************************
    int iOffSetHead = 0;
    for(int i = strSeq.length() - m_iKmerLen; i >= 0 ; i--)
    {
        string strKmer = strSeq.substr(i, m_iKmerLen);
        uiSeed = ConvertKmerToNum32(strKmer);
        itrDBGSeed = mpSDBG.find(uiSeed);
        if(itrDBGSeed != mpSDBG.end()) //证明找到了
        {
            bool bFindProbExon = false;
            for(vector<St_ExonInfo>::iterator itr = itrDBGSeed->second.vEI.begin();
                itr != itrDBGSeed->second.vEI.end(); itr++)
            {
                if(itr->stExonPos == stProbExonPos)
                {
                    bFindProbExon = true;
                    stCurLink = (*itr->vLink.begin());
                    break;
                }
            }
            if(bFindProbExon)
                break;
            else
                iOffSetHead++;
        }
        else
            iOffSetHead++;
    }

    if(itrDBGSeed == mpSDBG.end())
        return false;

    vector<unsigned int> vKmerNodeToHead;
    vector<char> vTagHead;
    iPathFullLen = 0;
    int iResultCodeHead = CheckDBGHead(stCurLink, iPathFullLen,
                                       vKmerNodeToHead, vTagHead, stProbExonPos, uiSeed, mpSDBG);

    if(iPathFullLen < CSDBGFULLLENMIN)
        return false;

    //*************************For Canclualte******************
    ///首先check个数，如果个数少的出奇，那么就是不正确的 (我们可以在最后的discussion进行讨论)
    unsigned int iMinNodeNum = (m_iReadLen - m_iKmerLen + 1) * .8 *.3; //注意我们这里是跳着取的 --> 这就是为什么要乘以0.3    
    if(vKmerNodeToTail.size() + vKmerNodeToHead.size() < iMinNodeNum)
    {
        //Check Premiun Node Number
        int iPremiunTailNum = 0;
        for(vector<char>::iterator itr = vTagTail.begin(); itr != vTagTail.end(); itr++)
        {
            if(*itr == 'T')
                iPremiunTailNum++;
        }

        int iPremiunHeadNum = 0;
        for(vector<char>::iterator itr = vTagHead.begin(); itr != vTagHead.end(); itr++)
        {
            if(*itr == 'H')
                iPremiunHeadNum++;
        }

        if(iPremiunTailNum >= 2 && iPremiunHeadNum >=2) // 我们认为能够成功hit到足够premiun的 (抓阄必会抓到一个，然后就是tail或head必须是一个)
        {}
        else
            return false;
    }

    ///For the case: To Tail
    bool bGoodToTail = false;
    int iHitNumber = 0;
    bool bGetTagE = false; //E: end
    bool bGetTagT = false; //T: tail
    //int iPremiunTailNum = 0;
    if(0 == iResultCodeTail) // 0: 说明已经遍历到了最后的tail节点
    {
        int i = 0;
        bool bDoLooseMatch = true;
        for(vector<unsigned int>::iterator itr = vKmerNodeToTail.begin(); itr != vKmerNodeToTail.end(); itr++, i++)
        {
            string strSeed = ConvertNum32ToKmer(*itr, m_iKmerLen);
            int iHitPos = strSeq.find(strSeed, iOffSetTail);
            if(iHitPos == -1 && vTagTail[i] == 'T' && bDoLooseMatch)
            {
                int iPos = iOffSetTail + 3;
                if(iPos + m_iKmerLen <= (int)strSeq.length()) // 这是往后所以 以length为定界
                {
                    string strCurSeqPattern = strSeq.substr(iPos, m_iKmerLen);
                    iHitPos = LooseMatch(strCurSeqPattern, strSeed, m_iKmerLen, iPos);
                }
                if(-1 == iHitPos)
                    bDoLooseMatch = false;
            }

            if(iHitPos != -1)
            {
                iHitNumber++;
                iOffSetTail = iHitPos;

                if(!bGetTagE)
                {
                    if(vTagTail[i] == 'E')
                        bGetTagE = true;
                }

                if(!bGetTagT)
                {
                    if(vTagTail[i] == 'T')
                        bGetTagT = true;
                }

//                if(vTagTail[i] == 'T')
//                    iPremiunTailNum++;
            }
        }
    }
    if(iHitNumber >= vKmerNodeToTail.size()*.8)  // To Tail Done!!
        bGoodToTail = true;
    else
    {
        //这个时候涉及到容错，既然我们有容错那么意味着，我们参与匹配的片段要够长，因此一定要跨越ET或者HS
        //同时用于匹配的部分也要长，并且hit的数量也不能太少
        if( ((int)vKmerNodeToTail.size() > m_iKmerLen) &&
            (iHitNumber >= ((int)vKmerNodeToTail.size() - m_iKmerLen)) &&
            bGetTagE && bGetTagT)
        {
            bGoodToTail = true;
        }
    }

    ///For the case: To Head
    bool bGoodToHead = false;
    iHitNumber = 0;
    int iHeadBackStartPos = m_iReadLen - iOffSetHead - 1;
    bool bGetTagS = false; //S: Start
    bool bGetTagH = false; //H: Head
    if(0 == iResultCodeHead) // 0: 说明已经遍历到了最后的Head节点
    {
        int i = 0;
        bool bDoLooseMatch = true;
        for(vector<unsigned int>::iterator itr = vKmerNodeToHead.begin(); itr != vKmerNodeToHead.end(); itr++, i++)
        {
            string strSeed = ConvertNum32ToKmer(*itr, m_iKmerLen);
            int iHitPos = strSeq.rfind(strSeed, iHeadBackStartPos);

            if(iHitPos == -1 && vTagHead[i] == 'H' && bDoLooseMatch)
            {
                int iPos = (iHeadBackStartPos - m_iKmerLen + 1)-3;
                if(iPos >=0) // 这是往前，所以以0为定界
                {
                    string strCurSeqPattern = strSeq.substr(iPos, m_iKmerLen);
                    iHitPos = LooseMatch(strCurSeqPattern, strSeed, m_iKmerLen, iPos);
                }
                if(-1 == iHitPos)
                    bDoLooseMatch = false;
            }

            if(iHitPos != -1)
            {
                iHitNumber++;
                iHeadBackStartPos = iHitPos + m_iKmerLen - 1;

                if(!bGetTagS)
                {
                    if(vTagHead[i] == 'S')
                        bGetTagS = true;
                }

                if(!bGetTagH)
                {
                    if(vTagHead[i] == 'H')
                        bGetTagH = true;
                }
            }
        }
    }
    if(iHitNumber >= vKmerNodeToHead.size()*.8)  // To Tail Done!!
        bGoodToHead = true;
    else
    {
//        //这个时候涉及到容错，既然我们有容错那么意味着，我们参与匹配的片段要够长，因此一定要跨越ET或者HS
//        //同时用于匹配的部分也要长，并且hit的数量也不能太少

//        int iMinGoodHitNumber = vKmerNodeToHead.size() > m_iKmerLen ? (vKmerNodeToHead.size() - m_iKmerLen):
//                                                                      (vKmerNodeToHead.size() - (m_iKmerLen*.6));
//        if( (vKmerNodeToHead.size() > m_iKmerLen * .6) &&
//            (iHitNumber >= iMinGoodHitNumber) &&
//            bGetTagS && bGetTagH)
//        {
//            bGoodToHead = true;
//        }

        if( ((int)vKmerNodeToHead.size() > m_iKmerLen) &&
            (iHitNumber >= ((int)vKmerNodeToHead.size() - m_iKmerLen)) &&
            bGetTagS && bGetTagH)
        {
            bGoodToHead = true;
        }
    }

    //Check 首尾是不是离的比较近
    int HeadTailDistance = abs((iOffSetTail + m_iKmerLen - 1) - (iHeadBackStartPos - m_iKmerLen));
    int iMinDistance = 2; //m_iReadLen * .40; //5;

    if(HeadTailDistance >= iMinDistance) // 我们需要前后的距离非常的近这样才符合需求
        return false;


    if(bGoodToTail && bGoodToHead)
    {
        //如果都不是, 那么我们要返回false
        if(!stProbCandi.pExon->GetIsCRNAHeadExon() && !stProbCandi.pExon->GetIsCRNATailExon())
        {
            return false;
        }

        ///********** This is the last step to save current candi *****************
        bool bFind = false;
        for(vector<St_Candidate>::iterator itr = vSelfCircCandi.begin(); itr != vSelfCircCandi.end(); itr++)
        {
            if(*itr == stProbCandi)
            {
                itr->iSupportNum++;
                bFind = true;
                break;
            }
            else if(stProbCandi.CheckSimilarity(*itr))
            {
               //选择短的
               if(stProbCandi.GetLength() < itr->GetLength()) // 这个新的短
               {
                   stProbCandi.iSupportNum = itr->iSupportNum + 1;
                   *itr = stProbCandi;
               }
               else
               {
                   itr->iSupportNum++;
               }

               bFind = true;
               break;
            }
        }
        if(!bFind) //找不到
        {
            vSelfCircCandi.push_back(stProbCandi);  //这里的理解方式应该为，从后往前，所以start是大的index，end是小的index
        }
        return true;
    }
    else
        return false;
}

bool ClsCircRNADetection::CheckCurRegCaseByDBG(St_Candidate& stProbCandi, string& strSeq,
                                               St_ExonPos& stProbExonPosDonor, St_ExonPos& stProbExonPosAcceptor,
                                               unordered_map<unsigned int, St_Node>& mpSDBG,
                                               St_Candidate& stCurCandi, vector<St_Candidate>& vRegCircCandi)
{
    // Donor:    means the tail of the exon hit with the reads
    // Acceptor: means the head of the exon hit with the reads


    //Do some normal filter --> something like self circ case
    /* 我们的策略如下:
     * (1) 从sequence的两端入手
     *     a) 左端: 从左到右，使用next进行: 边走边走快 --> 到最后Premiun的点的时候停下来  --> 记录两个数组，然后进行整合
     *        *有个要注意，如果多个没有搞定，那么我们就要后递增的走知道找到相应的好的点
     *     b) 右端: 从右到左，使用pre进行:
     */

    //from start to end --> 我们在这里期待的是: 能够有一半的到到达相应的尾部
    //--> 这三个变量之后会用到
    unordered_map<unsigned int, St_Node>::iterator itrDBGSeed;
    St_Link stCurLink;
    unsigned int uiSeed = 0;
    int iOffSetTail = 0;
    //<--

//    St_Raw_Exon& stExonDonor = (*m_pvChrom)[stProbExonPosDonor.ucChromIndex].vRG[stProbExonPosDonor.uiGeneIndex].\
//                               vRT[stProbExonPosDonor.ucTranscriptIndex].vRExon[stProbExonPosDonor.ucExonIndex];


    for(int i = 0; i < (int)strSeq.length() - m_iKmerLen + 1; i++)
    {
        string strKmer = strSeq.substr(i, m_iKmerLen);
        uiSeed = ConvertKmerToNum32(strKmer);
        itrDBGSeed = mpSDBG.find(uiSeed);
        if(itrDBGSeed != mpSDBG.end()) //证明找到了
        {
            bool bFindProbExon = false;
            for(vector<St_ExonInfo>::iterator itr = itrDBGSeed->second.vEI.begin();
                itr != itrDBGSeed->second.vEI.end(); itr++)
            {
                if(itr->stExonPos == stProbExonPosDonor)
                {
                    bFindProbExon = true;
                    stCurLink = (*itr->vLink.begin());
                    break;
                }
            }
            if(bFindProbExon)
                break;
            else
                iOffSetTail++;
        }
        else
            iOffSetTail++;
    }

    if(itrDBGSeed == mpSDBG.end())
        return false;

    //--> 然后我们来遍历
    ///先从左边开始遍历 ---> 目的是tail
    /// 干脆将节点都拿出来 -->

    // 2个条件
    // (1) 到了tail处
    // (2) walk的距离大于一个reads的长度
    ///To Tail Part

    vector<unsigned int> vKmerNodeToTail;
    vector<char> vTagTail;
    int iPathFullLen = 0;
     //这里一定要注意: 我们是跳着取的!! "/3"
    int iResultCodeTail = CheckDBGTail(stCurLink, iPathFullLen,
                                       vKmerNodeToTail, vTagTail, stProbExonPosDonor, uiSeed, mpSDBG);

//    //-->Cout all Head Node --> For debug
//    for(vector<unsigned int>::iterator itr = vKmerNodeToTail.begin(); itr != vKmerNodeToTail.end(); itr++)
//    {
//        cout << ConvertNum32ToKmer(*itr, m_iKmerLen) << endl;
//    }
//    //<--

    if(iPathFullLen < CSDBGFULLLENMIN)
        return false;

    ///From Begin Part --> Go! -->  Going to head part ************************
    int iOffSetHead = 0;
    for(int i = strSeq.length() - m_iKmerLen; i >= 0 ; i--)
    {
        string strKmer = strSeq.substr(i, m_iKmerLen);
        uiSeed = ConvertKmerToNum32(strKmer);
        itrDBGSeed = mpSDBG.find(uiSeed);
        if(itrDBGSeed != mpSDBG.end()) //证明找到了
        {
            bool bFindProbExon = false;
            for(vector<St_ExonInfo>::iterator itr = itrDBGSeed->second.vEI.begin();
                itr != itrDBGSeed->second.vEI.end(); itr++)
            {
                if(itr->stExonPos == stProbExonPosAcceptor)
                {
                    bFindProbExon = true;
                    stCurLink = (*itr->vLink.begin());
                    break;
                }
            }
            if(bFindProbExon)
                break;
            else
                iOffSetHead++;
        }
        else
            iOffSetHead++;
    }

    if(itrDBGSeed == mpSDBG.end())
        return false;

    vector<unsigned int> vKmerNodeToHead;
    vector<char> vTagHead;
    iPathFullLen = 0;
    int iResultCodeHead = CheckDBGHead(stCurLink, iPathFullLen,
                                       vKmerNodeToHead, vTagHead, stProbExonPosAcceptor, uiSeed, mpSDBG);

//    //-->Cout all Head Node --> For debug
//    for(vector<unsigned int>::iterator itr = vKmerNodeToHead.begin(); itr != vKmerNodeToHead.end(); itr++)
//    {
//        cout << ConvertNum32ToKmer(*itr, m_iKmerLen) << endl;
//    }
//    //<--

    if(iPathFullLen < CSDBGFULLLENMIN)
        return false;

    //*************************For Canclualte******************
    ///首先check个数，如果个数少的出奇，那么就是不正确的 (我们可以在最后的discussion进行讨论)
    unsigned int iMinNodeNum = (m_iReadLen - m_iKmerLen + 1) * .8 *.3; //注意我们这里是跳着取的 --> 这就是为什么要乘以0.3
    if(vKmerNodeToTail.size() + vKmerNodeToHead.size() < iMinNodeNum)
    {
        //Check Premiun Node Number
        int iPremiunTailNum = 0;
        for(vector<char>::iterator itr = vTagTail.begin(); itr != vTagTail.end(); itr++)
        {
            if(*itr == 'T')
                iPremiunTailNum++;
        }

        int iPremiunHeadNum = 0;
        for(vector<char>::iterator itr = vTagHead.begin(); itr != vTagHead.end(); itr++)
        {
            if(*itr == 'H')
                iPremiunHeadNum++;
        }

        if(iPremiunTailNum >= 2 && iPremiunHeadNum >=2) // 我们认为能够成功hit到足够premiun的 (抓阄必会抓到一个，然后就是tail或head必须是一个)
        {}
        else
            return false;
    }

    ///For the case: To Tail
    bool bGoodToTail = false;
    int iHitNumber = 0;
    bool bGetTagE = false; //E: end
    bool bGetTagT = false; //T: tail
    //int iPremiunTailNum = 0;
    if(0 == iResultCodeTail) // 0: 说明已经遍历到了最后的tail节点
    {
        int i = 0;
        bool bDoLooseMatch = true;
        for(vector<unsigned int>::iterator itr = vKmerNodeToTail.begin(); itr != vKmerNodeToTail.end(); itr++, i++)
        {
            string strSeed = ConvertNum32ToKmer(*itr, m_iKmerLen);
            int iHitPos = strSeq.find(strSeed, iOffSetTail);
            if(iHitPos == -1 && vTagTail[i] == 'T' && bDoLooseMatch)
            {
                int iPos = iOffSetTail + 3;
                if(iPos + m_iKmerLen <= (int)strSeq.length()) // 这是往后所以 以length为定界
                {
                    string strCurSeqPattern = strSeq.substr(iPos, m_iKmerLen);
                    iHitPos = LooseMatch(strCurSeqPattern, strSeed, m_iKmerLen, iPos);
                }
                if(-1 == iHitPos)
                {
                    //iOffSetTail += 3;
                    bDoLooseMatch = false;
                }
            }
            else if(iHitPos == -1 && i > 3 && vTagTail[i] == 'E' && bDoLooseMatch)
            {
               iOffSetTail += 3;
            }

            if(iHitPos != -1)
            {
                iHitNumber++;
                iOffSetTail = iHitPos;

                if(!bGetTagE)
                {
                    if(vTagTail[i] == 'E')
                        bGetTagE = true;
                }

                if(!bGetTagT)
                {
                    if(vTagTail[i] == 'T')
                        bGetTagT = true;
                }

//                if(vTagTail[i] == 'T')
//                    iPremiunTailNum++;
            }
        }
    }

    if(iHitNumber >= vKmerNodeToTail.size()*.8)  // To Tail Done!!
        bGoodToTail = true;
    else
    {
//        int iMinGoodHitNumber = vKmerNodeToTail.size() > m_iKmerLen ? (vKmerNodeToTail.size() - m_iKmerLen):
//                                                                      (vKmerNodeToTail.size() - (m_iKmerLen*.6));

//        if( (vKmerNodeToTail.size() > m_iKmerLen*.6) &&
//            (iHitNumber >= iMinGoodHitNumber) &&
//            bGetTagE && bGetTagT)
//        {
//            bGoodToTail = true;
//        }

        //这个时候涉及到容错，既然我们有容错那么意味着，我们参与匹配的片段要够长，因此一定要跨越ET或者HS
        //同时用于匹配的部分也要长，并且hit的数量也不能太少
        if( (vKmerNodeToTail.size() > m_iKmerLen*.6) &&
            (iHitNumber >= (vKmerNodeToTail.size() - (m_iKmerLen*.6))) &&
            bGetTagE && bGetTagT)
        {
            bGoodToTail = true;
        }
    }

    ///For the case: To Head
    bool bGoodToHead = false;
    iHitNumber = 0;
    int iHeadBackStartPos = m_iReadLen - iOffSetHead - 1;
    bool bGetTagS = false; //S: Start
    bool bGetTagH = false; //H: Head
    if(0 == iResultCodeHead) // 0: 说明已经遍历到了最后的Head节点
    {
        int i = 0;
        bool bDoLooseMatch = true;
        for(vector<unsigned int>::iterator itr = vKmerNodeToHead.begin(); itr != vKmerNodeToHead.end(); itr++, i++)
        {
            string strSeed = ConvertNum32ToKmer(*itr, m_iKmerLen);
            int iHitPos = strSeq.rfind(strSeed, iHeadBackStartPos);

            if(iHitPos == -1 && vTagHead[i] == 'H' && bDoLooseMatch)
            {
                int iPos = (iHeadBackStartPos - m_iKmerLen + 1)-3;
                if(iPos >=0) // 这是往前，所以以0为定界
                {
                    string strCurSeqPattern = strSeq.substr(iPos, m_iKmerLen);
                    iHitPos = LooseMatch(strCurSeqPattern, strSeed, m_iKmerLen, iPos);
                }
                if(-1 == iHitPos)
                {
                    //iHeadBackStartPos -= 3;
                    bDoLooseMatch = false;
                }
            }
            else if(iHitPos == -1 && i > 3 && vTagHead[i] == 'S' && bDoLooseMatch)
                iHeadBackStartPos -= 3;

            if(iHitPos != -1)
            {
                iHitNumber++;
                iHeadBackStartPos = iHitPos + m_iKmerLen - 1;

                if(!bGetTagS)
                {
                    if(vTagHead[i] == 'S')
                        bGetTagS = true;
                }

                if(!bGetTagH)
                {
                    if(vTagHead[i] == 'H')
                        bGetTagH = true;
                }
            }
        }
    }
    if(iHitNumber >= vKmerNodeToHead.size()*.8)  // To Tail Done!!
        bGoodToHead = true;
    else
    {
        //这个时候涉及到容错，既然我们有容错那么意味着，我们参与匹配的片段要够长，因此一定要跨越ET或者HS
        //同时用于匹配的部分也要长，并且hit的数量也不能太少
        int iMinGoodHitNumber = vKmerNodeToHead.size() > m_iKmerLen ? ((int)vKmerNodeToHead.size() - m_iKmerLen):
                                                                      ((int)vKmerNodeToHead.size() - (m_iKmerLen*.6));
        if( (vKmerNodeToHead.size() > m_iKmerLen * .6) &&
            (iHitNumber >= iMinGoodHitNumber) &&
            bGetTagS && bGetTagH)
        {
            bGoodToHead = true;
        }
    }

    //Check 首尾是不是离的比较近
    int HeadTailDistance = abs((iOffSetTail + m_iKmerLen - 1) - (iHeadBackStartPos - m_iKmerLen));
    int iMinDistance = 2; //m_iReadLen * .40; //5;

    if(HeadTailDistance >= iMinDistance) // 我们需要前后的距离非常的近这样才符合需求
       return false;

    if(bGoodToTail && bGoodToHead)
    {
#ifdef _DEBUG
        bool bDebugTarget = false;
        if(stProbCandi.iStartPos == 121116815 && stProbCandi.iEndPos == 121115832)
        {
            cout << "Find Target: " << "121116815, 121115832" << endl;
            bDebugTarget = true;
        }
#endif

        bool bFind = false;
        for(vector<St_Candidate>::iterator itr = vRegCircCandi.begin(); itr != vRegCircCandi.end(); itr++)
        {
            if(*itr == stProbCandi)
            {
#ifdef _DEBUG
                if(bDebugTarget)
                    cout << "*itr == stProbCandi" << endl;
#endif

                itr->iSupportNum++;
                bFind = true;
                stCurCandi = stProbCandi;
                break;
            }
            else if(stProbCandi.CheckSimilarity(*itr))
            {
               //选择短的
               if(stProbCandi.GetLength() < itr->GetLength()) // 这个新的短
               {
#ifdef _DEBUG
                   if(bDebugTarget)
                       cout << "stProbCandi.GetLength() < itr->GetLength()" << endl;
#endif

                   stProbCandi.iSupportNum = itr->iSupportNum + 1;
                   *itr = stProbCandi;
                   stCurCandi = stProbCandi;
               }
               else
               {
#ifdef _DEBUG
                   if(bDebugTarget)
                   {
                       cout << "stProbCandi.GetLength() > itr->GetLength()"
                            << IntToStr(itr->iStartPos) << ", "
                            << IntToStr(itr->iEndPos) << endl;
                   }
#endif

                   itr->iSupportNum++;
                   stCurCandi = *itr;
               }

               bFind = true;
               break;
            }
        }
        if(!bFind) //找不到
        {
#ifdef _DEBUG
            if(bDebugTarget)
                cout << "New Insert" << endl;
#endif
            vRegCircCandi.push_back(stProbCandi);  //这里的理解方式应该为，从后往前，所以start是大的index，end是小的index
            stCurCandi = stProbCandi;
        }
        return true;
    }
    else
        return false;
}

int ClsCircRNADetection::LooseMatch(string strSeq1, string strSeq2, int iLen, int iPos)
{
    if((int)strSeq1.length() != iLen || (int)strSeq2.length() != iLen)
        return -1;

    string iSubSeq11 = strSeq1.substr(0, iLen/3);
    string iSubSeq12 = strSeq1.substr(iLen/3, iLen/3);
    string iSubSeq13 = strSeq1.substr(iLen/3 * 2, iLen/3);

    string iSubSeq21 = strSeq2.substr(0, iLen/3);
    string iSubSeq22 = strSeq2.substr(iLen/3, iLen/3);
    string iSubSeq23 = strSeq2.substr(iLen/3 * 2, iLen/3);

    int iMatch = 0;
    if(iSubSeq11 == iSubSeq21)
        iMatch++;
    if(iSubSeq12 == iSubSeq22)
        iMatch++;
    if(iSubSeq13 == iSubSeq23)
        iMatch++;

    if(iMatch < 2)
        return -1;
    else
        return iPos;
}

int ClsCircRNADetection::CheckDBGTail(St_Link& stCurLink, int& iPathFullLen,
                                      vector<unsigned int>& vKmerNode, vector<char>& vTagTail,
                                      St_ExonPos& stProbExonPos, unsigned int uiPreSeed,
                                      unordered_map<unsigned int, St_Node>& mpSDBG)
{
    /*
     * 我现在知道怎么去找比对点了:
     * 其实比对点有两个性质
     * (1)保证它能够跟reads对上 --> 也就是这个东西确实是能够跟reads align上 --> 前后的分别10个premium的点
     * (2)保证它有区分性: 也就是能够跟那些share长度的exon有区分  --> 每个区分点的kmer要进行记录
     *
     * 在这里我们也没有必要保持所有的点都需要拿进去，拿一些就可以了，如果顺序能够弄好，并且能够全部拿到，那我觉得那么结果就一定是它了
     */
    bool bReachTail = false;
    int iResultCode = -1;
    int iMaxIteratorLen = m_iReadLen;
    int iWalkLen = 0;
    int iNumEI = 0;
    vKmerNode.clear();
    vTagTail.clear();

    while(!bReachTail && iWalkLen < iMaxIteratorLen)
    {
        iPathFullLen++;

        unsigned int uiCurSeed = stCurLink.uiNextKmer;
        //cout << ConvertNum64ToKmer(uiCurSeed, m_iKmerLen) << endl;
        unordered_map<unsigned int, St_Node>::iterator itrNextNode = mpSDBG.find(uiCurSeed);
        if(itrNextNode == mpSDBG.end())
        {
            bReachTail = true;
            iResultCode = 1;
        }
        else
        {
            bool bFind = false;
            vector<St_ExonInfo>::iterator itrExonInfo;
            for(vector<St_ExonInfo>::iterator itr = itrNextNode->second.vEI.begin();
                itr != itrNextNode->second.vEI.end(); itr++)
            {
                if(itr->stExonPos == stProbExonPos)
                {
                    itrExonInfo = itr;
                    bFind = true;
                    break;
                }
            }
            if(bFind)
            {
                St_Link stNextLink; // 实际上这个应该是CurLink，之前的CurLink应该是PreLink
                //找到这个时候的link  --> we call it as next link -->
                for(vector<St_Link>::iterator itrLink = itrExonInfo->vLink.begin();
                    itrLink != itrExonInfo->vLink.end(); itrLink++)
                {
                    if(itrLink->uiPreKmer == uiPreSeed)
                    {
                        //那么就是这个link
                        stNextLink = *itrLink;
                        break;
                    }
                }
                //<--

                int iCurNumEI = itrNextNode->second.vEI.size();
                bool bSaveCurSeed = false;
                //only save the node which comes from ending part (with tag: E, T and U)
                if(stCurLink.cNextTag == 'E' || stCurLink.cNextTag == 'T' || stCurLink.cNextTag == 'U')
                {
                    if(iNumEI == 0)
                    {
                        iNumEI = iCurNumEI;
                        vKmerNode.push_back(uiCurSeed); //这个还是当前的kmer node
                        vTagTail.push_back(stCurLink.cNextTag);
                        bSaveCurSeed = true;
                    }
                    else
                    {
//                        if(iNumEI != iCurNumEI)
//                        {
//                            vKmerNode.push_back(uiCurSeed); //这个还是当前的kmer node
//                            vTagTail.push_back(stCurLink.cNextTag);
//                            bSaveCurSeed = true;
//                        }
//                        else //即使相同, 3的倍数
                        {
                            if(iWalkLen < 3)
                            {
                                vKmerNode.push_back(uiCurSeed);
                                vTagTail.push_back(stCurLink.cNextTag);
                                bSaveCurSeed = true;
                            }
                            else if( iWalkLen % 3 == 0)
                            {
                                vKmerNode.push_back(uiCurSeed);
                                vTagTail.push_back(stCurLink.cNextTag);
                                bSaveCurSeed = true;
                            }
                            else
                            {}
                        }
                    }
                }
                //给完点了，现在更新link
                iWalkLen++;
                uiPreSeed = stCurLink.uiNextKmer;
                char cPreSeedTag = stCurLink.cNextTag;
                stCurLink = stNextLink;
                if(cPreSeedTag == 'T' && stCurLink.uiNextKmer == 0 && stCurLink.cNextTag == 'U')
                {
                    //we already reached to the end
                    if(!bSaveCurSeed)
                    {
                        vKmerNode.push_back(uiCurSeed);
                        vTagTail.push_back(stCurLink.cCurTag);  //这里要十分注意，已经变成了current tag们因为之前安的link已经变成了next link
                        //cout << endl << ConvertNum32ToKmer(uiCurSeed, m_iKmerLen) << endl;
                    }
                    bReachTail = true;
                    iResultCode = 0;
                }

            }
            else
            {
                bReachTail = true;
                iResultCode = 2;
            }
        }
    }
    return iResultCode;
}

//可以看看这两个可不可以合并: 但是这里最主要的是我们在这里需要遍历到head !!!  --> 搞起！！
int ClsCircRNADetection::CheckDBGHead(St_Link& stCurLink, int& iPathFullLen,
                                      vector<unsigned int>& vKmerNode, vector<char>& vTagHead,
                                      St_ExonPos& stProbExonPos, unsigned int uiNextSeed,
                                      unordered_map<unsigned int, St_Node>& mpSDBG)
{
    /*
     * 我现在知道怎么去找比对点了:
     * 其实比对点有两个性质
     * (1)保证它能够跟reads对上 --> 也就是这个东西确实是能够跟reads align上 --> 前后的分别10个premium的点
     * (2)保证它有区分性: 也就是能够跟那些share长度的exon有区分  --> 每个区分点的kmer要进行记录
     *
     * 在这里我们也没有必要保持所有的点都需要拿进去，拿一些就可以了，如果顺序能够弄好，并且能够全部拿到，那我觉得那么结果就一定是它了
     */
    bool bReachHead = false;
    int iResultCode = -1;
    int iMaxIteratorLen = m_iReadLen;
    int iWalkLen = 0;
    int iNumEI = 0;
    vKmerNode.clear();
    vTagHead.clear();

    while(!bReachHead && iWalkLen < iMaxIteratorLen)
    {
        iPathFullLen++;
        unsigned int uiCurSeed = stCurLink.uiPreKmer;
        unordered_map<unsigned int, St_Node>::iterator itrNextNode = mpSDBG.find(uiCurSeed);
        if(itrNextNode == mpSDBG.end())
        {
            bReachHead = true;
            iResultCode = 1;
        }
        else
        {
            bool bFind = false;
            vector<St_ExonInfo>::iterator itrExonInfo;
            for(vector<St_ExonInfo>::iterator itr = itrNextNode->second.vEI.begin();
                itr != itrNextNode->second.vEI.end(); itr++)
            {
                if(itr->stExonPos == stProbExonPos)
                {
                    itrExonInfo = itr;
                    bFind = true;
                    break;
                }
            }
            if(bFind)
            {
                St_Link stPreLink; // 实际上这个应该是CurLink，之前的CurLink应该是PreLink
                //找到这个时候的link  --> we call it as next link -->
                for(vector<St_Link>::iterator itrLink = itrExonInfo->vLink.begin();
                    itrLink != itrExonInfo->vLink.end(); itrLink++)
                {
                    if(itrLink->uiNextKmer == uiNextSeed)
                    {
                        //那么就是这个link
                        stPreLink = *itrLink;
                        break;
                    }
                }
                //<--

                int iCurNumEI = itrNextNode->second.vEI.size();
                bool bSaveCurSeed = false;

                //only save the node which comes from Head part (with tag: S, H and U)
                if(stCurLink.cNextTag == 'S' || stCurLink.cNextTag == 'H' || stCurLink.cNextTag == 'U')
                {
                    if(iNumEI == 0)
                    {
                        iNumEI = iCurNumEI;
                        vKmerNode.push_back(uiCurSeed); //这个还是当前的kmer node
                        vTagHead.push_back(stCurLink.cPreTag);
                        bSaveCurSeed = true;
                    }
                    else
                    {
//                        if(iNumEI != iCurNumEI)
//                        {
//                            vKmerNode.push_back(uiCurSeed); //这个还是当前的kmer node
//                            vTagHead.push_back(stCurLink.cPreTag);
//                            bSaveCurSeed = true;
//                        }
//                        else //即使相同, 3的倍数
                        {
                            if(iWalkLen < 3)
                            {
                                vKmerNode.push_back(uiCurSeed);
                                vTagHead.push_back(stCurLink.cPreTag);
                                bSaveCurSeed = true;
                            }
                            else if( iWalkLen % 3 == 0)
                            {
                                vKmerNode.push_back(uiCurSeed);
                                vTagHead.push_back(stCurLink.cPreTag);
                                bSaveCurSeed = true;
                            }
                            else
                            {}
                        }
                    }
                }
                //给完点了，现在更新link
                iWalkLen++;
                uiNextSeed = stCurLink.uiPreKmer;
                char cNextSeedTag = stCurLink.cPreTag;
                stCurLink = stPreLink;
                if(cNextSeedTag == 'H' && stCurLink.uiPreKmer == 0 && stCurLink.cPreTag == 'U')
                {
                    //we already reached to the end
                    if(!bSaveCurSeed)
                    {
                        vKmerNode.push_back(uiCurSeed);
                        vTagHead.push_back(stCurLink.cCurTag);
                    }
                    bReachHead = true;
                    iResultCode = 0;
                }

            }
            else
            {
                bReachHead = true;
                iResultCode = 2;
            }
        }
    }
    return iResultCode;
}
