#include "clsdebruijngraph.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <iostream>

const int PREMIUMLEN = 12; // 10 bps  --> (m_kmer + 1) --> for the exon with regular length
const int PREMIUMLENSHORT = 9; // 10 bps  --> (m_kmer + 1) --> For short Exon (we should different from regular length)
//#define SPECIEXON
const int ciSpeciExonStartPos = 153904679;
const int ciSpeciExonEndPos = 153904916;

ClsDeBruijnGraph::ClsDeBruijnGraph()
{
}

ClsDeBruijnGraph::~ClsDeBruijnGraph()
{
}

void ClsDeBruijnGraph::Init(St_Config& stConfig)
{
    m_strDNARefPath = stConfig.strRefPath; //Notice: this reference comes from DNA
    m_iKmerLen = stConfig.iKmerLen;
    m_iReadLen = stConfig.iReadLen;
    m_fKmerRatio = stConfig.fKmerRatio;
}

//Go later --> Go Now
void ClsDeBruijnGraph::BuildGraph(unordered_map<unsigned int, St_Node>& mpDBG,
                                  St_Row_Chrom* pChrom, St_Fasta* pCurRefFa,
                                  int iCurChromIndex)
{
    //Clear Old Data
    //m_vDBG.clear();
    //St_DBG stDBG;

    //2: Create DBG for each exon
    St_ExonPos stExonPos;
    ofstream ofsExon;
    ofsExon.open("./ExonRef.fa");

    cout << pChrom->strName << endl;

    stExonPos.ucChromIndex = iCurChromIndex;
    stExonPos.uiGeneIndex = 0;
    for(vector<St_Raw_Gene>::iterator itrGene = pChrom->vRG.begin();
        itrGene != pChrom->vRG.end(); itrGene++)
    {
        /*
        //Check if this gene is protein-coding gene
        if(itrGene->strBioType != GENEBIOTYPE[gbtCoding] &&
           itrGene->strBioType != GENEBIOTYPE[gbtPT] &&
           itrGene->strBioType != GENEBIOTYPE[gbtRI] &&
           itrGene->strBioType != GENEBIOTYPE[gbtAantisense] &&
           itrGene->strBioType != GENEBIOTYPE[gbtNMD] &&
           itrGene->strBioType != GENEBIOTYPE[gbtUP])  // only consider the coding gene --> check how about the result
        {
            stExonPos.uiGeneIndex++;
            continue;
        }*/

        stExonPos.bRC = itrGene->bRC; //To know it is in + strand or - strand
        stExonPos.ucTranscriptIndex = 0;
        for(vector<St_Raw_Transcript>::iterator itrTranscript = itrGene->vRT.begin();
            itrTranscript != itrGene->vRT.end(); itrTranscript++)
        {
            stExonPos.ucExonIndex = 0;
            //Iterator each exon
            for(vector<St_Raw_Exon>::iterator itrExon = itrTranscript->vRExon.begin();
                itrExon != itrTranscript->vRExon.end(); itrExon++)
            {
                //Build DBG
#ifdef SPECIEXON
                if( itrExon->iStart == ciSpeciExonStartPos &&
                    itrExon->iEnd == ciSpeciExonEndPos)  //  debug: 仅仅考虑当前debug的exon
#endif
                {
                    if(itrExon->GetIsSupportCircRNA())
                    {
                        BuildDBGForSingleExon(mpDBG, *itrExon, stExonPos,
                                              pCurRefFa->strSeq, ofsExon);
                    }
                }

                stExonPos.ucExonIndex++;
            }
            stExonPos.ucTranscriptIndex++;
        }
        stExonPos.uiGeneIndex++;
    }
    ofsExon.close();

//    // Output DBG of current exon for manually checking the accuracy of exon
//    ofstream ofs;
//    ofs.open("./tmp.seq");
//    // ofs << "Exon Sequence: " << endl;
//    // DisplayString(ofs, strExonSeq);
//    // ofs << endl << "DBG" << endl;

//    for(unordered_map<unsigned int, St_Node>::iterator itr = mpDBG.begin();
//        itr != mpDBG.end(); itr++)
//    {
//        ofs << ConvertNum32ToKmer(itr->first, iKmerLen) << ": ";

//        for(vector<St_Link>::iterator subItr = itr->second.vLink.begin();
//            subItr != itr->second.vLink.end(); subItr++)
//        {
//            ofs << ConvertNum32ToKmer(subItr->uiPreKmer, iKmerLen) << ", "
//                << ConvertNum32ToKmer(subItr->uiNextKmer, iKmerLen) << " | ";
//        }
//        ofs << endl;

//        for(vector<St_ExonPos>::iterator subItr = itr->second.vExonPos.begin();
//            subItr != itr->second.vExonPos.end(); subItr++)
//        {
//            ofs << IntToStr(subItr->ucChromIndex) << ", "
//                << IntToStr(subItr->uiGeneIndex) << ", "
//                << IntToStr(subItr->ucTranscriptIndex) << ", "
//                << IntToStr(subItr->ucExonIndex) << ", "
//                << " | ";
//        }
//        ofs << endl;

//        for(vector<unsigned char>::iterator subItr = itr->second.vNum.begin();
//            subItr != itr->second.vNum.end(); subItr++)
//        {
//            ofs << IntToStr(*subItr) << " | ";
//        }
//        ofs << endl;
//    }
}

//--->从这里来看如何进行移植
void ClsDeBruijnGraph::BuildDBGForSingleExon(unordered_map<unsigned int, St_Node>& mpDBG,
                                             St_Raw_Exon& stExon, St_ExonPos& stExonPos,
                                             string& strRef, ofstream& ofsExon)
{
    if(stExon.GetLength() < m_iKmerLen) // 如果太短了，我们就直接返回
        return;

    //我们一会儿再进行这个的测试
    //if( !pExon->GetIsSupportCircRNA() ) // 如果不support 我们就不建表
    //    return;

    //看来我们在这里还是只取了前后的一截，也就是不考虑中间的部分
    int iBoundaryLen = m_iReadLen - m_iKmerLen + 1 - 5; //m_iReadLen * m_fKmerRatio;
    int iTotalExtractLen = iBoundaryLen * 2; //+ 2 * (m_iKmerLen - 1); // 这里我们是需要分开考虑的

    unordered_map<unsigned int, St_ExonInfo> mpSEDBG;
    string strExonSeq = strRef.substr(stExon.iStart - 1, abs(stExon.iEnd - stExon.iStart) + 1);

    //--> Output Current Exon sequence
    string strName = "> " + IntToStr(stExon.iStart) + " : " +  IntToStr(stExon.iEnd);
//    ofsExon << strName << endl;
//    ofsExon << strExonSeq << endl;
    //<--

    ToUpper(strExonSeq);
    if(stExon.GetLength() <= iTotalExtractLen) ///Case 1: 对于短的exon构建相应的DBG
    {
        BuildByShortExon(mpSEDBG, strExonSeq, stExonPos);
    }
    else ///Case 2: Exon的长度足够去进行左右两边的采样,然后建立相应的DBG
    {
        string strHeadPartSeq = strExonSeq.substr(0, iBoundaryLen); //+ m_iKmerLen - 1);  //这样来定边界才是合理的
        string strTailPartSeq = strExonSeq.substr(strExonSeq.length() - iBoundaryLen, //- m_iKmerLen + 1,
                                           iBoundaryLen); //+ m_iKmerLen - 1);
        BuildByRegularExon(mpSEDBG, strHeadPartSeq, strTailPartSeq, stExonPos);
    }

//    //Output DBG of current exon  --> For check if the code works fine
//    ofs << "Exon Sequence: " << endl;
//    DisplayString(ofs, strExonSeq);
//    ofs << endl << "DBG" << endl;

//    for(unordered_map<unsigned int, St_SENode>::iterator itr = mpSingleExonDBG.begin();
//        itr != mpSingleExonDBG.end(); itr++)
//    {
//        ofs << ConvertNum32ToKmer(itr->first, m_iKmerLen) << ": ";

//        //traverse next node
//        for(vector<unsigned int>::iterator subItr = itr->second.vNextKmer.begin();
//            subItr != itr->second.vNextKmer.end(); subItr++)
//        {
//            ofs << ConvertNum32ToKmer(*subItr, m_iKmerLen) << " | ";
//        }
//        ofs << endl;

//        //traverse tags
//        for(vector<char>::iterator subItr = itr->second.vTag.begin();
//            subItr != itr->second.vTag.end(); subItr++)
//        {
//            ofs << *subItr << " | ";
//        }
//        ofs << endl;
//    }

    //Merge current exon DBG (unordered map) with the Main DBG (unordered map)
    for(unordered_map<unsigned int, St_ExonInfo>::iterator itr = mpSEDBG.begin();
        itr != mpSEDBG.end(); itr++)
    {
        if(mpDBG.find(itr->first) == mpDBG.end()) // This is a new node
        {
            St_Node stNode;
            stNode.vEI.push_back(itr->second);
            mpDBG.insert(std::pair<unsigned int, St_Node>(itr->first, stNode));
        }
        else //  The node is already existed
        {
            mpDBG[itr->first].vEI.push_back(itr->second);
        }
    }

    //release resources
    mpSEDBG.clear();
}

void ClsDeBruijnGraph::BuildByShortExon(unordered_map<unsigned int, St_ExonInfo>& mpSEDBG,
                      string& strExonSeq, St_ExonPos& stExonPos)
{
    mpSEDBG.clear();
    unsigned int uiPreKmer = 0;
    int iSplitPoint = (strExonSeq.length() - m_iKmerLen) / 2;  //从iSplitPoint这点开始已经就是“E”了

    //计算一下属于前面的premium node 有几个，以及属于后面的premium node有几个
    int iHeadNodeNum = PREMIUMLENSHORT; // 实际上这个就是head node的end pos
    if(iSplitPoint - PREMIUMLENSHORT < 0)
    {
        iHeadNodeNum = iSplitPoint + 1;
    }

    int iTailNodeNum = PREMIUMLENSHORT;
    if((int)strExonSeq.length() - (iSplitPoint + m_iKmerLen + PREMIUMLENSHORT) < 0)
    {
        iTailNodeNum =  strExonSeq.length() - (iSplitPoint + m_iKmerLen);
    }
    int iTailNodeStartPos = strExonSeq.length() - (m_iKmerLen + iTailNodeNum);

    for(int i = 0; i < (int)strExonSeq.length() - m_iKmerLen + 1; i++)
    {
        //1: Transfer string to integer
        string strKmer = strExonSeq.substr(i, m_iKmerLen);
        if(strKmer.find('N') == string::npos) //we discard the kmer which contain "N"
        {
            unsigned int uiKmer = ConvertKmerToNum32(strKmer);
            char cTag = (i <= iSplitPoint ? 'S' : 'E');
            if(i < iHeadNodeNum)
                cTag = 'H';
            else if(i > iTailNodeStartPos)
                cTag = 'T';

            UpdateSingleExonDBG(mpSEDBG, stExonPos, cTag, uiKmer, uiPreKmer, i);
            uiPreKmer = uiKmer;
        }
    }
}

void ClsDeBruijnGraph::UpdateSingleExonDBG(unordered_map<unsigned int, St_ExonInfo>& mpSEDBG,
                                           St_ExonPos& stExonPos, char cTag, unsigned int uiKmer,
                                           unsigned int uiPreKmer, int iItrIndex)
{
    St_Link stLink;
    stLink.cCurTag = cTag;

    if(mpSEDBG.find(uiKmer) == mpSEDBG.end()) // 是个新的节点 --> 新节点是一定需要插入的
    {
        St_ExonInfo stEI;
        stEI.stExonPos = stExonPos;

        if(iItrIndex == 0) //证明我是第一个，前面还啥都没有呢，所以不去进行任何操作
        {
            //Insert current Node
            stEI.vLink.push_back(stLink);
            mpSEDBG.insert(std::pair<unsigned int, St_ExonInfo>(uiKmer, stEI));
        }
        else // 当下的 uiPreKmer 已经是有意义的了
        {
            if(mpSEDBG.find(uiPreKmer) != mpSEDBG.end()) //The previous node is existed.
                                                         //this should be normal case!
            {
                //1: Update the existed Node (prevous one)
                unordered_map<unsigned int, St_ExonInfo>::iterator itrPreNode = mpSEDBG.find(uiPreKmer);
                (itrPreNode->second.vLink.end()-1)->uiNextKmer = uiKmer;
                (itrPreNode->second.vLink.end()-1)->cNextTag = cTag;

                //2: Insert current Node
                stLink.uiPreKmer = uiPreKmer;
                stLink.cPreTag = (itrPreNode->second.vLink.end()-1)->cCurTag;
                stEI.vLink.push_back(stLink);
                mpSEDBG.insert(std::pair<unsigned int, St_ExonInfo>(uiKmer, stEI));
            }
            else
            {
                cout << "There must be something wrong!" << endl;
                mpSEDBG.clear();
                return;
            }
        }
    }
    else //如果这是一个老的节点
    {
        if(iItrIndex == 0) // 在这里我们也成功避开了第一个节点 --> 不应该避开，每个节点都用该有对应的link
        {
            //这里不应该避开，而应该正确的去更新相应的link
            unordered_map<unsigned int, St_ExonInfo>::iterator itrCurNode = mpSEDBG.find(uiKmer);
            itrCurNode->second.vLink.push_back(stLink);
        }
        else
        {
            if(mpSEDBG.find(uiPreKmer) != mpSEDBG.end())
            {
                //Update the existed node (previous one)
                unordered_map<unsigned int, St_ExonInfo>::iterator itrPreNode = mpSEDBG.find(uiPreKmer);
                (itrPreNode->second.vLink.end()-1)->uiNextKmer = uiKmer;
                (itrPreNode->second.vLink.end()-1)->cNextTag = cTag;

                //Update the existed Node (current one)  --> 主要是更新link
                unordered_map<unsigned int, St_ExonInfo>::iterator itrCurNode = mpSEDBG.find(uiKmer);
                stLink.uiPreKmer = uiPreKmer;
                stLink.cPreTag = (itrPreNode->second.vLink.end()-1)->cCurTag;
                itrCurNode->second.vLink.push_back(stLink);
            }
            else
            {
                cout << "There must be something wrong!" << endl;
                mpSEDBG.clear();
                return;
            }
        }
    }
}

void ClsDeBruijnGraph::BuildByRegularExon(unordered_map<unsigned int, St_ExonInfo>& mpSEDBG,
                                          string& strHeadPartSeq, string& strTailPartSeq,
                                          St_ExonPos& stExonPos)
{
    mpSEDBG.clear();
    unsigned int uiPreKmer = 0;

    //For head part -->
    char cTag = 'S';
    for(int i = 0; i < (int)strHeadPartSeq.length() - m_iKmerLen + 1; i++)
    {
        if(i < PREMIUMLEN)
            cTag = 'H';
        else
            cTag = 'S';

        string strKmer = strHeadPartSeq.substr(i, m_iKmerLen);
        if(strKmer.find('N') == string::npos) //we discard the kmer which contain "N"
        {
            unsigned int uiKmer = ConvertKmerToNum32(strKmer);
            UpdateSingleExonDBG(mpSEDBG, stExonPos, cTag, uiKmer, uiPreKmer, i);
            uiPreKmer = uiKmer;
        }
    }

    //For tail part <--
    cTag = 'E';
    int iTailNodeStartPos = strTailPartSeq.length() - m_iKmerLen - PREMIUMLEN;
    uiPreKmer = 0;
    for(int i = 0; i < (int)strTailPartSeq.length() - m_iKmerLen + 1; i++)  // the 'i' does matter.
    {
        if(i > iTailNodeStartPos)
            cTag = 'T';
        else
            cTag = 'E';

        string strKmer = strTailPartSeq.substr(i, m_iKmerLen);
        if(strKmer.find('N') == string::npos) //we discard the kmer which contain "N"
        {
            unsigned int uiKmer = ConvertKmerToNum32(strKmer);
            UpdateSingleExonDBG(mpSEDBG, stExonPos, cTag, uiKmer, uiPreKmer, i);
            uiPreKmer = uiKmer;
        }
    }
}

//vector<St_DBG>& ClsDeBruijnGraph::GetDBG()
//{
//    return m_vDBG;
//}
