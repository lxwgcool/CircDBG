#include "clsclassificaiton.h"
#include "clsgtfparse.h"
#include "clsbasealgorithm.h"
#include "clsfastareader.h"
#include <unistd.h>
#include <map>
//#include <stdio.h>

const int IMBALANCELINE = 15; //2; //number of k-mer  --> 10 kmer = length 25 bps
const int ADDITIONALLINE = 8; //the real sequence length
const int MINOFFSET = 5;

string ALIGNSTATUS[msMax] = {"Imbalance", "LowQuality", "AdditionalPart", "AdditionalPart(Ref)",
                             "MultiSupport(Regular)",
                             "MultiSupport(Self)", "MultiSupport(R&S Diff)",
                             "MultiSupport(R&S Chimeric)", "Good", "Bad"};

ClsClassificaiton::ClsClassificaiton()
{
    m_strDBGResult = "";
    m_strGTF = "";
    m_strRef = "";
    m_strCircCandi = "";
    m_iReadLen = -1;
    m_iKmerLen = -1;
    m_bShowFlank = false;
}

ClsClassificaiton::~ClsClassificaiton()
{}

void ClsClassificaiton::Init(string strV1, string strV2,
                             string strV3, string strV4, string strV5,
                             int iV6, int iV7, bool bV8)
{
    m_strDBGResult = strV1;
    m_strGTF = strV2;
    m_strRef = strV3;
    m_strCircCandi = strV4;
    m_strBlastRootFolder = strV5;
    m_iReadLen = iV6;
    m_iKmerLen = iV7;
    m_bShowFlank = bV8;
}

//*********************************************

void ClsClassificaiton::ClassifyCirc()
{
    //Step 1: Parse DBG Result and save everything into the structure
    ParseDBGResult(); // Done

    //Step 2: Only Keep the CircCandi Comes From Configure File
    FilterCircCandi(); //Only Keep the candi recorded in m_strCircCandi

//    //Step 3: Parse GTF
//    vector<St_Row_Chrom> vChrom;
//    ClsGTFParse* pGTFParse = new ClsGTFParse();
//    pGTFParse->ReadGTF(vChrom, m_strGTF); //Done
//    pGTFParse->GetTagValue(vChrom, m_strRef);
//    delete pGTFParse;
//    pGTFParse = NULL;

    //Step 4: Get Blast Reference Sequence For each candidate
    SetBlastRefForCandi();//(vChrom);

    //Step 5: Use Blast to check Category for each circular candidate by each reads
    CheckCircCategory();

    //In one candi --> Check chimeric
    CheckChimeric();
}

void ClsClassificaiton::ParseDBGResult()
{
    if(access(m_strDBGResult.c_str(), 0) != 0)
    {
        cout << "File doesn't existed!" << endl;
        return;
    }

    m_vCandi.clear();

    ifstream infile;
    infile.open(m_strDBGResult.c_str(), ios::in);
    string strLine = "";
    St_Candidate stCandi;
    string strSeq = "";

    while(!infile.eof()) // check if reached the end of file
    {
        getline(infile, strLine);
        //Set the value for Donor Exon and Acceptor Exon
        if(strLine.find("Donor") != string::npos && strLine.find("Acceptor") != string::npos)
        {
            //Set value for donor and acceptor
            int iStartPos = 0;
            int iEndPos = 0;
            int iLen = 0;

            ///Get Donor Start
            iStartPos = strLine.find(":");
            iStartPos += 3;
            iEndPos = strLine.find(",", iStartPos);
            iLen = iEndPos - iStartPos;
            stCandi.stExonDonor.iStart = atoi(strLine.substr(iStartPos, iLen).c_str());

            ///Get Donor End
            iStartPos = iEndPos + 2;
            iEndPos = strLine.find(")", iStartPos);
            iLen = iEndPos - iStartPos;
            stCandi.stExonDonor.iEnd = atoi(strLine.substr(iStartPos, iLen).c_str());

            ///Get Acceptor Start
            iStartPos = strLine.find("(", iEndPos);
            iStartPos += 1;
            iEndPos = strLine.find(",", iStartPos);
            iLen = iEndPos - iStartPos;
            stCandi.stExonAcceptor.iStart = atoi(strLine.substr(iStartPos, iLen).c_str());

            ///Get Acceptor End
            iStartPos = iEndPos + 2;
            iEndPos = strLine.find(")", iStartPos);
            iLen = iEndPos - iStartPos;
            stCandi.stExonAcceptor.iEnd = atoi(strLine.substr(iStartPos, iLen).c_str());
        }
        else if(strLine.find("Candi") != string::npos) //Set Value for Candi
        {
            int iStartPos = 0;
            int iEndPos = 0;
            int iLen = 0;

            ///Get Candi End
            iStartPos = strLine.find(":");
            iStartPos += 2;
            iEndPos = strLine.find(" ", iStartPos);
            iLen = iEndPos - iStartPos;
            stCandi.iEndPos = atoi(strLine.substr(iStartPos, iLen).c_str());

            ///Get Candi Start
            iStartPos = iEndPos + 3;
            iEndPos = strLine.length();
            iLen = iEndPos - iStartPos;
            stCandi.iStartPos = atoi(strLine.substr(iStartPos, iLen).c_str());
        }
        else if(strLine.find("RC") != string::npos) // Set direction
        {
            bool bRC = true;
            if(strLine.find("No"))
            {
                bRC = false;
            }
            stCandi.bRC = bRC;
            stCandi.stExonDonor.bRC = bRC;
            stCandi.stExonAcceptor.bRC = bRC;
        }
        else if(strLine.find("Chrom") != string::npos) // Set Chrom Index
        {
            int iStartPos = strLine.find(" ") + 1;
            int iEndPos = strLine.length();
            int iLen = iEndPos - iStartPos;
            stCandi.strChromName = strLine.substr(iStartPos, iLen);
        }
        else if(strLine.find("Current Reads") != string::npos) // For Set the reads mapping direction
        {}
        else if(strLine.find("CircType") != string::npos)
        {
            int iStartPos = strLine.find(" ") + 1;
            int iEndPos = strLine.length();
            int iLen = iEndPos - iStartPos;
            switch(atoi(strLine.substr(iStartPos, iLen).c_str()))
            {
                case 0:
                    stCandi.enType = ctSelf;
                    break;
                case 1:
                    stCandi.enType = ctRegular;
                    break;
                default:
                    stCandi.enType = ctMax;
                    break;
            }
        }
        else if(strLine.find("***") != -1) // Section End, save current candi
        {
            bool bFind = false;
            for(vector<St_Candidate>::iterator itr = m_vCandi.begin(); itr != m_vCandi.end(); itr++)
            {
                if(*itr == stCandi)
                {
                    itr->vReadsSeq.push_back(strSeq);
                    stCandi.Clear();
                    bFind = true;
                    break;
                }
            }
            if(!bFind) // Is a new one
            {
                stCandi.vReadsSeq.push_back(strSeq);
                m_vCandi.push_back(stCandi);
                stCandi.Clear();
                strSeq = "";
            }
        }
        else //The reads sequence
        {
            strSeq = strLine;
        }
    }

    infile.close();
}

// Brief.txt (same directory as executable file)
void ClsClassificaiton::ParseBrief(string strCurPath, vector<St_Candidate>& vCandi)
{
    //Check if file existed
    if(access(strCurPath.c_str(), 0) != 0)
    {
        cout << "File doesn't exsited ! --> " << strCurPath << endl;
        return;
    }

    ifstream ifs;
    ifs.open(strCurPath.c_str());
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "")
            continue;

        //1: For Chrom Index
        int iStart = 0;
        int iEnd = strLine.find(" ");
        int iLen = iEnd - iStart;
        stCandi.strChromName = strLine.substr(iStart, iLen);

        //2: Start Pos
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        int iFirstValue = atoi(strLine.substr(iStart, iLen).c_str());

        //3: End Pos
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        int iSecondValue = atoi(strLine.substr(iStart, iLen).c_str());

        stCandi.iStartPos = iFirstValue <= iSecondValue ? iFirstValue : iSecondValue;
        stCandi.iEndPos = iFirstValue >= iSecondValue ? iFirstValue : iSecondValue;

        //4:Support Reads Count
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        stCandi.iSupportNum = atoi(strLine.substr(iStart, iLen).c_str());

        //5: Tag
        iStart = iEnd + 1;
        iEnd = strLine.find(" ", iStart);
        iLen = iEnd - iStart;
        stCandi.strTag = strLine.substr(iStart, iLen);

        //6: Circular Type
        iStart = iEnd + 1;
        iEnd = strLine.length();
        iLen = iEnd - iStart;
        string strType = strLine.substr(iStart, iLen);
        trim(strType);
        if(strType == "S")
            stCandi.enType = ctSelf;
        else if(strType == "R")
            stCandi.enType = ctRegular;
        else
            stCandi.enType = ctMax;

        //6: Save Current Candidate
        //if(stCandi.iSupportNum < 5) // half coverage
            vCandi.push_back(stCandi);
    }
    ifs.close();
}

void ClsClassificaiton::FilterCircCandi()
{
    if(access(m_strCircCandi.c_str(), 0) != 0)
    {
        cout << "Target Candi File doesn't exist!" << endl;
        return;
    }

    vector<St_Candidate> vBriefCandi;
    ParseBrief(m_strCircCandi, vBriefCandi);
    cout << "vBriefCandi Size: " << vBriefCandi.size() << endl;
    cout << "m_vCandi Size   : " << m_vCandi.size() << endl;

//    cout << endl << "Check If we have duplicate (vBriefCandi) -->" << endl;

//    int iDuplicateNum = 0;
//    for(vector<St_Candidate>::iterator itr = vBriefCandi.begin(); itr != vBriefCandi.end(); itr++)
//    {
//        bool bFind = false;
//        for(vector<St_Candidate>::iterator subItr = itr+1;
//            subItr != vBriefCandi.end(); subItr++)
//        {
//            if(*itr == *subItr)
//            {
//                bFind = true;
//                break;
//            }
//        }
//        if(bFind)
//        {
//            cout << itr->strChromName << " " << itr->iStartPos << " " << itr->iEndPos << endl;
//            iDuplicateNum++;
//        }
//    }

//    cout << endl;
//    cout << iDuplicateNum << endl << endl;



    for(vector<St_Candidate>::iterator itr = m_vCandi.end()-1; itr >= m_vCandi.begin(); itr--)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vBriefCandi.begin();
            subItr != vBriefCandi.end(); subItr++)
        {
            if(*itr == *subItr)
            {
                bFind = true;
                break;
            }
        }
        if(!bFind) // 没找到 --> 那么就删掉
        {
            //cout << itr->strChromName << " " << itr->iStartPos << " " << itr->iEndPos << endl;
            m_vCandi.erase(itr);
        }
    }
    cout << "Filtered m_vCandi Num: " << m_vCandi.size() << endl;

//    cout << endl;
//    int iMissNum = 0;
//    //Check Which do not find
//    for(vector<St_Candidate>::iterator itr = vBriefCandi.begin(); itr != vBriefCandi.end(); itr++)
//    {
//        bool bFind = false;
//        for(vector<St_Candidate>::iterator subItr = m_vCandi.begin(); subItr != m_vCandi.end(); subItr++)
//        {
//            if(*itr == *subItr)
//            {
//                bFind = true;
//                break;
//            }
//        }
//        if(!bFind)
//        {
//            cout << itr->strChromName << " " << itr->iStartPos << " " << itr->iEndPos << endl;
//            iMissNum++;
//        }
//    }

//    cout << endl;
//    cout << iMissNum << endl << endl;
}

void ClsClassificaiton::SetBlastRefForCandi()//(vector<St_Row_Chrom>& vChrom)
{
    if(m_vCandi.empty())
    {
        cout << "No Candidate !" << endl;
        return;
    }

    if(access(m_strRef.c_str(), 0) != 0)
    {
        cout << "Geno Reference File doesn't exist!" << endl;
        return;
    }

    //Read Genome Reference
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(m_strRef, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

    //If everything works fine, let's do it!!!
    for(vector<St_Candidate>::iterator itr = m_vCandi.begin(); itr != m_vCandi.end(); itr++)
    {
        //Step 1: Get the related Geno Rererence Sequence
        St_Fasta* pCurRefFa = NULL;
        for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
        {
            //cout << itrRef->strName << endl;
            string strRefName = "";
            if(itrRef->strName.find(' ') == string::npos)
                strRefName = itrRef->strName;
            else
                strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));
            //cout << "strRefName: " << strRefName << endl;

            if(strRefName == itr->strChromName)
            {
                pCurRefFa = &(*itrRef);
                break;
            }
        }
        if(pCurRefFa == NULL)
        {
            cout << "Do not find related reference chromosone: " << itr->strChromName << endl;
            continue;
        }
        ///<----

        //Step 2: Get Donor and Acceptor Sequence
        ///1: For Donor
        string strExonDonorSeq = pCurRefFa->strSeq.substr(itr->stExonDonor.iStart-1, itr->stExonDonor.GetLength());

        ///2: For Acceptor
        string strExonAcceptorSeq = pCurRefFa->strSeq.substr(itr->stExonAcceptor.iStart-1, itr->stExonAcceptor.GetLength());

        CreateNormalRef(strExonDonorSeq, strExonAcceptorSeq, *itr);
        CreateDonorSelfRef(strExonDonorSeq, *itr);
        CreateAcceptorSelfRef(strExonAcceptorSeq, *itr);
    }
}

//Do not worry about duplicate, since each segment comes from different exon
void ClsClassificaiton::CreateNormalRef(string& strExonDonorSeq, string& strExonAcceptorSeq,
                                        St_Candidate& stCandi)
{
    //if(strExonDonorSeq.length() + strExonAcceptorSeq.length() < m_iReadLen)
        //return;

    string strReadsHeadMappingPart = strExonDonorSeq;
    string strReadsTailMappingPart = strExonAcceptorSeq;

    switch(stCandi.enType)
    {
        case ctSelf:
//        {
//            int iLength = strExonDonorSeq.length();
//            strReadsHeadMappingPart = strExonDonorSeq.substr(strExonDonorSeq.length()-iLength/2, iLength/2);
//            strReadsTailMappingPart = strExonAcceptorSeq.substr(0, iLength/2);
//            break;
//        }
        case ctRegular:
        {
            if(strExonDonorSeq.length() > m_iReadLen) //取后面的部分
                strReadsHeadMappingPart = strExonDonorSeq.substr(strExonDonorSeq.length()-m_iReadLen, m_iReadLen);

            if(strExonAcceptorSeq.length() > m_iReadLen) //取后面的部分
                strReadsTailMappingPart = strExonAcceptorSeq.substr(0, m_iReadLen);
            break;
        }
        default:
            break;
    }

    //Step 3: Combine the valid part of Donor and the valid part of Acceptor
    stCandi.iBreakPoint =strReadsHeadMappingPart.length();
    //cout << stCandi.iBreakPoint << endl;
    stCandi.strBlastRef = strReadsHeadMappingPart + strReadsTailMappingPart;
}

void ClsClassificaiton::CreateDonorSelfRef(string& strExonDonorSeq, St_Candidate& stCandi)
{
//    //if(strExonDonorSeq.length() < m_iReadLen)
//    //    return;

//    int iRawLen = strExonDonorSeq.length() / 2;
//    string strReadsHeadMappingPart = strExonDonorSeq.substr(iRawLen, strExonDonorSeq.length() - iRawLen);
//    //cout << strReadsHeadMappingPart << endl << endl;
//    if(iRawLen > m_iReadLen) //取后面的部分
//        strReadsHeadMappingPart = strExonDonorSeq.substr(strExonDonorSeq.length()-m_iReadLen, m_iReadLen);

//    string strReadsTailMappingPart = strExonDonorSeq.substr(0, iRawLen); //strExonDonorSeq;
//    //cout << strReadsTailMappingPart << endl << endl;
//    if(iRawLen > m_iReadLen) //取后面的部分
//        strReadsTailMappingPart = strExonDonorSeq.substr(0, m_iReadLen);

//    //Step 3: Combine the valid part of Donor and the valid part of Acceptor
//    stCandi.iSelfDonorBP =strReadsHeadMappingPart.length();
//    stCandi.strSelfDonorRef = strReadsHeadMappingPart + strReadsTailMappingPart;

//    //stCandi.strSelfDonorRef = strExonDonorSeq;

      stCandi.strSelfDonorRef = strExonDonorSeq;
      stCandi.iSelfDonorBP = -1;
}

void ClsClassificaiton::CreateAcceptorSelfRef(string& strExonAcceptorSeq, St_Candidate& stCandi)
{
//    //if(strExonAcceptorSeq.length() < m_iReadLen)
//    //    return;

//    int iRawLen = strExonAcceptorSeq.length() / 2;

//    string strReadsHeadMappingPart = strExonAcceptorSeq.substr(iRawLen, strExonAcceptorSeq.length() - iRawLen);
//    if(iRawLen > m_iReadLen) //取后面的部分
//        strReadsHeadMappingPart = strExonAcceptorSeq.substr(strExonAcceptorSeq.length()-m_iReadLen, m_iReadLen);

//    string strReadsTailMappingPart = strExonAcceptorSeq.substr(0, iRawLen);
//    //cout << strReadsTailMappingPart << endl << endl;
//    if(iRawLen > m_iReadLen) //取后面的部分
//        strReadsTailMappingPart = strExonAcceptorSeq.substr(0, m_iReadLen);

//    //Step 3: Combine the valid part of Donor and the valid part of Acceptor
//    stCandi.iSelfAcceptorBP =strReadsHeadMappingPart.length();
//    stCandi.strSelfAcceptorRef = strReadsHeadMappingPart + strReadsTailMappingPart;

//    //stCandi.strSelfAcceptorRef = strExonAcceptorSeq;

    stCandi.strSelfAcceptorRef = strExonAcceptorSeq;
    stCandi.iSelfAcceptorBP = -1;
}

void ClsClassificaiton::CheckCircCategory() //Use blast to check mapping detail information
{
    if(m_vCandi.empty())
        return;

    string strCmd = (string)"mkdir -p " + get_current_dir_name() + "/Classify_Result";
    system(strCmd.c_str());
    strCmd = "rm ./Classify_Result/*";
    system(strCmd.c_str());

    ClsBlast* pBlast = new ClsBlast();
    pBlast->Init(m_strBlastRootFolder);

    vector<St_BlastResult> vBlastResult;

    ofstream ofs;
    ofs.open("./Classify_Result/AlignInfoSum.txt");


//    cout << "m_vCandi Size: " << m_vCandi.size() << endl;
//    map<string, int> mpCount;
//    for(vector<St_Candidate>::iterator itr = m_vCandi.begin(); itr != m_vCandi.end(); itr++)
//    {
//        mpCount[itr->strChromName]++;
//    }

//    for(map<string, int>::iterator itr = mpCount.begin(); itr != mpCount.end(); itr++)
//    {
//        cout << "Chrom " << itr->first << ": " << itr->second << endl;
//    }

//    return;

    for(vector<St_Candidate>::iterator itr = m_vCandi.begin(); itr != m_vCandi.end(); itr++)
    {
        if(itr->vReadsSeq.empty())
            continue;

        if(itr->strBlastRef == "")
            continue;

        //Step 1: -->Create Fasta File for Ref Seq (Blast Ref Seq)
        string strRefFa = pBlast->CreatFaFile("RefSeq", itr->strBlastRef, ftRef);

        for(vector<string>::iterator subItr = itr->vReadsSeq.begin();
            subItr != itr->vReadsSeq.end(); subItr++)
        {
            //Step 2: -->Create Fasta File for Query Seq (each reads)
            string strQueryFa = pBlast->CreatFaFile("QuerySeq", *subItr, ftQuery);
            vBlastResult.clear();
            pBlast->GetTwoSeqAlignResult(strQueryFa, strRefFa, vBlastResult, true, true);

            //Step 3:Analysis the blast alignment result and get the mapping Category  --> Go Later
            if(vBlastResult.size() == 1)
            {
                CheckBlastResult(*vBlastResult.begin(), *itr, *subItr, ofs);
            }
        }
    }
    ofs.close();

    delete pBlast;
    pBlast = NULL;

    //Output the statistic result
    ///1:取最多的作为最后的质量  --> 这个会覆盖掉一些特殊情况 --> 我们需要明天想想  --> Go Later
    /// 这里的质量可以这么进行就考虑
    /// 1:通过blast得出的质量包括以下4个: msImbalance, msLowQuality, msAdditionalPart, msGood
    ///   以上四个也是主质量
    /// 2:通过简单遍历得出的质量包括以下3个: msMultiSupportS, msMultiSupportSRDiff, msMultiSupportSRChimeric
    ///   以上三个是特殊的有意思的case (也可以说是复杂case)，这个需要单独进行相应的报告
    //所以在这里我们的统计逻辑应该是
    //1: 有good就是good
    //2: 没有good那么就看其他三种哪个大，如果一样, 那么 Additional Part > imbalance > PoorQuality

    int iSumGood = 0;
    int iSumBad = 0;
    int iSumAdditionalPart = 0;
    int iSumAdditionalPart_Ref = 0;
    int iSumImbalance = 0;
    int iSumLowQuality = 0;

    //Save all of good candidate
    ofstream ofsGood;
    ofsGood.open("./Classify_Result/Brief_sum_good.txt");

    for(vector<St_Candidate>::iterator itr = m_vCandi.begin(); itr != m_vCandi.end(); itr++)
    {
        int iAdditionalPart = 0;
        int iAdditionalPart_Ref = 0;
        int iImbalance = 0;
        int iLowQuality = 0;
        bool bFindGood = false;
        bool bFindBad = false;

        for(vector<En_MapStatus>::iterator subStatus = itr->vMapStatus.begin();
            subStatus != itr->vMapStatus.end(); subStatus++)
        {
            switch(*subStatus)
            {
                case msGood:
                    bFindGood = true;
                    break;
                case msBad:
                    bFindBad = true;
                    break;
                case msAdditionalPart:
                    iAdditionalPart++;
                    break;
                case msAdditionalPart_Ref:
                    iAdditionalPart_Ref++;
                    break;
                case msImbalance:
                    iImbalance++;
                    break;
                case msLowQuality:
                    iLowQuality++;
                    break;
                default:
                    break;
            }
            if(bFindGood)
                break;
        }

        if(bFindGood)
        {
            itr->enMajorMapStatus = msGood;
            iSumGood++;
            //---> Save good
            ofsGood << itr->strChromName << " " // << IntToStr(itr->ucChromIndex) << " "
                    << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                    << IntToStr(itr->iSupportNum) << " "
                    << itr->strTag << " "
                    << itr->GetTypeString() << endl;
            //<---
            continue;
        }        
        else
        {
            //找最大的 -->
            int iMax = (iAdditionalPart >= iImbalance ? iAdditionalPart : iImbalance);           
            iMax = (iMax >= iAdditionalPart_Ref ? iMax : iAdditionalPart_Ref);
            iMax = (iMax >= iLowQuality ? iMax : iLowQuality);

            if(iMax == 0)
            {
                if(bFindBad)
                {
                    itr->enMajorMapStatus = msBad;
                    //************ --> Output Bad
                    cout << "Bad: " << IntToStr(itr->iStartPos) << ", " << IntToStr(itr->iEndPos) << endl;
                    //<---------------
                    iSumBad++;
                }
            }
            else
            {
                if(iAdditionalPart == iMax)
                {
                    itr->enMajorMapStatus = msAdditionalPart;
                    iSumAdditionalPart++;
                }
                else if(iAdditionalPart_Ref == iMax)
                {
                    itr->enMajorMapStatus = msAdditionalPart_Ref;
                    iSumAdditionalPart_Ref++;
                }
                else if(iImbalance == iMax)
                {
                    itr->enMajorMapStatus = msImbalance;
                    iSumImbalance++;
                }
                else if(iLowQuality == iMax)
                {
                    itr->enMajorMapStatus = msLowQuality;
                    iSumLowQuality++;
                }
            }
        }
    }

    ofsGood.close();

    cout << "Total: " << IntToStr(m_vCandi.size()) << endl;
    cout << "Sum_Good              : " << IntToStr(iSumGood) << " --> "
                                  << GetRatio((float)iSumGood/m_vCandi.size()) << endl;
    cout << "Sum_Bad               : " << IntToStr(iSumBad) << " --> "
                                  << GetRatio((float)iSumBad/m_vCandi.size()) << endl;
    cout << "Sum_AdditionalPart     : " << IntToStr(iSumAdditionalPart) << " --> "
                                  << GetRatio((float)iSumAdditionalPart/m_vCandi.size()) << endl;
    cout << "Sum_AdditionalPart(Ref): " << IntToStr(iSumAdditionalPart_Ref) << " --> "
                                  << GetRatio((float)iSumAdditionalPart_Ref/m_vCandi.size()) << endl;
    cout << "Sum_Imbalance          : " << IntToStr(iSumImbalance) << " --> "
                                  << GetRatio((float)iSumImbalance/m_vCandi.size()) << endl;
    cout << "Sum_LowQuality         : " << IntToStr(iSumLowQuality) << " --> "
                                  << GetRatio((float)iSumLowQuality/m_vCandi.size()) << endl;      
    //Get Special Case:
    //GetSpecialMappingCase();
}

// Check the case of msMultiSupportS, msMultiSupportSRDiff, msMultiSupportSRChimeric
///--> We can do it later
void ClsClassificaiton::GetSpecialMappingCase()
{
    cout << endl << "-----------------------------" << endl << endl;

    //Get the following 4 special circular rna cases:---->
    /// 1: msMultiSupportR
    /// 2: msMultiSupportS
    /// 3: msMultiSupportSRDiff
    /// 4: msMultiSupportSRChimeric
    //Go <-----

    int iMultiSupportR = 0;
    int iMultiSupportS = 0;
    int iMultiSupportSRDiff = 0;
    int iMultiSupportSRChimeric = 0;

    int iSpecilGroupNum = 0;

    for(vector<St_Candidate>::iterator itrCandi = m_vCandi.begin(); itrCandi != m_vCandi.end(); itrCandi++)
    {
        int iSameNum = 0;
        vector<St_Candidate> vSelectCandi;
        for(vector<St_Candidate>::iterator itrLaterCandi = m_vCandi.end() - 1;
            itrLaterCandi > itrCandi; itrLaterCandi--)
        {
            bool bHitSame = false;
            for(vector<string>::iterator itrReads = itrCandi->vReadsSeq.begin();
                itrReads != itrCandi->vReadsSeq.end(); itrReads++)
            {
                for(vector<string>::iterator itrLaterReads = itrLaterCandi->vReadsSeq.begin();
                    itrLaterReads != itrLaterCandi->vReadsSeq.end(); itrLaterReads++)
                {
                    if(*itrLaterReads == *itrReads)
                    {
                        bHitSame = true;
                        break;
                    }
                }
                if(bHitSame)
                    break;
            }

            if(bHitSame)
            {
                vSelectCandi.push_back(*itrLaterCandi);
                m_vCandi.erase(itrLaterCandi);
                iSameNum++;
            }
        }

        if(iSameNum == 0)
            continue;
        else
        {
            iSpecilGroupNum++;
            //cout duplicate case
            cout << "Current Candi: " << "(" << IntToStr(itrCandi->iStartPos) << ", "
                 << IntToStr(itrCandi->iEndPos) << ")" << endl;
            cout << "Same Num: " << IntToStr(iSameNum) << endl;
            cout << "------------ Detail ----------" << endl;
            cout << IntToStr(vSelectCandi.size()) << endl;

            for(vector<St_Candidate>::iterator itrSelectCandi = vSelectCandi.begin();
                itrSelectCandi != vSelectCandi.end(); itrSelectCandi++)
            {
                switch(itrCandi->enType)
                {
                    case ctSelf:
                    {
                        if(itrSelectCandi->enType == ctSelf)
                        {
                            iMultiSupportS++;
                            cout << "S -- S" << endl;
                        }
                        else if(itrSelectCandi->enType == ctRegular)
                        {
                            if(itrCandi->CheckOneSideSimilarity(*itrSelectCandi))
                            {
                                //一边相等
                                iMultiSupportSRChimeric++;
                                cout << "S -- R (Chimeric)" << endl;
                            }
                            else
                            {
                                iMultiSupportSRDiff++;
                                cout << "S -- R (Diff)" << endl;
                            }
                        }
                        break;
                    }
                    case ctRegular:
                    {
                        if(itrSelectCandi->enType == ctSelf)
                        {
                            if(itrCandi->CheckOneSideSimilarity(*itrSelectCandi))
                            {
                                //一边相等
                                iMultiSupportSRChimeric++;
                                cout << "S -- R (Chimeric)" << endl;
                            }
                            else
                            {
                                iMultiSupportSRDiff++;
                                cout << "S -- R (Diff)" << endl;
                            }

                        }
                        else if(itrSelectCandi->enType == ctRegular)
                        {
                            iMultiSupportR++;
                            cout << "R -- R" << endl;
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
            cout << "***" << endl << endl;
        }
        vSelectCandi.clear();
    }

    //Statistic
    cout << "Special Group Num: " << IntToStr(iSpecilGroupNum) << endl;

    int iSum = iMultiSupportR + iMultiSupportS + iMultiSupportSRDiff + iMultiSupportSRChimeric;
    if(iSum == 0)
    {
        cout << "Nothing Strange" << endl;
    }
    else
    {
        cout << "Sum              : " << IntToStr(iSum) << endl;
        cout << "R -- R           : " << GetRatio((float)iMultiSupportR / iSum) << endl;
        cout << "S -- S           : " << GetRatio((float)iMultiSupportS / iSum) << endl;
        cout << "S -- R (Diff)    : " << GetRatio((float)iMultiSupportSRDiff / iSum) << endl;
        cout << "S -- R (Chimeric): " << GetRatio((float)iMultiSupportSRChimeric / iSum) << endl;
    }
}

void ClsClassificaiton::CheckBlastResult(St_BlastResult& stBlastResult, St_Candidate& stCandi,
                                         string& strCurReads, ofstream& ofs)
{
    if(stBlastResult.vBlastAlign.empty())
        return;

    //Get the best alignment
    vector<St_BlastAlignUnit>::iterator itrBestAlgn = stBlastResult.vBlastAlign.begin();

    for(vector<St_BlastAlignUnit>::iterator itr = stBlastResult.vBlastAlign.begin() + 1;
        itr != stBlastResult.vBlastAlign.end(); itr++)
    {
        if(itrBestAlgn->iAlignNum < itr->iAlignNum)
        {
            itrBestAlgn = itr;
        }
    }

    //Check Align Status
    St_BlastAlignUnit& stBA = *itrBestAlgn;  //Check later
    En_MapStatus enMapStatus = msMax;

    /// 1: Check Bad First: breakpoint do not dropped into the algined part
    if(stCandi.iBreakPoint <=  stBA.iStartPos || stCandi.iBreakPoint >=  stBA.iEndPos)
    {
        //--> do some further detections for this special case
        bool bSpecialAdditional = false;
        if(stCandi.iBreakPoint <= stBA.iStartPos) // BP is in the left
        {
            int iOffSetRight = stBA.iStartPos - stCandi.iBreakPoint;
            int iOffSetLeft = 10000;
            for(vector<St_BlastAlignUnit>::iterator itr = stBlastResult.vBlastAlign.begin();
                itr != stBlastResult.vBlastAlign.end(); itr++)
            {
                if(itr->iEndPos <= stCandi.iBreakPoint)
                {
                    if(iOffSetLeft > (stCandi.iBreakPoint - itr->iEndPos))
                        iOffSetLeft = stCandi.iBreakPoint - itr->iEndPos;
                }
            }
            if(iOffSetRight + iOffSetLeft < 21)
            {
                bSpecialAdditional = true;
            }
        }
        else if(stCandi.iBreakPoint >=  stBA.iEndPos) //BP is in the right
        {
            int iOffSetLeft = stCandi.iBreakPoint - stBA.iEndPos;
            int iOffSetRight = 10000;
            for(vector<St_BlastAlignUnit>::iterator itr = stBlastResult.vBlastAlign.begin();
                itr != stBlastResult.vBlastAlign.end(); itr++)
            {
                if(itr->iStartPos >= stCandi.iBreakPoint)
                {
                    if(iOffSetRight > (itr->iStartPos - stCandi.iBreakPoint))
                        iOffSetRight = itr->iStartPos - stCandi.iBreakPoint;
                }
            }
            if(iOffSetRight + iOffSetLeft < 25)
            {
                bSpecialAdditional = true;
            }
        }
        //<--
        if(bSpecialAdditional)
            enMapStatus = msAdditionalPart_Ref;
        else
            enMapStatus = msBad;
    }
    /// 2: Check the additional Part
    else if( abs(stBA.iAlignNum - m_iReadLen) >= ADDITIONALLINE)
        enMapStatus = msAdditionalPart;
    /// 3: Check if it is imbalance
    else if( (stCandi.iBreakPoint >= stBA.iStartPos && (stCandi.iBreakPoint - stBA.iStartPos <= m_iKmerLen)) ||
             (stCandi.iBreakPoint <= stBA.iEndPos && (stCandi.iEndPos - stCandi.iBreakPoint <= m_iKmerLen)))
        enMapStatus = msImbalance;
    else if(abs(stCandi.iBreakPoint - stBA.iStartPos - m_iKmerLen) < IMBALANCELINE ||
            abs(stBA.iEndPos - stCandi.iBreakPoint - m_iKmerLen) < IMBALANCELINE )
        enMapStatus = msImbalance;
    /// 4: Check the quality
    else if(stBA.fGapRatio > 0.1 || stBA.fIdentifyRatio < .99)
        enMapStatus = msLowQuality;    
    else
        /// 4: Check MultiSupport Cases (msMultiSupportR, msMultiSupportS,
        ///                              msMultiSupportSRDiff, msMultiSupportSRChimeric)
        /// 5: Check if it is Bad: msBad
        enMapStatus = msGood;

    stCandi.vMapStatus.push_back(enMapStatus);

    //Output Result -->
    ofs << "Donor(Start, End) & Acceptor(Start, End): "
         << "(" << IntToStr(stCandi.stExonDonor.iStart) << ", " << IntToStr(stCandi.stExonDonor.iEnd) << ")"
         << " -- "
         << "(" << IntToStr(stCandi.stExonAcceptor.iStart) << ", " << IntToStr(stCandi.stExonAcceptor.iEnd) << ")" << endl;

    ofs << "Candi Start & End: "<< IntToStr(stCandi.iStartPos) << " : "
         << IntToStr(stCandi.iEndPos) << endl;

    ofs << "BreakPoint: " << IntToStr(stCandi.iBreakPoint) << endl;
    ofs << "Align Status: " << ALIGNSTATUS[enMapStatus] << endl;

    ofs << "Current Reads: " << endl
         << strCurReads << endl;

    ofs << itrBestAlgn->strAlignInfo << endl;
    ofs << "***" << endl << endl;
    //<--
}

void ClsClassificaiton::CheckChimeric()  // we should re-identify the logic
{
    cout << endl << "****************Check Chimeric****************" << endl;
    //--> Check
    ClsBlast* pBlast = new ClsBlast();
    pBlast->Init(m_strBlastRootFolder);

    int iChimericNum = 0;
    for(vector<St_Candidate>::iterator itr = m_vCandi.begin(); itr != m_vCandi.end(); itr++)
    {
        if(itr->enType == ctSelf || itr->enMajorMapStatus == msBad)
            continue;

        if(itr->enMajorMapStatus == msAdditionalPart || itr->enMajorMapStatus == msLowQuality)
        {}
        else
            continue;

        //对于 Regular Circ 的例子
        if(itr->vReadsSeq.empty())
            continue;

        bool bFindDonorChimeric = false;
        bool bFindAcceptorChimeric = false;

        for(vector<string>::iterator subItr = itr->vReadsSeq.begin();
            subItr != itr->vReadsSeq.end(); subItr++)
        {
//            int iCutLen = m_iReadLen / 2;
//            string strHeadReads = subItr->substr(0, iCutLen);
//            string strTailReads = subItr->substr(iCutLen, subItr->length() - iCutLen);
//            string strQuery = strTailReads + strHeadReads;
//            //int iQueryBP = strTailReads.length();

//            //strQuery = *subItr;

            string strQueryFa = pBlast->CreatFaFile("QuerySeq", *subItr, ftQuery);

            ///Check Donor Chimeric
            bFindDonorChimeric = CheckDonorSelfAlign(strQueryFa, *itr, pBlast);
            if(bFindDonorChimeric)
            {
                cout << "Donor" << endl;
                cout << *subItr << endl;
                break;
            }

            ///Check Acceptor Chimeric
            bFindAcceptorChimeric = CheckAcceptorSelfAlign(strQueryFa, *itr, pBlast);
            if(bFindAcceptorChimeric)
            {
                cout << "Acceptor" << endl;
                cout << *subItr << endl;
                break;
            }
        }

        if(bFindDonorChimeric || bFindAcceptorChimeric)
        {
            iChimericNum++;
            cout << "Candi Start & End: "<< IntToStr(itr->iStartPos) << " : "
                 << IntToStr(itr->iEndPos) << endl;

            cout << endl << " ---" << endl << endl;
        }

    }

    cout << endl << ">>>>>>>>>>>>>>>>" << endl
         << "Chimeric Sum: " << IntToStr(iChimericNum) << " --> "
         << GetRatio((float)iChimericNum/m_vCandi.size()) << endl;

    delete pBlast;
    pBlast = NULL;
}

///************Do Self Donor Ref -->*************
bool ClsClassificaiton::CheckDonorSelfAlign(string strQueryFa, St_Candidate& stCandi, ClsBlast* pBlast)
{

    if(stCandi.strSelfDonorRef == "")
        return false;

//    cout << "Self Donor Ref:" << endl;
//    cout << stCandi.strSelfDonorRef << endl;
    string strRefFa = pBlast->CreatFaFile("RefSeq", stCandi.strSelfDonorRef, ftRef);

    vector<St_BlastResult> vBlastResult;
    pBlast->GetTwoSeqAlignResult(strQueryFa, strRefFa, vBlastResult, true, true);

    //Step 3:Analysis the blast alignment result and get the mapping Category  --> Go Later
    if(vBlastResult.size() != 1)
        return false;

    if(vBlastResult.begin()->vBlastAlign.size() < 2)
        return false;

    //Get two best alignment
    vector<St_BlastAlignUnit>::iterator itrBest1Algn = vBlastResult.begin()->vBlastAlign.begin();
    vector<St_BlastAlignUnit>::iterator itrBest2Algn = vBlastResult.begin()->vBlastAlign.begin() + 1;
    if(itrBest1Algn->iAlignNum < itrBest2Algn->iAlignNum)
    {
        vector<St_BlastAlignUnit>::iterator itrTemp = itrBest1Algn;
        itrBest1Algn = itrBest2Algn;
        itrBest2Algn = itrTemp;
    }

    for(vector<St_BlastAlignUnit>::iterator itrAlignUnit = vBlastResult.begin()->vBlastAlign.begin() + 2;
        itrAlignUnit !=  vBlastResult.begin()->vBlastAlign.end(); itrAlignUnit++)
    {
        if(itrBest1Algn->iAlignNum < itrAlignUnit->iAlignNum)
        {
            itrBest2Algn = itrBest1Algn;
            itrBest1Algn = itrAlignUnit;
        }
        else if(itrBest2Algn->iAlignNum < itrAlignUnit->iAlignNum)
        {
            itrBest2Algn = itrAlignUnit;
        }
    }

    //Check if one of those two best align is totally contained by another (we need to filt this case)
    if( (itrBest1Algn->iQueryStartPos >= itrBest2Algn->iQueryStartPos &&
         itrBest1Algn->iQueryEndPos <= itrBest2Algn->iQueryEndPos) ||
        (itrBest2Algn->iQueryStartPos >= itrBest1Algn->iQueryStartPos &&
         itrBest2Algn->iQueryEndPos <= itrBest1Algn->iQueryEndPos) )
        return false;

    //Check Align Status
    St_BlastAlignUnit& strBA1 = *itrBest1Algn;
    St_BlastAlignUnit& strBA2 = *itrBest2Algn;

    bool bHitDonorSelf = false;
    bool bHeadAlign = false;
    bool bTailAlign = false;

    if(strBA1.iStartPos < MINOFFSET || strBA2.iStartPos < MINOFFSET)
    {
        bHeadAlign = true;
    }
    if( abs(strBA1.iEndPos - (int)stCandi.strSelfDonorRef.length()) < MINOFFSET ||
        abs(strBA2.iEndPos - (int)stCandi.strSelfDonorRef.length()) < MINOFFSET)
    {
        bTailAlign = true;
    }

    if(bHeadAlign && bTailAlign)
        bHitDonorSelf = true;

    return bHitDonorSelf;
}

bool ClsClassificaiton::CheckAcceptorSelfAlign(string strQueryFa, St_Candidate& stCandi, ClsBlast* pBlast)
{
    ///************Do Self Acceptor Ref -->*************
    if(stCandi.strSelfAcceptorRef == "")
        return false;

    string strRefFa = pBlast->CreatFaFile("RefSeq", stCandi.strSelfAcceptorRef, ftRef);

    vector<St_BlastResult> vBlastResult;
    pBlast->GetTwoSeqAlignResult(strQueryFa, strRefFa, vBlastResult, true, true);


    //Step 3:Analysis the blast alignment result and get the mapping Category  --> Go Later
    if(vBlastResult.size() != 1)
        return false;

    if(vBlastResult.begin()->vBlastAlign.size() < 2)
        return false;

    //Get two best alignment
    vector<St_BlastAlignUnit>::iterator itrBest1Algn = vBlastResult.begin()->vBlastAlign.begin();
    vector<St_BlastAlignUnit>::iterator itrBest2Algn = vBlastResult.begin()->vBlastAlign.begin() + 1;
    if(itrBest1Algn->iAlignNum < itrBest2Algn->iAlignNum)
    {
        vector<St_BlastAlignUnit>::iterator itrTemp = itrBest1Algn;
        itrBest1Algn = itrBest2Algn;
        itrBest2Algn = itrTemp;
    }

    for(vector<St_BlastAlignUnit>::iterator itrAlignUnit = vBlastResult.begin()->vBlastAlign.begin() + 2;
        itrAlignUnit !=  vBlastResult.begin()->vBlastAlign.end(); itrAlignUnit++)
    {
        if(itrBest1Algn->iAlignNum < itrAlignUnit->iAlignNum)
        {
            itrBest2Algn = itrBest1Algn;
            itrBest1Algn = itrAlignUnit;
        }
        else if(itrBest2Algn->iAlignNum < itrAlignUnit->iAlignNum)
        {
            itrBest2Algn = itrAlignUnit;
        }
    }

    //Check if one of those two best align is totally contained by another (we need to filt this case)
    if( (itrBest1Algn->iQueryStartPos >= itrBest2Algn->iQueryStartPos &&
         itrBest1Algn->iQueryEndPos <= itrBest2Algn->iQueryEndPos) ||
        (itrBest2Algn->iQueryStartPos >= itrBest1Algn->iQueryStartPos &&
         itrBest2Algn->iQueryEndPos <= itrBest1Algn->iQueryEndPos) )
        return false;

    //Check Align Status
    St_BlastAlignUnit& strBA1 = *itrBest1Algn;
    St_BlastAlignUnit& strBA2 = *itrBest2Algn;

    bool bHitAcceptorSelf = false;
    bool bHeadAlign = false;
    bool bTailAlign = false;

    if(strBA1.iStartPos < MINOFFSET || strBA2.iStartPos < MINOFFSET)
    {
        bHeadAlign = true;
    }
    if( abs(strBA1.iEndPos - (int)stCandi.strSelfAcceptorRef.length()) < MINOFFSET ||
        abs(strBA2.iEndPos - (int)stCandi.strSelfAcceptorRef.length()) < MINOFFSET)
    {
        bTailAlign = true;
    }

    if(bHeadAlign && bTailAlign)
        bHitAcceptorSelf = true;

    return bHitAcceptorSelf;
}

