#include "clsparsestdcircrnainfo.h"
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include "../../ShareLibrary/clsbasealgorithm.h"
#include "../../ShareLibrary/clsfastareader.h"

const int OFFSET = 3;
const char* arryCellNameCopy1[ctMaxx] = {"H9", "glioblastoma", "brain"};

ClsParseStdCircRNAInfo::ClsParseStdCircRNAInfo()
{
}

ClsParseStdCircRNAInfo::~ClsParseStdCircRNAInfo()
{
}

void ClsParseStdCircRNAInfo::Init(string strTissue)
{
    //Initiate m_enCellType
    if(strTissue == arryCellNameCopy1[ctH9])
        m_enCellType = ctH9;
    else if(strTissue == arryCellNameCopy1[ctGlioblastoma])
        m_enCellType = ctGlioblastoma;
    else if(strTissue == arryCellNameCopy1[ctBrain])
        m_enCellType = ctBrain;
    else
        m_enCellType = ctMaxx;
}

void ClsParseStdCircRNAInfo::ParseCircRNADb(string strFilePath,
                                            vector<St_StdCircRNAInfo>& vInfo)
{
    //Parse CircBase Reuslt
    /* File Structure
     * 0: circRNA name
     * 1: chromosone name
     * 2: start
     * 3: end
     *    (在这里好像都是从小到大的)
     * 4: direction: + or -
     *-2: where is this circRNA come from
     */
    ifstream ifs;
    ifs.open(strFilePath.c_str());
    vInfo.clear();
    string strLine = "";
    St_StdCircRNAInfo stInfo; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "")
            continue;

        if(strLine.length() > 0 && strLine.substr(0, 1) == "#") //jump the comments
            continue;

        int iStartPos = 0;
        int iEndPos = strLine.find('\t', iStartPos);
        int iLen = iEndPos - iStartPos;

        //The 0 --> CircRNA Name
        string strCircRNAName = strLine.substr(iStartPos, iLen);
        stInfo.strName = strCircRNAName;

        //The 1 --> Chromosome Index
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strChromName = strLine.substr(iStartPos, iLen);
        iStartPos = 3;
        iLen = strChromName.length() - 3; //
        //注意，这里我们还是从0开始，这样我们可以直接通过下标去找相应的chromosone里面的exons的信息
        //stInfo.ucChromIndex = atoi(strChromName.substr(iStartPos, iLen).c_str()) - 1;
        stInfo.strChrom = strChromName;

        //The 2 --> Start position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stInfo.iStartPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 3 --> End position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stInfo.iEndPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 4 --> Direction
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strTag = strLine.substr(iStartPos, iLen);
        if(strTag == "-")
            stInfo.bRC = true;
        else
            stInfo.bRC = false;

        //The -2 --> where is this circRNA come from (cell types)
        iEndPos = strLine.rfind('\t');
        iStartPos = strLine.rfind('\t', iEndPos-1) + 1;
        iLen = iEndPos - iStartPos;
        string strCellList = strLine.substr(iStartPos, iLen);
        stInfo.vOrgCell.clear();
        //use " " to split different cells
        iStartPos = 0;
        iLen = 0;
        for(int i=0; i<(int)strCellList.length(); i++)
        {
            if(strCellList.substr(i, 1) != ",")
                iLen++;
            else
            {
                stInfo.vOrgCell.push_back(strCellList.substr(iStartPos, iLen));
                iStartPos = i + 1;
                iLen = 0;
            }

            if(i == (int)strCellList.length() - 1) // the last letter
            {
                stInfo.vOrgCell.push_back(strCellList.substr(iStartPos, iLen));
            }
        }

        vInfo.push_back(stInfo);
    }
    ifs.close();    
}

void ClsParseStdCircRNAInfo::FilterByCell(vector<St_StdCircRNAInfo>& vInfo, En_CellType enCellType)
{
    cout << "Org Size: " << IntToStr(vInfo.size()) << endl;

    if(vInfo.empty())
    {
        return;
    }

    for(vector<St_StdCircRNAInfo>::iterator itr = vInfo.end() - 1; itr >= vInfo.begin(); itr--)
    {
        bool bFind = false;
        for(vector<string>::iterator subItr = itr->vOrgCell.begin();
            subItr != itr->vOrgCell.end(); subItr++)
        {
            if(subItr->find(arryCellNameCopy1[enCellType]) != string::npos)
            {
                bFind = true;
                break;
            }
        }
        //Fail to be found
        if(!bFind)
        {
            vInfo.erase(itr);
        }
    }

    cout << "After Filter Size: " << IntToStr(vInfo.size()) << endl;
}

void ClsParseStdCircRNAInfo::ParseGTFWithCircRNADbStyle(string strGTFPath,
                                                        vector<St_Row_Chrom>& vChrom)
{
    ifstream ifs;
    ifs.open(strGTFPath.c_str());
    vChrom.clear();
    string strLine = "";

    St_Row_Chrom stChrom;
    St_Raw_Gene stGene;
    St_Raw_Transcript stTranscript;
    St_Raw_Exon stExon;

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "")
            continue;

        if(strLine.length() > 0 && strLine.substr(0, 1) == "#") //jump the comments
            continue;        

        /* Parse the line
         * 0 --> chrom index
         * 1 --> data type: only record exon will be fine
         * 2 --> start position
         * 3 --> end posiiton
         * 5 --> direction
         * 7 --> details of this data
         */
        int iStartPos = 0;
        int iEndPos = strLine.find('\t', iStartPos);
        int iLen = iEndPos - iStartPos;

        /// 0 --> chrom Index
        stChrom.strName = strLine.substr(iStartPos, iLen);

        /// 1 --> data type
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strDataType = strLine.substr(iStartPos, iLen);
        if(strDataType == "exon") // this is exon: we only keep exon here
        {
            //Clean the data
            stChrom.vRG.clear();
            stGene.vRT.clear();
            stTranscript.vRExon.clear();
        }
        else
            continue;

        /// 2 --> start position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stExon.iStart = atoi(strLine.substr(iStartPos, iLen).c_str());

        /// 3 --> End Positoin
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stExon.iEnd = atoi(strLine.substr(iStartPos, iLen).c_str());

        /// 5 --> Direction
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;

        string strTag = strLine.substr(iStartPos, iLen);
        if(strTag == "-")
            stExon.bRC = true;
        else
            stExon.bRC = false;

        /// 7 --> details of this data
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\n', iStartPos);
        iLen = iEndPos - iStartPos;

        string strDetail = strLine.substr(iStartPos, iLen);
        //0; gene id
        iStartPos = 0;
        iEndPos = strDetail.find(';', iStartPos);
        iLen = iEndPos - iStartPos;
        string strTemp = strDetail.substr(iStartPos, iLen);

        int iTempStart = strTemp.find("\"") + 1;
        stGene.strID = strTemp.substr(iTempStart, strTemp.length()-iTempStart-1);

        //1; gene name
        iStartPos = iEndPos + 1;
        iEndPos = strDetail.find(';', iStartPos);
        iLen = iEndPos - iStartPos;
        strTemp = strDetail.substr(iStartPos, iLen);
        iTempStart = strTemp.find("\"") + 1;
        stGene.strName = strTemp.substr(iTempStart, strTemp.length()-iTempStart-1);

        //2; transcript id
        iStartPos = iEndPos + 1;
        iEndPos = strDetail.find(';', iStartPos);
        iLen = iEndPos - iStartPos;
        strTemp = strDetail.substr(iStartPos, iLen);
        iTempStart = strTemp.find("\"") + 1;
        stTranscript.strID = strTemp.substr(iTempStart, strTemp.length()-iTempStart-1);

        //以上的信息都解析成功了，现在需要合并找到的信息
        bool bAttached = false;
        for(vector<St_Row_Chrom>::iterator itr = vChrom.begin(); itr != vChrom.end(); itr++)
        {
            if(itr->strName == stChrom.strName)
            {
                for(vector<St_Raw_Gene>::iterator itrGene = itr->vRG.begin();
                    itrGene != itr->vRG.end(); itrGene++)
                {
                    if(itrGene->strID == stGene.strID)
                    {
                        for(vector<St_Raw_Transcript>::iterator itrTrans = itrGene->vRT.begin();
                            itrTrans != itrGene->vRT.end(); itrTrans++)
                        {
                            if(itrTrans->strID == stTranscript.strID)
                            {
                                //We need to attach exon
                                itrTrans->vRExon.push_back(stExon);
                                bAttached = true;
                            }
                        }
                        if(!bAttached)
                        {
                            //attach this transcript
                            stTranscript.bRC = stExon.bRC;
                            stTranscript.vRExon.push_back(stExon);
                            itrGene->vRT.push_back(stTranscript);
                            bAttached = true;
                        }
                    }
                }
                if(!bAttached)
                {
                    //try to attach this gene
                    stTranscript.bRC = stExon.bRC;
                    stTranscript.vRExon.push_back(stExon);
                    stGene.bRC = stTranscript.bRC;
                    stGene.vRT.push_back(stTranscript);
                    itr->vRG.push_back(stGene);
                    bAttached = true;
                }
            }
        }
        if(!bAttached)
        {
            //Try to attach this chrom
            stTranscript.bRC = stExon.bRC;
            stTranscript.vRExon.push_back(stExon);
            stGene.bRC = stTranscript.bRC;
            stGene.vRT.push_back(stTranscript);
            stChrom.vRG.push_back(stGene);
            vChrom.push_back(stChrom);
            bAttached = true;
        }

        if(!bAttached)
        {
            cout << "Something Wrong here!" << endl;
        }
    }

    ifs.close();

    //Cout Result
    ofstream ofs;
    ofs.open("./gtf_brif.txt");
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
    {
        ofs << "Chromoson Name: " << itrChrom->strName << endl;
        for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); itrGene != itrChrom->vRG.end(); itrGene++)
        {
            ofs << "\t" << "Gene ID: " << itrGene->strID << " | " << "Name: " << itrGene->strName << " " << (itrGene->bRC ? "-" : "+") << endl;

            for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin(); itrRT != itrGene->vRT.end(); itrRT++)
            {
                ofs << "\t\t" << "Transcript ID: " << itrRT->strID << " --- " << IntToStr(itrRT->iStart)
                    << " | " << IntToStr(itrRT->iEnd) << " " << (itrRT->bRC ? "-" : "+") << endl;

                int iIndex = 0;
                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); itrExon != itrRT->vRExon.end(); itrExon++)
                {
                    ofs << "\t\t\t" << "Exon" << IntToStr(iIndex) << ": "<< itrExon->strHead2Bp
                        << " | " << itrExon->strTail2Bp << " " << (itrExon->bRC ? "-" : "+")
                        << "   <" << IntToStr(itrExon->iStart) << ", " << IntToStr(itrExon->iEnd) << ">"
                        << " Len: " << IntToStr(abs(itrExon->iEnd - itrExon->iStart)) << endl;

                    //-->Display Exon Sequence
//                    ofs << "Seq:" << endl;
//                    string strExonSeq = vFasta[atoi(itrChrom->strName.c_str()) - 1].strSeq.substr(itrExon->iStart-1, itrExon->GetLength());
//                    DisplayString(ofs, strExonSeq);
//                    ofs << endl;
                    //<--

                    iIndex++;
                }
            }
        }
    }

    //release reference sequence
    ofs.close();
}

enum En_UnHitStatus{uhsInExon, uhsInIntron, uhsNone, uhsMax};
const char* cUNHITSTATUS[uhsMax] = {"In Exon", "In Intron", "Unknown"};
struct St_UnHitStatus
{
    St_StdCircRNAInfo stCircRNAInfo;
    string strGeneName;
    string strTransName;

    //Left Status
    En_UnHitStatus enLeftStatus;
    ///Drop Exon case
    int iLSExonStart;
    int iLSExonEnd;
    int iLSExonLen;
    int iOffsetL1;
    int iOffsetL2;
    ///Drop Intron Case
    int iLSItronStart;
    int iLSItronEnd;
    int iOffsetLI1;
    int iOffsetLI2;

    //Right Status
    En_UnHitStatus enRightStatus;
    ///Drop Exon Case
    int iRSExonStart;
    int iRSExonEnd;
    int iRSExonLen;
    int iOffsetR1;
    int iOffsetR2;
    ///Drop Intron Case
    int iRSIntronStart;
    int iRSIntronEnd;
    int iOffsetRI1;
    int iOffsetRI2;

    St_UnHitStatus()
    {
        Clear();
    }

    void Clear()
    {
        strGeneName = "";
        strTransName = "";

        //Left Status
        enLeftStatus = uhsNone;
        ///Drop Exon case
        iLSExonStart = -1;
        iLSExonEnd = -1;
        iLSExonLen = -1;
        iOffsetL1 = -1;
        iOffsetL2 = -1;
        ///Drop Intron Case
        iLSItronStart = -1;
        iLSItronEnd = -1;
        iOffsetLI1 = -1;
        iOffsetLI2 = -1;

        //Right Status
        enRightStatus = uhsNone;
        ///Drop Exon Case
        iRSExonStart = -1;
        iRSExonEnd = -1;
        iRSExonLen = -1;
        iOffsetR1 = -1;
        iOffsetR2 = -1;
        ///Drop Intron Case
        iRSIntronStart = -1;
        iRSIntronEnd = -1;
        iOffsetRI1 = -1;
        iOffsetRI2 = -1;
    }
};

bool Sort_StdCircRNAInfo(St_StdCircRNAInfo stInfo1, St_StdCircRNAInfo stInfo2)
{
    if(stInfo1.strChrom <= stInfo2.strChrom)
        return true;
    else
        return false;
}

void ClsParseStdCircRNAInfo::CheckIrregularCircRNA(vector<St_StdCircRNAInfo>& vInfo,
                           vector<St_Row_Chrom>& vChrom) // 将所有不规则的CircRNA都输出
{
    //here, we need to create 1 type is enough
    //我们通过判断mapping的位置，来进行揣测
    //we have some estimate here
    //1: check the length of the boundary exons
    /// 1: if length is long enough, we use those wo directly  ->
    /// we need to append one more in the bounday
    /// 每一个circRNA Info 对应一个相应的realted reference

    cout << "Total Number of Std CircRNA: " << IntToStr(vInfo.size()) << endl;

    vector<St_StdCircRNAInfo> vHitInfo;
    vector<St_StdCircRNAInfo> vUnHitInfo;


    for(vector<St_StdCircRNAInfo>::iterator itr = vInfo.begin(); itr != vInfo.end(); itr++)
    {
        //1: 找到相应的transcript
        for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); //For each chromosome
            itrChrom != vChrom.end(); itrChrom++)
        {
            if( itr->strChrom == ("chr" + itrChrom->strName) )  // fins this chromosome
            {
                bool bFindRange = false;
                for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); //For each gene
                    itrGene != itrChrom->vRG.end(); itrGene++)
                {
                    for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin(); //For each trans
                        itrRT != itrGene->vRT.end(); itrRT++)
                    {
                        bool bHitLeft = false;
                        bool bHitRight = false;
                        for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); //For each exon
                            itrExon != itrRT->vRExon.end(); itrExon++)
                        {
                            if(!bHitLeft)
                            {
                                if(abs(itrExon->iStart - itr->iStartPos) < OFFSET || abs(itrExon->iEnd - itr->iStartPos) < OFFSET)
                                {
                                    bHitLeft = true;
                                }
                            }

                            if(!bHitRight)
                            {
                                if(abs(itrExon->iStart - itr->iEndPos) < OFFSET || abs(itrExon->iEnd - itr->iEndPos) < OFFSET)
                                {
                                    bHitRight = true;
                                }
                            }

                            if(bHitLeft && bHitRight)
                            {
                                bFindRange = true;
                                vHitInfo.push_back(*itr);
                                break;
                            }
                        }

                        if(bFindRange)
                        {
                            break;
                        }
                    }
                    if(bFindRange)
                    {
                        break;
                    }
                }
                if(!bFindRange)
                {
                    vUnHitInfo.push_back(*itr);
                }
            }
        }
    }

    //其实我们不关心对不到chrmosome的info
    cout << "Hit Number                      : " << IntToStr(vHitInfo.size()) << endl;
    cout << "Un-Hit Number (could find chrom): " << IntToStr(vUnHitInfo.size()) << endl;

    ///对Hit Info 进行处理
    const char* cArryHitCaseEnum[4] = {"S_I_H", "S_I_T", "B_I_H", "B_I_T"};  // S: Small, I: Index, H:Head, T:Tail.
    int iS_I_H_B_I_H = 0;
    int iS_I_H_B_I_T = 0;
    int iS_I_T_B_I_H = 0;
    int iS_I_T_B_I_T = 0;
    int iRCNum = 0;
    int iRegularNum = 0; // not reverse complimentary

    int iS_I_H_B_I_H_Regular = 0;
    int iS_I_H_B_I_H_RC = 0;

    int iS_I_H_B_I_T_Regular = 0;
    int iS_I_H_B_I_T_RC = 0;

    int iS_I_T_B_I_H_Regular = 0;
    int iS_I_T_B_I_H_RC = 0;

    int iS_I_T_B_I_T_Regular = 0;
    int iS_I_T_B_I_T_RC = 0;


    for(vector<St_StdCircRNAInfo>::iterator itr = vHitInfo.begin(); itr != vHitInfo.end(); itr++)
    {
        for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); //For each chromosome
            itrChrom != vChrom.end(); itrChrom++)
        {
            if( itr->strChrom == ("chr" + itrChrom->strName) )  // fins this chromosome
            {
                bool bFindRange = false;
                for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); //For each gene
                    itrGene != itrChrom->vRG.end(); itrGene++)
                {
                    for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin(); //For each trans
                        itrRT != itrGene->vRT.end(); itrRT++)
                    {
                        int iHitIndex1 = -1;
                        string strPart1 = "";
                        int iHitIndex2 = -1;
                        string strPart2 = "";
                        int iIndex = 0;

                        bool bHitLeft = false;
                        bool bHitRight = false;
                        for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); //For each exon
                            itrExon != itrRT->vRExon.end(); itrExon++)
                        {
                            if(!bHitLeft)
                            {
                                if(abs(itrExon->iStart - itr->iStartPos) < OFFSET)
                                {
                                    bHitLeft = true;
                                    iHitIndex1 = iIndex;
                                    strPart1 = "Head";
                                }

                                if(abs(itrExon->iEnd - itr->iStartPos) < OFFSET)
                                {
                                    bHitLeft = true;
                                    iHitIndex1 = iIndex;
                                    strPart1 = "Tail";
                                }
                            }

                            if(!bHitRight)
                            {
                                if(abs(itrExon->iStart - itr->iEndPos) < OFFSET)
                                {
                                    bHitRight = true;
                                    iHitIndex2 = iIndex;
                                    strPart2 = "Head";
                                }

                                if(abs(itrExon->iEnd - itr->iEndPos) < OFFSET)
                                {
                                    bHitRight = true;
                                    iHitIndex2 = iIndex;
                                    strPart2 = "Tail";
                                }
                            }

                            if(bHitLeft && bHitRight)
                            {
                                bFindRange = true;
                                break;
                            }
                            iIndex++;
                        }
                        if(bFindRange)
                        {
                            if(itrRT->bRC)
                                iRCNum++;
                            else
                                iRegularNum++;

                            if(iHitIndex1 <= iHitIndex2) // Index 1 is the smaller one
                            {
                                if(strPart1 == "Head") //S_I_H
                                {
                                    if(strPart2 == "Head") //B_I_H
                                    {
                                        iS_I_H_B_I_H++;

                                        if(!itrRT->bRC)
                                            iS_I_H_B_I_H_Regular++;
                                        else
                                            iS_I_H_B_I_H_RC++;

                                    }
                                    else if(strPart2 == "Tail")//B_I_T
                                    {
                                        iS_I_H_B_I_T++;

                                        if(!itrRT->bRC)
                                            iS_I_H_B_I_T_Regular++;
                                        else
                                            iS_I_H_B_I_T_RC++;
                                    }
                                }
                                else if(strPart1 == "Tail") //  S_I_T
                                {
                                    if(strPart2 == "Head") //B_I_H
                                    {
                                        iS_I_T_B_I_H++;

                                        if(!itrRT->bRC)
                                            iS_I_T_B_I_H_Regular++;
                                        else
                                            iS_I_T_B_I_H_RC++;

                                    }
                                    else if(strPart2 == "Tail")//B_I_T
                                    {
                                        iS_I_T_B_I_T++;

                                        if(!itrRT->bRC)
                                            iS_I_T_B_I_T_Regular++;
                                        else
                                            iS_I_T_B_I_T_RC++;

                                    }
                                }

                            }
                            else // Index 1 is the bigger one
                            {
                                if(strPart1 == "Head") //B_I_H
                                {
                                    if(strPart2 == "Head") //S_I_H
                                    {
                                        iS_I_H_B_I_H++;
                                        if(!itrRT->bRC)
                                            iS_I_H_B_I_H_Regular++;
                                        else
                                            iS_I_H_B_I_H_RC++;
                                    }
                                    else if(strPart2 == "Tail")//S_I_T
                                    {
                                        iS_I_T_B_I_H++;

                                        if(!itrRT->bRC)
                                            iS_I_T_B_I_H_Regular++;
                                        else
                                            iS_I_T_B_I_H_RC++;
                                    }
                                }
                                else if(strPart1 == "Tail") //  B_I_T
                                {
                                    if(strPart2 == "Head") //S_I_H
                                    {
                                        iS_I_H_B_I_T++;
                                        if(!itrRT->bRC)
                                            iS_I_H_B_I_T_Regular++;
                                        else
                                            iS_I_H_B_I_T_RC++;
                                    }
                                    else if(strPart2 == "Tail")//S_I_T
                                    {
                                        iS_I_T_B_I_T++;

                                        if(!itrRT->bRC)
                                            iS_I_T_B_I_T_Regular++;
                                        else
                                            iS_I_T_B_I_T_RC++;
                                    }
                                }
                            }

                            break;
                        }
                    }
                    if(bFindRange)
                    {
                        break;
                    }
                }
                if(bFindRange)
                    break;
            }
        }
    }

    cout << endl << "-------" << "Total Number: " << IntToStr(vHitInfo.size()) << endl;
    cout << "Small_Index_Head___Big_Index_Head: " << IntToStr(iS_I_H_B_I_H) << endl;
    cout << '\t' << "Regular: " << IntToStr(iS_I_H_B_I_H_Regular) << endl;
    cout << '\t' << "RC     : " << IntToStr(iS_I_H_B_I_H_RC) << endl;


    cout << "Small_Index_Head___Big_Index_Tail: " << IntToStr(iS_I_H_B_I_T) << endl;
    cout << '\t' << "Regular: " << IntToStr(iS_I_H_B_I_T_Regular) << endl;
    cout << '\t' << "RC     : " << IntToStr(iS_I_H_B_I_T_RC) << endl;

    cout << "Small_Index_Tail___Big_Index_Head: " << IntToStr(iS_I_T_B_I_H) << endl;
    cout << '\t' << "Regular: " << IntToStr(iS_I_T_B_I_H_Regular) << endl;
    cout << '\t' << "RC     : " << IntToStr(iS_I_T_B_I_H_RC) << endl;

    cout << "Small_Index_Tail___Big_Index_Tail: " << IntToStr(iS_I_T_B_I_T) << endl;
    cout << '\t' << "Regular: " << IntToStr(iS_I_T_B_I_T_Regular) << endl;
    cout << '\t' << "RC     : " << IntToStr(iS_I_T_B_I_T_RC) << endl;

    cout << endl;
    cout << "RC Direction Num    : " << IntToStr(iRCNum) << endl;
    cout << "Regular Diretion Num: " << IntToStr(iRegularNum) << endl << endl;

    ///对Un-Hit Info 进行处理    
    cout << endl << "Sort Un-Hit Info" << endl << endl;
    sort(vUnHitInfo.begin(), vUnHitInfo.end(), Sort_StdCircRNAInfo); //将一样的排到一起

    ///Check if it works fine
    cout << "Un-Hit Chroms" << endl << endl;
    for(vector<St_StdCircRNAInfo>::iterator itr = vUnHitInfo.begin(); itr != vUnHitInfo.end(); itr++)
    {
        cout << itr->strChrom << endl;
    }
    cout << endl << "----------" << endl;

    vector<St_UnHitStatus> vUnHitStatus;
    St_UnHitStatus stUnHitStatus;

    for(vector<St_StdCircRNAInfo>::iterator itr = vUnHitInfo.begin(); itr != vUnHitInfo.end(); itr++)
    {
        stUnHitStatus.Clear();
        stUnHitStatus.stCircRNAInfo = *itr;

        for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); //For each chromosome
            itrChrom != vChrom.end(); itrChrom++)
        {
            if( itr->strChrom == ("chr" + itrChrom->strName) )
            {
                bool bFindGene = false;

                for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin(); //For each gene
                    itrGene != itrChrom->vRG.end(); itrGene++)
                {
                    //Check if it is in this gene
                    if(itr->iStartPos >= (itrGene->iStart - OFFSET) &&
                       itr->iEndPos <= (itrGene->iEnd + OFFSET)) // 这就证明在这个Gene里面
                    {
                        bFindGene = true;
                        stUnHitStatus.strGeneName = itrGene->strID;
                        bool bFindTranscript = false;
                        for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin(); //For each trans
                            itrRT != itrGene->vRT.end(); itrRT++)
                        {
                            //Check if it is in this transcript
                            if(itr->iStartPos >= (itrRT->iStart - OFFSET) &&
                               itr->iEndPos <= (itrRT->iEnd + OFFSET)) // 这就证明在这个Transcript 里面
                            {
                                bFindTranscript = true;
                                stUnHitStatus.strTransName = itrRT->strID;
                                bool bHitLeft = false;
                                bool bHitRight = false;

                                int iLastEnd = -1;

                                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); //For each exon
                                    itrExon != itrRT->vRExon.end(); itrExon++)
                                {
                                    //Check Exon Part
                                    //For Circular RNA Left Side
                                    if(!bHitLeft)
                                    {
                                        if( itr->iStartPos >= (itrExon->iStart - OFFSET) &&  //找到了在Exon中
                                            itr->iStartPos <= (itrExon->iEnd + OFFSET))
                                        {
                                            stUnHitStatus.enLeftStatus = uhsInExon;
                                            stUnHitStatus.iLSExonStart = itrExon->iStart;
                                            stUnHitStatus.iLSExonEnd = itrExon->iEnd;
                                            stUnHitStatus.iLSExonLen = itrExon->iEnd - itrExon->iStart + 1;
                                            stUnHitStatus.iOffsetL1 = itr->iStartPos - itrExon->iStart;
                                            stUnHitStatus.iOffsetL2 = itrExon->iEnd - itr->iStartPos;
                                            bHitLeft = true;
                                        }
                                    }

                                    if(!bHitLeft)
                                    {
                                        if(iLastEnd > 0)
                                        {
                                            int iTempStart = -1;
                                            int iTempEnd = -1;
                                            if(!itrRT->bRC ) // 如果是正向的
                                            {
                                                iTempStart = iLastEnd;
                                                iTempEnd = itrExon->iStart;
                                            }
                                            else // 如果是负向的
                                            {
                                                iTempStart = itrExon->iEnd;
                                                iTempEnd = iLastEnd;

                                            }
                                            if( itr->iStartPos >= (iTempStart - OFFSET) &&  //找到了在Intron中
                                                itr->iStartPos <= (iTempEnd + OFFSET))
                                            {
                                                stUnHitStatus.enRightStatus = uhsInIntron;
                                                stUnHitStatus.iLSItronStart = iTempStart;
                                                stUnHitStatus.iLSItronEnd = iTempEnd;
                                                stUnHitStatus.iOffsetLI1 = itr->iStartPos - iTempStart;
                                                stUnHitStatus.iOffsetLI2 = iTempEnd - itr->iStartPos;
                                                bHitLeft = true;
                                            }
                                        }
                                    }

                                    //For Circular RNA Right Side
                                    if(!bHitRight)
                                    {
                                        if( itr->iEndPos >= (itrExon->iStart - OFFSET) &&  //找到了在Exon中
                                            itr->iEndPos <= (itrExon->iEnd + OFFSET))
                                        {
                                            stUnHitStatus.enRightStatus = uhsInExon;
                                            stUnHitStatus.iRSExonStart = itrExon->iStart;
                                            stUnHitStatus.iRSExonEnd = itrExon->iEnd;
                                            stUnHitStatus.iRSExonLen = itrExon->iEnd - itrExon->iStart + 1;
                                            stUnHitStatus.iOffsetR1 = itr->iEndPos - itrExon->iStart;
                                            stUnHitStatus.iOffsetR2 = itrExon->iEnd - itr->iEndPos;
                                            bHitRight = true;
                                        }
                                    }

                                    if(!bHitRight)
                                    {
                                        if(iLastEnd > 0)
                                        {
                                            int iTempStart = -1;
                                            int iTempEnd = -1;
                                            if(!itrRT->bRC ) // 如果是正向的
                                            {
                                                iTempStart = iLastEnd;
                                                iTempEnd = itrExon->iStart;
                                            }
                                            else // 如果是负向的
                                            {
                                                iTempStart = itrExon->iEnd;
                                                iTempEnd = iLastEnd;

                                            }
                                            if( itr->iEndPos >= (iTempStart - OFFSET) &&  //找到了在Intron中
                                                itr->iEndPos <= (iTempEnd + OFFSET))
                                            {
                                                stUnHitStatus.enRightStatus = uhsInIntron;
                                                stUnHitStatus.iRSIntronStart = iTempStart;
                                                stUnHitStatus.iRSIntronEnd = iTempEnd;
                                                stUnHitStatus.iOffsetRI1 = itr->iEndPos - iTempStart;
                                                stUnHitStatus.iOffsetRI2 = iTempEnd - itr->iEndPos;
                                                bHitRight = true;
                                            }
                                        }
                                    }

                                    if(bHitLeft && bHitRight)
                                    {
                                        break;
                                    }
                                    else
                                    {
                                        if(!itrRT->bRC ) // 如果是正向的
                                        {
                                            iLastEnd = itrExon->iEnd ;
                                        }
                                        else // 如果是负向的
                                        {
                                            iLastEnd = itrExon->iStart;
                                        }
                                    }
                                }
                                break;
                            }
                            if(bFindTranscript)
                                break;
                        }
                        break;
                    }
                    if(bFindGene)
                        break;
                }
                vUnHitStatus.push_back(stUnHitStatus);
                break;
            }
        }
    }

    //输出相应的结果
    cout << endl << "---------" << endl;
    cout << "Hit Chromosome Number: " << IntToStr(vUnHitStatus.size()) << endl;
    for(vector<St_UnHitStatus>::iterator itr = vUnHitStatus.begin(); itr != vUnHitStatus.end(); itr++)
    {
        cout << itr->stCircRNAInfo.strChrom << " -- "
             << itr->strGeneName << " -- "
             << itr->strTransName << endl;
        //Left Side
        cout << '\t' << "Left_Exon  : " << cUNHITSTATUS[itr->enLeftStatus] << "("
                                                        << IntToStr(itr->iLSExonStart) << ", "
                                                        << IntToStr(itr->iLSExonEnd) << ", "
                                                        << IntToStr(itr->iLSExonLen) << ", "
                                                        << IntToStr(itr->iOffsetL1) << ", "
                                                        << IntToStr(itr->iOffsetL2) << ")" << endl;
        cout << '\t' << "Left_Intron: " << cUNHITSTATUS[itr->enLeftStatus] << "("
                                                        << IntToStr(itr->iLSItronStart) << ", "
                                                        << IntToStr(itr->iLSItronEnd) << ", "
                                                        << IntToStr(itr->iOffsetLI1) << ", "
                                                        << IntToStr(itr->iOffsetLI2) << ")" << endl;

        //Right Side
        cout << '\t' << "Right_Exon  : " << cUNHITSTATUS[itr->enRightStatus] << "("
                                                        << IntToStr(itr->iRSExonStart) << ", "
                                                        << IntToStr(itr->iRSExonEnd) << ", "
                                                        << IntToStr(itr->iRSExonLen) << ", "
                                                        << IntToStr(itr->iOffsetR1) << ", "
                                                        << IntToStr(itr->iOffsetR2) << ")" << endl;
        cout << '\t' << "Right_Intron: " << cUNHITSTATUS[itr->enRightStatus] << "("
                                                        << IntToStr(itr->iRSIntronStart) << ", "
                                                        << IntToStr(itr->iRSIntronEnd) << ", "
                                                        << IntToStr(itr->iOffsetRI1) << ", "
                                                        << IntToStr(itr->iOffsetRI2) << ")" << endl;

    }

    //we could statistic the following result:
    //1: For case Exon & Exon
    int iExon_Exon_Num = 0;
    int iEE_HitEdge_NearEdge = 0;
    int iEE_TwoNearEdge = 0;
    int iEE_Others = 0;

    //2: For Exon & Intron
    int iExon_Intron_Num = 0;
    int iEI_TwoNearEdge = 0;
    int iEI_Others = 0;

    //3: For Intron & Intron
    int iIntron_Intron_Num = 0;
    int iII_TwoNearEdge = 0;
    int iII_Others = 0;

    //4: Contrain Unknow
    int iUnknown_Num = 0;

    for(vector<St_UnHitStatus>::iterator itr = vUnHitStatus.begin(); itr != vUnHitStatus.end(); itr++)
    {
        ///For case Exon & Exon
        if(itr->enLeftStatus == uhsInExon && itr->enRightStatus == uhsInExon)
        {
            if( ((abs(itr->iOffsetL1) <= OFFSET || abs(itr->iOffsetL2) <= OFFSET) &&
                 (abs(itr->iOffsetR1) <= 200 || abs(itr->iOffsetR2) <= 200)) ||
                ((abs(itr->iOffsetL1) <= OFFSET || abs(itr->iOffsetL2) <= OFFSET) &&
                 (abs(itr->iOffsetR1) <= 200 || abs(itr->iOffsetR2) <= 200))
              )
            {
                iEE_HitEdge_NearEdge++;
            }
            else if( (abs(itr->iOffsetL1) <= 200 || abs(itr->iOffsetL2) <= 200) &&
                     (abs(itr->iOffsetR1) <= 200 || abs(itr->iOffsetR2) <= 200) )
            {
                iEE_TwoNearEdge++;
            }
            else
                iEE_Others++;

            iExon_Exon_Num++;
        }
        ///For Exon & Intron
        else if(itr->enLeftStatus == uhsInExon && itr->enRightStatus == uhsInIntron)
        {
            if( (abs(itr->iOffsetL1) <= 200 || abs(itr->iOffsetL2) <= 200) &&
                (abs(itr->iOffsetRI1) <= 200 || abs(itr->iOffsetRI2) <= 200) )
            {
                iEI_TwoNearEdge++;
            }
            else
                 iEI_Others++;
            iExon_Intron_Num++;
        }
        else if(itr->enLeftStatus == uhsInIntron && itr->enRightStatus == uhsInExon)
        {
            if( (abs(itr->iOffsetLI1) <= 200 || abs(itr->iOffsetLI2) <= 200) &&
                (abs(itr->iOffsetR1) <= 200 || abs(itr->iOffsetR2) <= 200) )
            {
                iEI_TwoNearEdge++;
            }
            else
                 iEI_Others++;
            iExon_Intron_Num++;
        }
        //For Intron & Intron
        else if(itr->enLeftStatus == uhsInIntron && itr->enRightStatus == uhsInIntron)
        {
            if( (abs(itr->iOffsetLI1) <= 200 || abs(itr->iOffsetLI2) <= 200) &&
                (abs(itr->iOffsetRI1) <= 200 || abs(itr->iOffsetRI2) <= 200) )
            {
                iII_TwoNearEdge++;
            }
            else
                iII_Others++;
            iIntron_Intron_Num++;
        }
        else
            iUnknown_Num++;
    }
    //Output Result
    cout << endl << "==============" << endl;
    ///For Exon & Exon
    cout << "Exon & Exon: " << IntToStr(iExon_Exon_Num) << " ------- " << endl;
    cout << '\t' << "One_Hit_One_Near: " << IntToStr(iEE_HitEdge_NearEdge) << endl;
    cout << '\t' << "Two_Near_Edge   : " << IntToStr(iEE_TwoNearEdge) << endl;
    cout << '\t' << "Others          : " << IntToStr(iEE_Others) << endl;

    ///For Exon & Intron
    cout << "Exon & Intron: " << IntToStr(iExon_Intron_Num) << " ------- " << endl;
    cout << '\t' << "Two_Near_Edge   : " << IntToStr(iEI_TwoNearEdge) << endl;
    cout << '\t' << "Others          : " << IntToStr(iEI_Others) << endl;

    ///For Intron & Intron
    cout << "Intron & Intron: " << IntToStr(iIntron_Intron_Num) << " ------- " << endl;
    cout << '\t' << "Two_Near_Edge   : " << IntToStr(iII_TwoNearEdge) << endl;
    cout << '\t' << "Others          : " << IntToStr(iII_Others) << endl;

    ///For Unknown
    cout << "Contain Unknown: " << IntToStr(iUnknown_Num) << " ------- " << endl;

    vInfo.clear();
    vInfo = vHitInfo;
}

//what we need to do is that we use the raw annotation file
string ClsParseStdCircRNAInfo::GenerateCircularRelatedRef(vector<St_StdCircRNAInfo>& vInfo,//这里给的info应该是两边都能够比上exon边界的
                                                          vector<St_Row_Chrom>& vChrom,
                                                          string strRefPath)
{    
    /* 我们现在这里需要进行相应reference的组合。
     * 目的: 找到相应的支持这个junction point节点的reads, 看看有多少个
     * 关于我们如何制作相应的refernce:
     * (1) 将junciton point两端的exon进行组合: 得到ref1
     * (2) 将这个junction点中间的包含的所有的部分组合起来，得到ref1_comp
     * 通过获得ref1 mapping status + the mapping status of its pair in ref1_comp,
     * to check if current reads pair support this circula RNA
     */

    //我觉得我们一样要输出这些exon，来看看的到底考不靠谱

    //1: Parse Reference
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vRef;
    pFastaReader->ReadFastaRegular(strRefPath, vRef);
    delete pFastaReader;
    pFastaReader = NULL;

    ofstream ofs;
    string strCircRef = "../TempFile/CircRef.fa";
    ofs.open(strCircRef.c_str());

    ofstream ofsExon; //Related Exon
    string strRelatedExon = "../TempFile/RelatedExon.fa";
    ofsExon.open(strRelatedExon.c_str());

    int iSelfCircCase = 0;
    int iRegCircCase = 0;

    //generate circular rna related reference data
    int iPoorRefNum = 0;
    for(vector<St_StdCircRNAInfo>::iterator itr = vInfo.begin(); itr != vInfo.end(); itr++)
    {
        for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin();
            itrChrom != vChrom.end(); itrChrom++)
        {
            if(itr->strChrom == ("chr" + itrChrom->strName)) // 找到了这个chromosome
            {
                //Get the corresponding reference
                St_Fasta* pCurRef = NULL;
                for(vector<St_Fasta>::iterator itrRef = vRef.begin(); itrRef != vRef.end(); itrRef++)
                {
                    string strChromName = itrRef->strName.substr(0, itrRef->strName.find(" "));
                    if(strChromName == itrChrom->strName)
                    {
                        pCurRef = &(*itrRef);
                    }
                }

                if(pCurRef == NULL)
                    break;

                for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin();
                    itrGene != itrChrom->vRG.end(); itrGene++)
                {
                    //check if current circrna located in this range
                    if(itr->iStartPos >= itrGene->iStart && itr->iEndPos <= itrGene->iEnd)
                    {
                        for(vector<St_Raw_Transcript>::iterator itrRT = itrGene->vRT.begin();
                            itrRT != itrGene->vRT.end(); itrRT++)
                        {
                            if(itr->iStartPos >= itrRT->iStart && itr->iEndPos <= itrRT->iEnd)
                            {
                                vector<St_Raw_Exon> vTmpExon;
                                bool bFindStart = false;
                                bool bFindEnd = false;
                                for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin();
                                    itrExon != itrRT->vRExon.end(); itrExon++)
                                {
                                    //Check the corresponding exon  -->Go!!
                                    if(!bFindStart)
                                    {
                                        if( abs(itrExon->iStart - itr->iStartPos) < OFFSET ||
                                            abs(itrExon->iEnd - itr->iStartPos) < OFFSET )
                                        {
                                            bFindStart = true;
                                        }
                                    }

                                    if(!bFindEnd)
                                    {
                                        if( abs(itrExon->iStart - itr->iEndPos) < OFFSET ||
                                            abs(itrExon->iEnd - itr->iEndPos) < OFFSET )
                                        {
                                            bFindEnd = true;
                                        }
                                    }

                                    if(bFindStart || bFindEnd)
                                    {
                                        vTmpExon.push_back(*itrExon);
                                    }

                                    if(bFindStart && bFindEnd)
                                        break;
                                }
                                //我们现在已经得到了exon的list了，现在考虑如何构建相应的reference
                                if(!CreatCircRef(vTmpExon, pCurRef->strSeq,
                                                 itrRT->bRC, itr->strName, itr->enType,
                                                 ofs, ofsExon, iSelfCircCase, iRegCircCase))
                                    iPoorRefNum++;
                                break;
                            }
                        }
                        break;
                    }
                }
                break;
            }
        }
    }
    ofs.close();
    ofsExon.close();

    cout << endl << "---------------" << endl;
    cout << "Number of Poor Ref    : " << IntToStr(iPoorRefNum) << endl;
    cout << "Number of Self Circ   : " << IntToStr(iSelfCircCase) << endl;
    cout << "Number of Regular Circ: " << IntToStr(iRegCircCase) << endl;
    cout << "---------------" << endl << endl;

    return strCircRef;
}

//we need to guarantee the total lengthof circ ref should larger than 120 bps (1.2 times of reads length)
bool ClsParseStdCircRNAInfo::CreatCircRef(vector<St_Raw_Exon>& vExon, string& strRef,
                                          bool bRC, string strCircRNAName, En_CircType& enType,
                                          ofstream& ofs, ofstream& ofsExon, int& iSelfCircCase, int& iRegCircCase)
{
    if(vExon.empty())
        return false;

    static int siRefIndex = 1;

    string strCircRef = "";
    string strCircRefCompl = ""; // this is the complementory of current reference

    string strType = "";
    int iJunctionPos = 0;
    if(vExon.size() == 1) // case 1: self circ
    {
        string strCurExon = strRef.substr(vExon.begin()->iStart - 1,
                                          vExon.begin()->iEnd - vExon.begin()->iStart + 1);
        iJunctionPos = strCurExon.length();
        strCircRef = strCurExon + strCurExon;
        strCircRefCompl = strCurExon + strCurExon;
        strType = "Self Circ";
        enType = ctSelf;
        iSelfCircCase++;
    }
    else
    {
        //treat regular case (+)
        if(!bRC)  // only consider begin and end here
        {
            string strCurExonBegin = strRef.substr(vExon.begin()->iStart - 1,
                                                   vExon.begin()->iEnd - vExon.begin()->iStart + 1);
            string strCurExonEnd = strRef.substr((vExon.end()-1)->iStart - 1,
                                                 (vExon.end()-1)->iEnd - (vExon.end()-1)->iStart + 1);
            iJunctionPos = strCurExonEnd.length();
            strCircRef = strCurExonEnd + strCurExonBegin;

            //Create the complementary ref seq
            string strWholeSeq = strRef.substr(vExon.begin()->iStart - 1,
                                               (vExon.end()-1)->iEnd - vExon.begin()->iStart + 1);
            strCircRefCompl = strCurExonEnd + strWholeSeq;
        }
        else // reverse complemantory
        {
            vector<St_Raw_Exon>::iterator itrLeftExon = vExon.begin()->iEnd <= (vExon.end()-1)->iStart ?
                                                        vExon.begin():(vExon.end()-1);
            vector<St_Raw_Exon>::iterator itrRightExon = vExon.begin()->iEnd > (vExon.end()-1)->iStart ?
                                                        vExon.begin():(vExon.end()-1);
            string strCurExonBegin = strRef.substr(itrLeftExon->iStart - 1,
                                                   itrLeftExon->iEnd - itrLeftExon->iStart + 1);
            string strCurExonEnd = strRef.substr(itrRightExon->iStart - 1,
                                                 itrRightExon->iEnd - itrRightExon->iStart + 1);
            iJunctionPos = strCurExonEnd.length();
            strCircRef = strCurExonEnd + strCurExonBegin; //it's correct, still from big to small

            //Create the complementary ref seq
            string strWholeSeq = strRef.substr(itrLeftExon->iStart,
                                               itrRightExon->iEnd - itrLeftExon->iStart + 1);
            strCircRefCompl = strCurExonEnd + strWholeSeq;
        }
        strType = "Regular Case";
        enType = ctRegular;
        iRegCircCase++;
    }

    //For circ ref
    string strRefName = ">" + strCircRNAName + "_std" + " " + IntToStr(iJunctionPos) + " " + (bRC ? "-" : "+");
    ofs << strRefName << endl;
    ofs << strCircRef << endl;
    //For circ ref complementary
    string strRefComplName = ">" + strCircRNAName + "_compl" + " " + IntToStr(iJunctionPos) + " " +  (bRC ? "-" : "+");
    ofs << strRefComplName << endl;
    ofs << strCircRefCompl << endl;

    //For Related Exon  
    ofsExon << strRefName << endl;
    if(vExon.size() == 1)
    {
        string strCurExon = strRef.substr(vExon.begin()->iStart - 1,
                                          vExon.begin()->iEnd - vExon.begin()->iStart + 1);
        ofsExon << strCurExon << endl;
    }
    else
    {
        string strCurExonBegin = strRef.substr(vExon.begin()->iStart - 1,
                                               vExon.begin()->iEnd - vExon.begin()->iStart + 1);
        string strCurExonEnd = strRef.substr((vExon.end()-1)->iStart - 1,
                                             (vExon.end()-1)->iEnd - (vExon.end()-1)->iStart + 1);
        ofsExon << strCurExonBegin << endl;
        ofsExon << "-------" << endl;
        ofsExon << strCurExonEnd << endl;
    }


    siRefIndex++;

    if(strCircRef.length() < 101)
    {
        cout << "<" << IntToStr(vExon.begin()->iStart) << ", " << IntToStr(vExon.begin()->iEnd) << ">"
             << " " << IntToStr(abs(vExon.begin()->iStart - vExon.begin()->iEnd) + 1) << " --- "
             << "<" << IntToStr((vExon.end()-1)->iStart) << ", " << IntToStr((vExon.end()-1)->iEnd) << ">"
             << " " << IntToStr(abs((vExon.end()-1)->iStart - (vExon.end()-1)->iEnd) + 1)
             <<  " ------ " << strType << " " << IntToStr(vExon.size()) << " "
             << (bRC ? "-" : "+") << ": " << IntToStr(strCircRef.length()) << endl;
        return false;
    }
    else
        return true;
}

void ClsParseStdCircRNAInfo::GetValidStdCircRNAByChrom(vector<St_StdCircRNAInfo>& vTargetValidStdInfo,
                                                       vector<St_Row_Chrom>& vChrom,
                                                       string strDbPath, int iChromIndex)
{
    vTargetValidStdInfo.clear();

    //Read the whole info
    vector<St_StdCircRNAInfo> vInfo;
    this->ParseCircRNADb(strDbPath, vInfo);
    this->FilterByCell(vInfo, m_enCellType);

    if(vInfo.empty())
    {
        cout << "circRNADb is Empty!" << endl;
        return;
    }

    cout << "circRNADb Total Size: " << IntToStr(vInfo.size()) << endl;

    //Filter the un-related std circRNA  //--> 自己来写新的逻辑
    for(vector<St_StdCircRNAInfo>::iterator itr = vInfo.begin(); itr != vInfo.end(); itr++)
    {
        if(itr->strChrom != "chr" + IntToStr(iChromIndex))
        {
            continue;
        }

        //For current chromosome
        //Check if the boundary could be hit with the boundary of chromosome
        for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end(); itrChrom++)
        {
            if(itrChrom->strName != IntToStr(iChromIndex))
            {
                continue;
            }

            //cout << "Find Chrom!" << endl;

            bool bHitBoth = false;

            for(vector<St_Raw_Gene>::iterator itrGene = itrChrom->vRG.begin();
                itrGene != itrChrom->vRG.end(); itrGene++)
            {
                for(vector<St_Raw_Transcript>::iterator itrTranscript = itrGene->vRT.begin();
                    itrTranscript != itrGene->vRT.end(); itrTranscript++)
                {
                    bool bHitStart = false;
                    int iHitStartIndex = -1;
                    bool bHitEnd = false;
                    int iHitEndIndex = -1;
                    int iIndex = 0;
                    for(vector<St_Raw_Exon>::iterator itrExon = itrTranscript->vRExon.begin();
                        itrExon != itrTranscript->vRExon.end(); itrExon++)
                    {
                        //For Start
                        if(!bHitStart)
                        {
                            if(abs(itrExon->iStart - itr->iStartPos) <= OFFSET ||
                               abs(itrExon->iEnd - itr->iStartPos) <= OFFSET)
                            {
                                bHitStart = true;
                                iHitStartIndex = iIndex;
                            }
                        }

                        //For End
                        if(!bHitEnd)
                        {
                            if(abs(itrExon->iStart - itr->iEndPos) <= OFFSET ||
                               abs(itrExon->iEnd - itr->iEndPos) <= OFFSET)
                            {
                                bHitEnd = true;
                                iHitEndIndex = iIndex;
                            }
                        }

                        if(bHitStart && bHitEnd)
                        {
                            bHitBoth = true;
                            //Set Value
                            if(iHitStartIndex == iHitEndIndex)
                                itr->enType = ctSelf;
                            else
                                itr->enType = ctRegular;

                            vTargetValidStdInfo.push_back(*itr);
                            break;
                        }
                        iIndex++;
                    }
                    if(bHitBoth)
                        break;
                }
                if(bHitBoth)
                    break;
            }
            break;
        }
    }
}

void ClsParseStdCircRNAInfo::PrintStdHitInfoByBriefCandiStyle(vector<St_StdCircRNAInfo>& vTargetValidStdInfo,
                                                              vector<St_Row_Chrom>& vChrom,
                                                              string strDbPath, int iChromIndex)
{
    //Get the valid std target circRNA
    vTargetValidStdInfo.clear();
    GetValidStdCircRNAByChrom(vTargetValidStdInfo, vChrom, strDbPath, iChromIndex);
    cout << IntToStr(vTargetValidStdInfo.size()) << endl;

    ///Print Out separately: self, reg, and combine together
    ofstream ofsBriefSelf; // 这个是为了解析起来比较方便的文件
    ofsBriefSelf.open("./CircRNAdb_BriefSelf.txt");

    ofstream ofsBriefReg; // 这个是为了解析起来比较方便的文件
    ofsBriefReg.open("./CircRNAdb_BriefReg.txt");

    ofstream ofsBriefTotal; // 这个是为了解析起来比较方便的文件
    ofsBriefTotal.open("./CircRNAdb_BriefTotal.txt");

    for(vector<St_StdCircRNAInfo>::iterator itr = vTargetValidStdInfo.begin();
        itr != vTargetValidStdInfo.end(); itr++)
    {
        switch (itr->enType)
        {
            case ctSelf:
            {
                ofsBriefSelf << itr->strChrom << " "
                             << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                             << IntToStr(abs(itr->iStartPos - itr->iEndPos)) << " "
                             << "xx" << " "
                             << "S" << endl;

                ofsBriefTotal << itr->strChrom << " "
                              << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                              << IntToStr(abs(itr->iStartPos - itr->iEndPos)) << " "
                              << "xx" << " "
                              << "S" << endl;
                break;
            }
            case ctRegular:
            {
                ofsBriefReg << itr->strChrom << " "
                            << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                            << IntToStr(abs(itr->iStartPos - itr->iEndPos)) << " "
                            << "xx" << " "
                            << "R" << endl;

                ofsBriefTotal << itr->strChrom << " "
                              << IntToStr(itr->iStartPos) << " " << IntToStr(itr->iEndPos) << " "
                              << IntToStr(abs(itr->iStartPos - itr->iEndPos)) << " "
                              << "xx" << " "
                              << "R" << endl;
                break;
            }
            default:
                break;
        }
    }

    ofsBriefSelf.close();
    ofsBriefReg.close();
    ofsBriefTotal.close();

    ///Temporary: Only Keep the self circular case
    for(vector<St_StdCircRNAInfo>::iterator itr = vTargetValidStdInfo.end() - 1;
        itr >= vTargetValidStdInfo.begin(); itr--)
    {
        if(itr->enType == ctSelf)
        {} //Kept
        else
        {
            vTargetValidStdInfo.erase(itr);
        }
    }

    cout << "Final Valid Std CircRNA Size: " << IntToStr(vTargetValidStdInfo.size()) << endl;
}

