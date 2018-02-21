#include "clscomparebasefunction.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>

const int BOUNDARYOFFSET = 10; //我们发现有的offset差那么几个bp
const char* arryCellName[ctMaxx] = {"H9", "glioblastoma", "brain"};

ClsCompareBaseFunction::ClsCompareBaseFunction()
{
}

ClsCompareBaseFunction::~ClsCompareBaseFunction()
{}

void ClsCompareBaseFunction::Init(int iReadLen, string strTissue)
{
    m_iReadLen = iReadLen;

    //Initiate m_enCellType
    if(strTissue == arryCellName[ctH9])
        m_enCellType = ctH9;
    else if(strTissue == arryCellName[ctGlioblastoma])
        m_enCellType = ctGlioblastoma;
    else if(strTissue == arryCellName[ctBrain])
        m_enCellType = ctBrain;
    else
        m_enCellType = ctMaxx;
}

///---------------Base Parse Function-----------------
/// --------------------------------------------------
/// --------------------------------------------------
// Brief.txt (same directory as executable file)
void ClsCompareBaseFunction::ParseBrief(string strCurPath, vector<St_Candidate>& vCandi)
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

        if(strLine == "" || strLine == "\n")
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

// Standard Result  --------->Do it tomorrow
// we only do the testing based on chr1
void ClsCompareBaseFunction::ParseSimulationResult(string strSimulatePath,
                                                   vector<St_Candidate>& vCandi)
{
    ifstream ifs;
    ifs.open(strSimulatePath.c_str());
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "" || strLine == "\n")
            continue;

        //Get Chromosone Name
        int iStartPos = 0;
        int iEndPos = strLine.find(':', iStartPos);
        int iLen = iEndPos - iStartPos;
        string strChromName = strLine.substr(iStartPos, iLen);
        stCandi.strChromName = strChromName;

        //Get First Value
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('|', iStartPos);
        iLen = iEndPos - iStartPos;
        int iFirstValue = atoi(strLine.substr(iStartPos, iLen).c_str());

        //Get Second Value
        iStartPos = iEndPos + 1;
        iEndPos = strLine.length();
        iLen = iEndPos - iStartPos;
        int iSecondValue = atoi(strLine.substr(iStartPos, iLen).c_str());

        stCandi.iStartPos = iFirstValue <= iSecondValue ? iFirstValue : iSecondValue;
        stCandi.iEndPos = iFirstValue >= iSecondValue ? iFirstValue : iSecondValue;

        //在此处将重复的进行排除
        bool bExisted = false;
        for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
        {
            if( *itr == stCandi)
            {
                bExisted = true;
                break;
            }
        }

        if(!bExisted) // only record it when it is not existed
            vCandi.push_back(stCandi);
    }
    ifs.close();
}

//CircBase Standard Circular RNA
void ClsCompareBaseFunction::ParseCircBaseResult(string strCircBasePath, vector<St_Candidate>& vCandi)
{
    //Parse CircBase Reuslt
    /* File Structure
     * 0: chromosone name
     * 1: start
     * 2: end
     *    (在这里好像都是从小到大的)
     * 3: direction: + or -
     */
    ifstream ifs;
    ifs.open(strCircBasePath.c_str());
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "" || strLine == "\n")
            continue;

        if(strLine.length() > 0 && strLine.substr(0, 1) == "#")
            continue;

        int iStartPos = 0;
        int iEndPos = strLine.find('\t', iStartPos);
        int iLen = iEndPos - iStartPos;

        //The 0 --> chromosone name
        string strChrName = strLine.substr(iStartPos, iLen);
        iStartPos = 3;
        iLen = strChrName.length() - 3; //
        //注意，这里我们还是从0开始，这样我们可以直接通过下标去找相应的chromosone里面的exons的信息
        stCandi.strChromName = strChrName.substr(iStartPos, iLen); //- 1; //“减1”是为了到时候直接取值相应的exon更加方便

        //The 1 --> Start position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stCandi.iStartPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 2 --> End position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stCandi.iEndPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 3 --> Direction
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strTag = strLine.substr(iStartPos, iLen);
        if(strTag == "-")
            stCandi.bRC = true;
        else
            stCandi.bRC = false;

        //在此处将重复的进行排除
        bool bExisted = false;
        for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
        {
            if( *itr == stCandi)
            {
                bExisted = true;
                break;
            }
        }

        if(!bExisted) // only record it when it is not existed
            vCandi.push_back(stCandi);

    }
    ifs.close();
}

// CIRI Result
void ClsCompareBaseFunction::ParseCIRIResult(string strCIRIPath, vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(strCIRIPath.c_str(), ios::in);
    string strLine = "";
    St_Candidate stCandi;
    while(!ifs.eof()) // check if reached the end
    {
        getline(ifs, strLine);

        if(strLine == "" || strLine == "\n")
            continue;

        int iPos = strLine.find('\t');
        if(iPos < 0)
            continue;

        string strValue = strLine.substr(0, iPos);
        int iPosFirstSplit = strValue.find(':');
        int iPosSecondSplit = strValue.find('|');
        if(iPosFirstSplit < 0 || iPosSecondSplit < 0)
            continue;

        //Get Real value
        stCandi.strChromName = strValue.substr(0, iPosFirstSplit); //- 1;  //这里故意减去1，使得起始位置为0
        stCandi.iStartPos = atoi(strValue.substr(iPosFirstSplit+1,
                                                  iPosSecondSplit - iPosFirstSplit - 1).c_str());
        stCandi.iEndPos = atoi(strValue.substr(iPosSecondSplit + 1, strValue.length() -
                                                                iPosSecondSplit - 1).c_str());
        if(strLine.find('+') >= 0)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();

    //Delete Duplicate
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCandi.end() - 1; subItr > itr; subItr--)
        {
            if(*itr == *subItr)
                vCandi.erase(subItr);
        }
    }
}



//Find_Circ Result
void ClsCompareBaseFunction::ParseFindCircResult(string strFindCircPath, vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(strFindCircPath.c_str());
    string strLine = "";
    St_Candidate stCandi;
    while(!ifs.eof())
    {
        getline(ifs, strLine);

        if(strLine == "" || strLine == "\n")
            continue;

        int iTagPos = strLine.find('\t');
        if(iTagPos < 0)
            continue;

        //1: For chromosome
        int iStart = 0;
        int iLen = iTagPos - iStart;
        stCandi.strChromName = strLine.substr(iStart, iLen); //- 1; //注意这里我们没有-1 -->We need -1 to coincide with others

        //2: For Start Pos
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iStartPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For the second value
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iEndPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For Direction
        if(strLine.find('+', iStart) != string::npos)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();

    //Delete Duplicate
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCandi.end() - 1; subItr > itr; subItr--)
        {
            if(*itr == *subItr)
                vCandi.erase(subItr);
        }
    }
}

//CIRCexplorer Result
void ClsCompareBaseFunction::ParseCircExplorerResult(string strCIRCexplorerPath,
                                                     vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(strCIRCexplorerPath.c_str());
    string strLine = "";
    St_Candidate stCandi;

    while(!ifs.eof())
    {
        getline(ifs, strLine);

        if(strLine == "" || strLine == "\n")
            continue;

        int iTagPos = strLine.find('\t');
        if(iTagPos < 0)
            continue;

        //1: For chromosome
        int iStart = 0;
        int iLen = iTagPos - iStart;
        stCandi.strChromName = strLine.substr(iStart, iLen); // -1 //注意这里我们没有-1

        //2: For Start Pos
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iStartPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For the second value
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iEndPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For Direction
        if(strLine.find('+', iStart) != string::npos)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();

    //Delete Duplicate
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCandi.end() - 1; subItr > itr; subItr--)
        {
            if(*itr == *subItr)
                vCandi.erase(subItr);
        }
    }
}

void ClsCompareBaseFunction::ParseCircRNAFinderResult(string CircRNAFinderResult,
                                                      vector<St_Candidate>& vCandi)
{
    vCandi.clear();
    ifstream ifs;
    ifs.open(CircRNAFinderResult.c_str());
    string strLine = "";
    St_Candidate stCandi;

    while(!ifs.eof())
    {
        getline(ifs, strLine);

        if(strLine == "" || strLine == "\n")
            continue;

        int iTagPos = strLine.find('\t');
        if(iTagPos < 0)
            continue;

        //1: For chromosome
        int iStart = 0;
        int iLen = iTagPos - iStart;
        stCandi.strChromName = strLine.substr(iStart, iLen); //- 1; //注意这里我们没有-1

        //2: For Start Pos
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iStartPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For the second value
        iStart = iTagPos + 1;
        iTagPos = strLine.find('\t', iStart);
        iLen = iTagPos  - iStart;
        stCandi.iEndPos = atoi(strLine.substr(iStart, iLen).c_str());

        //For Direction
        if(strLine.find('+', iStart) != string::npos)
            stCandi.bRC = false;
        else
            stCandi.bRC = true;

        vCandi.push_back(stCandi);
    }
    ifs.close();

    //Delete Duplicate
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vCandi.end() - 1; subItr > itr; subItr--)
        {
            if(*itr == *subItr)
                vCandi.erase(subItr);
        }
    }
}

//Sub-Functions
//Return Value:
//   0: The chromosome does not existed
//   1: Both Hit
//  -1: Only one hit
//  -2: None Hit
int ClsCompareBaseFunction::CheckBothHitExonBoundary(vector<St_Row_Chrom>& vChrom, St_Candidate* pCandi)
{
    St_Row_Chrom *pRawChrom = NULL;
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin();
        itrChrom != vChrom.end(); itrChrom++)
    {
        if(pCandi->strChromName == itrChrom->strName)  // here, we just need to add 1 to coincide the rule of index
        {
            pRawChrom = &(*itrChrom);
            break;
        }
    }

    if(pRawChrom == NULL)
        return 0;

    bool bHitBothStartAndEnd = false;

    bool bFindStart = false;
    bool bFindEnd = false;

    for(vector<St_Raw_Gene>::iterator itrRG = pRawChrom->vRG.begin();
        itrRG != pRawChrom->vRG.end(); itrRG++)
    {
        for(vector<St_Raw_Transcript>::iterator itrRT = itrRG->vRT.begin();
            itrRT != itrRG->vRT.end(); itrRT++)
        {
            bFindStart = false;
            bFindEnd = false;
            for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); //
                itrExon != itrRT->vRExon.end(); itrExon++)
            {
                if(!bFindStart)
                {
                    if(abs(pCandi->iStartPos - itrExon->iStart) <= BOUNDARYOFFSET)
                        bFindStart = true;
                    else if(abs(pCandi->iStartPos - itrExon->iEnd) <= BOUNDARYOFFSET)
                        bFindStart = true;
                }

                if(!bFindEnd)
                {
                    if(abs(pCandi->iEndPos - itrExon->iStart) <= BOUNDARYOFFSET)
                        bFindEnd = true;
                    else if(abs(pCandi->iEndPos - itrExon->iEnd) <= BOUNDARYOFFSET)
                        bFindEnd = true;
                }

                if(bFindStart && bFindEnd)
                {
                    bHitBothStartAndEnd = true;
                    break;
                }
            }
            if(bHitBothStartAndEnd)
                break;
        }
        if(bHitBothStartAndEnd)
            break;
    }

    if(bHitBothStartAndEnd)
        return 1;
    else
    {
        if(!bFindStart && !bFindEnd)
            return -2; // Both of them are fail to be hit
        else
            return -1;
    }
}


//General function: compare two candidate
bool ClsCompareBaseFunction::CompareTwoCandi(St_Candidate* pCandi1, St_Candidate* pCandi2)
{
    bool bFind = false;
    int iStartCur = -1;
    int iEndCur = -1;
    if(pCandi1->iStartPos <= pCandi1->iEndPos)
    {
        iStartCur = pCandi1->iStartPos;
        iEndCur = pCandi1->iEndPos;
    }
    else
    {
        iStartCur = pCandi1->iEndPos;
        iEndCur = pCandi1->iStartPos;
    }

    int iStartStd = -1;
    int iEndStd = -1;

    if(pCandi2->iStartPos <= pCandi2->iEndPos)
    {
        iStartStd = pCandi2->iStartPos;
        iEndStd = pCandi2->iEndPos;
    }
    else
    {
        iStartStd = pCandi2->iEndPos;
        iEndStd = pCandi2->iStartPos;
    }

    //Make the comparison
    if(abs(iStartCur - iStartStd) <= BOUNDARYOFFSET &&
       abs(iEndCur - iEndStd) <= BOUNDARYOFFSET)
    {
        bFind = true;
    }
    return bFind;
}

//这里我们还需要判断一下，exon的长短分布问题
bool ClsCompareBaseFunction::CheckIfBoundaryFromSameExon(vector<St_Row_Chrom>& vChrom, St_Candidate* pCandi,
                                                         bool& bShortExon)
{
    St_Row_Chrom *pRawChrom = NULL;
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin();
        itrChrom != vChrom.end(); itrChrom++)
    {
        if(pCandi->strChromName == itrChrom->strName)  // here, we just need to add 1 to coincide the rule of index
        {
            pRawChrom = &(*itrChrom);
            break;
        }
    }

    if(pRawChrom == NULL)
        return false;

    bool bComeFromOneExon = false;

    for(vector<St_Raw_Gene>::iterator itrRG = pRawChrom->vRG.begin();
        itrRG != pRawChrom->vRG.end(); itrRG++)
    {
        for(vector<St_Raw_Transcript>::iterator itrRT = itrRG->vRT.begin();
            itrRT != itrRG->vRT.end(); itrRT++)
        {
            for(vector<St_Raw_Exon>::iterator itrExon = itrRT->vRExon.begin(); //
                itrExon != itrRT->vRExon.end(); itrExon++)
            {
                if( (abs(pCandi->iStartPos - itrExon->iStart) <= BOUNDARYOFFSET ||
                     abs(pCandi->iStartPos - itrExon->iEnd) <= BOUNDARYOFFSET ) &&
                    (abs(pCandi->iEndPos - itrExon->iStart) <= BOUNDARYOFFSET ||
                     abs(pCandi->iEndPos - itrExon->iEnd) <= BOUNDARYOFFSET )
                  )
                {
                    if(itrExon->GetLength() <= m_iReadLen)
                    {
                        bShortExon = true;
                    }
                    bComeFromOneExon = true;
                    break;
                }
            }
            if(bComeFromOneExon)
                break;
        }
        if(bComeFromOneExon)
            break;
    }

    if(bComeFromOneExon)
        return true;
    else
        return false;
}

//Set the tag for each candidate:
// 1) self circle
// 2) regular circle
// 3) other case (one side comes from the middle of exon or comes from itron)
//组合起来来实现这个功能，我们不去强行写成一个函数！！
void ClsCompareBaseFunction::SetCircTag(vector<St_Row_Chrom>& vChrom, vector<St_Candidate>& vCandi)
{
    int iBothHit_Self = 0;
    int iShortLenExon = 0;
    int iRegLenExon = 0;

    int iBothHit_Reg = 0;

    int iOneHit = 0;
    int iNoneHit = 0;

    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        int iHitStatus = CheckBothHitExonBoundary(vChrom, &(*itr));
        //cout << IntToStr(iHitStatus) << endl;
        switch(iHitStatus)
        {
            case 0:
                break;
            case 1:  // this means both side hit the real boundary of exon
            {
                bool bShorExon = false;
                bool bFromSameExon = CheckIfBoundaryFromSameExon(vChrom, &(*itr), bShorExon);
                if(bFromSameExon)
                {
                    itr->enType = ctSelf;
                    iBothHit_Self++;
                    if(bShorExon)
                        iShortLenExon++;
                    else
                        iRegLenExon++;
                }
                else
                {
                    itr->enType = ctRegular;
                    iBothHit_Reg++;
                }
                break;
            }
            case -1:
                iOneHit++;
                break;
            case -2:
                iNoneHit++;
                break;
        }
    }
    
/*
    cout << "**********************************" << endl;
    cout << "*******Hit Raw Chrom Result*******" << endl;
    cout << "**********************************" << endl;
    cout << "BothHit_Self   : " << IntToStr(iBothHit_Self) << "\t"
         << GetRatio((float)iBothHit_Self / vCandi.size()) << endl;

    cout << "   ShortLenExon: " << IntToStr(iShortLenExon) << "\t"
         << GetRatio((float)iShortLenExon / iBothHit_Self) << endl;

    cout << "   RegLenExon  : " << IntToStr(iRegLenExon) << "\t"
         << GetRatio((float)iRegLenExon / iBothHit_Self) << endl;

    cout << "BothHit_Reg    : " << IntToStr(iBothHit_Reg) << "\t"
         << GetRatio((float)iBothHit_Reg / vCandi.size()) << endl;

    cout << "OneHit         : " << IntToStr(iOneHit) << "\t"
         << GetRatio((float)iOneHit / vCandi.size()) << endl;

    cout << "NoneHit        : " << IntToStr(iNoneHit) << "\t"
         << GetRatio((float)iNoneHit / vCandi.size()) << endl << endl;
*/
}

void ClsCompareBaseFunction::GetCandiForSpeciChrom(vector<St_Row_Chrom>& vChrom,
                                                   vector<St_Candidate>& vCandiAll,
                                                   vector<St_Candidate>& vCandiForSpeciChrom,
                                                   string strChromName)
{
    for(vector<St_Candidate>::iterator itr = vCandiAll.begin(); itr != vCandiAll.end(); itr++)
    {
        if(itr->strChromName != strChromName) //我们在这里只考虑chromosome1的standard candidate
            continue;
        int iCheckResult = CheckBothHitExonBoundary(vChrom, &(*itr));

        if(iCheckResult == 1) // both map
        {
            vCandiForSpeciChrom.push_back(*itr);
        }
    }
}

void ClsCompareBaseFunction::SaveCandi(vector<St_Candidate>& vCandi, string strCandiPath)
{
    ofstream ofs;
    ofs.open(strCandiPath.c_str());
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        ofs << itr->strChromName << " "
            << IntToStr(itr->iStartPos) << " "
            << IntToStr(itr->iEndPos) << " "
            << IntToStr(itr->iSupportNum) << " "
            << "xx" << " "
            << itr->GetTypeString() << endl;
    }
    ofs.close();
}

void ClsCompareBaseFunction::SaveHitCandi(vector<St_Candidate>& vCandi, string strHitPath)
{
    ofstream ofs;
    ofs.open(strHitPath.c_str());
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        switch(itr->enType)
        {
            case ctSelf: //we only record the self circular case currently.
            {
                ofs << itr->strChromName << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << IntToStr(itr->iSupportNum) << " "
                    << itr->strTag << " "
                    << itr->GetTypeString() << endl;
                break;
            }
            case ctRegular: //we do not output this case currently, only focus on how improve the regular circ case !!!
            {
                ofs << itr->strChromName << " "
                    << IntToStr(itr->iStartPos) << " "
                    << IntToStr(itr->iEndPos) << " "
                    << IntToStr(itr->iSupportNum) << " "
                    << "xx" << " "
                    << itr->GetTypeString() << endl;
                break;
            }
            default:
            {
//                ofs << IntToStr(itr->ucChromIndex) << " "
//                    << IntToStr(itr->iStartPos) << " "
//                    << IntToStr(itr->iEndPos) << " "
//                    << IntToStr(itr->iSupportNum) << " "
//                    << "M" << " "
//                    << itr->GetTypeString() << endl;
                break;
            }
        }
    }
    ofs.close();
}

string ClsCompareBaseFunction::CreateConsensusBanchmark(string strMyProgram, string strCIRI, string strFind_Circ, string strCIRCExplorer)
{
    //1: Parse My Program
    vector<St_Candidate> vMyProgram;
    ParseBrief(strMyProgram, vMyProgram);
    cout << IntToStr(vMyProgram.size()) << endl;

    //2: Parse CIRI
    vector<St_Candidate> vCIRI;
    ParseBrief(strCIRI, vCIRI);
    cout << IntToStr(vCIRI.size()) << endl;

    //3: Parse Find_Circ
    vector<St_Candidate> vFindCirc;
    ParseBrief(strFind_Circ, vFindCirc);
    cout << IntToStr(vFindCirc.size()) << endl;

    //4: Parse CIRCExplorer
    vector<St_Candidate> vCIRCExplorer;
    ParseBrief(strCIRCExplorer, vCIRCExplorer);
    cout << IntToStr(vCIRCExplorer.size()) << endl;

    //5: Create Consensus Banchmark
    vector<St_CandiSoftwareSupport> vSum;

    CandiMerge(vMyProgram, vSum);
    CandiMerge(vCIRI, vSum);
    CandiMerge(vFindCirc, vSum);
    CandiMerge(vCIRCExplorer, vSum);

    cout << "Sum Size ---> " << IntToStr(vSum.size()) << endl;

    vector<St_Candidate> vConsensusCandi;
    for(vector<St_CandiSoftwareSupport>::iterator itr = vSum.begin(); itr != vSum.end(); itr++)
    {
        if(itr->iSoftwareCount >= 2)
        {
            vConsensusCandi.push_back(itr->stCandi);
        }
    }

    string strConsensusBanchmark = "./ConsensusBanchmark.txt";
    SaveHitCandi(vConsensusCandi, strConsensusBanchmark);

    return strConsensusBanchmark;
}

void ClsCompareBaseFunction::CandiMerge(vector<St_Candidate>& vCandi, vector<St_CandiSoftwareSupport>& vSum)
{
    for(vector<St_Candidate>::iterator itr = vCandi.begin(); itr != vCandi.end(); itr++)
    {
        bool bFind = false;
        if(!vSum.empty())
        {
            for(vector<St_CandiSoftwareSupport>::iterator subItr = vSum.begin(); subItr != vSum.end(); subItr++)
            {
                if( CompareTwoCandi(&(*itr), &(subItr->stCandi)) )
                {
                    subItr->stCandi.iSupportNum += itr->iSupportNum;
                    subItr->iSoftwareCount++;
                    bFind = true;
                    break;
                }
            }
        }

        if(!bFind)
        {
            vSum.push_back(St_CandiSoftwareSupport(*itr, 1));
        }
    }
}

void ClsCompareBaseFunction::ParseCircRNADbBrief(string strFilePath, vector<St_Candidate>& vCandi)
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
    vCandi.clear();
    string strLine = "";
    St_Candidate stCandi; //candidate atom

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
        //string strCircRNAName = strLine.substr(iStartPos, iLen);
        //stInfo.strName = strCircRNAName;

        //The 1 --> Chromosome Index
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strChromName = strLine.substr(iStartPos, iLen);
        iStartPos = 3;
        iLen = strChromName.length() - 3; //
        //注意，这里我们还是从0开始，这样我们可以直接通过下标去找相应的chromosone里面的exons的信息
        stCandi.strChromName = strChromName.substr(iStartPos, iLen); //- 1;
        ToUpper(stCandi.strChromName);
        //stInfo.strChrom = strChromName;

        //The 2 --> Start position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stCandi.iStartPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 3 --> End position
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        stCandi.iEndPos = atoi(strLine.substr(iStartPos, iLen).c_str());

        //The 4 --> Direction
        iStartPos = iEndPos + 1;
        iEndPos = strLine.find('\t', iStartPos);
        iLen = iEndPos - iStartPos;
        string strTag = strLine.substr(iStartPos, iLen);
        if(strTag == "-")
            stCandi.bRC = true;
        else
            stCandi.bRC = false;

        //The -2 --> where is this circRNA come from (cell types)
        iEndPos = strLine.rfind('\t');
        iStartPos = strLine.rfind('\t', iEndPos-1) + 1;
        iLen = iEndPos - iStartPos;
        string strCellList = strLine.substr(iStartPos, iLen);
        vector<string> vOrgCell;
        //use " " to split different cells
        iStartPos = 0;
        iLen = 0;
        for(int i=0; i<(int)strCellList.length(); i++)
        {
            if(strCellList.substr(i, 1) != ",")
                iLen++;
            else
            {
                vOrgCell.push_back(strCellList.substr(iStartPos, iLen));
                iStartPos = i + 1;
                iLen = 0;
            }

            if(i == (int)strCellList.length() - 1) // the last letter
            {
                vOrgCell.push_back(strCellList.substr(iStartPos, iLen));
            }
        }

        bool bFindTissue = false;
        for(vector<string>::iterator itr = vOrgCell.begin(); itr != vOrgCell.end(); itr++)
        {
            if(itr->find(arryCellName[m_enCellType]) != string::npos)
            {
                bFindTissue = true;
                break;
            }
        }

        if(bFindTissue)
        {
            vCandi.push_back(stCandi);
        }
    }
    ifs.close();
}
