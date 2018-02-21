#include "clsreadsmapping.h"
#include "../../ShareLibrary/clsbwa.h"
#include "unistd.h" // this is for: get_current_dir_name
#include "../../ShareLibrary/clsbasealgorithm.h"
#include "../../ShareLibrary/clsfastareader.h"
#include "../../ShareLibrary/bamtools/include/api/BamReader.h"

using namespace BamTools;

ClsReadsMapping::ClsReadsMapping()
{
}

struct St_CircRefMapStatus
{
    string strRefSeq;
    string strName;
    bool bCompl;
    int iRefID;
    bool bRC; // orientation
    int iKeyPos; // the position of junciton point
    vector<string> vReads;    
};

void ClsReadsMapping::MapReadsToCircRNADb(string strRefPath, string strReads1Path,
                                          string strReads2Path, En_ReadsType enReadsType)
{
    //Check if the sorted bam file exsited
    string strRootPath = GetHigherFolderPath(get_current_dir_name());
    string strSortedBamFile = strRootPath + "TempFile/Read.sorted.bam";
//    if(access(strSortedBamFile.c_str(), 0) == 0)
//    {
//        cout << "\"Read.sorted.bam\" already existed !" << endl;
//    }
//    else
    {
        //Remove everything:
        string strRemoveTempFiles= "rm " + strRootPath + "TempFile/CircRef.fa.amb "
                                         + strRootPath + "TempFile/CircRef.fa.ann "
                                         + strRootPath + "TempFile/CircRef.fa.bwt "
                                         + strRootPath + "TempFile/CircRef.fa.pac "
                                         + strRootPath + "TempFile/CircRef.fa.sa "
                                         + strRootPath + "TempFile/Read.bam "
                                         + strRootPath + "TempFile/Read.sam "
                                         + strRootPath + "TempFile/Read.sorted.bam "
                                         + strRootPath + "TempFile/Read.sorted.bam.bai";
        system(strRemoveTempFiles.c_str());

        if(enReadsType == rtPairEnd)
        {
            strSortedBamFile  = ClsBWA::GetInstance().CreateBamByPEReads(strRefPath, strReads1Path, strReads2Path,
                                                                     true, "", true, true, "", 12, true);
        }
        else if(enReadsType == rtSingle)
        {
            strSortedBamFile  = ClsBWA::GetInstance().CreateBamBySingleReads(strRefPath, strReads1Path, "", "",
                                                                             true, true, true, 12, true);
        }
//        string strReadsSumPath = "/home/lq/lxwg/WorkStudio/Prototype/CircRNA_DBG/TempFile/readsSum.fastq";
//        strSortedBamFile  = ClsBWA::GetInstance().CreateBamBySingleReads(strRefPath, strReadsSumPath);
    }

    //Use bam tools to check the mapping status
    ///1: Check all of CircInfo
    vector<St_CircRefMapStatus> vCircRefMapStatus;
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(strRefPath, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

    St_CircRefMapStatus stMS;
    int iRefID = 1; // Notice: the Ref is start from 1!!!!
    for(vector<St_Fasta>::iterator itr = vFasta.begin(); itr != vFasta.end(); itr++)
    {
        //Parse the file name
        ///1. For CircRNA Name
        int iStartPos = 0;
        int iEndPos = itr->strName.find(" ", iStartPos);
        int iLen = iEndPos - iStartPos;
        string strRawName = itr->strName.substr(iStartPos, iLen);
        int iSplitPos = strRawName.find("_");
        stMS.strName = strRawName.substr(0, iSplitPos); // set the name
        if(strRawName.find("compl") != string::npos) // find it
            stMS.bCompl = true;
        else
            stMS.bCompl = false;

        ///2. For Key Pos
        iStartPos = iEndPos + 1;
        iEndPos = itr->strName.find(" ", iStartPos);
        iLen = iEndPos - iStartPos;
        stMS.iKeyPos = atoi(itr->strName.substr(iStartPos, iLen).c_str());

        ///3: For Orientation
        iStartPos = iEndPos + 1;
        stMS.bRC = (itr->strName.substr(iStartPos, 1) == "+" ? false: true);

        stMS.strRefSeq = itr->strSeq;

        stMS.iRefID = iRefID;

        vCircRefMapStatus.push_back(stMS);
        iRefID++;
    }

    //Go to check sorted bam file
    //1: parse bam file & output the expected reads
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strSortedBamFile);
    pBamReader->OpenIndex(strSortedBamFile + ".bai");
    BamAlignment al;
    ofstream ofsTemp;
    ofsTemp.open("TempDirect.fa");
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        //Case 1: do not map  ------->
        if(!al.IsMapped() || !al.IsMateMapped() || al.QueryBases == "") // need all of them be mapped
            continue;

        ///Get the corresponding circ-ref  (Notice here: the ref id is start from 0 !!!!!)
        St_CircRefMapStatus& stTempMS = vCircRefMapStatus[al.RefID];

        if(stTempMS.bCompl) // we only consider the reads in reference related
            continue;

        //Case 2: If the reads maps fine
        ///1: Check if there is soft clip
        int iStart = al.Position;
        int iEnd = iStart + al.AlignedBases.length();
        /*
        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin();
            itr != al.CigarData.end(); itr++)
        {
            switch(itr->Type)
            {
                case 'M': // alignment match (can be a sequence match or mismatch)
                    iStart += itr->Length;
                    break;
                case 'I': // insertion to the reference
                    iStart += itr->Length;
                    break;
                case 'D': // deletion from the reference
                case 'N':  // skipped region from the reference
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    bFindSoftClip = true;
                    break;
                case 'P': // padding (silent deletion from padded reference)
                case '=': // sequence match
                case 'X': // sequence mismatch
                    iStart += itr->Length;
                    break;
            }
            if(bFindSoftClip)
                break;
        }*/

        //check its mats ---> Go!!!
        if(stTempMS.iKeyPos > iStart && stTempMS.iKeyPos < iEnd) //at least need to be contained
        {
            //check it's complementary
//            St_CircRefMapStatus& stTempMSMate = vCircRefMapStatus[al.MateRefID];
//            if(stTempMSMate.strName == stTempMS.strName) // both of current and mat come from one original circular rna
            {
                string strAlignSeq = al.AlignedBases + (al.IsReverseStrand() ? " -" : " +");
                stTempMS.vReads.push_back(strAlignSeq);//(al.QueryBases);
            }
        }

        //Cout the mapping relationship directly
        ofsTemp << ">" << vFasta[al.RefID].strName << " " << IntToStr(al.RefID) << endl;
        ofsTemp << vFasta[al.RefID].strSeq << endl;
        ofsTemp << "----------- " << al.AlignedBases << endl;
    }
    ofsTemp.close();
    delete pBamReader;
    pBamReader = NULL;

    //Cout the file
    ofstream ofs;
    ofs.open("CircRefMappingResult.sudo.fa");
    int iValidHitNum = 0;
    int iTotalNum = 0;
    int iRCValidNum = 0;
    int iRCSum = 0;
    int iRegValidNum = 0;
    int iRegSum = 0;

    for(vector<St_CircRefMapStatus>::iterator itr = vCircRefMapStatus.begin();
        itr != vCircRefMapStatus.end(); itr++)
    {
        if(itr->bCompl)
            continue;

        //We only output the case which supported by some reads
        if(itr->vReads.empty())
            continue;

        iTotalNum++;

        ofs << ">" << itr->strName << " " << IntToStr(itr->iKeyPos) << " "
            << (itr->bRC ? "-" : "+") << endl;
        ofs << itr->strRefSeq << endl;
        for(vector<string>::iterator subItr = itr->vReads.begin(); subItr != itr->vReads.end(); subItr++)
        {
            ofs << " ------ " << *subItr << endl;
        }

        if(!itr->vReads.empty())
        {
            iValidHitNum++;
            if(itr->bRC)
                iRCValidNum++;
            else
                iRegValidNum++;
        }

        if(itr->bRC)
            iRCSum++;
        else
            iRegSum++;
    }
    ofs.close();

    cout << endl;
    cout << "Valid Hit Num: " << IntToStr(iValidHitNum) << endl;
    cout << "Total Num    : " << IntToStr(iTotalNum) << endl;
    cout << "RC Valid Num : " << IntToStr(iRCValidNum) << endl;
    cout << "RC Sum       : " << IntToStr(iRCSum) << endl;
    cout << "Reg Valid Num: " << IntToStr(iRegValidNum) << endl;
    cout << "Reg Sum      : " << IntToStr(iRegSum) << endl;
}

