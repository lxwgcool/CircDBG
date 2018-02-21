#include "clspreevaluation.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <iostream>

const char* arryCellNameCopy2[ctMaxx] = {"H9", "glioblastoma", "brain"};

ClsPreEvaluation::ClsPreEvaluation(St_Config& stConfig)
{
    this->m_stConfig = stConfig;

    //Initiate m_enCellType
    if(m_stConfig.strTissue == arryCellNameCopy2[ctH9])
        m_enCellType = ctH9;
    else if(m_stConfig.strTissue == arryCellNameCopy2[ctGlioblastoma])
        m_enCellType = ctGlioblastoma;
    else if(m_stConfig.strTissue == arryCellNameCopy2[ctBrain])
        m_enCellType = ctBrain;
    else
        m_enCellType = ctMaxx;
}

ClsPreEvaluation::~ClsPreEvaluation()
{
}

void ClsPreEvaluation::GeneralEvaluate()
{
    //Parse Std CircRNA Info
    ClsParseStdCircRNAInfo* pParseStdCircRNAInfo =  new ClsParseStdCircRNAInfo();
    pParseStdCircRNAInfo->Init(m_stConfig.strTissue);

    ///Get CircRNA Info
    vector<St_StdCircRNAInfo> vStdCircRNAInfo;
    pParseStdCircRNAInfo->ParseCircRNADb(m_stConfig.strCircRNADb, vStdCircRNAInfo);
    pParseStdCircRNAInfo->FilterByCell(vStdCircRNAInfo, m_enCellType);

    ////Cout which chromosome get involved in this event
    vector<string> vChromMark;
    for(vector<St_StdCircRNAInfo>::iterator itr = vStdCircRNAInfo.begin();
        itr != vStdCircRNAInfo.end(); itr++)
    {
        bool bFind = false;
        for(vector<string>::iterator subItr = vChromMark.begin();
            subItr != vChromMark.end(); subItr++)
        {
            if(itr->strChrom == *subItr)
            {
                bFind = true;
                break;
            }
        }
        if(!bFind)
            vChromMark.push_back(itr->strChrom);
    }

    /// ---> output all of chromosones get involved in this circular events.
    cout << endl << "-------------" << endl;
    for(vector<string>::iterator subItr = vChromMark.begin();
        subItr != vChromMark.end(); subItr++)
    {
        cout << *subItr << endl;
    }
    cout << "-------------" << endl << endl;
    /// <---

    ///Get GTF file
    //vector<St_Row_Chrom> vChrom;
    //pParseStdCircRNAInfo->ParseGTFWithCircRNADbStyle(m_stConfig.strGtfPath, vChrom);

    /////Get the full set annotation file (Full Set GTF)
    vector<St_Row_Chrom> vChrom;
    ClsGTFParse* pGtfParse = new ClsGTFParse();   
    pGtfParse->ReadGTF(vChrom, m_stConfig.strGtfPath);

    //Cout那些没能hit到exon两端的circular rna
    pParseStdCircRNAInfo->CheckIrregularCircRNA(vStdCircRNAInfo, vChrom); //after this: vStdCircRNAInfo --> HitCircRNAInfo

//    string strCircRefPath = pParseStdCircRNAInfo->GenerateCircularRelatedRef(vStdCircRNAInfo, vChrom,
//                                                                             m_stConfig.strRefPath);

    delete pParseStdCircRNAInfo;
    pParseStdCircRNAInfo = NULL;

    //Map reads to Circ-Reference
    ClsReadsMapping* pReadsMapping = new ClsReadsMapping();
//    pReadsMapping->MapReadsToCircRNADb(strCircRefPath, m_stConfig.strReads1Path, m_stConfig.strReads2Path,
//                                       rtPairEnd);
    delete pReadsMapping;
    pReadsMapping = NULL;
}

void ClsPreEvaluation::GetValidStdCandiByChormoson(int iChromIndex)
{
    //Get the full set annotation file (Full Set GTF)
    vector<St_Row_Chrom> vChrom;
    ClsGTFParse* pGtfParse = new ClsGTFParse();
    pGtfParse->ReadGTF(vChrom, m_stConfig.strGtfPath);
    delete pGtfParse;
    pGtfParse = NULL;

    //Get & Print valid std circRNA
    cout << "Print the Customized Candidate" << endl;
    ClsParseStdCircRNAInfo* pParseStdCircRNAInfo =  new ClsParseStdCircRNAInfo();
    pParseStdCircRNAInfo->Init(m_stConfig.strTissue);
    vector<St_StdCircRNAInfo> vTargetValidStdInfo;
    pParseStdCircRNAInfo->PrintStdHitInfoByBriefCandiStyle(vTargetValidStdInfo, vChrom,
                                                           m_stConfig.strCircRNADb, iChromIndex);

//    //Create reference related to
//    string strCircRefPath = pParseStdCircRNAInfo->GenerateCircularRelatedRef(vTargetValidStdInfo, vChrom,
//                                                                             m_stConfig.strRefPath);
//    //Do mapping
//    ClsReadsMapping* pReadsMapping = new ClsReadsMapping();
//    pReadsMapping->MapReadsToCircRNADb(strCircRefPath, m_stConfig.strReads1Path, m_stConfig.strReads2Path, rtSingle);
//    delete pReadsMapping;
//    pReadsMapping = NULL;


    delete pParseStdCircRNAInfo;
    pParseStdCircRNAInfo = NULL;
}



