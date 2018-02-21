#include "clsresultcomparison.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <fstream>
#include <cstring>
#include <unistd.h>

const string TOOLNAME[swMax] = {"My_Program", "CIRI", "Find_circ",
                                "CIRCExplorer", "CircRNAFinder", "CircMarker"};

ClsResultComparison::ClsResultComparison(): m_pCBF(NULL), m_iReadLen(0)
{
}

ClsResultComparison::~ClsResultComparison()
{
    if(m_pCBF != NULL)
    {
        delete m_pCBF;
        m_pCBF = NULL;
    }
}

void ClsResultComparison::Init(St_Config& stConfig)
{
    m_pCBF = new ClsCompareBaseFunction();
    m_pCBF->Init(stConfig.iReadLen, stConfig.strTissue);

    m_iReadLen = stConfig.iReadLen;
    m_iCIRIHitNum = 0;
    m_fCIRITP = 0;

    m_iFindCircHitNUm = 0;
    m_fFindCircTP = 0;

    m_iMaxSupportNum = stConfig.iMaxSupportNum;
}

void ClsResultComparison::CompareMyProgramWithSimulator(vector<St_Row_Chrom>& vChrom,
                                                        string strMyProgramPath, string strSimulatorBenchmark)
{
    if(access(strMyProgramPath.c_str(), 0) != 0 || access(strSimulatorBenchmark.c_str(), 0) != 0)
    {
        cout << "Compare (MyProgram) With (Simulator): One File Missing" << endl;
        return;
    }
    //Parse Simulator
    vector<St_Candidate> vCandiSimulator;
    m_pCBF->ParseSimulationResult(strSimulatorBenchmark, vCandiSimulator);
    m_pCBF->SetCircTag(vChrom, vCandiSimulator);
    string strSimulator_HitPath = "./Comparison_Result/SimulatorHit.txt";
    m_pCBF->SaveCandi(vCandiSimulator, strSimulator_HitPath);

    //Compare
    CompareThirdPartyHitAndStdDBHit(vChrom, strMyProgramPath, strSimulator_HitPath, swMyProgram);
}

///CIRI with Simulator
void ClsResultComparison::CompareCIRIWithSimulator(vector<St_Row_Chrom>& vChrom,
                                                   string strCIRIPath, string strSimulatorBenchmark)
{
    if(access(strCIRIPath.c_str(), 0) != 0 || access(strSimulatorBenchmark.c_str(), 0) != 0)
    {
        cout << "Compare (CIRI) With (Simulator): One File Missing" << endl;
        return;
    }

    //Parse Simulator
    vector<St_Candidate> vCandiSimulator;
    m_pCBF->ParseSimulationResult(strSimulatorBenchmark, vCandiSimulator);
    m_pCBF->SetCircTag(vChrom, vCandiSimulator);
    string strSimulator_HitPath = "./Comparison_Result/SimulatorHit.txt";
    m_pCBF->SaveCandi(vCandiSimulator, strSimulator_HitPath);

    //Parse CIRI
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRIPath, vCandiCIRI);
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Org_Hit.txt";
    m_pCBF->SaveCandi(vCandiCIRI, strCIRIHitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCIRIHitPath, strSimulator_HitPath, swCIRI);
}

///CIRCExplorer with Simulator
void ClsResultComparison::CompareCIRCExplorerWithSimulator(vector<St_Row_Chrom>& vChrom,
                                                           string strCircExplorerPath, string strSimulatorBenchmark)
{
    if(access(strCircExplorerPath.c_str(), 0) != 0 || access(strSimulatorBenchmark.c_str(), 0) != 0)
    {
        cout << "Compare (CircExplorer) With (Simulator): One File Missing" << endl;
        return;
    }

    //Parse Simulator
    vector<St_Candidate> vCandiSimulator;
    m_pCBF->ParseSimulationResult(strSimulatorBenchmark, vCandiSimulator);
    m_pCBF->SetCircTag(vChrom, vCandiSimulator);
    string strSimulator_HitPath = "./Comparison_Result/SimulatorHit.txt";
    m_pCBF->SaveCandi(vCandiSimulator, strSimulator_HitPath);

    //Parse CIRCExplorer
    vector<St_Candidate> vCandiCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorerPath, vCandiCircExplorer);
    m_pCBF->SetCircTag(vChrom, vCandiCircExplorer);
    string strCirc_Explorer_HitPath = "./Comparison_Result/Circ_Explorer_Org_Hit.txt";
    m_pCBF->SaveCandi(vCandiCircExplorer, strCirc_Explorer_HitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCirc_Explorer_HitPath, strSimulator_HitPath, swCIRCExplorer);
}

///Find-Circ with Simulator
void ClsResultComparison::CompareFindCircWithSimulator(vector<St_Row_Chrom>& vChrom,
                                                       string strFindCircPath, string strSimulatorBenchmark)
{
    if(access(strFindCircPath.c_str(), 0) != 0 || access(strSimulatorBenchmark.c_str(), 0) != 0)
    {
        cout << "Compare (FindCirc) With (Simulator): One File Missing" << endl;
        return;
    }

    //Parse Simulator
    vector<St_Candidate> vCandiSimulator;
    m_pCBF->ParseSimulationResult(strSimulatorBenchmark, vCandiSimulator);
    m_pCBF->SetCircTag(vChrom, vCandiSimulator);
    string strSimulator_HitPath = "./Comparison_Result/SimulatorHit.txt";
    m_pCBF->SaveCandi(vCandiSimulator, strSimulator_HitPath);

    //Parse FindCirc
    vector<St_Candidate> vCandiFindCirc;
    m_pCBF->ParseFindCircResult(strFindCircPath, vCandiFindCirc);
    m_pCBF->SetCircTag(vChrom, vCandiFindCirc);
    string strFind_Circ_HitPath = "./Comparison_Result/Find_Circ_Org_Hit.txt";
    m_pCBF->SaveCandi(vCandiFindCirc, strFind_Circ_HitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strFind_Circ_HitPath, strSimulator_HitPath, swFind_circ);
}

///CircRNAFinder with Simulator
void ClsResultComparison::CompareCircRNAFinderWithSimulator(vector<St_Row_Chrom>& vChrom,
                                                            string strCircRNAFinder, string strSimulatorBenchmark)
{
    if(access(strCircRNAFinder.c_str(), 0) != 0 || access(strSimulatorBenchmark.c_str(), 0) != 0)
    {
        cout << "Compare (CircRNAFinder) With (Simulator): One File Missing" << endl;
        return;
    }

    //Parse Simulator
    vector<St_Candidate> vCandiSimulator;
    m_pCBF->ParseSimulationResult(strSimulatorBenchmark, vCandiSimulator);
    m_pCBF->SetCircTag(vChrom, vCandiSimulator);
    string strSimulator_HitPath = "./Comparison_Result/SimulatorHit.txt";
    m_pCBF->SaveCandi(vCandiSimulator, strSimulator_HitPath);

    //Parse CircRNAFinder
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Org_Hit.txt";
    m_pCBF->SaveCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCircRNAFinderHit, strSimulator_HitPath, swCircRNAFinder);
}

///CircMarker with Simulator
void ClsResultComparison::CompareCircMarkerWithSimulator(vector<St_Row_Chrom>& vChrom,
                                                         string strCircMarkerPath, string strSimulatorBenchmark)
{
    if(access(strCircMarkerPath.c_str(), 0) != 0 || access(strSimulatorBenchmark.c_str(), 0) != 0)
    {
        cout << "Compare (CircMarker) With (Simulator): One File Missing" << endl;
        return;
    }

    //Parse Simulator
    vector<St_Candidate> vCandiSimulator;
    m_pCBF->ParseSimulationResult(strSimulatorBenchmark, vCandiSimulator);
    m_pCBF->SetCircTag(vChrom, vCandiSimulator);
    string strSimulator_HitPath = "./Comparison_Result/SimulatorHit.txt";
    m_pCBF->SaveCandi(vCandiSimulator, strSimulator_HitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCircMarkerPath, strSimulator_HitPath, swCircMarker);
}
//<---


void ClsResultComparison::CompareMyProgramAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                                       string strMyProgramPath, string strCircRNAdbPath,
                                                       string strChromName, ofstream* pOfs)
{
    if(access(strMyProgramPath.c_str(), 0) != 0 || access(strCircRNAdbPath.c_str(), 0) != 0)
    {
        cout << "Warning: CompareMyProgramAndCircRNAdb: one path do not exist!" << endl;
        return;
    }

    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strMyProgramPath, strCircRNADb_HitPath, swMyProgram,
                                    strChromName, pOfs);
}

void ClsResultComparison::CompareCIRIAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                                  string strCIRIPath, string strCircRNAdbPath,
                                                  string strChromName, ofstream* pOfs)
{
    if(access(strCIRIPath.c_str(), 0) != 0 || access(strCircRNAdbPath.c_str(), 0) != 0)
    {
        cout << "Warning: CompareCIRIAndCircRNAdb: one path do not exist!" << endl;
        return;
    }

    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);

    //Parse CIRI
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRIPath, vCandiCIRI);
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Org_Hit.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCIRIHitPath, strCircRNADb_HitPath, swCIRI,
                                    strChromName, pOfs);
}

void ClsResultComparison::CompareFindCircAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                                  string strFindCircPath, string strCircRNAdbPath,
                                                  string strChromName, ofstream* pOfs)
{
    if(access(strFindCircPath.c_str(), 0) != 0 || access(strCircRNAdbPath.c_str(), 0) != 0)
    {
        cout << "Warning: CompareFindCircAndCircRNAdb: one path do not exist!" << endl;
        return;
    }

    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);

    //Parse FindCirc
    vector<St_Candidate> vCandiFindCirc;
    m_pCBF->ParseFindCircResult(strFindCircPath, vCandiFindCirc);
    m_pCBF->SetCircTag(vChrom, vCandiFindCirc);
    string strFind_Circ_HitPath = "./Comparison_Result/Find_Circ_Org_Hit.txt";
    m_pCBF->SaveHitCandi(vCandiFindCirc, strFind_Circ_HitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strFind_Circ_HitPath, strCircRNADb_HitPath, swFind_circ,
                                    strChromName, pOfs);
}

void ClsResultComparison::CompareCircExplorerAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                                          string strCircExplorerPath, string strCircRNAdbPath,
                                                          string strChromName, ofstream* pOfs)
{
    if(access(strCircExplorerPath.c_str(), 0) != 0 || access(strCircRNAdbPath.c_str(), 0) != 0)
    {
        cout << "Warning: CompareCircExplorerAndCircRNAdb: one path do not exist!" << endl;
        return;
    }

    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);

    vector<St_Candidate> vCandiCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorerPath, vCandiCircExplorer);
    m_pCBF->SetCircTag(vChrom, vCandiCircExplorer);
    string strCirc_Explorer_HitPath = "./Comparison_Result/Circ_Explorer_Org_Hit.txt";
    m_pCBF->SaveHitCandi(vCandiCircExplorer, strCirc_Explorer_HitPath);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCirc_Explorer_HitPath, strCircRNADb_HitPath, swCIRCExplorer,
                                    strChromName, pOfs);
}

void ClsResultComparison::CompareCircRNAFinderAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                                          string strCircRNAFinder, string strCircRNAdbPath,
                                                          string strChromName, ofstream* pOfs)
{
    if(access(strCircRNAFinder.c_str(), 0) != 0 || access(strCircRNAdbPath.c_str(), 0) != 0)
    {
        cout << "Warning: CompareCircRNAFinderAndCircRNAdb: one path do not exist!" << endl;
        return;
    }

    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);

    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Org_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCircRNAFinderHit, strCircRNADb_HitPath, swCircRNAFinder,
                                    strChromName, pOfs);
}

void ClsResultComparison::CompareCircMarkerAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                                       string strCircMarker, string strCircRNAdbPath,
                                                       string strChromName, ofstream* pOfs)
{
    if(access(strCircMarker.c_str(), 0) != 0 || access(strCircRNAdbPath.c_str(), 0) != 0)
    {
        cout << "Warning: CompareCircMarkerAndCircRNAdb: one path do not exist!" << endl;
        return;
    }

    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);

    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Org_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareThirdPartyHitAndStdDBHit(vChrom, strCircMarkerHit, strCircRNADb_HitPath, swCircMarker,
                                    strChromName, pOfs);
}

void ClsResultComparison::CompareThirdPartyHitAndStdDBHit(vector<St_Row_Chrom>& vChrom,
                                                          string strThirdpartyHitPath, string strStdDBHitPath,
                                                          En_Software enSoftware,
                                                          string strChromName, ofstream* pOfs)
{
    // 在这里，我们只跟两边都能对上exon boundary的circ base的result进行比较
    ///Step1: Get my Result
    vector<St_Candidate> v3rdPartyToolHit;
    m_pCBF->ParseBrief(strThirdpartyHitPath, v3rdPartyToolHit);

    if(enSoftware == swMyProgram)
    {
        for(vector<St_Candidate>::iterator itr = v3rdPartyToolHit.end() - 1;
            itr >= v3rdPartyToolHit.begin(); itr--)
        {
            if(itr->iSupportNum > m_iMaxSupportNum)
                v3rdPartyToolHit.erase(itr);
        }
    }

    ///Step2: Get the valid candi in CircBase: (need to combine both vChrom and CircBase)
    vector<St_Candidate> vDBHit; //DB is circBase
    m_pCBF->ParseBrief(strStdDBHitPath, vDBHit); //这里已经是所有能够hit上的candidate了
    //我们在这里 all of chromosomes could be considered
    vector<St_Candidate> vHitDBForSpeciChrom;
    if(strChromName != "")
    {
        m_pCBF->GetCandiForSpeciChrom(vChrom, vDBHit, vHitDBForSpeciChrom, strChromName);
    }
    else
        vHitDBForSpeciChrom = vDBHit;

    ///Step3: Make the comparison: here "start pos" is smaller than "end pos"
    ofstream ofsUnHit3rdPartyTool;
    ofstream ofsHit3rdPartyTool;
    ofstream ofsUnHitStd;
    switch(enSoftware)
    {
        case swMyProgram:
            ofsUnHit3rdPartyTool.open("./Comparison_Result/My_Program_UnHit_DB.txt");
            ofsHit3rdPartyTool.open("./Comparison_Result/My_Program_Hit_DB.txt");
            ofsUnHitStd.open("./Comparison_Result/My_Program_Std_From_DB_UnHit.txt");
            break;
        case swCIRI:
            ofsUnHit3rdPartyTool.open("./Comparison_Result/CIRI_UnHit_DB.txt");
            ofsHit3rdPartyTool.open("./Comparison_Result/CIRI_Hit_DB.txt");
            ofsUnHitStd.open("./Comparison_Result/CIRI_Std_From_DB_UnHit.txt");
            break;
        case swFind_circ:
            ofsUnHit3rdPartyTool.open("./Comparison_Result/Find_circ_UnHit_DB.txt");
            ofsHit3rdPartyTool.open("./Comparison_Result/Find_circ_Hit_DB.txt");
            ofsUnHitStd.open("./Comparison_Result/Find_circ_Std_From_DB_UnHit.txt");
            break;
        case swCIRCExplorer:
            ofsUnHit3rdPartyTool.open("./Comparison_Result/CIRCExplorer_UnHit_DB.txt");
            ofsHit3rdPartyTool.open("./Comparison_Result/CIRCExplorer_Hit_DB.txt");
            ofsUnHitStd.open("./Comparison_Result/CIRCExplorer_Std_From_DB_UnHit.txt");
            break;
        case swCircRNAFinder:
            ofsUnHit3rdPartyTool.open("./Comparison_Result/CircRNAFinder_UnHit_DB.txt");
            ofsHit3rdPartyTool.open("./Comparison_Result/CircRNAFinder_Hit_DB.txt");
            ofsUnHitStd.open("./Comparison_Result/CircRNAFinder_Std_From_DB_UnHit.txt");
            break;
        case swCircMarker:
            ofsUnHit3rdPartyTool.open("./Comparison_Result/CircMarker_UnHit_DB.txt");
            ofsHit3rdPartyTool.open("./Comparison_Result/CircMarker_Hit_DB.txt");
            ofsUnHitStd.open("./Comparison_Result/CircMarker_Std_From_DB_UnHit.txt");
            break;
        default:
            break;
    }

    //-->
    int iSelfCircNum = 0;
    int iRegCircNum = 0;

    int iHitSum = 0;
    int iSelfHit = 0;
    int iRegHit = 0;
    //<--

    for(vector<St_Candidate>::iterator itr = v3rdPartyToolHit.begin(); itr != v3rdPartyToolHit.end(); itr++)
    {
        if(itr->enType == ctSelf)
            iSelfCircNum++;
        else if(itr->enType == ctRegular)
            iRegCircNum++;        

        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vHitDBForSpeciChrom.begin();
            subItr != vHitDBForSpeciChrom.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName) //check if they are in same chromosome
                continue;

            bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                iHitSum++;
                if(itr->enType == ctSelf)
                    iSelfHit++;
                else if(itr->enType == ctRegular)
                    iRegHit++;

                subItr->bHit = true;
                break;
            }
        }

        if(!bFind)
        {
            ofsUnHit3rdPartyTool << itr->strChromName << " "
                                 << IntToStr(itr->iStartPos) << " "
                                 << IntToStr(itr->iEndPos) << " "
                                 << IntToStr(itr->iSupportNum) << " "
                                 << itr->strTag << " "
                                 << itr->GetTypeString() << endl;
        }
        else //Find it
        {
            //Only Record the self circle --> ONLY FOR NOW!
            if(itr->enType == ctSelf || itr->enType == ctRegular)
            {
                ofsHit3rdPartyTool << itr->strChromName << " "
                                   << IntToStr(itr->iStartPos) << " "
                                   << IntToStr(itr->iEndPos) << " "
                                   << IntToStr(itr->iSupportNum) << " "
                                   << itr->strTag << " "
                                   << itr->GetTypeString() << endl;
            }
        }
    }
    ofsUnHit3rdPartyTool.close();
    ofsHit3rdPartyTool.close();

    //Check the unhit candidate --> this is very important, sicne those candidate are the real one
    //                              which we want to find them back
    for(vector<St_Candidate>::iterator itr = vHitDBForSpeciChrom.begin();
        itr != vHitDBForSpeciChrom.end(); itr++)
    {
        if(itr->bHit)
            continue;
        else  //If does not hit by other
        {
            ofsUnHitStd << itr->strChromName << " "
                        << IntToStr(itr->iStartPos) << " "
                        << IntToStr(itr->iEndPos) << " "
                        << IntToStr(itr->iSupportNum) << " "
                        << itr->strTag << " "
                        << itr->GetTypeString() << endl;
        }
    }
    ofsUnHitStd.close();

    //Calculate F1 Score -->
    int iStdSum = vHitDBForSpeciChrom.size();
    int iToolResultSum = v3rdPartyToolHit.size();
    int iToolHitSum = iSelfHit + iRegHit;
    float fPrecision = (float)iToolHitSum / iToolResultSum;
    float fRecall = (float)iToolHitSum / iStdSum;
    float fF1 = 2 * fPrecision * fRecall / (fPrecision + fRecall);
    //<--

    //Print Result
    if(pOfs != NULL)
    {
        if(strChromName != "")
        {
            (*pOfs) << "Chrom " << strChromName;
        }
        else
            (*pOfs) << "Chrom All";

        cout << " DB CircRNA Size: " << vHitDBForSpeciChrom.size() << endl;

        (*pOfs) << "**********ThirdPartyTools: "<< TOOLNAME[enSoftware] << "**********" << endl;
        (*pOfs) << "Self Circ Num :" << IntToStr(iSelfCircNum) << " | "
             << GetRatio((float)iSelfCircNum / v3rdPartyToolHit.size()) << endl;
        (*pOfs) << "Reg Circ Num  :" << IntToStr(iRegCircNum) << " | "
             << GetRatio((float)iRegCircNum / v3rdPartyToolHit.size()) << endl;
        (*pOfs) << "Find Sum      : " << IntToStr(iSelfCircNum + iRegCircNum)
                << " | " << GetRatio((float)(iSelfCircNum + iRegCircNum) / v3rdPartyToolHit.size())
                << endl;
        (*pOfs) << "Org Detect Sum:" << IntToStr(v3rdPartyToolHit.size()) << endl << endl;

        (*pOfs) << "----Hit Case----" << endl;
        (*pOfs) << "Hit Self Num  :" << IntToStr(iSelfHit) << " | "
             << GetRatio((float)iSelfHit / iSelfCircNum) << endl;
        (*pOfs) << "Hit Reg Num   :" << IntToStr(iRegHit) << " | "
             << GetRatio((float)iRegHit / iHitSum) << endl;
        (*pOfs) << "Hit Sum       : " << IntToStr(iSelfHit + iRegHit) << endl;
        (*pOfs) << "-->Precision: " << GetRatio(fPrecision) << " "
                << "-->Recall: " << GetRatio(fRecall) << " "
                << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;
        (*pOfs) << endl;
    }
    else
    {
        cout << "Chrom " << strChromName << " DB CircRNA Size: " << vHitDBForSpeciChrom.size() << endl;

        cout << "**********ThirdPartyTools: "<< TOOLNAME[enSoftware] << "**********" << endl;
        cout << "Self Circ Num :" << IntToStr(iSelfCircNum) << " | "
             << GetRatio((float)iSelfCircNum / v3rdPartyToolHit.size()) << endl;
        cout << "Reg Circ Num  :" << IntToStr(iRegCircNum) << " | "
             << GetRatio((float)iRegCircNum / v3rdPartyToolHit.size()) << endl;
        cout << "Find Sum      : " << IntToStr(iSelfCircNum + iRegCircNum)
                << " | " << GetRatio((float)(iSelfCircNum + iRegCircNum) / v3rdPartyToolHit.size())
                << endl;
        cout << "Org Detect Sum:" << IntToStr(v3rdPartyToolHit.size()) << endl << endl;

        cout << "----Hit Case----" << endl;
        cout << "Hit Self Num  :" << IntToStr(iSelfHit) << " | "
             << GetRatio((float)iSelfHit / iSelfCircNum) << endl;
        cout << "Hit Reg Num   :" << IntToStr(iRegHit) << " | "
             << GetRatio((float)iRegHit / iHitSum) << endl;
        cout << "Hit Sum       : " << IntToStr(iSelfHit + iRegHit) << endl;
        cout << "-->Precision: " << GetRatio(fPrecision) << " "
             << "-->Recall: " << GetRatio(fRecall) << " "
             << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;
        cout << endl;
    }
}

//Compare My Program with others
void ClsResultComparison::CheckIntersetBTMyProgramAndCIRI(vector<St_Row_Chrom>& vChrom,
                                     string strMyProgram, string strCIRI, ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareOrgResultWithOtherTool(strMyProgram, strCIRIHitPath, swMyProgram, swCIRI, pOfs);
}

void ClsResultComparison::CheckIntersetBTMyProgramAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                         string strMyProgram, string strFindCirc, ofstream* pOfs)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And Find-circ<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And Find-circ<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);

    CompareOrgResultWithOtherTool(strMyProgram, strFindCircHit, swMyProgram, swFind_circ, pOfs);
}

void ClsResultComparison::CheckIntersetBTMyProgramAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                             string strMyProgram, string strCircExplorer,
                                             ofstream* pOfs)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strCircExplorer.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    CompareOrgResultWithOtherTool(strMyProgram, strCircExplorerHit, swMyProgram, swCIRCExplorer);
}

void ClsResultComparison::CheckIntersetBTMyProgramAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                             string strMyProgram, string strCircRNAFinder,
                                             ofstream* pOfs)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strCircRNAFinder.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strMyProgram, strCircRNAFinderHit, swMyProgram, swCircRNAFinder);
}

void ClsResultComparison::CheckIntersetBTMyProgramAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                                                string strMyProgram, string strCircMarker,
                                                                ofstream* pOfs)
{
    if(access(strMyProgram.c_str(), 0) != 0 || access(strCircMarker.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT My Program And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareOrgResultWithOtherTool(strMyProgram, strCircMarkerHit, swMyProgram, swCircMarker);
}

//*********Compare CIRI with Others: MyProgram, Find_Circ, Circ_Explorer*********************
// 比对 My program VS CIRI
//                   Find_Circ
//                   Circ_Explorer
void ClsResultComparison::CheckIntersetBTCIRIAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                     string strCIRI, string strMyProgram, ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And My Program<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And My Program<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareOrgResultWithOtherTool(strCIRIHitPath, strMyProgram, swCIRI, swMyProgram, pOfs);
}

void ClsResultComparison::CheckIntersetBTCIRIAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                     string strCIRI, string strFindCirc, ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);

    CompareOrgResultWithOtherTool(strCIRIHitPath, strFindCircHit, swCIRI, swFind_circ, pOfs);
}

void ClsResultComparison::CheckIntersetBTCIRIAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                     string strCIRI, string strCircExplorer, ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strCircExplorer.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And CircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And CircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    //2: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    CompareOrgResultWithOtherTool(strCIRIHitPath, strCircExplorerHit, swCIRI, swCIRCExplorer);
}

void ClsResultComparison::CheckIntersetBTCIRIAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                             string strCIRI, string strCircRNAFinder,
                                             ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strCircRNAFinder.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    //2: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strCIRIHitPath, strCircRNAFinderHit, swCIRI, swCircRNAFinder);
}

void ClsResultComparison::CheckIntersetBTCIRIAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                                           string strCIRI, string strCircMarker,
                                                           ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strCircMarker.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CIRI And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareOrgResultWithOtherTool(strCIRIHitPath, strCircMarkerHit, swCIRI, swCircMarker);
}

//*********Compare Find_Circ with Others: Ciri, MyProgram, Circ_Explorer*********************
// 比对 My program VS CIRI
//                   Find_Circ
//                   Circ_Explorer
void ClsResultComparison::CheckIntersetBTFindCircAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                     string strFindCirc, string strMyProgram, ofstream* pOfs)
{
    if(access(strFindCirc.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset Find_Circ And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset Find_Circ And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);

    CompareOrgResultWithOtherTool(strFindCircHit, strMyProgram, swFind_circ, swMyProgram, pOfs);
}

void ClsResultComparison::CheckIntersetBTFindCircAndCIRI(vector<St_Row_Chrom>& vChrom,
                                     string strFindCirc, string strCIRI, ofstream* pOfs)
{
    if(access(strCIRI.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset Find_Circ And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset Find_Circ And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareOrgResultWithOtherTool(strFindCircHit, strCIRIHitPath, swFind_circ, swCIRI, pOfs);
}

void ClsResultComparison::CheckIntersetBTFindCircAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                     string strFindCirc, string strCircExplorer, ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset Find_Circ And CIRCExplorer<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset Find_Circ And CIRCExplorer<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);

    //2: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    CompareOrgResultWithOtherTool(strFindCircHit, strCircExplorerHit, swFind_circ, swCIRCExplorer, pOfs);
}

void ClsResultComparison::CheckIntersetBTFindCircAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                                                 string strFindCirc, string strCircRNAFinder,
                                                                 ofstream* pOfs)
{
    if(access(strFindCirc.c_str(), 0) != 0 || access(strCircRNAFinder.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT strFindCirc And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT strFindCirc And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);


    //2: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strFindCircHit, strCircRNAFinderHit, swFind_circ, swCircRNAFinder);
}

void ClsResultComparison::CheckIntersetBTFindCircAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                                                string strFindCirc, string strCircMarker,
                                                                ofstream* pOfs)
{
    if(access(strFindCirc.c_str(), 0) != 0 || access(strCircMarker.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT FindCirc And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             << endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT FindCirc And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             << endl;
    }

    //1: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);

    //2: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareOrgResultWithOtherTool(strFindCircHit, strCircMarkerHit, swFind_circ, swCircMarker);
}

//*********Compare Circ_Explorer with Others: MyProgram, Ciri, Find_Circ*********************
// 比对 My program VS CIRI
//                   Find_Circ
//                   Circ_Explorer
void ClsResultComparison::CheckIntersetBTCircExplorerAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                     string strCircExplorer, string strMyProgram, ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CIRCExplorer And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CIRCExplorer And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    CompareOrgResultWithOtherTool(strCircExplorerHit, strMyProgram, swCIRCExplorer, swMyProgram, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndCIRI(vector<St_Row_Chrom>& vChrom,
                                     string strCircExplorer, string strCIRI, ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strCIRI.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CIRCExplorer And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CIRCExplorer And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareOrgResultWithOtherTool(strCircExplorerHit, strCIRIHitPath, swCIRCExplorer, swCIRI, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                     string strCircExplorer, string strFindCirc, ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CIRCExplorer And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CIRCExplorer And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);


    CompareOrgResultWithOtherTool(strCircExplorerHit, strFindCircHit, swCIRCExplorer, swFind_circ, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                             string strCircExplorer, string strCircRNAFinder,
                                             ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strCircRNAFinder.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircExplorer And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircExplorer And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    //2: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strCircExplorerHit, strCircRNAFinderHit, swCIRCExplorer, swCircRNAFinder);
}

void ClsResultComparison::CheckIntersetBTCircExplorerAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                                                string strCircExplorer, string strCircMarker,
                                                                ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strCircMarker.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircExplorer And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircExplorer And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    //2: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareOrgResultWithOtherTool(strCircExplorerHit, strCircMarkerHit, swCIRCExplorer, swCircMarker);
}
//*****************************For CircRNAFinder Compare with others****************************
void ClsResultComparison::CheckIntersetBTCircRNAFinderAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                     string strCircRNAFinder, string strMyProgram, ofstream* pOfs)
{
    if(access(strCircRNAFinder.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CircRNAFinder And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CircRNAFinder And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strCircRNAFinderHit, strMyProgram, swCircRNAFinder, swMyProgram, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircRNAFinderAndCIRI(vector<St_Row_Chrom>& vChrom,
                                     string strCircRNAFinder, string strCIRI, ofstream* pOfs)
{
    if(access(strCircRNAFinder.c_str(), 0) != 0 || access(strCIRI.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircRNAFinder And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircRNAFinder And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareOrgResultWithOtherTool(strCircRNAFinderHit, strCIRIHitPath, swCircRNAFinder, swCIRI, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircRNAFinderAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                     string strCircRNAFinder, string strFindCirc, ofstream* pOfs)
{
    if(access(strCircRNAFinder.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircRNAFinder And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircRNAFinder And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);


    CompareOrgResultWithOtherTool(strCircRNAFinderHit, strFindCircHit, swCircRNAFinder, swFind_circ, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircRNAFinderAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                             string strCircRNAFinder, string strCircExplorer,
                                             ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strCircRNAFinder.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircRNAFinder And CircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircRNAFinder And CircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    //2: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strCircRNAFinderHit, strCircExplorerHit, swCircRNAFinder, swCIRCExplorer);
}

void ClsResultComparison::CheckIntersetBTCircRNAFinderAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                           string strCircRNAFinder, string strCircMarker,
                                           ofstream* pOfs)
{
    if(access(strCircRNAFinder.c_str(), 0) != 0 || access(strCircMarker.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircRNAFinder And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT CircRNAFinder And CircMarker<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //2: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    //2: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareOrgResultWithOtherTool(strCircRNAFinderHit, strCircMarkerHit, swCircRNAFinder, swCircMarker);
}

//**************CircMarker compare with other tools***********************************************
void ClsResultComparison::CheckIntersetBTCircMarkerAndMyPrgram(string strCircMarker, string strMyProgram,
                                                               ofstream* pOfs)
{
    if(access(strCircMarker.c_str(), 0) != 0 || access(strMyProgram.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CircMarker And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset CircMarker And MyProgram<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    CompareOrgResultWithOtherTool(strCircMarkerHit, strMyProgram, swCircMarker, swMyProgram, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircMarkerAndCIRI(vector<St_Row_Chrom>& vChrom,
                                     string strCircMarker, string strCIRI, ofstream* pOfs)
{
    if(access(strCircMarker.c_str(), 0) != 0 || access(strCIRI.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircMarker And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircMarker And CIRI<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    //2: Parse CIRI Result
    vector<St_Candidate> vCandiCIRI;
    m_pCBF->ParseCIRIResult(strCIRI, vCandiCIRI);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCandiCIRI);
    string strCIRIHitPath = "./Comparison_Result/CIRI_Full_Result.txt";
    m_pCBF->SaveHitCandi(vCandiCIRI, strCIRIHitPath);

    CompareOrgResultWithOtherTool(strCircMarkerHit, strCIRIHitPath, swCircMarker, swCIRI, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircMarkerAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                     string strCircMarker, string strFindCirc, ofstream* pOfs)
{
    if(access(strCircMarker.c_str(), 0) != 0 || access(strFindCirc.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircMarker And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset strCircMarker And Find_Circ<<<<<<<<<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    //2: Parse Find_Circ Result
    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseFindCircResult(strFindCirc, vFindCirc);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vFindCirc);
    string strFindCircHit = "./Comparison_Result/Find_Circ_Full_Result.txt";
    m_pCBF->SaveHitCandi(vFindCirc, strFindCircHit);


    CompareOrgResultWithOtherTool(strCircMarkerHit, strFindCircHit, swCircMarker, swFind_circ, pOfs);
}

void ClsResultComparison::CheckIntersetBTCircMarkerAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                             string strCircMarker, string strCircExplorer,
                                             ofstream* pOfs)
{
    if(access(strCircExplorer.c_str(), 0) != 0 || access(strCircMarker.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT strCircMarker And strCircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT strCircMarker And strCircExplorer<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);

    //2: Parse CircExplorer Result
    vector<St_Candidate> vCircExplorer;
    m_pCBF->ParseCircExplorerResult(strCircExplorer, vCircExplorer);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircExplorer);
    string strCircExplorerHit = "./Comparison_Result/CircExplorer_Hit.txt";
    m_pCBF->SaveHitCandi(vCircExplorer, strCircExplorerHit);

    CompareOrgResultWithOtherTool(strCircMarkerHit, strCircExplorerHit, swCircMarker, swCIRCExplorer);
}

void ClsResultComparison::CheckIntersetBTCircMarkerAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                           string strCircMarker, string strCircRNAFinder,
                                           ofstream* pOfs)
{
    if(access(strCircMarker.c_str(), 0) != 0 || access(strCircRNAFinder.c_str(), 0) != 0)
    {
        cout << "One of file do not existed!" << endl;
        return;
    }

    if(pOfs != NULL)
    {
        (*pOfs) << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT strCircMarker And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }
    else
    {
        cout << endl
             << "******************************************************************" << endl
             << ">>>>>>>>>>>>Check Interset BT strCircMarker And CircRNAFinder<<<<<<<<<" << endl
             << "******************************************************************"
             <<endl;
    }

    //1: Parse CircMarker Result
    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);
    string strCircMarkerHit = "./Comparison_Result/CircMarker_Hit.txt";
    m_pCBF->SaveHitCandi(vCircMarker, strCircMarkerHit);


    //2: Parse CircRNAFinder Result
    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseCircRNAFinderResult(strCircRNAFinder, vCircRNAFinder);
    ///3.1: Set Tag Value for each CIRI result
    m_pCBF->SetCircTag(vChrom, vCircRNAFinder);
    string strCircRNAFinderHit = "./Comparison_Result/CircRNAFinder_Hit.txt";
    m_pCBF->SaveHitCandi(vCircRNAFinder, strCircRNAFinderHit);

    CompareOrgResultWithOtherTool(strCircMarkerHit, strCircRNAFinderHit, swCircMarker, swCircRNAFinder);
}


//*********************************************************************************************

void ClsResultComparison::CompareOrgResultWithOtherTool(string strOrg, string strOtherTool,
                                                        En_Software enOrg, En_Software enOther,
                                                        ofstream* pOfs)
{
    //1: Parse Org
    vector<St_Candidate> vOrg;
    m_pCBF->ParseBrief(strOrg, vOrg);

    //2: Parse OtherTool
    vector<St_Candidate> vOtherTool;
    m_pCBF->ParseBrief(strOtherTool, vOtherTool);

    ofstream ofsUnHitOtherVSOrg;
    string strFileName = string("./Comparison_Result/") + "Un_Hit_" + TOOLNAME[enOther] + "_VS_" + TOOLNAME[enOrg];
    ofsUnHitOtherVSOrg.open(strFileName.c_str());

    //3: collect the number of self and number of regular
    int iOrgReg = 0;
    int iOrgSelf = 0;
    int iOrgSum = vOrg.size();
    for(vector<St_Candidate>::iterator itr = vOrg.begin(); itr != vOrg.end(); itr++)
    {
        if(itr->enType == ctSelf)
            iOrgSelf++;
        else if(itr->enType == ctRegular)
            iOrgReg++;
    }

    //4: get the hit status between other tool and org result
    int iOtherToolReg = 0;
    int iOtherToolSelf = 0;

    int iOtherToolRegHit = 0;
    int iOtherToolSelfHit = 0;
    int iOtherHitSum = 0;
    int iOtherSum = vOtherTool.size();
    for(vector<St_Candidate>::iterator itr = vOtherTool.begin(); itr != vOtherTool.end(); itr++)
    {
        //update the self and regular circular rna of other tool
        if(itr->enType == ctSelf)
            iOtherToolSelf++;
        else if(itr->enType == ctRegular)
            iOtherToolReg++;

        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vOrg.begin(); subItr != vOrg.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName)
                continue;

            bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                iOtherHitSum++;

                if(itr->enType == ctSelf)
                    iOtherToolSelfHit++;
                else if(itr->enType == ctRegular)
                    iOtherToolRegHit++;
                break;
            }
        }
        //record the unhit case
        if(!bFind)
        {
            ofsUnHitOtherVSOrg << itr->strChromName << " "
                               << IntToStr(itr->iStartPos) << " "
                               << IntToStr(itr->iEndPos) << " "
                               << IntToStr(itr->iSupportNum) << " "
                               << itr->strTag << " "
                               << itr->GetTypeString() << endl;
        }
    }

    ofsUnHitOtherVSOrg.close();

    if(pOfs != NULL)
    {
        (*pOfs) << ">>>>>>>>>>Self Case<<<<<<<<<<" << endl << endl;
        (*pOfs) << TOOLNAME[enOrg] << " Size: " << IntToStr(iOrgSelf) << endl;

        (*pOfs) << "\t size \t hit(%) \t un-hit(%)" << endl;
        (*pOfs) << TOOLNAME[enOther] << "\t" << IntToStr(iOtherToolSelf) << "\t"
             << IntToStr(iOtherToolSelfHit) << " | "
             << ((iOtherToolSelf != 0) ? GetRatio((float)iOtherToolSelfHit / iOtherToolSelf) :"NAN")
             << "\t"
             << IntToStr(iOtherToolSelf - iOtherToolSelfHit) << " | "
             << ((iOtherToolSelf != 0) ? GetRatio(1 - (float)iOtherToolSelfHit / iOtherToolSelf) : "NAN")
             << endl;

        (*pOfs) << endl << ">>>>>>>>>>Reg Case<<<<<<<<<<" << endl << endl;
        (*pOfs) << TOOLNAME[enOrg] << " Size: " << IntToStr(iOrgReg) << endl;

        (*pOfs) << "\t size \t hit(%) \t un-hit(%)" << endl;
        (*pOfs) << TOOLNAME[enOther] << "\t" << IntToStr(iOtherToolReg) << "\t"
             << IntToStr(iOtherToolRegHit) << " | "
             << ((iOtherToolReg != 0) ? GetRatio((float)iOtherToolRegHit / iOtherToolReg) : "NAN")
             << "\t"
             << IntToStr(iOtherToolReg - iOtherToolRegHit) << " | "
             << ((iOtherToolReg != 0) ? GetRatio(1 - (float)iOtherToolRegHit / iOtherToolReg) : "NAN")
             << endl;

        (*pOfs) << endl << ">>>>>>>>>>Summarize<<<<<<<<<<" << endl << endl;
        (*pOfs) << "\t Org_Sum \t Other_Sum \t Hit_Total" << endl;
        (*pOfs) << IntToStr(iOrgSum) << "\t" << IntToStr(iOtherSum) << "\t" << IntToStr(iOtherHitSum)
                << "\t" << GetRatio((float)iOtherHitSum/iOtherSum) << endl;
    }
    else
    {
        cout << ">>>>>>>>>>Self Case<<<<<<<<<<" << endl << endl;
        cout << TOOLNAME[enOrg] << " Size: " << IntToStr(iOrgSelf) << endl;

        cout << "\t size \t hit(%) \t un-hit(%)" << endl;
        cout << TOOLNAME[enOther] << "\t" << IntToStr(iOtherToolSelf) << "\t"
             << IntToStr(iOtherToolSelfHit) << " | "
             << ((iOtherToolSelf != 0) ? GetRatio((float)iOtherToolSelfHit / iOtherToolSelf) :"NAN")
             << "\t"
             << IntToStr(iOtherToolSelf - iOtherToolSelfHit) << " | "
             << ((iOtherToolSelf != 0) ? GetRatio(1 - (float)iOtherToolSelfHit / iOtherToolSelf) : "NAN")
             << endl;

        cout << endl << ">>>>>>>>>>Reg Case<<<<<<<<<<" << endl << endl;
        cout << TOOLNAME[enOrg] << " Size: " << IntToStr(iOrgReg) << endl;

        cout << "\t size \t hit(%) \t un-hit(%)" << endl;
        cout << TOOLNAME[enOther] << "\t" << IntToStr(iOtherToolReg) << "\t"
             << IntToStr(iOtherToolRegHit) << " | "
             << ((iOtherToolReg != 0) ? GetRatio((float)iOtherToolRegHit / iOtherToolReg) : "NAN")
             << "\t"
             << IntToStr(iOtherToolReg - iOtherToolRegHit) << " | "
             << ((iOtherToolReg != 0) ? GetRatio(1 - (float)iOtherToolRegHit / iOtherToolReg) : "NAN")
             << endl;

        cout << endl << ">>>>>>>>>>Summarize<<<<<<<<<<" << endl << endl;
        cout << "\t Org_Sum \t Other_Sum \t Hit_Total" << endl;
        cout << IntToStr(iOrgSum) << "\t" << IntToStr(iOtherSum) << "\t" << IntToStr(iOtherHitSum)
             << "\t" << GetRatio((float)iOtherHitSum/iOtherSum)<< endl;
    }
}


//**************Get Intersection between two tools' result*****************************************
void ClsResultComparison::GetIntersectionBT2Results(vector<St_Row_Chrom>& vChrom,
                                                    string strResult1Path, string strResult2Path,
                                                    En_Software enTool1, En_Software enTool2)
{
    //1: Parse Result1
    vector<St_Candidate> vResult1;
    GetHitExonBoundaryResult(vChrom, vResult1, strResult1Path, enTool1);

    //2: Parse Result2
    vector<St_Candidate> vResult2;
    GetHitExonBoundaryResult(vChrom, vResult2, strResult2Path, enTool2);

    //Get the intersection
    string strFile = "./Comparison_Result/Intersect_BT_" + TOOLNAME[enTool1] + "_" + TOOLNAME[enTool2] + ".txt";
    ofstream ofsIntsectBT;
    ofsIntsectBT.open(strFile.c_str());

    for(vector<St_Candidate>::iterator itr = vResult1.begin(); itr != vResult1.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vResult2.begin(); subItr != vResult2.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName)
                continue;

            bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                ofsIntsectBT<< itr->strChromName << " "
                            << IntToStr(itr->iStartPos) << " "
                            << IntToStr(itr->iEndPos) << " "
                            << IntToStr(itr->iSupportNum) << " "
                            << itr->strTag << " "
                            << itr->GetTypeString() << endl;
            }
        }
    }
    ofsIntsectBT.close();
}

void ClsResultComparison::GetBenchmark(string strFileSumPath, string strBenchmarkPath)
{
    vector<St_Candidate> vResultSum;
    m_pCBF->ParseBrief(strFileSumPath, vResultSum);

    ofstream ofsBenchmark;
    ofsBenchmark.open(strBenchmarkPath.c_str());

    for(vector<St_Candidate>::iterator itr = vResultSum.begin(); itr < vResultSum.end(); itr++)
    {
        int iTimes = 1;
        for(vector<St_Candidate>::iterator subItr = vResultSum.end() - 1; subItr > itr; subItr--)
        {
            if(itr->strChromName != subItr->strChromName)
                continue;

            bool bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                iTimes++;
                vResultSum.erase(subItr);
            }
        }

        if(iTimes >= 2)
        {
            ofsBenchmark<< itr->strChromName << " "
                        << IntToStr(itr->iStartPos) << " "
                        << IntToStr(itr->iEndPos) << " "
                        << IntToStr(itr->iSupportNum) << " "
                        << itr->strTag << " "
                        << itr->GetTypeString() << endl;
        }
    }

    ofsBenchmark.close();
}

float ClsResultComparison::CheckConsensusBasedAccuracy(string strBenchmarkPath,
                                                       En_Software enTool1, En_Software enTool2)
{
    //For Benchmark
    vector<St_Candidate> vBenchmark;
    m_pCBF->ParseBrief(strBenchmarkPath, vBenchmark);

    //For Current Tool
    string strFile = "./Comparison_Result/Intersect_BT_" + TOOLNAME[enTool1] + "_" + TOOLNAME[enTool2] + ".txt";
    vector<St_Candidate> vResult;
    m_pCBF->ParseBrief(strFile, vResult);

    int iHitNum = 0;
    for(vector<St_Candidate>::iterator itr = vResult.begin(); itr != vResult.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vBenchmark.begin(); subItr != vBenchmark.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName)
                continue;

            bool bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                iHitNum++;
            }
        }
    }

    float fAccuracy = float(iHitNum) / vResult.size();
    cout << TOOLNAME[enTool1] << " " << "Org: " <<  IntToStr(vResult.size())
                              << " " << "Hit: " <<  IntToStr(iHitNum)
                              << " " << "Accuracy Ratio: " << GetRatio(fAccuracy)
                              << endl;
    return fAccuracy;
}

float ClsResultComparison::CheckConsensusBasedSensitivity(string strBenchmarkPath,
                                                         En_Software enTool1, En_Software enTool2)
{
    //For Benchmark
    vector<St_Candidate> vBenchmark;
    m_pCBF->ParseBrief(strBenchmarkPath, vBenchmark);

    //For Current Tool
    string strFile = "./Comparison_Result/Intersect_BT_" + TOOLNAME[enTool1] + "_" + TOOLNAME[enTool2] + ".txt";
    vector<St_Candidate> vResult;
    m_pCBF->ParseBrief(strFile, vResult);

    int iHitNum = 0;
    for(vector<St_Candidate>::iterator itr = vBenchmark.begin(); itr != vBenchmark.end(); itr++)
    {
        for(vector<St_Candidate>::iterator subItr = vResult.begin(); subItr != vResult.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName)
                continue;

            bool bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                iHitNum++;
            }
        }
    }

    float fSensitivity = float(iHitNum) / vBenchmark.size();

    cout << TOOLNAME[enTool1] << " " << "Benchmark: " <<  IntToStr(vBenchmark.size())
                              << " " << "Findback: " <<  IntToStr(iHitNum)
                              << " " << "Sensitivity Ratio: " << GetRatio(fSensitivity)
                              << endl;
    return fSensitivity;
}

void ClsResultComparison::CheckIntersectionResultHitDB(vector<St_Row_Chrom>& vChrom, string strCircRNAdbPath,
                                                       En_Software enTool1, En_Software enTool2)
{
    //For DB
    vector<St_Candidate> vCandiCircRNAdb;
    m_pCBF->ParseCircRNADbBrief(strCircRNAdbPath, vCandiCircRNAdb);
    m_pCBF->SetCircTag(vChrom, vCandiCircRNAdb);
    ///cout << IntToStr(vCandiCircRNAdb.size()) << endl;

    string strCircRNADb_HitPath = "./Comparison_Result/CircRNAdbHit.txt";
    m_pCBF->SaveHitCandi(vCandiCircRNAdb, strCircRNADb_HitPath);
    m_pCBF->ParseBrief(strCircRNADb_HitPath, vCandiCircRNAdb); //这里已经是所有能够hit上的candidate了

    ///cout << IntToStr(vCandiCircRNAdb.size()) << endl;
    //For Current Tool
    string strFile = "./Comparison_Result/Intersect_BT_" + TOOLNAME[enTool1] + "_" + TOOLNAME[enTool2] + ".txt";
    vector<St_Candidate> vResult;
    m_pCBF->ParseBrief(strFile, vResult);
    ///cout << IntToStr(vResult.size()) << endl;

    int iHitNum = 0;
    for(vector<St_Candidate>::iterator itr = vCandiCircRNAdb.begin(); itr != vCandiCircRNAdb.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vResult.begin(); subItr != vResult.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName) //check if they are in same chromosome
                continue;

            bFind = m_pCBF->CompareTwoCandi(&(*itr), &(*subItr));

            if(bFind)
            {
                iHitNum++;
                break;
            }
        }
    }

    cout << TOOLNAME[enTool1] << " " << "Benchmark: " <<  IntToStr(vCandiCircRNAdb.size())
                              << " " << "Findback: " <<  IntToStr(iHitNum)
                              << " " << "Hit Ratio: " << GetRatio(float(iHitNum) / vCandiCircRNAdb.size())
                              << endl;
}

//-------------------> The function for Debug
void ClsResultComparison::CompareMyResultAndOtherSum(string strMyResult, string strCIRI, string strFind_Circ,
                                string strCIRCexplorer, string strCircRNAFinder, string strCircMarker,
                                bool bFindMyUnique)
{
    //Read those 6 dataset:
    vector<St_Candidate> vMyResult;
    m_pCBF->ParseBrief(strMyResult, vMyResult);

    vector<St_Candidate> vCIRI;
    m_pCBF->ParseBrief(strCIRI, vCIRI);

    vector<St_Candidate> vFindCirc;
    m_pCBF->ParseBrief(strFind_Circ, vFindCirc);

    vector<St_Candidate> vCIRCexplorer;
    m_pCBF->ParseBrief(strCIRCexplorer, vCIRCexplorer);

    vector<St_Candidate> vCircRNAFinder;
    m_pCBF->ParseBrief(strCircRNAFinder, vCircRNAFinder);

    vector<St_Candidate> vCircMarker;
    m_pCBF->ParseBrief(strCircMarker, vCircMarker);

    vector<St_Candidate> vSum;
//    vSum.insert(vSum.end(), vCIRI.begin(), vCIRI.end());
//    vSum.insert(vSum.end(), vFindCirc.begin(), vFindCirc.end());
//    vSum.insert(vSum.end(), vCIRCexplorer.begin(), vCIRCexplorer.end());
//    vSum.insert(vSum.end(), vCircRNAFinder.begin(), vCircRNAFinder.end());
    vSum.insert(vSum.end(), vCircMarker.begin(), vCircMarker.end());

    for(vector<St_Candidate>::iterator itr = vSum.end() - 1; itr >= vSum.begin(); itr--)
    {
        for(vector<St_Candidate>::iterator subItr = itr - 1; subItr >= vSum.begin(); subItr--)
        {
            if(itr->CheckSimilarity(*subItr))
            {
                vSum.erase(itr);
                break;
            }
        }
    }

    //Check Intersection BT my program and Sum --> 主要是看sum 有的和没有的！！
    ofstream ofsDebugShare("./Comparison_Result/__debug_shared");
    ofstream ofsDebugUnhit("./Comparison_Result/__debug_unhit");
    int iSharedNum = 0;
    int iUnhitNum = 0;

    vector<St_Candidate> vTargetSet;
    vector<St_Candidate> vPattenSet;
    if(bFindMyUnique)
    {
        vTargetSet = vMyResult;
        vPattenSet = vSum;
    }
    else
    {
       vTargetSet = vSum;
       vPattenSet = vMyResult;
    }

    for(vector<St_Candidate>::iterator itr = vTargetSet.begin(); itr != vTargetSet.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_Candidate>::iterator subItr = vPattenSet.begin(); subItr != vPattenSet.end(); subItr++)
        {
            if(itr->strChromName != subItr->strChromName)
                continue;

            if(m_pCBF->CompareTwoCandi(&(*itr), &(*subItr)))
            {
                bFind = true;
                break;
            }
        }

        if(bFind)
        {
            ofsDebugShare << itr->strChromName << " "
                          << IntToStr(itr->iStartPos) << " "
                          << IntToStr(itr->iEndPos) << " "
                          << IntToStr(itr->iSupportNum) << " "
                          << itr->strTag << " "
                          << itr->GetTypeString() << endl;
            iSharedNum++;
        }
        else //Do not find it
        {
            ofsDebugUnhit << itr->strChromName << " "
                          << IntToStr(itr->iStartPos) << " "
                          << IntToStr(itr->iEndPos) << " "
                          << IntToStr(itr->iSupportNum) << " "
                          << itr->strTag << " "
                          << itr->GetTypeString() << endl;
            iUnhitNum++;
        }
    }

    cout << "Sum       : " << IntToStr(vSum.size()) << endl;
    cout << "Shared Num: " << IntToStr(iSharedNum) << endl;
    cout << "Unhit Num : " << IntToStr(iUnhitNum) << endl;

    cout << endl << "vMyResult Size: " << IntToStr(vMyResult.size()) << endl;

    ofsDebugShare.close();
    ofsDebugUnhit.close();
}

void ClsResultComparison::GetHitExonBoundaryResult(vector<St_Row_Chrom>& vChrom, vector<St_Candidate>& vResult,
                                                   string strResultPath, En_Software enTool)
{
    vResult.clear();
    switch (enTool)
    {
        case swMyProgram:
        {
            m_pCBF->ParseBrief(strResultPath, vResult);
            break;
        }
        case swCIRI:
        {
            m_pCBF->ParseCIRIResult(strResultPath, vResult);
            m_pCBF->SetCircTag(vChrom, vResult);
            string strHitPath = "./Comparison_Result/Tmp_Hit.txt";
            m_pCBF->SaveHitCandi(vResult, strHitPath);
            m_pCBF->ParseBrief(strHitPath, vResult);
            break;
        }
        case swFind_circ:
        {
            m_pCBF->ParseFindCircResult(strResultPath, vResult);
            m_pCBF->SetCircTag(vChrom, vResult);
            string strHitPath = "./Comparison_Result/Tmp_Hit.txt";
            m_pCBF->SaveHitCandi(vResult, strHitPath);
            m_pCBF->ParseBrief(strHitPath, vResult);
            break;
        }
        case swCIRCExplorer:
        {
            m_pCBF->ParseCircExplorerResult(strResultPath, vResult);
            m_pCBF->SetCircTag(vChrom, vResult);
            string strHitPath = "./Comparison_Result/Tmp_Hit.txt";
            m_pCBF->SaveHitCandi(vResult, strHitPath);
            m_pCBF->ParseBrief(strHitPath, vResult);
            break;
        }
        case swCircRNAFinder:
        {
            m_pCBF->ParseCircExplorerResult(strResultPath, vResult);
            m_pCBF->SetCircTag(vChrom, vResult);
            string strHitPath = "./Comparison_Result/Tmp_Hit.txt";
            m_pCBF->SaveHitCandi(vResult, strHitPath);
            m_pCBF->ParseBrief(strHitPath, vResult);
            break;
        }
        case swCircMarker:
        {
            m_pCBF->ParseBrief(strResultPath, vResult);
            break;
        }
        default:
            break;
    }
}
//<-------------------

