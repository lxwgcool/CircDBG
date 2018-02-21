#ifndef CLSRESULTCOMPARISON_H
#define CLSRESULTCOMPARISON_H

#include "clscomparebasefunction.h"

enum En_Software{swMyProgram=0, swCIRI, swFind_circ, swCIRCExplorer, swCircRNAFinder, swCircMarker, swMax};

class ClsResultComparison
{
public:
    ClsResultComparison();
    ~ClsResultComparison();

public:
    void Init(St_Config& stConfig);

    //-->Compare Tools with simulator
    ///MyProgram with Simulator
    void CompareMyProgramWithSimulator(vector<St_Row_Chrom>& vChrom,
                                       string strMyProgramPath, string strSimulatorBenchmark);

    ///CIRI with Simulator
    void CompareCIRIWithSimulator(vector<St_Row_Chrom>& vChrom,
                                  string strCIRIPath, string strSimulatorBenchmark);

    ///CIRCExplorer with Simulator
    void CompareCIRCExplorerWithSimulator(vector<St_Row_Chrom>& vChrom,
                                          string strCircExplorerPath, string strSimulatorBenchmark);

    ///Find-Circ with Simulator
    void CompareFindCircWithSimulator(vector<St_Row_Chrom>& vChrom,
                                      string strFindCircPath, string strSimulatorBenchmark);

    ///CircRNAFinder with Simulator
    void CompareCircRNAFinderWithSimulator(vector<St_Row_Chrom>& vChrom,
                                           string strCircRNAFinder, string strSimulatorBenchmark);

    ///eCircMarker with Simulator
    void CompareCircMarkerWithSimulator(vector<St_Row_Chrom>& vChrom,
                                        string strCircMarkerPath, string strSimulatorBenchmark);
    //<---

    //--->
    ///Compare with the standard database: circRNAdb
    void CompareMyProgramAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                      string strMyProgramPath, string strCircRNAdbPath,
                                      string strChromName="", ofstream* pOfs=NULL);

    void CompareCIRIAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                 string strCIRIPath, string strCircRNAdbPath,
                                 string strChromName="", ofstream* pOfs=NULL);

    void CompareFindCircAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                     string strFindCircPath, string strCircRNAdbPath,
                                     string strChromName="", ofstream* pOfs=NULL);

    void CompareCircExplorerAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                         string strCircExplorerPath, string strCircRNAdbPath,
                                         string strChromName="", ofstream* pOfs=NULL);

    void CompareCircRNAFinderAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                          string strCircRNAFinder, string strCircRNAdbPath,
                                          string strChromName="", ofstream* pOfs=NULL);

    void CompareCircMarkerAndCircRNAdb(vector<St_Row_Chrom>& vChrom,
                                       string strCircMarker, string strCircRNAdbPath,
                                       string strChromName="", ofstream* pOfs=NULL);
    //<---

    ///Compare third part result
    //*********Compare My Program with Others: Ciri, Find_Circ, Circ_Explorer*********************
    // 比对 My program VS CIRI
    //                   Find_Circ
    //                   Circ_Explorer
    void CheckIntersetBTMyProgramAndCIRI(vector<St_Row_Chrom>& vChrom,
                                         string strMyProgram, string strCIRI, ofstream* pOfs=NULL);

    void CheckIntersetBTMyProgramAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                             string strMyProgram, string strFindCirc, ofstream* pOfs=NULL);

    void CheckIntersetBTMyProgramAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                                 string strMyProgram, string strCircExplorer,
                                                 ofstream* pOfs=NULL);

    void CheckIntersetBTMyProgramAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                                 string strCIRI, string strCircRNAFinder,
                                                 ofstream* pOfs=NULL);

    void CheckIntersetBTMyProgramAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                               string strCIRI, string strCircMarker,
                                               ofstream* pOfs=NULL);


    //*********Compare CIRI with Others: MyProgram, Find_Circ, Circ_Explorer*********************
    // 比对 My program VS CIRI
    //                   Find_Circ
    //                   Circ_Explorer
    void CheckIntersetBTCIRIAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                         string strCIRI, string strMyProgram, ofstream* pOfs=NULL);
    void CheckIntersetBTCIRIAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                         string strCIRI, string strFindCirc, ofstream* pOfs=NULL);
    void CheckIntersetBTCIRIAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                         string strCIRI, string strCircExplorer, ofstream* pOfs=NULL);
    void CheckIntersetBTCIRIAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                             string strCIRI, string strCircRNAFinder,
                                             ofstream* pOfs=NULL);

    void CheckIntersetBTCIRIAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                          string strCIRI, string strCircMarker,
                                          ofstream* pOfs=NULL);

    //*********Compare Find_Circ with Others: Ciri, MyProgram, Circ_Explorer*********************
    // 比对 My program VS CIRI
    //                   Find_Circ
    //                   Circ_Explorer
    void CheckIntersetBTFindCircAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                         string strFindCirc, string strMyProgram, ofstream* pOfs=NULL);
    void CheckIntersetBTFindCircAndCIRI(vector<St_Row_Chrom>& vChrom,
                                         string strFindCirc, string strCIRI, ofstream* pOfs=NULL);
    void CheckIntersetBTFindCircAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                         string strFindCirc, string strCircExplorer, ofstream* pOfs=NULL);
    void CheckIntersetBTFindCircAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                                 string strFindCirc, string strCircRNAFinder,
                                                 ofstream* pOfs=NULL);
    void CheckIntersetBTFindCircAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                              string strFindCirc, string strCircMarker,
                                              ofstream* pOfs=NULL);

    //*********Compare Circ_Explorer with Others: MyProgram, Ciri, Find_Circ*********************
    // 比对 My program VS CIRI
    //                   Find_Circ
    //                   Circ_Explorer
    void CheckIntersetBTCircExplorerAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                         string strCircExplorer, string strMyProgram, ofstream* pOfs=NULL);
    void CheckIntersetBTCircExplorerAndCIRI(vector<St_Row_Chrom>& vChrom,
                                         string strCircExplorer, string strCIRI, ofstream* pOfs=NULL);
    void CheckIntersetBTCircExplorerAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                         string strCircExplorer, string strFindCirc, ofstream* pOfs=NULL);
    void CheckIntersetBTCircExplorerAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                                 string strCircExplorer, string strCircRNAFinder,
                                                 ofstream* pOfs=NULL);

    void CheckIntersetBTCircExplorerAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                               string strCircExplorer, string strCircMarker,
                                               ofstream* pOfs=NULL);

    //*********Compare CircRNAFinder with Others: MyProgram, Ciri, Find_Circ, Circ_Explorer, CircMarker*********************
    // 比对 My program VS CIRI
    //                   Find_Circ
    //                   Circ_Explorer
    void CheckIntersetBTCircRNAFinderAndMyPrgram(vector<St_Row_Chrom>& vChrom,
                                         string strCircRNAFinder, string strMyProgram, ofstream* pOfs=NULL);
    void CheckIntersetBTCircRNAFinderAndCIRI(vector<St_Row_Chrom>& vChrom,
                                         string strCircRNAFinder, string strCIRI, ofstream* pOfs=NULL);
    void CheckIntersetBTCircRNAFinderAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                         string strCircRNAFinder, string strFindCirc, ofstream* pOfs=NULL);
    void CheckIntersetBTCircRNAFinderAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                                 string strCircRNAFinder, string strCircExplorer,
                                                 ofstream* pOfs=NULL);
    void CheckIntersetBTCircRNAFinderAndCircMarker(vector<St_Row_Chrom>& vChrom,
                                               string strCircRNAFinder, string strCircMarker,
                                               ofstream* pOfs=NULL);

    //**************CircMarker compare with other tools***********************************************
    void CheckIntersetBTCircMarkerAndMyPrgram(string strCircMarker, string strMyProgram,
                                              ofstream* pOfs=NULL);

    void CheckIntersetBTCircMarkerAndCIRI(vector<St_Row_Chrom>& vChrom,
                                         string strCircMarker, string strCIRI, ofstream* pOfs=NULL);

    void CheckIntersetBTCircMarkerAndFindCirc(vector<St_Row_Chrom>& vChrom,
                                         string strCircMarker, string strFindCirc, ofstream* pOfs=NULL);

    void CheckIntersetBTCircMarkerAndCircExplorer(vector<St_Row_Chrom>& vChrom,
                                                 string strCircMarker, string strCircExplorer,
                                                 ofstream* pOfs=NULL);

    void CheckIntersetBTCircMarkerAndCircRNAFinder(vector<St_Row_Chrom>& vChrom,
                                               string strCircMarker, string strCircRNAFinder,
                                               ofstream* pOfs=NULL);

    //**************Get Intersection between two tools' result*****************************************
    void GetIntersectionBT2Results(vector<St_Row_Chrom>& vChrom,
                                   string strResult1Path, string strResult2Path,
                                   En_Software enTool1, En_Software enTool2);
    void GetBenchmark(string strFileSumPath, string strBenchmarkPath);

    float CheckConsensusBasedAccuracy(string strBenchmarkPath, En_Software enTool1, En_Software enTool2);

    float CheckConsensusBasedSensitivity(string strBenchmarkPath, En_Software enTool1, En_Software enTool2);

    void CheckIntersectionResultHitDB(vector<St_Row_Chrom>& vChrom, string strCircRNAdbPath,
                                      En_Software enTool1, En_Software enTool2);

    //********************For Debug ***********************
    void CompareMyResultAndOtherSum(string strMyResult, string strCIRI, string strFind_Circ,
                                    string strCIRCexplorer, string strCircRNAFinder, string strCircMarker,
                                    bool bFindMyUnique = false);

private:
    //--> Compare Third party tools (include my result with standard circRNA database "circRNAdb")
    void CompareThirdPartyHitAndStdDBHit(vector<St_Row_Chrom>& vChrom,
                                         string strThirdpartyHitPath, string strStdDBHitPath,
                                         En_Software enSoftware,
                                         string strChromName="", ofstream* pOfs=NULL);

    void CompareOrgResultWithOtherTool(string strOrg, string strOtherTool,
                                       En_Software enOrg, En_Software enOther,
                                       ofstream* pOfs=NULL);

    void GetHitExonBoundaryResult(vector<St_Row_Chrom>& vChrom, vector<St_Candidate>& vResult,
                                  string strResultPath, En_Software enTool);

private:
    ClsCompareBaseFunction* m_pCBF;
    int m_iReadLen;
    int m_iCIRIHitNum;
    float m_fCIRITP;

    int m_iFindCircHitNUm;
    float m_fFindCircTP;

    int m_iMaxSupportNum;
};

#endif // CLSRESULTCOMPARISON_H
