#include <iostream>
#include "clsconfig.h"
//#include "clspreevaluation.h"
#include "clsdebruijngraph.h"
#include "clscircrnadetection.h"
#include "clsresultcomparison.h"
#include "../../ShareLibrary/clsfastqreader.h"
#include "../../ShareLibrary/clsfastareader.h"
#include "../../ShareLibrary/clsbasealgorithm.h"
#include <unistd.h>
#include <pthread.h>
#include <string.h>

using namespace std;

//#define NUM_THREADS 24
//#define HUMAN
//#define HEK293
#define MOUSE

void CircRNADetection(St_Config& stConfig);
void ResultComparison(St_Config& stConfig);
void FindCircByMultiThreads(ClsCircRNADetection* pCircRNADetection, St_Config& stConfig,
                            vector<St_Fastq>& vReads, ClsDeBruijnGraph* vDBG,
                            vector<St_Row_Chrom>& vChrom, vector<St_Fasta>& vFasta);

int main(int argc, char **argv)
{
    if(argc == 1 ||
       strcmp(argv[1], "-h") == 0 ||
       strcmp(argv[1], "--help") == 0)
    {
        cout << "**********************************" << endl;
        cout << "How to use CircDBG" << endl;
        cout << "./CircRNADBG ./config.ini" << endl;
        cout << "**********************************" << endl;
        return 0;
    }

    //1. Read Config File
    ClsConfig* pConfig = new ClsConfig();
    St_Config stConfig;
    pConfig->ReadConfig(stConfig, argv[1]);
    delete pConfig;
    pConfig = NULL;

//    //2. Make evaluation (Check if we could arvhieve the idea result theoretically)
//    ClsPreEvaluation* pPE = new ClsPreEvaluation(stConfig);
//    //pPE->GetValidStdCandiByChormoson(1);  // get the first chromosome
//    delete pPE;
//    pPE = NULL;

    //3. The main body of circular rna detection
    if(stConfig.bDoDetection)
    {
        cout << "*********CircRNADetection*********" << endl;
        CircRNADetection(stConfig);
    }

    //4: Make result comparison
    if(stConfig.bDoComparison)
    {
        cout << "*********ResultComparison*********" << endl;
        ResultComparison(stConfig);
    }

    return 0;
}

void CircRNADetection(St_Config& stConfig)
{
    //---> Do First by some unspoken reason <---
    string strCmd = (string)"mkdir -p " + get_current_dir_name() + "/Detection_Result";
    system(strCmd.c_str());

    //1. Parse GTF File
    cout << "Step 1: Parse GTF File" << endl;
    vector<St_Row_Chrom> vChrom;
    ClsGTFParse* pGTFParse = new ClsGTFParse();
    pGTFParse->Init(stConfig.strGtfPath, stConfig.strRefPath,
                    stConfig.iKmerLen, stConfig.iReadLen, stConfig.fKmerRatio);
    pGTFParse->ReadGTF(vChrom, stConfig.strGtfPath); //Done
    pGTFParse->GetTagValue(vChrom);
    delete pGTFParse;
    pGTFParse = NULL;

    //2: try to create De Bruijn Graph based on Reference and annotation file
    cout << "Step 2: Create DBG --> Move create DBG to the thread" << endl;
    ClsDeBruijnGraph* pDBG = new ClsDeBruijnGraph();
    pDBG->Init(stConfig); //Move create DBG to the thread

    //3: Read Reads
    cout << "Read Reads"  << endl;
    ClsFastqReader* pFqReader = new ClsFastqReader();
    vector<St_Fastq> vReads;
    pFqReader->ReadFastqFile(stConfig.strReads1Path, vReads, true);
    pFqReader->ReadFastqFile(stConfig.strReads2Path, vReads, true);
    delete pFqReader;
    pFqReader = NULL;

    //4: Read Ref
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(stConfig.strRefPath, vFasta);
    delete pFastaReader;
    pFastaReader = NULL;

    //4. Circular RNA Detection  --> Find Circ
    cout << "Step 3: Detect Circular RNA" << endl;
    strCmd = (string)"mkdir -p " + get_current_dir_name() + "/Detection_Result";
    system(strCmd.c_str());
    strCmd = "rm ./Detection_Result/*";
    system(strCmd.c_str());
    ClsCircRNADetection* pCircRNADetection = new ClsCircRNADetection();
    pCircRNADetection->Init(stConfig, vChrom);
    FindCircByMultiThreads(pCircRNADetection, stConfig, vReads, pDBG,
                           vChrom, vFasta); //Detect Circular RNA by Miltiple Threads

    delete pDBG;
    pDBG = NULL;

    delete pCircRNADetection;
    pCircRNADetection = NULL;

    vChrom.clear();
}

struct St_FindCircDBG
{
    vector<St_Fastq>* pvReads;
    ClsDeBruijnGraph* pDBG;
    ClsCircRNADetection* pCircRNADetection;
    St_Row_Chrom* pChrom;
    St_Fasta* pRef;
    int iChromIndex;
    string strChromName;

    St_FindCircDBG()
    {
        Reset();
    }

    void Init(vector<St_Fastq>* pV1, ClsDeBruijnGraph* pV2, ClsCircRNADetection* pV3,
              St_Row_Chrom* pV4, St_Fasta* pV5, int iV6, string strV7)
    {
        pvReads = pV1;
        pDBG = pV2;
        pCircRNADetection = pV3;
        pChrom = pV4;
        pRef = pV5;
        iChromIndex = iV6;
        strChromName = strV7;
    }

    void Reset()
    {
        pvReads = NULL;
        pDBG = NULL;
        pCircRNADetection = NULL;
        pChrom = NULL;
        pRef = NULL;
        iChromIndex = -1;
        strChromName = "";
    }
};

void* FindCircForSingleChrom(void* pValue)
{
    St_FindCircDBG* pFindCircDBG = (St_FindCircDBG*)pValue;

    //1: Create Graph
    unordered_map<unsigned int, St_Node> mpDBG;
    pFindCircDBG->pDBG->BuildGraph(mpDBG, pFindCircDBG->pChrom,
                                   pFindCircDBG->pRef, pFindCircDBG->iChromIndex);

    //2: Find Circular RNA
    pFindCircDBG->pCircRNADetection->FindCirc(mpDBG, pFindCircDBG->pRef,
                                              *(pFindCircDBG->pvReads),
                                              pFindCircDBG->strChromName);
    mpDBG.clear();
    pthread_exit(NULL);
}

void FindCircByMultiThreads(ClsCircRNADetection* pCircRNADetection, St_Config& stConfig,
                            vector<St_Fastq>& vReads, ClsDeBruijnGraph* pDBG,
                            vector<St_Row_Chrom>& vChrom, vector<St_Fasta>& vFasta)
{
    //1: Get related reference and chromosome.
    int iCurChromIndex = 0;
    int iThreadsNum = 0;
    vector<St_FindCircDBG> vFindCircDBG;
    St_FindCircDBG stFindCircDBG;
    for(vector<St_Row_Chrom>::iterator itrChrom = vChrom.begin(); itrChrom != vChrom.end();
        itrChrom++, iCurChromIndex++)
    {
        if(iThreadsNum >= (int)vFasta.size())
            break;

        for(vector<St_Fasta>::iterator itrRef = vFasta.begin(); itrRef != vFasta.end(); itrRef++)
        {
            string strRefName = "";
            if(itrRef->strName.find(' ') == string::npos)
                strRefName = itrRef->strName;
            else
                strRefName = itrRef->strName.substr(0, itrRef->strName.find(' '));

            if(strRefName == itrChrom->strName)
            {
                stFindCircDBG.Reset();
                stFindCircDBG.Init(&vReads, pDBG, pCircRNADetection,
                                   &(*itrChrom), &(*itrRef), iCurChromIndex, strRefName);
                vFindCircDBG.push_back(stFindCircDBG);
                iThreadsNum++;
                break;
            }
        }
    }

    cout << iThreadsNum << " Number of thread will be created!" << endl;
    if(iThreadsNum <= stConfig.iThreadsNum)
        cout << "Good: The number of real thread is less than configuration!" << endl;
    else
        cout << "Warnning: The number of real thread is more than configuration! The speed may be slow down" << endl;


    //Step 1: Define how may threads you want to use
    pthread_t threads[iThreadsNum];
    pthread_attr_t attr;

    //Step 2: Initialize and set thread joinable
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    //Step 3: Do the main body of multiple threads function
    for(int i = 0; i < iThreadsNum; i++)
    {
        //cout << "main() : creating thread, " << i << endl;

        int rc = pthread_create(&threads[i], NULL, FindCircForSingleChrom, (void *)&vFindCircDBG[i]);

        if(rc)
        {
            cout << "Error:unable to create thread," << rc << endl;
            exit(-1);
        }
    }

    //4: free attribute and wait for the other threads
    pthread_attr_destroy(&attr);

    void *status;
    //5: Join those threads together --> to make sure the remaining part of the code will be run only after those threads been finished
    for(int i = 0; i < iThreadsNum; i++)
    {
        int rc = pthread_join(threads[i], &status);
        if (rc)
        {
            cout << "Error:unable to join," << rc << endl;
            exit(-1);
        }

        //cout << "Main: completed thread id :" << i ;
        //cout << "  exiting with status :" << status << endl;
    }

    //6: Do the remaning thing
    cout << endl << "ALL SET!!!" << endl;

    cout << "Combine the result together" << endl;
    string strCmd = "rm ./Detection_Result/Brief_sum.txt";
    system(strCmd.c_str());
    strCmd = "cat ./Detection_Result/Brief* > ./Detection_Result/Brief_sum.txt";
    system(strCmd.c_str());
}

void ComparisonSimulator(St_Config& stConfig, vector<St_Row_Chrom>& vChrom,
                         ClsResultComparison* pResultCompare);

void ComparisonNormal(St_Config& stConfig, vector<St_Row_Chrom>& vChrom,
                      ClsResultComparison* pResultCompare);

void ComparisonIntersection(St_Config& stConfig, vector<St_Row_Chrom>& vChrom,
                            ClsResultComparison* pResultCompare);

void ResultComparison(St_Config& stConfig)
{
    string strCmd = (string)"mkdir -p " + get_current_dir_name() + "/Comparison_Result";
    system(strCmd.c_str());
    strCmd = "rm ./Comparison_Result/*";
    system(strCmd.c_str());

    //1. Parse GTF File
    vector<St_Row_Chrom> vChrom;
    ClsGTFParse* pGTFParse = new ClsGTFParse();
    pGTFParse->Init(stConfig.strGtfPath, stConfig.strRefPath, stConfig.iKmerLen,
                    stConfig.iReadLen, stConfig.fKmerRatio);
    pGTFParse->ReadGTF(vChrom, stConfig.strGtfPath); //Done
    cout << "GetTagValue(vChrom)" << endl;
    pGTFParse->GetTagValue(vChrom);

//    for(vector<St_Row_Chrom>::iterator itr = vChrom.begin(); itr != vChrom.end(); itr++)
//    {
//        string strRefPath = "RNA_Ref_" + itr->strName + ".txt";
//        cout << strRefPath << ": " << IntToStr(itr->vRG.size()) << endl;
//        pGTFParse->GetRNARef(strRefPath, itr->strName, itr->vRG);
//    }

    delete pGTFParse;
    pGTFParse = NULL;

    ClsResultComparison* pResultCompare = new ClsResultComparison();
    pResultCompare->Init(stConfig);

    ///Compare with simulation Benchmark
    cout << endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> Compare with simulation Benchmark <<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    //ComparisonSimulator(stConfig, vChrom, pResultCompare);

    ///Compare with CircRNAdb
    cout << endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> Compare with CircRNAdb <<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    //ComparisonNormal(stConfig, vChrom, pResultCompare);

    ///Normal Comparison (intersection)
    cout << endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> Normal Comparison (intersection) <<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    ComparisonIntersection(stConfig, vChrom, pResultCompare);

    delete pResultCompare;
    pResultCompare = NULL;
}

void ComparisonSimulator(St_Config& stConfig, vector<St_Row_Chrom>& vChrom,
                         ClsResultComparison* pResultCompare)
{
    //MyProgram with Simulator
    pResultCompare->CompareMyProgramWithSimulator(vChrom, stConfig.strMyResult,
                                                  stConfig.strSimulatorBenchmark);

    //CIRI with Simulator
    pResultCompare->CompareCIRIWithSimulator(vChrom, stConfig.strCiriResult,
                                             stConfig.strSimulatorBenchmark);

    //CIRCExplorer with Simulator
    pResultCompare->CompareCIRCExplorerWithSimulator(vChrom, stConfig.strCIRCExplorerResult,
                                                     stConfig.strSimulatorBenchmark);

    //Find-Circ with Simulator
    pResultCompare->CompareFindCircWithSimulator(vChrom, stConfig.strFindCircResult,
                                                 stConfig.strSimulatorBenchmark);

    //CircRNAFinder with Simulator
    pResultCompare->CompareCircRNAFinderWithSimulator(vChrom, stConfig.strCircRNAFinderResult,
                                                      stConfig.strSimulatorBenchmark);

    //CircMarker with Simulator
    pResultCompare->CompareCircMarkerWithSimulator(vChrom, stConfig.strCircMarkerResult,
                                                   stConfig.strSimulatorBenchmark);
}


void ComparisonIntersection(St_Config& stConfig, vector<St_Row_Chrom>& vChrom,
                            ClsResultComparison* pResultCompare)
{
    //Clear the old file
    string strCmd = (string)"mkdir -p " + get_current_dir_name() + "/Comparison_Result";
    system(strCmd.c_str());
    strCmd = "rm ./Comparison_Result/Intersect_*";
    system(strCmd.c_str());

#ifdef HUMAN
    ///----> Human
    //My Program
    string strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/glioblastoma/SRR5133906/Detection_Result/Brief_sum.txt";
    string strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/glioblastoma/SRR4095542/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swMyProgram, swMyProgram);

    //CIRI
    strResult1 = "/scratch/xil14026/circRNA/data/CIRI/glioblastoma/SRR5133906/whole_chrom/CIRI_Result.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRI/glioblastoma/SRR4095542/CIRI_Result.txt";;
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRI, swCIRI);

    //Find-Circ
    strResult1 = "/scratch/xil14026/circRNA/data/find_circ/glioblastoma/whole_chrom/SRR5133906/GRCh37.75.chromosome_test_out/circ_candidates.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/find_circ/glioblastoma/whole_chrom/SRR4095542/GRCh37.75.chromosome_test_out/circ_candidates.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swFind_circ, swFind_circ);

    //CIRCExplorer
    strResult1 = "/scratch/xil14026/circRNA/data/CIRCexplorer/glioblastoma/SRR5133906/CIRCexplorer_circ.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRCexplorer/glioblastoma/SRR4095542/CIRCexplorer_circ.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRCExplorer, swCIRCExplorer);

    //CircRNAFinder
    strResult1 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/glioblastoma/SRR5133906/resultoutputs_filteredJunctions.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/glioblastoma/SRR4095542/resultoutputs_filteredJunctions.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircRNAFinder, swCircRNAFinder);

    //CircMarker
    strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/glioblastoma/SRR5133906/Detection_Result/Brief_sum.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/glioblastoma/SRR4095542/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircMarker, swCircMarker);
    ///<-----
#else
  #ifdef HEK293
    //My Program
    string strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/HEK293/treated/Detection_Result/Brief_sum.txt";
    string strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/HEK293/un-treated/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swMyProgram, swMyProgram);

    //CIRI
    strResult1 = "/scratch/xil14026/circRNA/data/CIRI/HEK293/treated/CIRI_Result.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRI/HEK293/un-treated/CIRI_Result.txt";;
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRI, swCIRI);

    //Find-Circ
    strResult1 = "/scratch/xil14026/circRNA/data/find_circ/HEK293/treated/GRCh37.75.chromosome_test_out/circ_candidates.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/find_circ/HEK293/un-treated/GRCh37.75.chromosome_test_out/circ_candidates.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swFind_circ, swFind_circ);

    //CIRCExplorer
    strResult1 = "/scratch/xil14026/circRNA/data/CIRCexplorer/HEK293/treated/CIRCexplorer_circ.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRCexplorer/HEK293/un-treated/CIRCexplorer_circ.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRCExplorer, swCIRCExplorer);

    //CircRNAFinder
    strResult1 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/HEK293/treated/resultoutputs_filteredJunctions.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/HEK293/un-treated/resultoutputs_filteredJunctions.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircRNAFinder, swCircRNAFinder);

    //CircMarker
    strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/HEK293/treated/Detection_Result/Brief_sum.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/HEK293/un-treated/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircMarker, swCircMarker);
  #else
    #ifndef MOUSE
    //My Program
    string strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/mouse_brain/SRR2185851_untreated/Detection_Result/Brief_sum.txt";
    string strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/mouse_brain/SRR2219951_Treated/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swMyProgram, swMyProgram);

    //CIRI
    strResult1 = "/scratch/xil14026/circRNA/data/CIRI/mouse_brain/whole_genome/ribominus_SRR2185851/CIRI_Result.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRI/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/CIRI_Result.txt";;
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRI, swCIRI);

    //Find-Circ
    strResult1 = "/scratch/xil14026/circRNA/data/find_circ/mouse_brain/whole_genome/ribominus_SRR2185851/GRCm38.chromosome_test_out/circ_candidates.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/find_circ/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/GRCm38.chromosome_test_out/circ_candidates.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swFind_circ, swFind_circ);

    //CIRCExplorer
    strResult1 = "/scratch/xil14026/circRNA/data/CIRCexplorer/mouse_brain/whole_genome/ribominus_SRR2185851/CIRCexplorer_circ.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRCexplorer/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/CIRCexplorer_circ.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRCExplorer, swCIRCExplorer);

    //CircRNAFinder
    strResult1 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/mouse_brain/ribominus_SRR2185851/resultoutputs_filteredJunctions.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/mouse_brain/RNaseR_PolyA_SRR2219951/resultoutputs_filteredJunctions.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircRNAFinder, swCircRNAFinder);

    //CircMarker
    strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/mouse_brain/whole_genome/ribominus_SRR2185851/Detection_Result/Brief_sum.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircMarker, swCircMarker);
    #endif
  #endif
#endif

#ifdef MOUSE
    //My Program
    string strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/mouse_brain/SRR2219951_Treated/Detection_Result/Brief_sum.txt";
    string strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircRNADBG/CircRNA_DBG/mouse_brain/SRR2185851_untreated/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swMyProgram, swMyProgram);

    //CIRI
    strResult1 = "/scratch/xil14026/circRNA/data/CIRI/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/CIRI_Result.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRI/mouse_brain/whole_genome/ribominus_SRR2185851/CIRI_Result.txt";;
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRI, swCIRI);

    //Find-Circ
    strResult1 = "/scratch/xil14026/circRNA/data/find_circ/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/GRCm38.chromosome_test_out/circ_candidates.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/find_circ/mouse_brain/whole_genome/ribominus_SRR2185851/GRCm38.chromosome_test_out/circ_candidates.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swFind_circ, swFind_circ);

    //CIRCExplorer
    strResult1 = "/scratch/xil14026/circRNA/data/CIRCexplorer/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/CIRCexplorer_circ.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/CIRCexplorer/mouse_brain/whole_genome/ribominus_SRR2185851/CIRCexplorer_circ.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCIRCExplorer, swCIRCExplorer);

    //CircRNAFinder
    strResult1 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/mouse_brain/RNaseR_PolyA_SRR2219951/resultoutputs_filteredJunctions.bed";
    strResult2 = "/scratch/xil14026/circRNA/data/OtherTools/circRNA_finder/jobs/mouse_brain/ribominus_SRR2185851/resultoutputs_filteredJunctions.bed";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircRNAFinder, swCircRNAFinder);

    //CircMarker
    strResult1 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/mouse_brain/whole_genome/RNaseR_PolyA_SRR2219951/Detection_Result/Brief_sum.txt";
    strResult2 = "/scratch/xil14026/circRNA/data/Test_RealData/2017_autumn/CircMarker/CircRNA_Detect_Draft/Jobs/mouse_brain/whole_genome/ribominus_SRR2185851/Detection_Result/Brief_sum.txt";
    pResultCompare->GetIntersectionBT2Results(vChrom, strResult1, strResult2, swCircMarker, swCircMarker);
#endif

    ///Combine the intersection result above
    string strIntersectionSum = "./Comparison_Result/Intersect_BT_Sum.txt";
    strCmd = "cat ./Comparison_Result/Intersect_BT* > " + strIntersectionSum;
    system(strCmd.c_str());

    string strBenchmarkPath = "./Comparison_Result/Intersect_Benchmark.txt";
    pResultCompare->GetBenchmark(strIntersectionSum, strBenchmarkPath);

    //Consensus based accuracy
    cout<< endl << "********** Consensus based Accuracy **********" << endl;
    float fPreciseMyProgram = pResultCompare->CheckConsensusBasedAccuracy(strBenchmarkPath, swMyProgram, swMyProgram);
    float fPreciseCIRI = pResultCompare->CheckConsensusBasedAccuracy(strBenchmarkPath, swCIRI, swCIRI);
    float fPreciseFindCirc = pResultCompare->CheckConsensusBasedAccuracy(strBenchmarkPath, swFind_circ, swFind_circ);
    float fPreciseCIRCExplorer = pResultCompare->CheckConsensusBasedAccuracy(strBenchmarkPath, swCIRCExplorer, swCIRCExplorer);
    float fPreciseCircRNAFinder = pResultCompare->CheckConsensusBasedAccuracy(strBenchmarkPath, swCircRNAFinder, swCircRNAFinder);
    float fPreciseCircMarker = pResultCompare->CheckConsensusBasedAccuracy(strBenchmarkPath, swCircMarker, swCircMarker);
    cout << endl;

    //Consensus based sensitivity
    cout<< endl << "********** Consensus based Sensitivity **********" << endl;
    float fRecallMyProgram = pResultCompare->CheckConsensusBasedSensitivity(strBenchmarkPath, swMyProgram, swMyProgram);
    float fRecallCIRI = pResultCompare->CheckConsensusBasedSensitivity(strBenchmarkPath, swCIRI, swCIRI);
    float fRecallFindCirc = pResultCompare->CheckConsensusBasedSensitivity(strBenchmarkPath, swFind_circ, swFind_circ);
    float fRecallCIRCExplorer = pResultCompare->CheckConsensusBasedSensitivity(strBenchmarkPath, swCIRCExplorer, swCIRCExplorer);
    float fRecallCircRNAFinder = pResultCompare->CheckConsensusBasedSensitivity(strBenchmarkPath, swCircRNAFinder, swCircRNAFinder);
    float fRecallCircMarker = pResultCompare->CheckConsensusBasedSensitivity(strBenchmarkPath, swCircMarker, swCircMarker);
    cout << endl;

    cout<< endl << "********** F1 Scoure **********" << endl;
    float fF1 = (float)2 * fPreciseMyProgram * fRecallMyProgram / (fPreciseMyProgram + fRecallMyProgram);
    cout << "MyProgram    : " << "-->Precision: " << GetRatio(fPreciseMyProgram) << " "
                              << "-->Recall: " << GetRatio(fRecallMyProgram) << " "
                              << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;

    fF1 = (float)2 * fPreciseCIRI * fRecallCIRI / (fPreciseCIRI + fRecallCIRI);
    cout << "CIRI         : " << "-->Precision: " << GetRatio(fPreciseCIRI) << " "
                              << "-->Recall: " << GetRatio(fRecallCIRI) << " "
                              << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;

    fF1 = (float)2 * fPreciseFindCirc * fRecallFindCirc / (fPreciseFindCirc + fRecallFindCirc);
    cout << "Find-Circ    : " << "-->Precision: " << GetRatio(fPreciseFindCirc) << " "
                              << "-->Recall: " << GetRatio(fRecallFindCirc) << " "
                              << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;


    fF1 = (float)2 * fPreciseCIRCExplorer * fRecallCIRCExplorer / (fPreciseCIRCExplorer + fRecallCIRCExplorer);
    cout << "CIRCExplorer : " << "-->Precision: " << GetRatio(fPreciseCIRCExplorer) << " "
                              << "-->Recall: " << GetRatio(fRecallCIRCExplorer) << " "
                              << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;

    fF1 = (float)2 * fPreciseCircRNAFinder * fRecallCircRNAFinder / (fPreciseCircRNAFinder + fRecallCircRNAFinder);
    cout << "CircRNAFinder: " << "-->Precision: " << GetRatio(fPreciseCircRNAFinder) << " "
                              << "-->Recall: " << GetRatio(fRecallCircRNAFinder) << " "
                              << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;

    fF1 = (float)2 * fPreciseCircMarker * fRecallCircMarker / (fPreciseCircMarker + fRecallCircMarker);
    cout << "CircMarker   : " << "-->Precision: " << GetRatio(fPreciseCircMarker) << " "
                              << "-->Recall: " << GetRatio(fRecallCircMarker) << " "
                              << "-->F1 Score: " << GetRatio(fF1, 2, false) << endl;

    //Check Intersection between after intersection result and standard DB --> circRNAdb
    cout<< endl << "********** Intersection between after intersection result and standard DB **********" << endl;
    pResultCompare->CheckIntersectionResultHitDB(vChrom, stConfig.strCircRNAdbResult, swMyProgram, swMyProgram);
    pResultCompare->CheckIntersectionResultHitDB(vChrom, stConfig.strCircRNAdbResult, swCIRI, swCIRI);
    pResultCompare->CheckIntersectionResultHitDB(vChrom, stConfig.strCircRNAdbResult, swFind_circ, swFind_circ);
    pResultCompare->CheckIntersectionResultHitDB(vChrom, stConfig.strCircRNAdbResult, swCIRCExplorer, swCIRCExplorer);
    pResultCompare->CheckIntersectionResultHitDB(vChrom, stConfig.strCircRNAdbResult, swCircRNAFinder, swCircRNAFinder);
    pResultCompare->CheckIntersectionResultHitDB(vChrom, stConfig.strCircRNAdbResult, swCircMarker, swCircMarker);
}

void ComparisonNormal(St_Config& stConfig, vector<St_Row_Chrom>& vChrom, ClsResultComparison* pResultCompare)
{
    ///Compare my result and circRNAdb
    cout << "Compare my result and circRNAdb" << endl;
    pResultCompare->CompareMyProgramAndCircRNAdb(vChrom, stConfig.strMyResult,
                                                 stConfig.strCircRNAdbResult);

    ///Compare CIRI and circRNAdb
    cout << "Compare CIRI and circRNAdb" << endl;
    pResultCompare->CompareCIRIAndCircRNAdb(vChrom, stConfig.strCiriResult,
                                            stConfig.strCircRNAdbResult);

    ///Compare Find Circ and circRNAdb
    cout << "Compare Find Circ and circRNAdbs" << endl;
    pResultCompare->CompareFindCircAndCircRNAdb(vChrom, stConfig.strFindCircResult,
                                                stConfig.strCircRNAdbResult);

    ///Compare CIRCexplorer and CircRNAdb
    cout << "Compare CIRCexplorer and CircRNAdb" << endl;
    pResultCompare->CompareCircExplorerAndCircRNAdb(vChrom, stConfig.strCIRCExplorerResult,
                                                    stConfig.strCircRNAdbResult);

    ///Compare CircRNAFinder and CircRNAdb
    cout << "Compare CircRNAFinder and CircRNAdb" << endl;
    pResultCompare->CompareCircRNAFinderAndCircRNAdb(vChrom, stConfig.strCircRNAFinderResult,
                                                    stConfig.strCircRNAdbResult);

    ///Compare CircMarker and CircRNAdb
    cout << "Compare CircMarke and CircRNAdb" << endl;
    pResultCompare->CompareCircMarkerAndCircRNAdb(vChrom, stConfig.strCircMarkerResult,
                                                  stConfig.strCircRNAdbResult);


//    // Make the comparison -->
//    pResultCompare->CompareMyResultAndOtherSum("./Comparison_Result/My_Program_Hit.txt",
//                                               "./Comparison_Result/CIRI_Org_Hit.txt",
//                                               "./Comparison_Result/Find_Circ_Org_Hit.txt",
//                                               "./Comparison_Result/Circ_Explorer_Org_Hit.txt",
//                                               "./Comparison_Result/CircRNAFinder_Org_Hit.txt",
//                                               "./Comparison_Result/CircMarker_Hit.txt", true);


    //****Compare the intersection of the result from two different third party tools*****
    ///My Program
    pResultCompare->CheckIntersetBTMyProgramAndCIRI(vChrom,
                                                    stConfig.strMyResult, stConfig.strCiriResult);
    pResultCompare->CheckIntersetBTMyProgramAndFindCirc(vChrom,
                                                    stConfig.strMyResult, stConfig.strFindCircResult);
    pResultCompare->CheckIntersetBTMyProgramAndCircExplorer(vChrom,
                                                    stConfig.strMyResult, stConfig.strCIRCExplorerResult);
    pResultCompare->CheckIntersetBTMyProgramAndCircRNAFinder(vChrom,
                                                    stConfig.strMyResult, stConfig.strCircRNAFinderResult);
    pResultCompare->CheckIntersetBTMyProgramAndCircMarker(vChrom,
                                                    stConfig.strMyResult, stConfig.strCircMarkerResult);


    ///CIRI
    pResultCompare->CheckIntersetBTCIRIAndMyPrgram(vChrom,
                                                   stConfig.strCiriResult, stConfig.strMyResult);
    pResultCompare->CheckIntersetBTCIRIAndFindCirc(vChrom,
                                                   stConfig.strCiriResult, stConfig.strFindCircResult);
    pResultCompare->CheckIntersetBTCIRIAndCircExplorer(vChrom,
                                                       stConfig.strCiriResult, stConfig.strCIRCExplorerResult);
    pResultCompare->CheckIntersetBTCIRIAndCircRNAFinder(vChrom,
                                                       stConfig.strCiriResult, stConfig.strCircRNAFinderResult);
    pResultCompare->CheckIntersetBTCIRIAndCircMarker(vChrom,
                                                       stConfig.strCiriResult, stConfig.strCircMarkerResult);

    ///FindCirc
    pResultCompare->CheckIntersetBTFindCircAndMyPrgram(vChrom,
                                                       stConfig.strFindCircResult, stConfig.strMyResult);
    pResultCompare->CheckIntersetBTFindCircAndCIRI(vChrom,
                                                   stConfig.strFindCircResult, stConfig.strCiriResult);
    pResultCompare->CheckIntersetBTFindCircAndCircExplorer(vChrom,
                                                           stConfig.strFindCircResult, stConfig.strCIRCExplorerResult);
    pResultCompare->CheckIntersetBTFindCircAndCircRNAFinder(vChrom,
                                                           stConfig.strFindCircResult, stConfig.strCircRNAFinderResult);
    pResultCompare->CheckIntersetBTFindCircAndCircMarker(vChrom,
                                                           stConfig.strFindCircResult, stConfig.strCircMarkerResult);

    ///CircExplorer
    pResultCompare->CheckIntersetBTCircExplorerAndMyPrgram(vChrom,
                                                           stConfig.strCIRCExplorerResult, stConfig.strMyResult);
    pResultCompare->CheckIntersetBTCircExplorerAndCIRI(vChrom,
                                                       stConfig.strCIRCExplorerResult, stConfig.strCiriResult);
    pResultCompare->CheckIntersetBTCircExplorerAndFindCirc(vChrom,
                                                           stConfig.strCIRCExplorerResult, stConfig.strFindCircResult);
    pResultCompare->CheckIntersetBTCircExplorerAndCircRNAFinder(vChrom,
                                                           stConfig.strCIRCExplorerResult, stConfig.strCircRNAFinderResult);
    pResultCompare->CheckIntersetBTCircExplorerAndCircMarker(vChrom,
                                                           stConfig.strCIRCExplorerResult, stConfig.strCircMarkerResult);

    ///CircRNAFinder
    pResultCompare->CheckIntersetBTCircRNAFinderAndMyPrgram(vChrom,
                                                           stConfig.strCircRNAFinderResult, stConfig.strMyResult);
    pResultCompare->CheckIntersetBTCircRNAFinderAndCIRI(vChrom,
                                                       stConfig.strCircRNAFinderResult, stConfig.strCiriResult);
    pResultCompare->CheckIntersetBTCircRNAFinderAndFindCirc(vChrom,
                                                           stConfig.strCircRNAFinderResult, stConfig.strFindCircResult);
    pResultCompare->CheckIntersetBTCircRNAFinderAndCircExplorer(vChrom,
                                                           stConfig.strCircRNAFinderResult, stConfig.strCIRCExplorerResult);
    pResultCompare->CheckIntersetBTCircRNAFinderAndCircMarker(vChrom,
                                                           stConfig.strCircRNAFinderResult, stConfig.strCircMarkerResult);

    ///CircMarker
    pResultCompare->CheckIntersetBTCircMarkerAndMyPrgram(stConfig.strCircMarkerResult, stConfig.strMyResult);
    pResultCompare->CheckIntersetBTCircMarkerAndCIRI(vChrom,
                                                       stConfig.strCircMarkerResult, stConfig.strCiriResult);
    pResultCompare->CheckIntersetBTCircMarkerAndFindCirc(vChrom,
                                                           stConfig.strCircMarkerResult, stConfig.strFindCircResult);
    pResultCompare->CheckIntersetBTCircMarkerAndCircExplorer(vChrom,
                                                           stConfig.strCircMarkerResult, stConfig.strCIRCExplorerResult);
    pResultCompare->CheckIntersetBTCircMarkerAndCircRNAFinder(vChrom,
                                                           stConfig.strCircMarkerResult, stConfig.strCircRNAFinderResult);
}

