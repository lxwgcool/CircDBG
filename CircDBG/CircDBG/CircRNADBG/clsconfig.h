#ifndef CLSCONFIG_H
#define CLSCONFIG_H

#include <string>
using namespace std;

struct St_Config
{
    //General
    string strRefPath;
    string strGtfPath;
    string strReads1Path;
    string strReads2Path;    

    //Parameters
    int iKmerLen;
    int iReadLen;
    float fKmerRatio;
    int iThreadsNum;
    bool bDoDetection;
    bool bDoComparison;
    int iMaxSupportNum;

    //Std Circular RNA Info
    string strCircRNADb;
    string strCircRNAdbResult;
    string strTissue;    

    //Result Comparsion
    string strMyResult;
    string strCiriResult;
    string strCIRCExplorerResult;
    string strFindCircResult;
    string strCircRNAFinderResult;
    string strCircMarkerResult;
    string strSimulatorBenchmark;

    St_Config():strRefPath(""), strGtfPath(""), strReads1Path(""), strReads2Path(""),
                iKmerLen(-1), iReadLen(-1), fKmerRatio(0), iThreadsNum(1),
                bDoDetection(true), bDoComparison(true), iMaxSupportNum(999),
                strCircRNADb(""), strCircRNAdbResult(""), strTissue(""),
                strMyResult(""), strCiriResult(""), strCIRCExplorerResult(""),
                strFindCircResult(""), strCircRNAFinderResult(""),
                strCircMarkerResult(""), strSimulatorBenchmark("")
    {}
};

class ClsConfig
{
public:
    ClsConfig();

public:
    void ReadConfig(St_Config& stConfig, char* cpIniPath);
    bool CheckConfig(St_Config& stConfig);
};

#endif // CLSCONFIG_H
