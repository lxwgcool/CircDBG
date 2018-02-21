#include "clsconfig.h"
#include "../../ShareLibrary/clsreadconfigini.h"
#include "stdlib.h"
#include <iostream>
using namespace std;

ClsConfig::ClsConfig()
{
}

void ClsConfig::ReadConfig(St_Config& stConfig, char* cpIniPath)
{
    ClsReadConfigIni* pIni = new ClsReadConfigIni();
    pIni->ReadIni(cpIniPath); // The first valid parameters

    for(vector<St_Section>::iterator itr = pIni->GetConfigInfo().vSect.begin();
        itr != pIni->GetConfigInfo().vSect.end(); itr++)
    {
        if(itr->strName == "General")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "Reference")
                    stConfig.strRefPath = itrmp->second;
                else if(itrmp->first == "GTF")
                    stConfig.strGtfPath = itrmp->second;
                else if(itrmp->first == "Reads1")
                    stConfig.strReads1Path = itrmp->second;
                else if(itrmp->first == "Reads2")
                    stConfig.strReads2Path = itrmp->second;
            }
        }
        else if(itr->strName == "StdCircularRNAInfo")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "CircRNADb")
                    stConfig.strCircRNADb = itrmp->second;
                else if(itrmp->first == "CircRNAdbResult")
                    stConfig.strCircRNAdbResult = itrmp->second.c_str();
                else if(itrmp->first == "Tissue")
                    stConfig.strTissue = itrmp->second.c_str();
                else{}
            }
        }
        else if(itr->strName == "Parameter")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "KmerLen")
                    stConfig.iKmerLen = atoi(itrmp->second.c_str());
                else if(itrmp->first == "ReadLen")
                    stConfig.iReadLen = atoi(itrmp->second.c_str());
                else if(itrmp->first == "KmerRatio")
                    stConfig.fKmerRatio = (float)atoi(itrmp->second.c_str()) / 100;
                else if(itrmp->first == "ThreadsNum")
                    stConfig.iThreadsNum = atoi(itrmp->second.c_str());
                else if(itrmp->first == "DoDetection")
                    stConfig.bDoDetection = (atoi(itrmp->second.c_str()) != 0 ? true : false);
                else if(itrmp->first == "DoComparison")
                    stConfig.bDoComparison = (atoi(itrmp->second.c_str()) != 0 ? true : false);
                else if(itrmp->first == "MaxSupportNum")
                    stConfig.iMaxSupportNum = atoi(itrmp->second.c_str());
                else
                {}
            }
        }
        else if(itr->strName == "Comparison")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "MyResult")
                    stConfig.strMyResult = itrmp->second.c_str();                
                else if(itrmp->first == "CiriResult")
                    stConfig.strCiriResult = itrmp->second.c_str();
                else if(itrmp->first == "CIRCExplorerResult")
                    stConfig.strCIRCExplorerResult = itrmp->second.c_str();
                else if(itrmp->first == "FindCircResult")
                    stConfig.strFindCircResult = itrmp->second.c_str();
                else if(itrmp->first == "CircRNAFinderResult")
                    stConfig.strCircRNAFinderResult = itrmp->second.c_str();
                else if(itrmp->first == "CircMarkerResult")
                    stConfig.strCircMarkerResult = itrmp->second.c_str();
                else if(itrmp->first == "SimulatorBenchmark")
                    stConfig.strSimulatorBenchmark = itrmp->second.c_str();
                else
                {}
            }
        }
        else{}
    }

    delete pIni;
    pIni = NULL;
}

bool ClsConfig::CheckConfig(St_Config& stConfig)
{
    return true;
}
