#include <iostream>
#include "clsblast.h"
#include "clsreadconfigini.h"
#include "clsfastareader.h"
#include "cstdlib"  // This library is special for "atoi"
#include <map>
#include "clsclassificaiton.h"
#include <string.h>

using namespace std;

struct St_Config
{
    string strDBGResult;
    string strGTF;
    string strRef;
    string strCircCandi;
    string strBlastRootFolder;
    int iReadLen;
    int iKmerLen;
    bool bShowFlank;

    St_Config():strDBGResult(""), strGTF(""), strRef(""), strCircCandi(""), strBlastRootFolder(""),
                iReadLen(-1),iKmerLen(-1), bShowFlank(true)
    {}
};

int CheckFromIni(int argc, char** argv);
int CheckFromGroupRef(int argc, char** argv);

enum En_ReadsType{rtFasta=0, rtFastq, rtMax};

int main(int argc, char** argv)
{
    if(argc == 1 ||
       strcmp(argv[1], "-h") == 0 ||
       strcmp(argv[1], "--help") == 0)
    {
        cout << "**********************************" << endl;
        cout << "How to use CircAssistant" << endl;
        cout << "./CircAssistant ./config.ini" << endl;
        cout << "**********************************" << endl;
        return 0;
    }

    //Read Config File --->
    St_Config stConfig;
    ClsReadConfigIni* pIni = new ClsReadConfigIni();
    pIni->ReadIni(argv[1]); // The first valid parameters
    for(vector<St_Section>::iterator itr = pIni->GetConfigInfo().vSect.begin();
        itr != pIni->GetConfigInfo().vSect.end(); itr++)
    {
        if(itr->strName == "BlastChecking")
        {
            for(map<string, string>::iterator itrmp = itr->m_mpKeyValue.begin();
                itrmp != itr->m_mpKeyValue.end(); itrmp++)
            {
                if(itrmp->first == "DBGResult")
                    stConfig.strDBGResult = itrmp->second;
                else if(itrmp->first == "GTF")
                    stConfig.strGTF = itrmp->second;
                else if(itrmp->first == "Ref")
                    stConfig.strRef = itrmp->second;
                else if(itrmp->first == "CircCandi")
                    stConfig.strCircCandi = itrmp->second;
                else if(itrmp->first == "ReadLen")
                    stConfig.iReadLen = atoi(itrmp->second.c_str());
                else if(itrmp->first == "KmerLen")
                    stConfig.iKmerLen = atoi(itrmp->second.c_str());
                else if(itrmp->first == "BlastRootFolder")
                    stConfig.strBlastRootFolder = itrmp->second;
                else if(itrmp->first == "ShowFlank")
                {
                    int iType = atoi(itrmp->second.c_str());
                    if(iType == 0)
                        stConfig.bShowFlank = false;
                    else
                        stConfig.bShowFlank = true;
                }
            }
        }
    }
    delete pIni;
    pIni = NULL;
    //<---
    //CheckFromGroupRef(argc, argv);

    ClsClassificaiton *pClassify = new ClsClassificaiton();
    pClassify->Init(stConfig.strDBGResult, stConfig.strGTF, stConfig.strRef,
                    stConfig.strCircCandi, stConfig.strBlastRootFolder,
                    stConfig.iReadLen, stConfig.iKmerLen,
                    stConfig.bShowFlank);
    pClassify->ClassifyCirc();
    delete pClassify;
    pClassify = NULL;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////



