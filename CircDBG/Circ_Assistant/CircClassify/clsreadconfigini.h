#ifndef CLSREADCONFIGINI_H
#define CLSREADCONFIGINI_H
#include "string"
#include <map>
#include <vector>
using namespace std;

struct St_Section
{
    string strName;
    map<string, string> m_mpKeyValue;

    St_Section():strName("")
    {}

    void Init()
    {
        strName = "";
        m_mpKeyValue.clear();
    }
};

struct St_ConfigIni
{
    vector<St_Section> vSect;
};

class ClsReadConfigIni
{
public:
    ClsReadConfigIni();
    ~ClsReadConfigIni();

public:
    void ReadIni(const char* czIniPath);

    St_ConfigIni& GetConfigInfo()
    {
        return m_stIni;
    }

private:
    St_ConfigIni m_stIni;
};

#endif // CLSREADCONFIGINI_H
