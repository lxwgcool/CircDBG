TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lz
LIBS += -pthread
CONFIG += c++11

SOURCES += main.cpp \
    ../../ShareLibrary/clsbasealgorithm.cpp \
    ../../ShareLibrary/clsbwa.cpp \
    ../../ShareLibrary/clsfastareader.cpp \
    ../../ShareLibrary/clsfastqreader.cpp \
    ../../ShareLibrary/clsgtfparse.cpp \
    ../../ShareLibrary/clsjellyfish.cpp \
    ../../ShareLibrary/clskmeralgorithm.cpp \
    ../../ShareLibrary/clsreadconfigini.cpp \
    ../../ShareLibrary/clsvelvet.cpp \
    ../../ShareLibrary/FastqFileParse/kseq_test.c \
    clsparsestdcircrnainfo.cpp \
    clsconfig.cpp \
    clsreadsmapping.cpp \
    clspreevaluation.cpp \
    clsdebruijngraph.cpp \
    clscircrnadetection.cpp \
    clsresultcomparison.cpp \
    clscomparebasefunction.cpp

HEADERS += \
    ../../ShareLibrary/clsbasealgorithm.h \
    ../../ShareLibrary/clsbwa.h \
    ../../ShareLibrary/clsfastareader.h \
    ../../ShareLibrary/clsfastqreader.h \
    ../../ShareLibrary/clsgtfparse.h \
    ../../ShareLibrary/clsjellyfish.h \
    ../../ShareLibrary/clskmeralgorithm.h \
    ../../ShareLibrary/clsreadconfigini.h \
    ../../ShareLibrary/clsvelvet.h \
    ../../ShareLibrary/FastqFileParse/kseq.h \
    clsparsestdcircrnainfo.h \
    clsconfig.h \
    clsreadsmapping.h \
    clspreevaluation.h \
    clsdebruijngraph.h \
    clscircrnadetection.h \
    clsresultcomparison.h \
    clscomparebasefunction.h


unix:!macx: LIBS += -L$$PWD/../../ShareLibrary/bamtools/lib/ -lbamtools

INCLUDEPATH += $$PWD/../../ShareLibrary/bamtools/include
DEPENDPATH += $$PWD/../../ShareLibrary/bamtools/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../../ShareLibrary/bamtools/lib/libbamtools.a

OTHER_FILES += \
    ../Document/Readme.txt \
    ../Document/DBG_Implementation_Logic.txt \
    ../Document/Debug_Schedule.txt \
    ../Document/Data_Result.txt \
    ../Document/Data_Finish_Status.txt
