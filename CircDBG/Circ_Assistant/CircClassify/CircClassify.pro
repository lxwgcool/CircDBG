TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lz

SOURCES += main.cpp \
    clsreadconfigini.cpp \
    clsfastareader.cpp \
    clsbasealgorithm.cpp \
    clsblast.cpp \
    clsclassificaiton.cpp \
    clsgtfparse.cpp

HEADERS += \
    clsreadconfigini.h \
    clsfastareader.h \
    clsbasealgorithm.h \
    clsblast.h \
    clsclassificaiton.h \
    clsgtfparse.h

