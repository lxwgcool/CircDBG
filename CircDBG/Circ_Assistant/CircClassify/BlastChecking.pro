TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lz

SOURCES += main.cpp \
    clsbasealgorithm.cpp \
    clsreadconfigini.cpp \
    clsmethod.cpp \
    clsbasefunction.cpp \
    clsfastareader.cpp

HEADERS += \
    clsbasealgorithm.h \
    clsreadconfigini.h \
    clsmethod.h \
    clsbasefunction.h \
    clsfastareader.h

