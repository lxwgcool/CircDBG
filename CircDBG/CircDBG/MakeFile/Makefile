SharedLibraryPath=../../ShareLibrary/

#BAMTOOLS=$(SharedLibraryPath)bamtools/include
#BAMTOOLS_LD=$(SharedLibraryPath)bamtools/lib

BAMTOOLS=/apps/bamtools/2.3.0/include
BAMTOOLS_LD=/apps/bamtools/2.3.0/lib

#RootPath = ../ScaffoldGapFillerTest/
#NeedleManPath = ../Depends/needleman_wunsch/
#NormDependsPath = ../Depends/
#INC=-I$(SharedLibraryPath)
#INC_PARAMS=$(foreach d, $(INC), -I$d)
#Notice: "-c" will cause the error

CC = g++
INC_DIR = $(SharedLibraryPath) #Identify include path
#CFLAGS = -O3 -Wall -std=c++0x -static -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)
CFLAGS = -O3 -Wall -static -std=c++0x
LIBS += -lz
LIBS += -lpthread
#INCLUDEPATH += tr1/memory

OBJ = ./obj/

all: CircRNADBG

clsbasealgorithm.o : $(SharedLibraryPath)clsbasealgorithm.cpp $(SharedLibraryPath)clsbasealgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbasealgorithm.cpp

clsfastareader.o: $(SharedLibraryPath)clsfastareader.cpp $(SharedLibraryPath)clsfastareader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastareader.cpp

clsfastqreader.o: $(SharedLibraryPath)clsfastqreader.cpp $(SharedLibraryPath)clsfastqreader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastqreader.cpp

clsbwa.o : $(SharedLibraryPath)clsbwa.cpp $(SharedLibraryPath)clsbwa.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbwa.cpp

clsgtfparse.o: $(SharedLibraryPath)clsgtfparse.cpp $(SharedLibraryPath)clsgtfparse.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsgtfparse.cpp

clsreadconfigini.o: $(SharedLibraryPath)clsreadconfigini.cpp $(SharedLibraryPath)clsreadconfigini.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsreadconfigini.cpp

clskmeralgorithm.o: $(SharedLibraryPath)clskmeralgorithm.cpp $(SharedLibraryPath)clskmeralgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clskmeralgorithm.cpp

clsconfig.o: ../CircRNADBG/clsconfig.cpp ../CircRNADBG/clsconfig.h
	$(CC) $(CFLAGS) -c ../CircRNADBG/clsconfig.cpp

clsparsestdcircrnainfo.o: ../CircRNADBG/clsparsestdcircrnainfo.cpp ../CircRNADBG/clsparsestdcircrnainfo.h
	$(CC) $(CFLAGS) -c ../CircRNADBG/clsparsestdcircrnainfo.cpp

#clsreadsmapping.o: ../CircRNADBG/clsreadsmapping.cpp ../CircRNADBG/clsreadsmapping.h
#	$(CC) $(CFLAGS) -c ../CircRNADBG/clsreadsmapping.cpp

clscircrnadetection.o: ../CircRNADBG/clscircrnadetection.cpp ../CircRNADBG/clscircrnadetection.h
	$(CC) $(CFLAGS) -c ../CircRNADBG/clscircrnadetection.cpp

clsdebruijngraph.o: ../CircRNADBG/clsdebruijngraph.cpp ../CircRNADBG/clsdebruijngraph.h
	$(CC) $(CFLAGS) -c ../CircRNADBG/clsdebruijngraph.cpp

#clspreevaluation.o: ../CircRNADBG/clspreevaluation.cpp ../CircRNADBG/clspreevaluation.h
#	$(CC) $(CFLAGS) -c ../CircRNADBG/clspreevaluation.cpp

clscomparebasefunction.o: ../CircRNADBG/clscomparebasefunction.cpp ../CircRNADBG/clscomparebasefunction.h
	$(CC) $(CFLAGS) -c ../CircRNADBG/clscomparebasefunction.cpp

clsresultcomparison.o: ../CircRNADBG/clsresultcomparison.cpp ../CircRNADBG/clsresultcomparison.h
	$(CC) $(CFLAGS) -c ../CircRNADBG/clsresultcomparison.cpp

main.o: ../CircRNADBG/main.cpp 
	$(CC) $(CFLAGS) -c ../CircRNADBG/main.cpp

#CircRNADBG: clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsparsestdcircrnainfo.o clsreadsmapping.o clscircrnadetection.o clsdebruijngraph.o clspreevaluation.o clscomparebasefunction.o clsresultcomparison.o main.o 
#	$(CC) $(CFLAGS) -o CircRNADBG clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsparsestdcircrnainfo.o clsreadsmapping.o clscircrnadetection.o clsdebruijngraph.o clspreevaluation.o clscomparebasefunction.o clsresultcomparison.o main.o -lz -lbamtools -lpthread

CircRNADBG: clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsparsestdcircrnainfo.o clscircrnadetection.o clsdebruijngraph.o clscomparebasefunction.o clsresultcomparison.o main.o 
	$(CC) $(CFLAGS) -o CircRNADBG clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsgtfparse.o clsreadconfigini.o clskmeralgorithm.o clsconfig.o clsparsestdcircrnainfo.o clscircrnadetection.o clsdebruijngraph.o clscomparebasefunction.o clsresultcomparison.o main.o -lz -lpthread

RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) CheckScafInfo


