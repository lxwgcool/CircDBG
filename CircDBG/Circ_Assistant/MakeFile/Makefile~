CODE=../CircClassify/

CC = g++
CFLAGS = -O3 -Wall -std=c++0x -static -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)
LIBS += -lz
OBJ = ./obj/

all: CircAssistant

clsreadconfigini.o : $(CODE)clsreadconfigini.cpp $(CODE)clsreadconfigini.h
	$(CC) $(CFLAGS) -c $(CODE)clsreadconfigini.cpp

clsgtfparse.o : $(CODE)clsgtfparse.cpp $(CODE)clsgtfparse.h
	$(CC) $(CFLAGS) -c $(CODE)clsgtfparse.cpp

clsfastareader.o : $(CODE)clsfastareader.cpp $(CODE)clsfastareader.h
	$(CC) $(CFLAGS) -c $(CODE)clsfastareader.cpp

clsclassificaiton.o : $(CODE)clsclassificaiton.cpp $(CODE)clsclassificaiton.h
	$(CC) $(CFLAGS) -c $(CODE)clsclassificaiton.cpp

clsblast.o : $(CODE)clsblast.cpp $(CODE)clsblast.h
	$(CC) $(CFLAGS) -c $(CODE)clsblast.cpp

clsbasealgorithm.o : $(CODE)clsbasealgorithm.cpp $(CODE)clsbasealgorithm.h
	$(CC) $(CFLAGS) -c $(CODE)clsbasealgorithm.cpp

main.o:  $(CODE)main.cpp 
	$(CC) $(CFLAGS) -c $(CODE)main.cpp


CircAssistant: clsreadconfigini.o clsgtfparse.o clsfastareader.o clsclassificaiton.o clsblast.o clsbasealgorithm.o main.o 
	$(CC) $(CFLAGS) -o CircAssistant clsreadconfigini.o clsgtfparse.o clsfastareader.o clsclassificaiton.o clsblast.o clsbasealgorithm.o main.o -lz


RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) CheckScafInfo


