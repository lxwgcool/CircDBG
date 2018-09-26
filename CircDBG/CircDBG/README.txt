*************************
How to compile CircDBG
*************************
$ git clone https://github.com/lxwgcool/CircDBG.git
$ cd ./CircDBG/CircDBG/CircDBG/MakeFile
$ make clean
$ make

Notice: 
1: Since CircDBG depends on zlib, please install zlib first (or load the zlib module first if you compile it in HPC).
2: Since CircDBG uses some c++11 features, please use gcc 5.4.0 or the higher version to compile the code.

The executable file is named "CircRNADBG"

*********************
How to use CircDBG
*********************
1: The "config.ini" associated with CircDBG is located in "./CircDBG/CircDBG/CircDBG/config.ini"
2: Edit "config.ini" file to figure out the location of reads, reference, annotation file and Reads length. 
   The corresponding option name as below: 
   (1) Sequence Reads  : Reads1 and Reads2
   (2) Reference       : Reference
   (3) Annotation file : GTF
   (4) Reads Length    : ReadsLen

   Here is the example:

	[General]
	Reference=/scratch/hpc-xin/circRNA/customer/ref/hg19.fa
	GTF=/scratch/hpc-xin/circRNA/customer/gtf/gencode.v19.annotation.gtf
	Reads1=/scratch/hpc-xin/circRNA/customer/reads/H1input_A_R1.fq
	Reads2=
	[StdCircularRNAInfo]
	CircRNADb=
	CircRNAdbResult=
	Tissue=
	[Parameter]
	KmerLen=15
	ReadLen=101
	KmerRatio=60
	ThreadsNum=4
	DoDetection=1
	DoComparison=0
	MaxSupportNum=999
	[Comparison]
	MyResult=
	CiriResult=
	CIRCExplorerResult=
	FindCircResult=
	CircRNAFinderResult=
	CircMarkerResult=
	SimulatorBenchmark=     

3: How to run circMarker:
   ./CircRNADBG ./config.ini

4: Result: 
  (1) The results locate in the folder "./Detection_Result"
  (2) ./Candidate_*.txt: Recorded the circular splicing junciton of each potential circular RNA for each chromosome respectively. 
  (3) ./Brief_*.txt    : Breif format of "./Candidate_*.txt" for each chromosome respectively.
  (4) ./Brief_sum.txt  : The summary of /Brief_*.txt.

5: Notice: 
  (1) The chromosome should be listed in gtf file sequencially, which means that the order should be 1,2,3,4,5 and etc. The chromosome should start at the first one.  
  (2) CircDBG only support gtf format. 

