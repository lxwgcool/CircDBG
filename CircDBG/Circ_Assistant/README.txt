****************************
How to compile CircAssistant
****************************
$ git clone https://github.com/lxwgcool/CircDBG.git
$ cd ./CircDBG/CircDBG/Circ_Assistant/MakeFile
$ make clean
$ make

Notice: 
1: Since CircAssistant depends on zlib, please install zlib first (or load the zlib module first if you compile it in HPC).
2: Since CircAssistant uses some c++11 features, please use gcc 5.4.0 or the higher version to compile the code.

The executable file is named "CircAssistant"

************************
How to use CircAssistant
************************
1: The "config.ini" associated with CircDBG is located in "./CircDBG/CircDBG/Circ_Assistant/config.ini"
2: Edit "config.ini" file to figure out the location of detail result of CircDBG, brief result of CircDBG, reference, annotation file, Blast folder and . 
   The corresponding option name as below: 
   (1) Detail result of CircDBG: DBGResult
   (2) Brief result of CircDBG : CircCandi
   (3) Reference       	       : Reference
   (4) Annotation file 	       : GTF
   (5) Reads Length    	       : ReadsLen

   Here is the example:

	[BlastChecking]
	DBGResult=/home/lq/lxwg/WorkStudio/Prototype/Circ_Assistant/Data/DBG_Result.txt
	GTF=/home/lq/lxwg/WorkStudio/Prototype/Circ_Assistant/Data/chr13.gtf
	Ref=/home/lq/lxwg/WorkStudio/Prototype/Circ_Assistant/Data/Homo_sapiens.GRCh37.75.dna.chromosome.13.fa
	CircCandi=/home/lq/lxwg/WorkStudio/Prototype/Circ_Assistant/Data/Brief.txt
	BlastRootFolder=/home/lq/lxwg/WorkStudio/Prototype/Circ_Assistant/ThirdPartyTools/blast
	ReadLen=101
	KmerLen=15
	ShowFlank=1  

3: How to run circMarker:
   ./CircAssistant ./config.ini

4: Result: 
  (1) The results locate in the folder "./Classify_Result"
  (2) ./AlignInfoSum.txt  : Alignment status of each circular RNA
  (3) ./Brief_sum_good.txt: Summary of good circRNA

5: Notice: 
  (1) The chromosome should be listed in gtf file sequencially, which means that the order should be 1,2,3,4,5 and etc. The chromosome should start at the first one.  
  (2) CircDBG only support gtf format. 

