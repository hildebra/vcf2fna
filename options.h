#pragma once
//#include "IO.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
//#include <cstring>
#include <map>
//#include <list>
#include <stdlib.h>
//#include <algorithm>
#include <math.h>
//#include <cmath>
#include <time.h>
//#include <random>
#include <assert.h>
//#include <unordered_map>
//#include <numeric>
//#include <future>
//#include <mutex>
//#include <chrono>
//#include <random>
#include "include/robin_hood.h"
#include <algorithm>    
#include <sys/stat.h>



void helpMsg();



#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#pragma warning(disable:4996)
#else
#define _gzipread
#endif
#define notRpackage


#ifdef _gzipread
#include "gzstream.h"
#endif


using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;


//Version history
//0.1.0: 10.5.25: first version, to replace vcf2cons_mpi.pl in MG-TK
//0.2:11.5.25: first working version, complete with contig, vcf, gff and depth file read, vcf filtering, write contigs, genes and AA seqs out, 
const string vcf2fnaVERSION = "0.2";




bool isGZfile(std::string fi);
istream* openGZUZ(const string& inF);
ostream* writeGZUZ(const string& outF);



/*
inline bool file_exists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}
*/




struct options
{
public:
	options(int argc, char** argv);
	~options() {}

	//vars
	std::string refFasta = "";
	std::string outfna = "";
	string outGeneNT, outGeneAA;
	std::string tmp = "";
	std::string inVCF = "";
	std::string depthF = "";
	string gffFile = "";
	string outputTypes = "ANC";//each letter is one additional output
	uint threads = 1;
	int minDepthPar;
	int minCallQual;
	float minFS, minMQ0F, minBQBZ, minSP;//filtering of vcf
	bool addHDTags;//add info to contig header needed by MGTK
	bool skipEmptyContigs;//skip empty contigs in printing contig file



//functions within options
	void announce();
	void print_details();

};


//usage: 
// -depthF C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/m21BR347s3-smd.bam.coverage -ref C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/scaffolds.fasta.filt -gff C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/genes.gff -inVCF C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/test.vcf -t 1 -minCallDepth 1 -minCallQual 20 -oCtg C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.fna -oGeneNT C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.gene.fna -oGeneAA C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.gene.faa
