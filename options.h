#pragma once
//#include "IO.h"
#include <stdio.h>
#include <cstring>
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
#include <algorithm>    
#include <sys/stat.h>
#include <list>
#include <memory>
#include "include/robin_hood.h"



void helpMsg();



#if defined(WIN32) || defined(_WIN32) || (defined(__WIN32) && !defined(__CYGWIN__))
#pragma warning(disable:4996)
#endif
#define notRpackage

// gzstream is an optional, externally supplied dependency.  Defining
// VCF2FNA_USE_GZSTREAM enables transparent .gz input/output; builds without
// that dependency remain portable and report a clear error for .gz paths.
#ifdef VCF2FNA_USE_GZSTREAM
#include "gzstream.h"
#define VCF2FNA_HAS_GZIP 1
#endif


using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;


//Version history
//0.1.0: 10.5.25: first version, to replace vcf2cons_mpi.pl in MG-TK
//0.2:11.5.25: first working version, complete with contig, vcf, gff and depth file read, vcf filtering, write contigs, genes and AA seqs out, 
//0.21: 18.12.25:  more reports on depth statistics, possibility to not output any files added (or only specific subfiles)
//0.22: 12.1.26: INDELS now also handled in gene sequences
//0.23: 5.2.26: fix indel options that were too strict
//0.24: 16.2.26: bug removed that removed contigs with 0 non-Ns from gene NT & AA files. Also added stats on genes & NTs & AA written 
//0.25: 17.2.26: added depth stat per gene, more precise reporting on actually implemented SNPs
//0.26: 20.2.26: further bugs addressed
//0.27: 26.2.26: added filter stats reports for SNPs
//0.28: 27.2.26: added ts/tv and mutMatrix
//0.29: dynamic depth-dependent filtering, indel prox
//0.30: 4.3.26: made options accessible
//0.31: 6.3.26: high freq SNP acceptance, masking only for depth
//0.32: 6.3.26: multi-allelic integration
//0.40: 21.7.26: correctness, VCF/GFF validation, bacterial translation, and CLI hardening
const string vcf2fnaVERSION = "0.40";




bool isGZfile(std::string fi);
istream* openGZUZ(const string& inF);
ostream* writeGZUZ(const string& outF);



/*
inline bool file_exists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}
*/

vector<string> splitByComma(const string& str);
vector<int> splitByCommaI(const string& str);
string combCommas(const vector<string>& vec);
string combCommas(const vector<int>& vec);


struct options
{
public:
	options(int argc, char** argv);
	~options() {}

	//functions within options
	void announce();
	void print_details();

	//vars
	std::string refFasta = "";
	std::string outfna = "";
	string outGeneNT, outGeneAA;
	std::string tmp = "";
	vector<string> inVCF;
	vector<string> depthF;
	string gffFile = "";
	string outputTypes = "";//each letter is one additional output
	vector<string> seqPlatform;//informational labels: ill, PB, ONT, or unspecified
	string vcfFilterPolicy = "technical";//technical, all, ignore
	uint threads = 1;
	uint indelRange = 5;
	int maxRefMismatches = 10;
	vector<int> minDepthPar;
	int minAltReads = 4;
	float minAltFreq = 0.05f;
	int minCallQual;
	float minMQ0F, minBQBZ, minSP;//filtering of vcf
	bool addHDTags;//add info to contig header needed by MGTK
	bool skipEmptyContigs;//skip empty contigs in printing contig file
	bool skipEmptyGenes;//skip empty genes in printing gene file
	bool debug1;//prints vcf lines..
	bool reportINDELs;
	bool maskMinorAllele;
	float minCallQualAdaptive;//if >0, adds a mean-depth-scaled QUAL threshold
	float depthFilterScale; // if DP < mean contig depth *x, filter. Default: 0.25
	float maxDepthFilterScale; // if DP > mean contig depth *x, filter. 0 disables. Default: 3

};


//usage example 1: 
// -depthF C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/m21BR347s3-smd.bam.coverage -ref C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/scaffolds.fasta.filt -gff C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/genes.gff -inVCF C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/test.vcf -t 1 -minCallDepth 1 -minCallQual 20 -oCtg C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.fna -oGeneNT C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.gene.fna -oGeneAA C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.gene.faa
//usage example 2:
// -seqPlatform ill,PB -minCallDepth 2,1 -depthF C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/m21BR347s3-smd.bam.coverage,C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/m21BR347s3-smd.bam.coverage -ref C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/scaffolds.fasta.filt -gff C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/genes.gff -inVCF C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/test.vcf,C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/test.vcf -t 1 -minCallQual 20 -oCtg C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.fna -oGeneNT C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.gene.fna -oGeneAA C:/Users/hildebra/OneDrive/science/data/test/vcf2fasta_muco/cons.new.gene.faa


//trial data:
// /ei/projects/1/115b210e-aa45-4469-9eb1-d0d85d879fc0/data/mucosodom_runs/Mucosodom_Run5/m22BR624s6
