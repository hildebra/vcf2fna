#include "fasta.h"


class VCFReader {
public:
	VCFReader(options* opts, refAssembly* R);
	~VCFReader() {}
private:
	//functions
	void read_vcf_file(std::istream* fp);
	void parse_single_vcf(string line);
	void parseFields(const string& s);
	bool splitXtra(const string& s);
	void parseINFO(const string& s);


	//variables
	refAssembly* refG;
	options* opts;
	string header;
	//string curSeq;
	fasta* cF;
	string vcfFile;
	string curChrom;
	int cntAvContigs;
	int lnCnt, snpCNT, indelCNT,snpFILT,indelFILT,unsrSNP,unsrINDEL;

	//parameters for filtering
	int minQual, minDep;
	float minFS, minMQ0F; //minimum values for filtering FS and MQ0F
	float minBQBZ, minSP; //minimum values for filtering BQBZ and SP

	robin_hood::unordered_map  <string, int> seenCtgs;

	//for vcf entry:
	string chrom, pos, id, ref, alt, qual, filter, info, format,xtra;
	//vcfstats
	int snpRepl, snpKept;
	bool fieldsSet;//set if GT:DP:SP:ADF:ADR:AD has been deparsed..
	int GT, DP, SP, ADF, ADR, AD, BQBZ, IDV, IMF, MQ0F, FS;
	int numFields;
	vector<string> fields;
	vector<int> DP4;//Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases
	float AF1val, AF2val, FSval, MQ0Fval, BQBZval, SPval, IDVval, IMFval,DPval, RPBZval, SCBZval, MQBZval;
};