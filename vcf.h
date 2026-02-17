#include "fasta.h"
#include <unordered_map>

class vcfFields {
public:
	vcfFields(options* o, int vv) :opts(o),
		numFields(0),vcfFileNum(vv),
		GT(-1), DP(-1), SP(-1), ADF(-1), ADR(-1), AD(-1), BQBZ(-1), IDV(-1), IMF(-1), MQ0F(-1), FS(-1),
		fieldsSet(false){}
	~vcfFields() {}
	void parseFields(const string& s);
	bool isSet() { return fieldsSet; }
	int getNumFields() { return numFields; }
	bool filter(vector<string>& fields) {
		bool isFiltered = false;	
		if (false && DP > -1) {
			int DPv = stoi(fields[DP]);
			if (SP != -1 && (DPv - stoi(fields[SP]) < opts->minDepthPar[vcfFileNum])) { isFiltered = true; }
			if (BQBZ != -1 && (stof(fields[BQBZ]) < opts->minBQBZ)) { isFiltered = true; }
			if (SP != -1 && (stof(fields[SP]) > opts->minSP)) { isFiltered = true; }
			//these are in info field..
			//if (MQ0F != -1 && (DPv - stoi(fields[MQ0F]) < minDep)) { return; }
			//if (FS != -1 && (stof(fields[FS]) < minFS)) { return; }
			//if (MQ0F != -1 && (stof(fields[MQ0F]) > minMQ0F)) { return; }
			//DPval = stoi(fields[DP]);
		}
		return isFiltered;
	}

	//variables
	options* opts;
	int numFields;
	int vcfFileNum;//in which of several vcf files are you??
	int GT, DP, SP, ADF, ADR, AD, BQBZ, IDV, IMF, MQ0F, FS;
	bool fieldsSet;//set if GT:DP:SP:ADF:ADR:AD has been deparsed..

};


//vcf entry, but saved in memory
class VCFmem {
public:
	VCFmem(const string& line, vcfFields* VF);
	~VCFmem() {}
	bool evalVCFentry(options* opts, vcfFields* VF);
	bool evalVCFentry(options* opts, vcfFields* VF, VCFmem* v2);
	const string& getChrom() { return chrom; }
	//const string& getChrom() { return chrom; }
	bool isINDEL() { return isIndel; }
	bool isSNP() { return isSnp; }
	bool conflicted() { return conflict; }
	bool filtered() { return isFiltered; }
	bool unsure() { return isUnsure; }
	bool majorAllele() { return (AF1val > 0.51f || altFreq > 0.51f); }
	float getFreq() { return altFreq; }
	int getPos() { return posN; }
	string getRef() { return ref; }
	string getAlt() { return alt; }
	//void getFields(vcfFields* VF){ VF->filter(fields); }
private:
	bool splitXtra(string& s, vcfFields* VF);
	void parseINFO(const string& s);



	string chrom, id, ref, alt, filter;// , info, format, xtra;

	bool fieldsSet;//set if GT:DP:SP:ADF:ADR:AD has been deparsed..
	float QUALval;
	int posN;// , GT, DP, SP, ADF, ADR, AD, BQBZ, IDV, IMF, MQ0F, FS;
	vector<string> fields;
	vector<int> DP4;//Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases
	float altFreq, AF1val, AF2val, FSval, MQ0Fval, BQBZval, SPval, IDVval, IMFval, DPval, RPBZval, SCBZval, MQBZval;
	bool isIndel, isSnp;// , isMNP, isComplex, isMixed, isBND, isSNPfiltered, isIndelfiltered;
	bool isFiltered, isUnsure;
	bool conflict;
};

typedef unordered_map<int, VCFmem*> map2vcfm;

class VCFcollection {
public:
	VCFcollection():curContig(""), curMems(nullptr), VMems(0) {}
	~VCFcollection();

	VCFmem* recoverVCF(const string& cont, int pos);
	void addVCF(VCFmem*);
	size_t size();

private:
	string curContig;
	map2vcfm* curMems;
	unordered_map <string, map2vcfm* > VMems;
};


struct VariantStats {
	int snpCNT, indelCNT, snpFILT, indelFILT, unsrSNP, unsrINDEL, indelUsed, SNPused, conflictCnt;
	VariantStats() :snpCNT(0), indelCNT(0), snpFILT(0), indelFILT(0), unsrSNP(0), unsrINDEL(0), indelUsed(0), SNPused(0), conflictCnt(0) {}
};

class VCFReader {
public:
	VCFReader(options* opts, refAssembly* R);
	~VCFReader() { delete Vstats; }
private:
	//functions
	void read_vcf_file(std::istream* fp, int VCFnum, VCFcollection* Vcol);
	void collectVCFstats(VCFmem* vcf);
	//void parse_single_vcf(string line);
	//bool evalVCFentry();


	//variables
	refAssembly* refG;
	options* opts;
	string header;
	//string curSeq;
	fasta* cF;
	vector<string> vcfFile;
	string curChrom;
	int cntAvContigs;
	int lnCnt;
	//int snpCNT, indelCNT, snpFILT, indelFILT, unsrSNP, unsrINDEL, indelUsed, SNPused, conflictCnt;

	//parameters for filtering
	int minQual, minDep;
	float minFS, minMQ0F; //minimum values for filtering FS and MQ0F
	float minBQBZ, minSP; //minimum values for filtering BQBZ and SP

	robin_hood::unordered_map  <string, int> seenCtgs;

	//for vcf entry:
	//string chrom, id, ref, alt, filter, info, format,xtra;
	//vcfstats
	int snpRepl, snpKept;

	VCFcollection* VCF2ptr;
	VariantStats* Vstats;
};
