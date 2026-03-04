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
/*	bool filter(vector<string>& fields) {
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
	*/

	//variables
	options* opts;
	int numFields;
	int vcfFileNum;//in which of several vcf files are you??
	int GT, DP, SP, ADF, ADR, AD, BQBZ, IDV, IMF, MQ0F, FS;
	bool fieldsSet;//set if GT:DP:SP:ADF:ADR:AD has been deparsed..

};

class mutMatrix {
public:
	mutMatrix() { 
		mutmatrix=vector<vector<int>>(4, vector<int>(4,0)); 
	}
	void addMut(string ref, string alt) {
		if (ref == alt) { return; }
		if (ref == "." || ref == "*" || ref == "<" || alt == "." || alt == "*" || alt == "<") { return; }
		if (ref == "I" || ref == "D" || alt == "I" || alt == "D") { return; }
		if (ref == "N" || alt == "N") { return; }
		int i1(0), i2(0);
		if (ref == "C") { i1 = 1; }else if (ref == "T") { i1 = 2; }else if (ref == "G") { i1 = 3; }
		if (alt == "C") { i2 = 1; }else if (alt == "T") { i2 = 2; }else if (alt == "G") { i2 = 3; }
		mutmatrix[i1][i2]++;
	}
	int getTS() {
		//		//A:0, C:1, T:2, G:3
		//transitions are A<->G and C<->T, so 0<->3 and 1<->2
		int ts = mutmatrix[0][3] + mutmatrix[3][0] + mutmatrix[1][2] + mutmatrix[2][1];
		return ts;
	}
	int getTV() {
		//transversions are the other 8 possible changes
		int tv = mutmatrix[0][2] + mutmatrix[0][1] + 
			mutmatrix[1][0] + mutmatrix[1][3] + 
			mutmatrix[2][0] + mutmatrix[2][3] + 
			mutmatrix[3][1] + mutmatrix[3][2];
		return tv;
	}
	float getTsTv() {
		int ts = getTS(); int tv = getTV();
		//if (tv == 0) { return 1.f; }
		if (tv == 0 && ts == 0) { return -1.f; }
		return float(ts) / float(tv);
	}
private:
	vector<vector<int>> mutmatrix;
};

struct VCFfilterStats {
public:
	void printStats() {
		cout << "Filtering stats: ";
		cout << "QUAL: " << QUAL << ", ";
		cout << "Adapt_QUAL: " << QUALadaptive << ", ";
		cout << "indelProx: " << indelProx << ", ";
		cout << "MQ0F: " << MQ0Ffilt <<", ";
		cout << "FS: " << FSfilt <<", ";
		cout << "MQBZ: " << MQBZ <<", ";
		cout << "RPBZ: " << RPBZ <<", ";
		cout << "SCBZ: " << SCBZ <<", ";		
		cout << "BQBZ: " << BQBZ <<", ";
		cout << "SP: " << SP << ", ";
		cout << "DP: " << DP << endl;

		//mutmatrix stats
		cout << "TsTv ratio: Major Allele: " << muma->getTsTv() << " (" << muma->getTS() << "/ " << muma->getTV() << ")" << "; ";
		cout << " major filtered: " << mumaFilt->getTsTv() << " (" << mumaFilt->getTS() << "/" << mumaFilt->getTV() << ")" << "; ";
		cout << " minor Allele: " << mumaLowFreq->getTsTv() << " (" << mumaLowFreq->getTS() << "/" << mumaLowFreq->getTV() << ")" << "; ";
		cout << " minor filtered: " << mumaLowFreqFilt->getTsTv() << " (" << mumaLowFreqFilt->getTS() << "/" << mumaLowFreqFilt->getTV() << ")" << endl;
	}
	int MQ0Ffilt, FSfilt, MQBZ, RPBZ, SCBZ, QUAL, BQBZ, SP, DP, 
		QUALadaptive, indelProx;
	//capture the kind of mutation allowed through
	mutMatrix* muma;
	//mutation matrix for confirmed, filtered and unsure variants
	mutMatrix* mumaFilt;mutMatrix* mumaLowFreq; 
	mutMatrix* mumaLowFreqFilt;
	
	//, DP, SP, ADF, ADR, AD, BQBZ, IDV, IMF, MQ0F, FS;
	VCFfilterStats() :MQ0Ffilt(0), FSfilt(0), MQBZ(0), RPBZ(0), SCBZ(0), 
		QUAL(0), BQBZ(0),SP(0), DP(0), QUALadaptive(0), indelProx(0){
		muma = new mutMatrix();
		mumaFilt = new mutMatrix(); 
		mumaLowFreq = new mutMatrix();
		mumaLowFreqFilt = new mutMatrix();
	}
	~VCFfilterStats() {
		delete muma, mumaFilt, mumaLowFreq, mumaLowFreqFilt;
	}
	//, //DP(0), SP(0), ADF(0), ADR(0), AD(0), BQBZ(0), IDV(0), IMF(0), MQ0F(0), FS(0) {}
};


//vcf entry, but saved in memory
class VCFmem {
public:
	VCFmem(const string& line, vcfFields* VF, VCFfilterStats*);
	~VCFmem() {  }
	bool evalVCFentry(options* opts, float mD,int,int);
	bool evalVCFentry(options* opts,  VCFmem* v2, float mD,int,int);
	const string& getChrom() { return chrom; }
	//const string& getChrom() { return chrom; }
	bool isINDEL() { return isIndel; }
	bool isSNP() { return isSnp; }
	bool conflicted() { return conflict; }
	bool filtered() { return isFiltered; }
	bool unsure() { return isUnsure; }
	bool majorAllele() { return (getFreq() >= 0.51f); }
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
	VCFfilterStats* filterStats;
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

class VCFmulti {
public:
	VCFmulti():mems(0){}
	void evalVCFs(fasta* cF, VCFcollection* VCF2ptr, VariantStats*, VCFfilterStats* ,options*);
	void addVCFmem(VCFmem* v) { mems.push_back(v); }
private:
	list<VCFmem*> mems;
};




class VCFReader {
public:
	VCFReader(options* opts, refAssembly* R);
	~VCFReader() { delete Vstats; delete filterStats;}
private:
	//functions
	void read_vcf_file(std::istream* fp, int VCFnum, VCFcollection* Vcol);
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
	VCFfilterStats* filterStats;

};
