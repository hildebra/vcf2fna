#pragma once

#include "options.h"



using namespace std;


class VCFmem; // forward declaration
class fasta; // forward declaration
class VCFfilterStats; // forward declaration
void ini_AA();


struct VariantStats {
	int snpCNT, indelCNT, snpFILT, indelFILT, unsrSNP, 
		unsrINDEL, indelUsed, SNPused, conflictCnt,
		majorAllele, majorAlleleFilt, minorAllele, minorAlleleFilt;
	VariantStats() :snpCNT(0), indelCNT(0), snpFILT(0), indelFILT(0), unsrSNP(0), 
		unsrINDEL(0), indelUsed(0), SNPused(0), conflictCnt(0),
		majorAllele(0), majorAlleleFilt(0), minorAllele(0), minorAlleleFilt(0){}
	void printSNPstats() {
		cout << "SNP stats: " << endl;
		cout << "  - Found " << snpCNT << " SNPs and " << indelCNT << " INDELS." << endl;
		cout << "  - Filtered " << majorAlleleFilt << ", "<< minorAlleleFilt << "; " << indelFILT << " entries (major, minor SNP; INDEL)." << endl;
		cout << "  - Passing Filters: " << majorAllele << ", " << minorAllele << "; " << indelUsed << " entries (major, minor SNPs; INDELS). Conflicts resolved: " << conflictCnt << endl;
		cout << "  - Unsure: " << unsrSNP << ";" << unsrINDEL << " (SNPs; INDELS - replaced with N)" << endl;
	}
};

struct OutputStats {
	int totalContigs, totalCtgNTs, totalGenes, totalGeneNTs, totalGeneAAs;
	OutputStats() :totalContigs(0), totalCtgNTs(0), totalGenes(0), totalGeneNTs(0), totalGeneAAs(0) {}
	void printStats() {
		cout << "In total, wrote " << totalContigs << " contigs with "<< totalCtgNTs<< " valid NTs, ";
		cout << totalGenes << " genes, ";
		cout << totalGeneNTs << " gene NTs, ";
		cout << totalGeneAAs << " gene AAs." << endl;
	}
};

class gene
{
public:
	gene() : geneID(""), geneStart(0), geneEnd(0), 
			geneLength(0), geneStrand(true), numOnContig(-1), accumDepth(0){}
	gene(string id, int sta, int end);// : geneID(id), geneStart(sta), geneEnd(end), geneLength(end - sta), geneStrand(true), numOnContig(-1) {}
	gene(gene* GG);
	~gene() {}
	void setStrand(string s) { if (s == "-") geneStrand = false; else geneStrand = true; }
	void setType(string t) { type = t; }
	void setTranslationTable(int t) { translationTable = t; }
	void setNumOnContig(int n) { numOnContig = n; }
	void setPartial(string p) { partial = p; }

	void addAccumDepth(int d) { accumDepth += d; }
	float getAvgDepth() { return (float)accumDepth / geneLength; }

	string geneNT(const string& seq);
	string geneAA(const string& seq);
	int getIdx() { return numOnContig; }
	string createHDtag(const string& seq, fasta* fa, int& nonNs );

private:
	//variables
	string geneID;	// gene name
	int geneStart;	// start position of the gene
	int geneEnd;	// end position of the gene
	int geneLength;	// length of the gene
	bool geneStrand;	// strand of the gene: + = true, - = false
	string type;	// type of the gene: CDS, gene, exon, intron, pseudogene, ncRNA_gene, tRNA, rRNA, miRNA, mRNA, etc.
	int translationTable;	// translation table for the gene
	int numOnContig;	// number of genes on the contig
	string partial;	// partial gene: 00:no, 01:5', 10:3', 11:both
	int accumDepth; //accumulated depth in the gene, used for calculating average depth

	//functions
	void reverseComplement( string& seq);

	friend class geneCollection;
};


class fasta;//fwd declaration

class geneCollection
{
public:
	geneCollection():genes(0, nullptr), genesMut(0, nullptr), mutsPrepared(false) {}
	~geneCollection();// { for (size_t i = 0; i < genes.size(); i++) { if (genes[i] != nullptr) { delete genes[i]; } } }
	void addGene(gene* g) { genes.push_back(g); }
	void writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, 
		bool doAA, fasta* fa, OutputStats* Ostats);
	size_t size() { return genes.size(); }
	void push_back(gene* g) { genes.push_back(g); }
	void correctCoords(int pos, int altL, int refL);
	void prepMuts();
	void depthInGenes(int sta, int sto, int depth);
private:
	vector<gene*> genes;
	vector<gene*> genesMut;
	bool mutsPrepared;

};



class fasta
{
public:
	fasta(string s, string h) :seq(s), mutSeq(""),mutSeqDone(false),
		header(h), seqUse(s.length(), false), length(s.length()),
		SNPsCnt(0), UnctCnt(0),  SNPsPos(0), SNPfreqs(0),conflictCnt(0),
		INDELcnt(0),INDELpos(0),INDELfreq(0), depthAccum(0), geneCol(nullptr)
	{
		geneCol = new geneCollection();
	}
	~fasta() { delete geneCol; }
//functions
	string getSeq() const { return seq; }
	string getHeader() const { return header; }
	int getLength() const { return (int)seq.length(); }
	void setSeq(string s) { seq = s; }
	void resetCnts() { length = seq.length(); }
	void ntVariant(VCFmem* vx, VariantStats* Vstats, VCFfilterStats*,options*);
	string write(options* opts, OutputStats* Ostats);
	void addGene(string id, int sta, int end,string strand,string type,int transTab,string partial);
	//void writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, bool doAA);
	int maskSeq(int start, int end, bool repl=false);
	int unmaskSeq(int start, int end);
	long getSNPcount() { return SNPsCnt; }
	long getNonNcount() { return length - getNumNs(); }
	void writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, bool doAA, OutputStats* Ostats) {
		geneCol->writeAllGenes(opts, NTs, AAs, doNT, doAA, this, Ostats);
	}
	void prepMuts() { geneCol->prepMuts();  }
	geneCollection* getGeneCollection() { return geneCol; } 

	list<int>& getSNPsPos() { return SNPsPos; }
	list<float>& getSNPfreqs() { return SNPfreqs; }
	void addDepth(int sta, int sto, int depth) { depthAccum += (long)(sto- sta) * depth; }
	float getAvgDepth() { return (float)depthAccum / float(seq.length()); }


private:
	//functions
	void SNP(int pos, string r, string a, float freq);
	void INDEL(int pos, string r, string a, float freq);
	string getMutatedHeader(bool);
	//int getNumNs() { return count(seq.begin(), seq.end(), 'N'); }
	int getNumNs() { return count(seqUse.begin(), seqUse.end(), false); }
	string getMutatedSeq(options* opts);


	//variables
	string seq, mutSeq;
	bool mutSeqDone;
	string header;
	vector<bool> seqUse; //use this position? true = use, false = mask with N
	int length;
	//stats
	int SNPsCnt, UnctCnt;
	list<int> SNPsPos;
	list<float> SNPfreqs;
	int conflictCnt;
	int INDELcnt;
	vector<int> INDELpos;
	list<float> INDELfreq;
	vector<string> INDELref, INDELalt;
	long depthAccum; //accumulated depth in the sequence, used for calculating average depth

	// list of genes in the fasta file, from gff
	geneCollection* geneCol;

	friend class geneCollection; // Allows geneCollection to see vars
	//friend class gene; // Allows geneCollection to see vars

};



class refAssembly
{
public:
	refAssembly(options* opt);
	~refAssembly(void);
	fasta* getFasta(string id);
	int getNSeqs() { return NSeqs; }
	void setFasta(string id, fasta*);
	bool isSequence(string id);
	//string getHeader() const { return headers[0]; }
	void readDepth();
	void writeOutputs();
	void readGFF();

private:
//functions
	void writeContigs(OutputStats* Ostats);
	void writeGenes(OutputStats* Ostats);
	void processDepth(const string&, int minDepth);



	int NSeqs;	// number of sequences in the fasta file
	string refFasta;	// reference fasta file name

	//vector<int> seqLength;	// length of the sequence
	//vector<int> seqNonNs;	// ID of the sequence
	//vector<string> sequences;	// sequence of the fasta file
	vector<fasta*> fastas;
	//vector<string> headers;
	//void readFastaFile(FILE* fp);
	robin_hood::unordered_map  <string, int> hd2ID;
	options* opts;


	//functions
	int replaceWithNs(std::string& seq, const int start, const int end, char replaceWith = 'N');

};
