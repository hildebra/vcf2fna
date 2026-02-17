#pragma once

#include "options.h"



using namespace std;


class VCFmem; // forward declaration

void ini_AA();


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
	gene() : geneID(""), geneStart(0), geneEnd(0), geneLength(0), geneStrand(true), numOnContig(-1){}
	gene(string id, int sta, int end);// : geneID(id), geneStart(sta), geneEnd(end), geneLength(end - sta), geneStrand(true), numOnContig(-1) {}
	gene(gene* GG);
	~gene() {}
	void setStrand(string s) { if (s == "-") geneStrand = false; else geneStrand = true; }
	void setType(string t) { type = t; }
	void setTranslationTable(int t) { translationTable = t; }
	void setNumOnContig(int n) { numOnContig = n; }
	void setPartial(string p) { partial = p; }

	string geneNT(const string& seq);
	string geneAA(const string& seq);
	int getIdx() { return numOnContig; }
	string createHDtag(const string& seq, list<int>& SNPsPos, list<float>& SNPfreqs,
		int& nonNs );

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
		INDELcnt(0),INDELpos(0),INDELfreq(0), geneCol(nullptr)
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
	void ntVariant(VCFmem* vx);
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

	// list of genes in the fasta file, from gff
	geneCollection* geneCol;

	friend class geneCollection; // Allows geneCollection to see vars

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
	void maskAllSeqs(const string&, int minDepth);



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
