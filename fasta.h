#pragma once

#include "options.h"
#include <list>
#include <algorithm>



using namespace std;


void ini_AA();



class gene
{
public:
	gene() : geneID(""), geneStart(0), geneEnd(0), geneLength(0), geneStrand(true), numOnContig(-1){}
	gene(string id, int sta, int end) : geneID(id), geneStart(sta), geneEnd(end), geneLength(end-sta), geneStrand(true), numOnContig(-1) {}
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
};


class fasta
{
public:
	fasta(string s, string h) :seq(s),header(h),length(s.length()),
		SNPsCnt(0), UnctCnt(0), SNPsPos(0), SNPfreqs(0),conflictCnt(0){}
	~fasta() { for (size_t i = 0; i < genes.size(); i++) { if (genes[i] != nullptr) { delete genes[i]; } } }
	string getSeq() const { return seq; }
	string getHeader() const { return header; }
	int getLength() const { return (int)seq.length(); }
	void setSeq(string s) { seq = s; }
	void resetCnts() { length = seq.length(); }
	void SNP(int pos, string r, string a, float freq);
	string write(bool hdTags,bool skipEmpty);
	void addGene(string id, int sta, int end,string strand,string type,int transTab,string partial);
	void writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, bool doAA);
private:
	//functions
	string getMutatedSeq() { return seq; }  //TODO
	string getMutatedHeader(bool);
	int getNumNs() { return count(seq.begin(), seq.end(), 'N'); }


	//variables
	string seq;
	string header;
	int length;
	//stats
	int SNPsCnt, UnctCnt;
	list<int> SNPsPos;
	list<float> SNPfreqs;
	int conflictCnt;

	vector<gene*> genes;	// list of genes in the fasta file
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
	void writeContigs();
	void writeGenes();



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
