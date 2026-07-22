#pragma once

#include "options.h"
#include <cstdint>



using namespace std;


class VCFmem; // forward declaration
class fasta; // forward declaration
class VCFfilterStats; // forward declaration


struct VariantStats {
	int64_t snpCNT, indelCNT, snpFILT, indelFILT, unsrSNP,
		unsrINDEL, indelUsed, SNPused, conflictCnt,
		majorAllele, majorAlleleFilt, minorAllele, minorAlleleFilt,
		SNPlowCov, refMismatch;
	VariantStats() :snpCNT(0), indelCNT(0), snpFILT(0), indelFILT(0), unsrSNP(0), 
		unsrINDEL(0), indelUsed(0), SNPused(0), conflictCnt(0),
		majorAllele(0), majorAlleleFilt(0), minorAllele(0), minorAlleleFilt(0),
		SNPlowCov(0), refMismatch(0){}
	void printSNPstats() {
		cout << "SNP stats: " << endl;
		cout << "  - Found " << snpCNT << " SNPs and " << indelCNT << " INDELS." << endl;
		cout << "  - Filtered " << majorAlleleFilt << ", "<< minorAlleleFilt << "; " << indelFILT << " entries (major, minor SNP; INDEL)." << endl;
		cout << "  - Passing Filters: " << majorAllele << ", " << minorAllele << "; " << indelUsed << " entries (major, minor SNPs; INDELS). Conflicts resolved: " << conflictCnt << endl;
		cout << "  - Unsure: " << unsrSNP << ";" << unsrINDEL << " (SNPs; INDELS - replaced with N)" << endl;
		cout << "  - low Coverage: " << SNPlowCov << endl;
		cout << "  - Reference/coordinate mismatches skipped: " << refMismatch << endl;
	}
};

struct OutputStats {
	int64_t totalContigs, totalCtgNTs, totalGenes, totalGeneNTs, totalGeneAAs;
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
			geneLength(0), depthReferenceLength(0), geneStrand(true), translationTable(11), phase(0),
			numOnContig(-1), accumDepth(0){}
	gene(string id, int sta, int end);// : geneID(id), geneStart(sta), geneEnd(end), geneLength(end - sta), geneStrand(true), numOnContig(-1) {}
	gene(gene* GG);
	~gene() {}
	void setStrand(string s) { if (s == "-") geneStrand = false; else geneStrand = true; }
	void setType(string t) { type = t; }
	void setTranslationTable(int t) { translationTable = t; }
	void setPhase(int p) { phase = p; }
	void setSegments(const vector<pair<int, int> >& s);
	void setNumOnContig(int n) { numOnContig = n; }
	void setPartial(string p) { partial = p; }

	void addAccumDepth(int64_t d) { accumDepth += d; }
	float getAvgDepth() { return depthReferenceLength > 0 ?
		static_cast<float>(accumDepth) / static_cast<float>(depthReferenceLength) : 0.f; }

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
	int depthReferenceLength; // original annotated length represented by accumDepth
	vector<pair<int, int> > segments; // zero-based, inclusive genomic CDS intervals
	bool geneStrand;	// strand of the gene: + = true, - = false
	string type;	// type of the gene: CDS, gene, exon, intron, pseudogene, ncRNA_gene, tRNA, rRNA, miRNA, mRNA, etc.
	int translationTable;	// translation table for the gene
	int phase; // GFF3 phase: bases before the first complete codon in the oriented CDS
	int numOnContig;	// number of genes on the contig
	string partial;	// partial at left/right genomic boundaries (Prodigal partial=XY semantics)
	int64_t accumDepth; // accumulated depth in the gene, used for calculating average depth

	//functions
	void reverseComplement( string& seq);
	void recalculateGeometry();
	bool fivePrimePartial() const;
	int positionOffset(int genomicPos) const;

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
	void correctCoords(int pos, const string& ref, const string& alt);
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
	fasta(string s, string h);
	~fasta() { delete geneCol; }
//functions
	string getSeq() const { return seq; }
	string getHeader() const { return header; }
	string getSequenceId() const { return sequenceId; }
	int getLength() const { return (int)seq.length(); }
	void setSeq(string s) { seq = s; }
	bool validateVariantReference(VCFmem* vx, VariantStats* stats, const options* opts) const;
	void ntVariant(VCFmem* vx, VariantStats* Vstats, VCFfilterStats*,options*);
	string write(options* opts, OutputStats* Ostats);
	void addGene(string id, int sta, int end,string strand,string type,int transTab,
		string partial, int phase = 0, const vector<pair<int, int> >& segments = vector<pair<int, int> >());
	//void writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, bool doAA);
	int maskSeq(int start, int end, bool repl=false);
	int unmaskSeq(int start, int end);
	int64_t getSNPcount() { return SNPsCnt; }
	int64_t getNonNcount(options* opts);
	void writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, bool doAA, OutputStats* Ostats) {
		prepMuts();
		geneCol->writeAllGenes(opts, NTs, AAs, doNT, doAA, this, Ostats);
	}
	void prepMuts() { geneCol->prepMuts();  }
	geneCollection* getGeneCollection() { return geneCol; } 

	list<int>& getSNPsPos() { return SNPsPos; }
	list<float>& getSNPfreqs() { return SNPfreqs; }
	void addDepth(int sta, int sto, int depth, size_t source);
	bool isDepthResolved(int pos) const {
		return pos >= 0 && static_cast<size_t>(pos) < seqUse.size() &&
			seqUse[static_cast<size_t>(pos)];
	}
	float getAvgDepth() { return seq.empty() ? 0.f : static_cast<float>(depthAccum) / static_cast<float>(seq.length()); }
	float getAvgDepth(size_t source) const {
		return seq.empty() || source >= depthAccumBySource.size() ? 0.f :
			static_cast<float>(depthAccumBySource[source]) / static_cast<float>(seq.length());
	}


private:
	//functions
	void SNP(int pos, string r, string a, float freq);
	void INDEL(int pos, string r, string a, float freq);
	void maskUncertainSpan(int pos, size_t refLength);
	string getMutatedHeader(bool, const string& renderedSeq);
	bool referenceAlleleMatches(int pos, const string& ref, string& reason) const;
	string getMutatedSeq(options* opts);


	//variables
	string seq, mutSeq;
	bool mutSeqDone;
	string header, sequenceId;
	vector<bool> seqUse; // depth-resolution mask: true = emit the base, false = emit N
	//stats
	int64_t SNPsCnt, UnctCnt;
	list<int> SNPsPos;
	list<float> SNPfreqs;
	int64_t conflictCnt;
	int64_t INDELcnt;
	vector<int> INDELpos;
	list<float> INDELfreq;
	vector<string> INDELref, INDELalt;
	int64_t depthAccum; // accumulated depth in the sequence, used for calculating average depth
	vector<int64_t> depthAccumBySource; // per depth/VCF source; seqUse remains their union

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
	void processDepth(const string&, int minDepth, size_t source);



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

};
