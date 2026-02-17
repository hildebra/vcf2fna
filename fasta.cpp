#include "fasta.h"
#include "vcf.h"


//****************************************************************


char DNA_AA[256][256][256];

void ini_AA() {
	DNA_AA['T']['C']['A'] = 'S'; DNA_AA['T']['C']['C'] = 'S'; DNA_AA['T']['C']['G'] = 'S'; DNA_AA['T']['C']['T'] = 'S'; // Serine
	DNA_AA['T']['T']['C'] = 'F'; DNA_AA['T']['T']['T'] = 'F'; // Phenylalanine
	DNA_AA['T']['T']['A'] = 'L'; DNA_AA['T']['T']['G'] = 'L'; // Leucine
	DNA_AA['T']['A']['C'] = 'Y'; DNA_AA['T']['A']['T'] = 'Y'; // Tyrosine
	DNA_AA['T']['A']['A'] = '*'; DNA_AA['T']['A']['G'] = '*'; DNA_AA['T']['G']['A'] = '*'; // Stop
	DNA_AA['T']['G']['C'] = 'C'; DNA_AA['T']['G']['T'] = 'C'; // Cysteine
	DNA_AA['T']['G']['G'] = 'W'; // Tryptophan
	DNA_AA['C']['T']['A'] = 'L'; DNA_AA['C']['T']['C'] = 'L'; DNA_AA['C']['T']['G'] = 'L'; DNA_AA['C']['T']['T'] = 'L'; // Leucine
	DNA_AA['C']['C']['A'] = 'P'; DNA_AA['C']['C']['C'] = 'P'; DNA_AA['C']['C']['G'] = 'P'; DNA_AA['C']['C']['T'] = 'P'; // Proline
	DNA_AA['C']['A']['C'] = 'H'; DNA_AA['C']['A']['T'] = 'H'; // Histidine
	DNA_AA['C']['A']['A'] = 'Q'; DNA_AA['C']['A']['G'] = 'Q'; // Glutamine
	DNA_AA['C']['G']['A'] = 'R'; DNA_AA['C']['G']['C'] = 'R'; DNA_AA['C']['G']['G'] = 'R'; DNA_AA['C']['G']['T'] = 'R'; // Arginine
	DNA_AA['A']['T']['A'] = 'I'; DNA_AA['A']['T']['C'] = 'I'; DNA_AA['A']['T']['T'] = 'I'; // Isoleucine
	DNA_AA['A']['T']['G'] = 'M'; // Methionine
	DNA_AA['A']['C']['A'] = 'T'; DNA_AA['A']['C']['C'] = 'T'; DNA_AA['A']['C']['G'] = 'T'; DNA_AA['A']['C']['T'] = 'T'; // Threonine
	DNA_AA['A']['A']['C'] = 'N'; DNA_AA['A']['A']['T'] = 'N'; // Asparagine
	DNA_AA['A']['A']['A'] = 'K'; DNA_AA['A']['A']['G'] = 'K'; // Lysine
	DNA_AA['A']['G']['C'] = 'S'; DNA_AA['A']['G']['T'] = 'S'; // Serine
	DNA_AA['A']['G']['A'] = 'R'; DNA_AA['A']['G']['G'] = 'R'; // Arginine
	DNA_AA['G']['T']['A'] = 'V'; DNA_AA['G']['T']['C'] = 'V'; DNA_AA['G']['T']['G'] = 'V'; DNA_AA['G']['T']['T'] = 'V'; // Valine
	DNA_AA['G']['C']['A'] = 'A'; DNA_AA['G']['C']['C'] = 'A'; DNA_AA['G']['C']['G'] = 'A'; DNA_AA['G']['C']['T'] = 'A'; // Alanine
	DNA_AA['G']['A']['C'] = 'D'; DNA_AA['G']['A']['T'] = 'D'; // Aspartic Acid
	DNA_AA['G']['A']['A'] = 'E'; DNA_AA['G']['A']['G'] = 'E'; // Glutamic Acid
	DNA_AA['G']['G']['A'] = 'G'; DNA_AA['G']['G']['C'] = 'G'; DNA_AA['G']['G']['G'] = 'G'; DNA_AA['G']['G']['T'] = 'G'; // Glycine
}


static char complement(char b)
{
	switch (b)
	{
	case 'A': return 'T'; case 'T': return 'A'; case 'G': return 'C'; case 'C': return 'G';
	case 'a': return 't'; case 't': return 'a'; case 'g': return 'c'; case 'c': return 'g';
	case 'w': return 'w'; case 'W': return 'W';
	case 's': return 's'; case 'S': return 'S';
	case 'y': return 'r'; case 'Y': return 'R';
	case 'r': return 'y'; case 'R': return 'Y';
	case 'k': return 'm'; case 'K': return 'M';
	case 'm': return 'k'; case 'M': return 'K';
	case 'b': return 'v'; case 'd': return 'h'; case 'h': return 'd'; case 'v': return 'b';
	case 'B': return 'V'; case 'D': return 'H'; case 'H': return 'D'; case 'V': return 'B';
	case 'N': return 'N'; case 'n': return 'n';

	}
	return '?';
}
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << a_value;
	return std::move(out).str();
}



//class gene

gene::gene(string id, int sta, int end) : 
	geneID(id), geneStart(sta), geneEnd(end), geneLength(end - sta), geneStrand(true), 
	numOnContig(-1) , accumDepth(0)
{
	assert(geneEnd >= geneStart);
}

gene::gene(gene* GG) : geneID(GG->geneID), geneStart(GG->geneStart), geneEnd(GG->geneEnd), geneLength(GG->geneLength),
geneStrand(GG->geneStrand), type(GG->type), translationTable(GG->translationTable),
numOnContig(GG->numOnContig), partial(GG->partial), accumDepth(GG->accumDepth)
{
}


void gene::reverseComplement (string& seq) {
	transform(
		begin(seq),
		end(seq),
		begin(seq),
		complement);
	reverse(seq.begin(), seq.end());

}



string gene::geneNT(const string& seq) {
	string ret("");
	if (seq.length() < geneEnd) {
		throw std::runtime_error("Error: Gene sequence length is less than expected: " + to_string(seq.length()) + " < " + to_string(geneEnd));
	}
	if (seq.length() < geneStart) {
		throw std::runtime_error("Error: Gene sequence length is less than expected: " + to_string(seq.length()) + " < " + to_string(geneStart));
	}
	ret = seq.substr(geneStart, geneEnd - geneStart + 1); 
	if (!geneStrand){
		reverseComplement(ret);
	} 

	return ret;
}
//assumes that the correct reading frame is given..
string gene::geneAA(const string& seq) {
	string ret("");
	if (float(seq.length()) != float(geneEnd - geneStart+1) ) {
		float tmp= float(geneEnd - geneStart) / 3.f;
		throw std::runtime_error("Error: Gene sequence length is less than expected: " + to_string(seq.length()) + " < " + to_string(geneEnd) + "-"+ to_string(geneStart)+"/3\n");
	}



	for (unsigned i = 0; i < seq.length(); i = i + 3) {
		char P1 = toupper(seq[i]); char P2 = toupper(seq[i + 1]); char P3 = toupper(seq[i + 2]);
		//P1=toupper(P1); P2=toupper(P2); P3=toupper(P3);
		char protein = '?';
		if (P1 == 'N' || P2 == 'N' || P3 == 'N') {
			protein = 'X'; //unknown
		} else if (P1 == 'X' || P2 == 'X' || P3 == 'X') {
			protein = 'X'; //unknown
		} else if (P1 == '-' || P2 == '-' || P3 == '-') {
			protein = '-'; //gap
		} else if (P1 == '.' || P2 == '.' || P3 == '.') {
			protein = '.'; //gap
		} else{
			protein = DNA_AA[P1][P2][P3];
		}
		
		ret += protein ;
		if (protein == '*' ) break;
	}

	return ret;
}



//***************************************************************

//cF->ntVariant(posN, ref, "N", -1.f); 
void fasta::ntVariant(VCFmem* vx, VariantStats* Vstats) {
	//don't include if position is not deep enough in first place..
	if (vx->getPos() < seqUse.size() && !seqUse[vx->getPos()]) { return; }
	if (vx->isINDEL()) {
		Vstats->indelCNT++;
		if (vx->majorAllele()) {
			if (vx->filtered()) {
				Vstats->unsrINDEL++; //cF->ntVariant(posN, ref, "N", -1.f); 
			}
			else {
				INDEL(vx->getPos(), vx->getRef(), vx->getAlt(), vx->getFreq());
				Vstats->indelUsed++;
			}
		}
		if (vx->filtered()) {
			Vstats->indelFILT++;
		}

		return; //skip indels for now
	}

	if (!vx->isSNP()) {
		return;//not sure what this is.. skip for now
	}

	Vstats->snpCNT++;

	//implement SNPs (maybe)
	if (vx->majorAllele()) {
		if (vx->filtered()) {
			Vstats->unsrSNP++;
			SNP(vx->getPos(), vx->getRef(), "N", vx->getFreq());
		} else {
			Vstats->SNPused++;
			SNP(vx->getPos(), vx->getRef(), vx->getAlt(), vx->getFreq());
			if (vx->conflicted()) { conflictCnt++; Vstats->conflictCnt++;
			}
		}
	} 
	if (vx->filtered()) {
		Vstats->snpFILT++;
	}


}

//structural variant (indel) of some sort
void fasta::INDEL(int pos, string r, string a, float freq)
{
	INDELcnt++;
	INDELpos.push_back(pos);
	INDELfreq.push_back(freq);
	INDELref.push_back(r);
	INDELalt.push_back(a);//for now only catalogue.. needs to be applied later..
}

//single nucleotide polymorphism
void fasta::SNP(int pos, string r, string a, float freq)
{ 
	if (freq < 0.f || a == "N") { 
		maskSeq(pos, pos + 1); //replace with N
		UnctCnt++; return; 
	}
	if (seq[pos] != r[0] && toupper(seq[pos]) != 'N') {
		char tmp = seq[pos];
		throw std::runtime_error("vcf2fna:: Error: SNP position does not match reference sequence " + header + " at position " + std::to_string(pos) + ": " + seq[pos] + " != " + r);
	}
	seq[pos] = a[0];
	SNPsCnt++; 
	SNPsPos.push_back(pos);
	SNPfreqs.push_back(freq);
}


string fasta::getMutatedHeader(bool hdTags) {
	string ret = ">" + header;
	if (!hdTags) {
		return ret;
	}
	//add tags for MGTK
	int cntNs = getNumNs();
	ret += " COV=" + to_string(seq.length()- cntNs) + " REPL=" + to_string(SNPsCnt);
	ret += " POS=";

	std::list<int>::iterator it;
	for (it = SNPsPos.begin(); it != SNPsPos.end(); ++it) {
		if (it == SNPsPos.begin()) {ret+=to_string(*it);
		}else {ret += "," + to_string(*it);
		}
	}

	ret += " FR="; 
	std::list<float>::iterator it2;
	for (it2 = SNPfreqs.begin(); it2 != SNPfreqs.end(); ++it2) {
		if (it2 == SNPfreqs.begin()) {ret+= to_string_with_precision(*it2,2);
		}else {ret+= "," + to_string_with_precision(*it2,2);
		}
	}
	ret += " FREQT=";
	//TODO create histogram of frequencies:
	

	ret += " CONFL=" + to_string(conflictCnt);


	return ret;
}  

string gene::createHDtag(const string& seq, fasta* fa, int& nonNs) {
	list<int>& SNPsPos = fa->getSNPsPos();
	list<float>& SNPfreqs = fa->getSNPfreqs();
	//create a tag emulating scores from contig2fasta.py
	std::list<int>::iterator it;
	std::list<float>::iterator it2;
	string hd("");
	hd += " D=" + to_string_with_precision(this->getAvgDepth(),2);
	hd += " P=[";
	string freqV;
	it2 = SNPfreqs.begin();
	vector<int> sumFreqs(11, 0);
	for (it = SNPsPos.begin(); it != SNPsPos.end(); ++it) {
		if (*it < geneStart || *it > geneEnd) { continue; }
		//fequence vector
		if (it == SNPsPos.begin()) {
			freqV += to_string_with_precision(*it2, 2);
		}
		else {
			freqV += "," + to_string_with_precision(*it2, 2);
		}
		sumFreqs[int(((*it2) * 10.f) + 0.5f)]++; //convert to integer
		//position vector
		if (it == SNPsPos.begin()) {
			hd += to_string(*it - geneStart);
		}
		else {
			hd += "," + to_string(*it - geneStart);
		}
		it2++;
		if (it2 == SNPfreqs.end()) {
			break;
		}
	}

	int sumMidFreqs(sumFreqs[2] + sumFreqs[3] + sumFreqs[4] + sumFreqs[5] + sumFreqs[6] + sumFreqs[7]);
	int sumFixFreqs(sumFreqs[0] + sumFreqs[1] + sumFreqs[8] + sumFreqs[9] + sumFreqs[10]);
	float oCSP(0.f);
	if (sumMidFreqs + sumFixFreqs > 0) {
		oCSP = float(sumMidFreqs) / float(sumMidFreqs + sumFixFreqs);
	}
	float CSP(0.f);
	if (nonNs == -1) {
		nonNs = seq.length() - count(seq.begin(), seq.end(), 'N');
	}
	if (nonNs>0) {
		CSP = float(sumMidFreqs) / float(nonNs);
	}
	
	hd += "] F = [" + freqV + "]";
	hd+= " oCSP=" + to_string_with_precision(CSP,3) +" CSP="+to_string_with_precision(oCSP, 3); 
	return hd;
}

string fasta::getMutatedSeq(options* opts) {
	if (mutSeqDone) { return mutSeq; }
	mutSeq = seq;
	//SNPs
	for (size_t i = 0; i < mutSeq.length(); ++i) {
		if (!seqUse[i]) {
			mutSeq[i] = 'N'; //replace with N
		}
	}
	if (INDELpos.size() == 0 || opts->reportINDELs) {
		return mutSeq;
	}
	//INDELs
	assert(INDELpos.size() == INDELalt.size());
	assert(INDELpos.size() == INDELref.size());
	assert(INDELpos.size() == INDELfreq.size());

	//std::list<int>::iterator itPos; int i = (int)(INDELpos.size() );
	//for (itPos = INDELpos.end(); itPos != INDELpos.begin(); itPos-- ) {
	//	i--;
	for (int i = (int(INDELpos.size()) - 1); i >= 0; i--) { // start from end of read that positions are not affected by earlier indels
		int pos =  INDELpos[i];//(*itPos);
		string ref = INDELref[i];
		string alt = INDELalt[i];
		if (alt == "N" || alt == ".") {
			continue; //skip uncertain indels
		}
		if (alt == "<DEL>" || alt == "<INS>") {
			cerr << "symbolic indels found:" << alt << " at position " << pos << ". Skipping these variants." << endl;
			continue; //skip symbolic indels
		}
		//apply indel
		//deletion
		int refL = (int)(ref.length());
		int altL = (int)(alt.length());
		int difLen = refL - altL;
		string obsRef = mutSeq.substr(pos, refL);
		//int delLen = refL;
		//first delete ref length from position
		if (refL > 0) {
			mutSeq.erase(pos , refL);
		}
		//and replace with alt sequence if any
		if (alt.length() > 0) {
			mutSeq.insert(pos+1 , alt);
		}
		//now apply to every gene coordinate
		geneCol->correctCoords(pos, altL, refL);
	}
	mutSeqDone = true;



	return mutSeq; 
}  //


int fasta::maskSeq(int start, int end,  bool repl) {
	//unmask sequence from start to end (exclusive)
	assert(seq.length() == seqUse.size());

	//int replCnt(0);
	if (start < 0 || end > seq.length() || start >= end) {
		throw std::out_of_range("Invalid range for unmasking: " + to_string(start) + ", " + to_string(end) + ", " + to_string(seq.length()));
	}
	for (int i = start; i < end; ++i) {
		//if (seqUse[i] == repl) { continue; }// Skip if already 'N'
		
		seqUse[i] = repl; //unmask position
		//replCnt++;
	}
	int replCnt = end - start;
	return replCnt;
}
int fasta::unmaskSeq(int start, int end) {
	//unmask sequence from start to end (exclusive)
	return maskSeq(start,end,true);
}


geneCollection::~geneCollection() { 
	for (size_t i = 0; i < genes.size(); i++) { if (genes[i] != nullptr) { delete genes[i]; } }
	for (size_t i = 0; i < genesMut.size(); i++) { if (genesMut[i] != nullptr) { delete genesMut[i]; } }
}

//gets info from depth file, apply to genes potentially in frame
void geneCollection::depthInGenes(int sta, int sto, int depth) {
	for (size_t i = 0; i < genes.size(); ++i) {
		if (genes[i]->geneEnd < sta || genes[i]->geneStart > sto) {
			continue; //not affected
		}
		//overlapping depth region
		int sta2 = max(sta, genes[i]->geneStart);
		int sto2 = min(sto, genes[i]->geneEnd);
		int overlap = sto2 - sta2 ;
		genes[i]->addAccumDepth(overlap * depth);
	}
}
void geneCollection::prepMuts() {
	if (mutsPrepared) { return; }
	mutsPrepared = true; 

	genesMut.clear(); genesMut.resize(genes.size(), nullptr);
	for (size_t i = 0; i < genes.size(); i++) {
		genesMut[i] = new gene(*genes[i]); //copy gene
	}
}


void geneCollection::correctCoords(int pos, int altL, int refL) {
	//needs to go through all genes and correct their coordinates based on indel at pos
	//only genes after the indel position are affected
	//alt:20 -ref:10 = 10   ; alt:5 - ref:10 = -5
	int geneDiff = altL - refL;
	for (size_t i = 0; i < genesMut.size(); ++i) {
		if (genesMut[i]->geneEnd < pos && genesMut[i]->geneStart < pos) {
			continue; //not affected
		}
		if (genesMut[i]->geneStart > pos) {
			//fully after indel
			genesMut[i]->geneStart += geneDiff;
			genesMut[i]->geneEnd += geneDiff;
		}
		else {
			//overlapping indel
			genesMut[i]->geneEnd += geneDiff;
			//start remains the same
		}
	}
}


void geneCollection::writeAllGenes(options* opts, string& NTs, string& AAs, 
			bool doNT, bool doAA, fasta* fa, OutputStats* Ostats) {
	//bool HDtags(opts->addHDTags); //bool skipE(opts->skipEmptyContigs);
	int NonNcnt = fa->getNumNs();// (seq.length() - count(seq.begin(), seq.end(), 'N'));
	if (NonNcnt == fa->getLength() && opts->skipEmptyContigs) {
		return; //skip empty sequences
	}
	


	string ctg = fa->getMutatedSeq(opts);
	if (ctg.empty()) {
		return; //skip empty sequences
	}
	//string ret("");
	for (size_t i = 0; i < genesMut.size(); ++i) {
		string geneSeq = genesMut[i]->geneNT(ctg);
		int geneNs = count(geneSeq.begin(), geneSeq.end(), 'N');
		int geneNonNs = (int)geneSeq.length() - geneNs;
		if (opts->skipEmptyGenes && (geneSeq.length() == 0 || 
			!(geneSeq.length() - geneNs ) )) {
			continue; //skip empty sequences
		}
		string hd = ">"+fa->header + "_" + to_string(genes[i]->getIdx()+1);
		if (opts->addHDTags) {
			hd+=genes[i]->createHDtag(geneSeq,fa, geneNonNs);
		}
		if (doNT) {
			NTs += hd + "\n";
			NTs += geneSeq + "\n";
			Ostats->totalGeneNTs += geneNonNs;
		}
		if (doAA){
			AAs += hd + "\n";
			string AAseq = genesMut[i]->geneAA(geneSeq);
			AAs += AAseq+"\n";
			int geneXs = count(AAseq.begin(), AAseq.end(), 'X');
			int geneNonXs = (int)AAseq.length() - geneXs;
			Ostats->totalGeneAAs += geneNonXs;
		}
		Ostats->totalGenes++;
	}
	
}


string fasta::write(options* opts, OutputStats* Ostats){//
	bool hdTags = opts->addHDTags;
	bool skipEmpty = opts->skipEmptyContigs;
	string ret("");
	int numNs = getNumNs();
	if (skipEmpty && numNs==getLength()){
		return ret; //skip empty sequences
	}
	prepMuts();

	Ostats->totalContigs++;
	Ostats->totalCtgNTs += (getLength() - numNs);

	ret += getMutatedHeader(hdTags) + "\n";
	ret += getMutatedSeq(opts) + "\n";
	return ret;
}

void fasta::addGene(string id, int sta, int end, string strand, string type, int transTab,string partial) {
	gene* AG = new gene(id,sta,end);
	AG->setStrand(strand);
	AG->setType(type);
	AG->setTranslationTable(transTab);
	AG->setNumOnContig(geneCol->size());
	AG->setPartial(partial);
	geneCol->push_back(AG);
}





refAssembly::refAssembly(options* opt):
	NSeqs(0), refFasta(opt->refFasta),  
	fastas(0),
	hd2ID(0), opts(opt)
{
	const std::string& filename = opt->refFasta;
	cout << "Reading reference assembly from: " << filename << endl;
	istream* file = openGZUZ(filename);
	long totalBp(0);
	std::string line; string sequence(""); string header("");
	while (std::getline(*file, line)) {
		if (line.empty() ) {
			continue; // Skip empty lines and headers
		}
		if (line[0] == '>') {					
			// Store the previous sequence before starting a new one
			if (!sequence.empty()) {
				fastas.push_back(new fasta(sequence, header));
				totalBp += sequence.length();
				sequence.clear();
				NSeqs++;
			}
			// Optionally, you can store the header or do something with it
			header = line.substr(1); // Remove '>' character
			continue;
		}
		sequence += line;
	}
	delete file;
	//add last readin batch:
	if (!sequence.empty()) {
		fastas.push_back(new fasta(sequence, header));
		totalBp += sequence.length();
		sequence.clear();
		NSeqs++;
	}


	//build index for sequences from headers
	for (size_t i = 0; i < NSeqs; ++i) {
		hd2ID[fastas[i]->getHeader()] = i;
	}

	cout << "Number of sequences: " << NSeqs << " containing "<< totalBp<< "bp."<< endl;
}

refAssembly::~refAssembly() {
	for (size_t i = 0; i < NSeqs; ++i) {
		delete fastas[i]; // Free the memory allocated for each fasta object
	}
}

void refAssembly::setFasta(string id, fasta* f){
	auto it = hd2ID.find(id);
	if (it != hd2ID.end()) {
		fastas[it->second] = f; 
	} else {
		throw std::runtime_error("Header not found: " + id);
	}
}

fasta* refAssembly::getFasta(string id)
{ 
	auto it = hd2ID.find(id);
	if (it != hd2ID.end()) {
		return fastas[it->second];//->getSeq(); 
	} else {
		throw std::runtime_error("Header not found: " + id);
		return nullptr; // This line will never be reached due to the exception above
	}
}
bool refAssembly::isSequence(string id) {
	auto it = hd2ID.find(id);
	if (it != hd2ID.end()) {
		return true;
	} else {
		return false;
	}
	return false; // This line will never be reached due to the exception above
}

int refAssembly::replaceWithNs(std::string& seq, const int start, const int end, char replaceWith) {
	int repl(0);
	
	if (start < 0 || end > seq.length() || start >= end) {
		throw std::out_of_range("Invalid range for replacement with Ns:" + to_string(start) + ", " + to_string(end) + ", "+ to_string(seq.length()));
	}
	for (int i = start; i < end; ++i) {
		if (seq[i] == 'N' || seq[i] == 'n') {
			continue; // Skip if already 'N'
		}
		seq[i] = replaceWith;
		repl++;
	}
	return repl;
}

void refAssembly::readGFF() {
	istream* file = openGZUZ(opts->gffFile);
	

	cout << "Reading GFF file from: " << opts->gffFile << endl;
	string line; string curChrom(""); fasta* curFasta(nullptr);
	int lastPos(0); string INFO("");
	while (std::getline(*file, line)) {
		if (line.empty() ) {
			continue; // Skip empty lines
		}
		if (line[0] == '#') {
			INFO = line.substr(1); // Store the INFO line
			continue; // Skip empty lines
		}
		std::istringstream iss(line);
		string chrom, ver,type, strand,other,add;
		int sta,sto;
		int TT = 11; //default type for GFF3: gene
		float score;
		string partial = "";
		if (!(iss >> chrom >> ver>> type >> sta>> sto >> score >> strand >> other >> add)) {
			throw std::runtime_error("Error reading GFF file: " + line +"\n" + opts->gffFile);
		}
		//corrections to get right coordinates in C++ string
		sta--; sto--;

		if (curChrom != chrom) {//switch to new chromosome
			curFasta = getFasta(chrom);
		}
		if (curFasta == nullptr) {
			throw std::runtime_error("Contig not found in GFF file: " + chrom);
		}
		//parse TT
		std::string segment; std::stringstream test(INFO);
		while (getline(test, segment, ';')) {
			if (segment.substr(0, 13) == "transl_table=") {
				TT = stoi(segment.substr(13));
				break;
			}
		}
		test = stringstream(add);
		while (getline(test, segment, ';')) {
			if (segment.substr(0, 8) == "partial=") {
				partial = segment.substr(8);
				break;
			}
		}

		curFasta->addGene("",sta,sto,strand,type,TT, partial);
		lastPos = sto;
	}
	//cout << "GFF filter: Number of positions replaced with N: " << lastPos << endl;
	delete file;
}


void refAssembly::processDepth(const string& inF, int minDepth) {
	
	//open connection to file
	istream* in = openGZUZ(inF);

	cout << "Reading depth file from: " << inF << endl;

	string line; string curChrom(""); //int curChromIdx(-1); string curChromSeq("");
	fasta* curFasta(nullptr); int curChromL(0);
	int lastPos(0);
	long totalChromL(0);
	long cntPosKept(0), cntPosRm(0);//counting how many positions are kept and removed based on depth profile
	int curPosKept(0), curPosRm(0);
	vector<long> posKeptPerDepth(10 , 0); //up to depth 1000


	geneCollection* genes(nullptr);
	int curGeneDepthAccum(0); int geneIdx(0);
	
	//start reading depth file line by line
	//example: "A.14.1640M2__C15_L=476= 25      104     2"
	while (std::getline(*in, line)) {
		if (line.empty()) {
			continue; // Skip empty lines
		}
		std::istringstream iss(line);
		string header;
		int sta, sto, depth;
		if (!(iss >> header >> sta >> sto >> depth)) {
			throw std::runtime_error("Error reading depth file: " + line + "\n" + inF);
		}
		if (curChrom != header) {//switch to new chromosome
			curPosRm = curChromL - curPosKept;
			cntPosRm += curPosRm;
			cntPosKept += curPosKept;
			curPosKept = curPosRm = 0;
			
			curChrom = header;
			curFasta = getFasta(curChrom);
			curChromL = curFasta->getLength();
			totalChromL += curChromL;
			genes = curFasta->getGeneCollection();
		}
		
		//some basic checks:
		if (sto > curChromL) {
			throw std::runtime_error("Error: Stop post seq length: " + line + "\n" + inF);
		}
		genes->depthInGenes(sta, sto, depth);
		//replace  seq with N where < depthThreshold
		if (depth >= minDepth) {
			curPosKept += curFasta->unmaskSeq(sta, sto);
		}
		lastPos = sto;
		//for stats only, check at which depth positions where kept
		for (int depthIdx = 0; depthIdx < posKeptPerDepth.size(); depthIdx++) {
			if (depth >= (int)depthIdx) {
				posKeptPerDepth[depthIdx] += (sto - sta);
			} else {
				break;
			}
		}
	}
	delete in;
	
	//last iteration
	curPosRm = curChromL - curPosKept;
	cntPosRm += curPosRm;
	cntPosKept += curPosKept;
	//general stats on depth filter
	cout << "Depth filter (" << minDepth << "): Number of positions kept: " << cntPosKept << " removed: " << cntPosRm << endl;

	//stats on depth filter per depth bin
	cout << "Theoretical depth filter stats per depth bin, increasing from 0 to "<< posKeptPerDepth.size() << " (of max "<< totalChromL<<"): " ;
	for (size_t depthIdx = 0; depthIdx < posKeptPerDepth.size(); depthIdx++) {
		cout << posKeptPerDepth[depthIdx];
		if (depthIdx+1 != posKeptPerDepth.size()) {cout<< ", "; }
	}
	cout << endl;


}


void refAssembly::readDepth() {

	for (size_t i = 0; i < opts->depthF.size(); i++) {
		int minDepth = opts->minDepthPar[i];
		string inF(opts->depthF[i]);
		this->processDepth(inF, minDepth);
	}
}


void refAssembly::writeOutputs() {
	if (opts->outputTypes.length() > 0) {
		cout << "Writing consensus contigs and genes to output files.." << endl;
	}

	OutputStats* Ostats = new OutputStats();

	if (opts->outputTypes.find("C") != string::npos) {
		writeContigs(Ostats);
	}
	if (opts->outputTypes.find("N") != string::npos || opts->outputTypes.find("A") != string::npos) {
		writeGenes(Ostats);
	}
	if (true || opts->outputTypes.length() == 0) {//simulate writing contigs to get bp written message
		long totalBp(0);
		for (size_t i = 0; i < fastas.size(); ++i) {
			totalBp += fastas[i]->getNonNcount();
		}
		cout << "Total bp that can be determined: " << totalBp << " in " << fastas.size() << " entries." << endl;
	}

	Ostats->printStats();
	delete Ostats;

}
void refAssembly::writeGenes(OutputStats* Ostats) {
	string outFastaNT(opts->outGeneNT);
	string outFastaAA (opts->outGeneAA); 
	bool wrAA(false), wrNT(false);
	ostream* outFileAA(nullptr); ostream* outFileNT(nullptr);

	//setting up what needs to be done..
	if (opts->outputTypes.find("A") != string::npos) {
		wrAA = true;
		if (outFastaAA.empty()) {throw std::runtime_error("No output gene file given, but requested (-oGeneAA )");}
		outFileAA = writeGZUZ(outFastaAA);
		//ostream* outFile = writeGZUZ(outFasta);
	}
	if (opts->outputTypes.find("N") != string::npos) {
		wrNT = true;
		if (outFastaNT.empty()) { throw std::runtime_error("No output gene file given, but requested (-oGeneNT )"); }
		outFileNT = writeGZUZ(outFastaNT);
	}
	
	//main translation work happens here:
	for (size_t i = 0; i < fastas.size(); ++i) {
		string NT(""), AA("");
		fastas[i]->prepMuts();
		fastas[i]->writeAllGenes(opts, NT, AA, wrNT, wrAA, Ostats);
		if (wrAA) {
			*outFileAA << AA;
		}
		if (wrNT) {
			*outFileNT << NT;
		}
	}

	//close open file connections..
	if (wrAA) {delete outFileAA;}
	if (wrNT) {delete outFileNT;}
	
}

void refAssembly::writeContigs(OutputStats* Ostats) {
	const string& outFasta(opts->outfna);
	if (outFasta.empty()) {
		cout<<"No output contig file requested..";
		return;
	}
	//could be normal or gzipped output requested via .gz
	ostream* outFile = writeGZUZ(outFasta);

	for (size_t i = 0; i < fastas.size(); ++i) {
		(*outFile) << fastas[i]->write(opts, Ostats);
	}
	//close open file connections..
	delete outFile;
}