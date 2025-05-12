#include "fasta.h"



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

void fasta::SNP(int pos, string r, string a, float freq)
{ 
	if (seq[pos] != r[0] && toupper(seq[pos]) != 'N') {
		char tmp = seq[pos];
		throw std::runtime_error("Error: SNP position does not match reference sequence "+ header +" at position " + std::to_string(pos) + ": " + seq[pos] + " != " + r);
	}
	seq[pos] = a[0]; 
	if (freq < 0.f || a == "N") { UnctCnt++; return; }
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

string gene::createHDtag(const string& seq, list<int>& SNPsPos, list<float>& SNPfreqs,
	int& nonNs) {
	//create a tag emulating scores from contig2fasta.py
	std::list<int>::iterator it;
	std::list<float>::iterator it2;
	string hd("");
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

void fasta::writeAllGenes(options* opts, string& NTs, string& AAs, bool doNT, bool doAA) {
	bool HDtags(opts->addHDTags); bool skipE(opts->skipEmptyContigs);
	string ctg = getMutatedSeq();
	if (ctg.empty()) {
		return; //skip empty sequences
	}
	int NonNcnt(seq.length() - count(seq.begin(), seq.end(), 'N'));
	if (NonNcnt == 0 && skipE) {
		return; //skip empty sequences
	}
	//string ret("");
	for (size_t i = 0; i < genes.size(); ++i) {
		string geneSeq = genes[i]->geneNT(ctg);
		string hd = ">"+header + "_" + to_string(genes[i]->getIdx()+1);
		if (HDtags) {
			hd+=genes[i]->createHDtag(geneSeq, SNPsPos, SNPfreqs, NonNcnt);
		}
		if (doNT) {
			NTs += hd + "\n";
			NTs += geneSeq + "\n";
		} 
		if (doAA){
			AAs += hd + "\n";
			AAs += genes[i]->geneAA(geneSeq)+"\n";
		}
	}
	
}


string fasta::write(bool hdTags,bool skipEmpty) {
	string ret("");
	if (skipEmpty && getNumNs()==getLength()){
		return ret; //skip empty sequences
	}
	ret += getMutatedHeader(hdTags) + "\n";
	ret += getMutatedSeq() + "\n";
	return ret;
}

void fasta::addGene(string id, int sta, int end, string strand, string type, int transTab,string partial) {
	gene* AG = new gene(id,sta,end);
	AG->setStrand(strand);
	AG->setType(type);
	AG->setTranslationTable(transTab);
	AG->setNumOnContig(genes.size());
	AG->setPartial(partial);
	genes.push_back(AG);
}





refAssembly::refAssembly(options* opt):
	NSeqs(0), refFasta(opt->refFasta),  
	fastas(0),
	hd2ID(0), opts(opt)
{
	const std::string& filename = opt->refFasta;
	cout << "Reading reference assembly from: " << filename << endl;
	istream* file = openGZUZ(filename);
	std::string line; string sequence(""); string header("");
	while (std::getline(*file, line)) {
		if (line.empty() ) {
			continue; // Skip empty lines and headers
		}
		if (line[0] == '>') {					
			// Store the previous sequence before starting a new one
			if (!sequence.empty()) {
				fastas.push_back(new fasta(sequence, header));
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
		sequence.clear();
		NSeqs++;
	}


	//build index for sequences from headers
	for (size_t i = 0; i < NSeqs; ++i) {
		hd2ID[fastas[i]->getHeader()] = i;
	}

	cout << "Number of sequences: " << NSeqs << endl;
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
	cout << "GFF filter: Number of positions replaced with N: " << lastPos << endl;
	delete file;
}

void refAssembly::readDepth() {
	
	//open connection to file
	istream* in = openGZUZ(opts->depthF);
	

	cout << "Reading depth file from: " << opts->depthF << endl;
	int minDepth = opts->minDepthPar;

	string line; string curChrom(""); int curChromIdx(-1); string curChromSeq("");
	int lastPos(0);
	int cntPosKept(0), cntPosRm(0);//counting how many positions are kept and removed based on depth profile
	int curPosKept(0), curPosRm(0);
	while (std::getline(*in, line)) {
		if (line.empty()) {
			continue; // Skip empty lines
		}
		std::istringstream iss(line);
		string header;
		int sta,sto, depth;
		if (!(iss >> header >> sta>> sto>>depth)) {
			throw std::runtime_error("Error reading depth file: " + line +"\n" + opts->depthF);
		}
		if (curChrom != header) {//switch to new chromosome
			if (curChromIdx >= 0) {
				//next chromosome: copy old results over..
				if (lastPos < fastas[curChromIdx]->getLength()) {//any missing bases at the end of the sequence?:
					int repl = replaceWithNs(curChromSeq, lastPos, fastas[curChromIdx]->getLength());
					curPosRm += repl; curPosKept += (fastas[curChromIdx]->getLength() - lastPos - repl);
					//cerr << curChrom << endl;//DEBUG
				}
				fastas[curChromIdx]->setSeq(curChromSeq);
				//collect stats..
				cntPosKept += curPosKept; cntPosRm += curPosRm;
			}
			curPosKept = 0; curPosRm = 0;
			curChrom = header;
			auto it = hd2ID.find(header);
			if (it != hd2ID.end()) {
				curChromIdx = it->second;
			} else {
				throw std::runtime_error("Header not found in depth file: " + header);
			}
			curChromSeq = fastas[curChromIdx]->getSeq();
		}
		//some baesic checks:
		if (sto > curChromSeq.length()) {
			throw std::runtime_error("Error reading depth file: " + line + "\n" + opts->depthF);
		}
		//replace  seq with N where < depthThreshold
		if (lastPos < sta) {
			int repl = replaceWithNs(curChromSeq, 0, sta);
			curPosRm += repl; curPosKept += (sta - repl);
		}
		if (depth < minDepth) {
			int repl = replaceWithNs(curChromSeq, sta, sto);
			curPosRm += repl; curPosKept += (sto - sta - repl);
		} else {
			;//nothing...
			curPosKept += (sto - sta);
		}
		lastPos = sto;


	}

	cout << "Depth filter ("<< minDepth << "): Number of positions kept: " << cntPosKept << " removed: "<< cntPosRm << endl;
	//if (file.is_open()) {file.close();}
	delete in;
}


void refAssembly::writeOutputs() {
	cout << "Writing consensus contigs and genes to output files.." << endl;
	if (opts->outputTypes.find("C") != string::npos) {
		writeContigs();
	}
	if (opts->outputTypes.find("N") != string::npos || opts->outputTypes.find("A") != string::npos) {
		writeGenes();
	}

}
void refAssembly::writeGenes() {
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
		fastas[i]->writeAllGenes(opts, NT, AA, wrNT, wrAA);
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

void refAssembly::writeContigs() {
	const string& outFasta(opts->outfna);
	if (outFasta.empty()) {
		cout<<"No output contig file requested..";
		return;
	}
	//could be normal or gzipped output requested via .gz
	ostream* outFile = writeGZUZ(outFasta);

	for (size_t i = 0; i < fastas.size(); ++i) {
		(*outFile) << fastas[i]->write(opts->addHDTags,opts->skipEmptyContigs);
	}
	//close open file connections..
	delete outFile;
}