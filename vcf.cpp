#include "vcf.h"


VCFReader::VCFReader(options* opt, refAssembly* R) : 
	refG(R), opts(opt),
	vcfFile(opts->inVCF),cF(nullptr),header(""),// curSeq(""),
	curChrom(""), cntAvContigs(0),
	lnCnt(0), snpCNT(0), indelCNT(0), snpFILT(0),indelFILT(0), unsrSNP(0), unsrINDEL(0),
	minQual(opts->minCallQual), minDep(opts->minDepthPar),
	minFS(opts->minFS), minMQ0F(opts->minMQ0F),
	minBQBZ(opts->minBQBZ), minSP(opts->minSP),
	seenCtgs(0),
	chrom(""), pos(""), id(""), ref(""), alt(""), qual(""), filter(""), info(""), format(""),xtra(""),
	snpRepl(0), snpKept(0),
	fieldsSet(false),
	GT(-1), DP(-1), SP(-1), ADF(-1), ADR(-1), AD(-1), BQBZ(-1), IDV(-1), IMF(-1), MQ0F(-1), FS(-1),
	numFields(0),
	fields(0), DP4(4),
	AF1val(-1.f),FSval(0.f), MQ0Fval(0.f), BQBZval(0.f), SPval(0.f), IDVval(0.f), IMFval(0.f)
{
	// Open the VCF file

	/*
	std::ifstream fp;
	std::istream& in = (vcfFile != "-")
		? [&]() -> std::istream& {
		fp.open(vcfFile);
		if (!fp) abort();
		return fp;
		}()
		: std::cin;
*/
	istream* in = openGZUZ(vcfFile);
	read_vcf_file(in);
	delete in;

	//some stat message at end..
	//cout << "Parsed " << snpCNT << " SNPs and "<< indelCNT <<" INDELS from " << lnCnt << " VCF lines." << endl;

	/*
	if (fp.is_open()) {
		fp.close();
	}
	*/
}


void VCFReader::parseFields(const string& s) {
	
	if (fieldsSet) { return; }

	// Parse the fields in the VCF file
	string segment; stringstream test(s);
	std::vector<std::string> seglist;

	while (getline(test, segment, ':')){
		seglist.push_back(segment);
	}
	for (size_t i=0;i<seglist.size();i++){
		if (seglist[i] == "GT") GT = i;
		else if (seglist[i] == "FS") FS = i;// Fisher's exact test P-value to detect strand bias 
		else if (seglist[i] == "MQ0F") MQ0F = i;//Fraction of reads with zero mapping quality
		else if (seglist[i] == "IDV") IDV = i;//Maximum number of raw reads supporting an indel 
		else if (seglist[i] == "IMF") IMF = i;//Maximum fraction of raw reads supporting an indel 
		else if (seglist[i] == "BQBZ") BQBZ = i;//Mann-Whitney U test of Base Quality Bias 
		else if (seglist[i] == "DP") DP = i;//Number of high-quality bases
		else if (seglist[i] == "SP") SP = i;//Phred-scaled strand bias P-value
		else if (seglist[i] == "ADF") ADF = i;//Allelic depths on the forward strand 
		else if (seglist[i] == "ADR") ADR = i;//Allelic depths on the reverse strand
		else if (seglist[i] == "AD") AD = i;//Allelic depth
	}
	numFields = seglist.size(); fields.resize(numFields);
	fieldsSet = true;
}

bool VCFReader::splitXtra(const string& s) {
	// Split the string into segments using ':' as the delimiter
	string segment; stringstream test(s);
	int i = 0;
	while (getline(test, segment, ':')) {
		fields[i] = segment;
		i++;
	}
	if (i != numFields) {
		cerr << "Warning: number of fields (" << i << ") in VCF xtra field does not match expected number (" << numFields << ")." << " :: " << s << endl;
	}
	return true;
}


void VCFReader::read_vcf_file(std::istream* fp) {
	std::string line;
	while (getline(*fp, line)) {
		
		//std::cout << line << std::endl;;
		if (line.empty()) continue;
		if (line[0] == '#') {
			if (line[1] == '#') {
				// meta line
				string test = line.substr(2, 7);
				if (line.substr(2, 7) == "contig=") {
					string tmpContig(line.substr(13, (line.find(",len") - 13)));
					cntAvContigs++;// .push_back(line.substr(13, (line.find(",len") - 13)));
					if (!refG->isSequence(tmpContig)) {
						cerr << "Warning: contig " << tmpContig << " not found in reference assembly." << endl;
					}
				}
				//meta.push_back(line);
			} else {
				// header line
				header = line;
				string test = header.substr(0,45);
				if (header.substr(0,45) != "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT") {
					cerr << "Warning: VCF header does not match expected format." << endl<< header<< endl;
				}
			}
		} else {
			lnCnt++;
			// data line
			//data.push_back(line);
			parse_single_vcf(line);
		}
	}
	cout << endl<<"Read " << cntAvContigs << " contigs and "<< lnCnt<<" entries from VCF file." << endl;
	cout << "Found " << snpFILT << "/" << snpCNT << " SNPs and " << indelFILT << "/" << indelCNT << " INDELS."<<endl;
	cout << unsrSNP<<","<< unsrINDEL << " SNPs and INDELS replaced with N." << endl;
}


void VCFReader::parse_single_vcf(string line) {
	// Parse a single VCF line
//std::cout << line << std::endl;
	std::istringstream ss(line);

	ss >> chrom >> pos >> id >> ref >> alt;
	int posN = stoi(pos);
	bool isFILT(false);
	posN--;//counted differently in vcf to C++
	if (alt == "." || alt =="*" || alt=="<*>") { isFILT = false; return; } //nothing needed here...
	ss >> qual;
	int QUALval(stoi(qual));
	

	//take care of loading correct fasta..
	if (chrom != curChrom) {
		//not needed: working with pointers..
		//if (cF != nullptr) { refG->setFasta(curChrom, cF); }
		curChrom = chrom;
		if (seenCtgs.find(chrom) != seenCtgs.end()) {
			throw std::runtime_error("Chromosome occurring multiple times. Is the vcf sorted?: " + chrom);
		}
		else {
			seenCtgs[chrom] = 1;
		}
		cF = refG->getFasta(chrom);
	}


	//first check if the filter is set to "PASS" or not, if not, skip this line
	ss >> filter >> info >> format >> xtra;
	if (!fieldsSet) { parseFields(format); }
	splitXtra(xtra);
	//currently deactivated: "," format ??
	if (false && DP > -1) {
		int DPv = stoi(fields[DP]);
		if (SP != -1 && (DPv - stoi(fields[SP]) < minDep)) { isFILT = true; }
		if (BQBZ != -1 && (stof(fields[BQBZ]) < minBQBZ)) { isFILT = true; }
		if (SP != -1 && (stof(fields[SP]) > minSP)) { isFILT = true; }
		//these are in info field..
		//if (MQ0F != -1 && (DPv - stoi(fields[MQ0F]) < minDep)) { return; }
		//if (FS != -1 && (stof(fields[FS]) < minFS)) { return; }
		//if (MQ0F != -1 && (stof(fields[MQ0F]) > minMQ0F)) { return; }
		//DPval = stoi(fields[DP]);
	}

	parseINFO(info);
	//filter2 based on the INFO field
	if (SPval > (30.f + DPval / 2.f)  ) { isFILT = true; } //minSP
	if (abs(BQBZval) > (3.1f + DPval / 40.f) ) { isFILT = true; } //prev minBQBZ, now following htslib rec
	
	//custom Filter
	if (MQ0Fval > minMQ0F) { isFILT = true; }
	if (FSval < minFS) { isFILT = true; }
	
	float altFreq = 0.f;
	//assert(DP4[2] + DP4[3] + DP4[1] + DP4[0] == (int)DPval); //is actually not always the same..
	//two ways of calculation altFreq
	if (AF1val > 0.f) { 
		altFreq = AF1val; 
	}else {
		altFreq = float(DP4[2] + DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0]);
	}
	assert(altFreq > 0);




	if (info.substr(0, 5) == "INDEL") {
		indelCNT++;
		//https://www.htslib.org/workflow/filter.html
		//DV < 2 ||     IMF < 0.02+(($qual+1)/($qual+31))*(($qual+1)/($qual+31))/4 || \
    DP > ($DP/2) * (1.7 + 12/($qual+20)) || MQBZ < -(5+DP/20) || RPBZ+SCBZ > 9"
		if (IMFval > 0.1f) { return; } //skip indels with high IMF
		if (RPBZval > 6.f || SCBZval > 6.f) { return; } //skip indels with high RPBZ
		if (RPBZval + SCBZval > 9.f) { return; }
		if (MQBZval < -(5.f + DPval / 20.f)) { return; } //skip indels with high MQBZ
		if (QUALval < minQual) { isFILT = true; } 
		if (isFILT) {
			indelFILT++;
			if (altFreq > 0.51f || AF2val > 0.51f) { unsrINDEL++; cF->SNP(posN, ref, "N",-1.f); } //replace with N
			return;
		}


		return;//skip indels
	} else {
		assert(ref.length() == 1); assert(alt.length() == 1);
		snpCNT++;
		//ok this seems to be a SNP 
		//bcftools rec filter: MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 32 || SCBZ > 3
		//or:
		//DP>2*$DP || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)
		if (MQBZval < -3.f) { return; } //skip snps with high MQBZ
		if (abs(RPBZval) > 3.5f) { return; } //skip snps with high RPBZ
		if (SCBZval > 2.f + DPval/30.f) { return; } //skip snps with high SCBZ
		if (QUALval < minQual) { isFILT = true; }


		if (isFILT) {
			snpFILT++;
			//high freq but unsure? replace with N..
			if (altFreq > 0.51f || AF2val>0.51f) { unsrSNP++; cF->SNP(posN, ref, "N",-1.f); }
			return;
		}


		if (altFreq > 0.51f) {//only replace consensus..
			cF->SNP(posN, ref, alt,altFreq);
			snpCNT++;
		}
	}


}


void VCFReader::parseINFO(const string& info) {
	// Parse the DP4 field from the INFO field
	//DP4 = {0, 0, 0, 0};
	DP4[0] = DP4[1] = DP4[2] = DP4[3] = 0;
	MQ0Fval = BQBZval = SPval = IDVval = IMFval = DPval = RPBZval = SCBZval = 0.f;
	FSval = 1.f; AF1val = 0.f; AF2val = 0.f;
	std::string segment; std::stringstream test(info);
	while (getline(test, segment, ';')) {
		if (segment.substr(0, 4) == "DP4=") {
			std::string dp4 = segment.substr(4);
			std::stringstream dp(dp4);
			std::string value;int i = 0;
			while (getline(dp, value, ',')) {
				DP4[i] = stoi(value);
				i++;
			}
			//break;
		}else if (segment.substr(0, 3) == "FS=") {
			FSval = stof(segment.substr(3));
		}else if (segment.substr(0, 4) == "AF1=") {
			AF1val = stof(segment.substr(4));
		}else if (segment.substr(0, 4) == "AF2=") {
			AF2val = stof(segment.substr(4));
		}else if (segment.substr(0, 3) == "DP=") {
			DPval = stof(segment.substr(3));
		} else if (segment.substr(0, 5) == "MQ0F=") {
			MQ0Fval = stof(segment.substr(5));
		} else if (segment.substr(0, 5) == "BQBZ=") {
			BQBZval = stof(segment.substr(5));
		} else if (segment.substr(0, 3) == "SP=") {
			SPval = stof(segment.substr(3));
		}else if (segment.substr(0, 4) == "IDV=") {
			IDVval = stoi(segment.substr(4));
		}
		else if (segment.substr(0, 4) == "IMF=") {
			IMFval = stof(segment.substr(4));
		}else if (segment.substr(0, 5) == "RPBZ=") {
			RPBZval = stof(segment.substr(5));
		}else if (segment.substr(0, 5) == "SCBZ=") {
			SCBZval = stof(segment.substr(5));
		}else if (segment.substr(0, 5) == "MQBZ=") {
			MQBZval = stof(segment.substr(5));
		}
	}
}