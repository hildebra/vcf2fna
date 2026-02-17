#include "vcf.h"


VCFReader::VCFReader(options* opt, refAssembly* R): 
	refG(R), opts(opt),
	vcfFile(opts->inVCF),cF(nullptr),header(""),// curSeq(""),
	curChrom(""), cntAvContigs(0),
	lnCnt(0), //snpCNT(0), indelCNT(0), snpFILT(0),indelFILT(0), 
	//unsrSNP(0), unsrINDEL(0), SNPused(0), indelUsed(0), conflictCnt(0),
	seenCtgs(0),
	snpRepl(0), snpKept(0),
	VCF2ptr(nullptr), Vstats(nullptr)
	//VF(nullptr)
{


	Vstats = new VariantStats();

	// Open the VCF file

	// could be two files, separated by a comma. assumes second file is from PacBio
	if (vcfFile.size()==2) {
		VCFcollection * VCFcolTmp = new VCFcollection();
		istream* in2 = openGZUZ(vcfFile[1]);
		read_vcf_file(in2, 1, VCFcolTmp);
		delete in2;
		VCF2ptr = VCFcolTmp;
		istream* in1 = openGZUZ(vcfFile[0]);
		read_vcf_file(in1,0,nullptr);
		delete in1; 
		delete VCF2ptr;

	} else if (vcfFile.size()==1){//just one file
		istream* in = openGZUZ(vcfFile[0]);
		read_vcf_file(in,0,nullptr);
		delete in;
	}
	else {
		cerr << "Unknown number of vcffiles, currently not supported:\n" << combCommas(vcfFile) << endl;
	}

}


void VCFReader::read_vcf_file(istream* fp, int VCFnum, VCFcollection* Vcol) {
	vcfFields* VF = new vcfFields(opts, VCFnum);
	std::string line;
	VCFmem* vcf2nd(nullptr);
	bool DE1 = opts->debug1;
	while (getline(*fp, line)) {		
		//std::cout << line << std::endl;;
		if (line.empty()) continue;
		if (line[0] == '#') {
			if (line[1] == '#') {
				// meta line
				string test = line.substr(2, 7);
				if (line.substr(2, 7) == "contig=") {
					string tmpContig(line.substr(13, (line.find(",len") - 13)));
					cntAvContigs++;//
					
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
			if (DE1) { cerr << lnCnt << ":: " << line << endl; }
			
			//reading vcf entry..
			VCFmem* vcf = new VCFmem(line,VF);
			if (DE1) { cerr << "X"; }

			if (Vcol != nullptr) {
				Vcol->addVCF(vcf);
				continue;
			}
			//take care of loading correct fasta..
			if (vcf->getChrom() != curChrom) {
				curChrom = vcf->getChrom();
				if (seenCtgs.find(curChrom) != seenCtgs.end()) {
					throw std::runtime_error("Chromosome occurring multiple times. Is the vcf sorted?: " + curChrom);
				}
				else {
					seenCtgs[curChrom] = 1;
				}
				cF = refG->getFasta(curChrom);
			}

			//check if double vcf entry exists for position..
			if (VCF2ptr != nullptr) {
				vcf2nd = VCF2ptr->recoverVCF(curChrom, vcf->getPos());
				if (vcf2nd != nullptr) {//double vcf call!!
					vcf->evalVCFentry(opts, VF, vcf2nd);
				} else {//normal single vcf call..
					vcf->evalVCFentry(opts, VF);
				}
			} else {//normal "single" vcf call
				vcf->evalVCFentry(opts, VF);//eval if actual ntVariant/INDEL should be filtered or not
			}
			if (DE1) { cerr << "Y "; }
			
			
			//now apply variant to contig
			cF->ntVariant(vcf);

			//collect stats for this vcf entry
			collectVCFstats(vcf); 

			delete vcf;
		}
	}
	delete VF;
	if (Vcol == nullptr) {
		cout << endl << "Read " << cntAvContigs << " contigs and " << lnCnt << " entries from VCF file." << endl;
		cout << "  - Found " << Vstats->snpCNT << " SNPs and " << Vstats->indelCNT << " INDELS." << endl;
		cout << "  - Filtered " << Vstats->snpFILT << ";" << Vstats->indelFILT << " entries (SNP; INDEL)." << endl;
		cout << "  - Passed " << Vstats->SNPused << ";" << Vstats->indelUsed << " SNPs and INDELS. Conflicts resolved: "<< Vstats->conflictCnt << endl;
		cout << "  - Unsure " << Vstats->unsrSNP << ";" << Vstats->unsrINDEL << " SNPs and INDELS replaced with N." << endl;
	}
	else {
		cout << "Read and stored " << Vcol->size() << " vcf entries\n";
		cntAvContigs = 0; lnCnt = 0; //reset vars
	}
}

void VCFReader::collectVCFstats(VCFmem* vcf)
{

	//count stats up...
	if (vcf->isINDEL()) {
		Vstats->indelCNT++;
		if (vcf->filtered()) {
			Vstats->indelFILT++;
			if (vcf->majorAllele()) {
				Vstats->unsrINDEL++; //cF->ntVariant(posN, ref, "N", -1.f); 
			} //replace with N
		}
		else if (vcf->majorAllele()) {//only replace consensus..
			//currently not used
			Vstats->indelUsed++;
		}
	}

	//if (DE1) { cerr << "Z "; }

	if (vcf->isSNP()) {
		Vstats->snpCNT++;
		if (vcf->filtered()) {
			Vstats->snpFILT++;
			//high freq but unsure? replace with N..
			if (vcf->majorAllele()) {
				Vstats->unsrSNP++;
			}
		}
		else if (vcf->majorAllele()) {//only replace consensus..
			Vstats->SNPused++;
			if (vcf->conflicted()) { Vstats->conflictCnt++; }

		}
	}
}

/*
minQual(opts->minCallQual), minDep(opts->minDepthPar),
minFS(opts->minFS), minMQ0F(opts->minMQ0F),
minBQBZ(opts->minBQBZ), minSP(opts->minSP),
*/

bool VCFmem::evalVCFentry(options* opts, vcfFields* VF) {
	isFiltered = false; altFreq = 0.f;
	if (alt == "." || alt == "*" || alt == "<*>") { //nothing needed here...
		isIndel = false; isSnp = false; isFiltered = false; return false;
	}
	if (VF->filter(fields) ) {
		isFiltered = true;
	}
	//filter2 based on the INFO field
	if (SPval > (30.f + DPval / 2.f)) { isFiltered = true; } //minSP
	if (abs(BQBZval) > (3.1f + DPval / 40.f)) { isFiltered = true; } //prev minBQBZ, now following htslib rec

	//custom Filter
	if (MQ0Fval > opts->minMQ0F) { isFiltered = true; }
	if (FSval < opts->minFS) { isFiltered = true; }


	//assert(DP4[2] + DP4[3] + DP4[1] + DP4[0] == (int)DPval); //is actually not always the same..
	//two ways of calculation altFreq
	if (AF1val > 0.f) {
		altFreq = AF1val;
	}
	else if ((DP4[2] + DP4[3] + DP4[1] + DP4[0]) > 0) {
		altFreq = float(DP4[2] + DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0]);
	}
	assert(altFreq > 0);

	if (isIndel) {
		//https://www.htslib.org/workflow/filter.html
		//DV < 2 ||     IMF < 0.02+(($qual+1)/($qual+31))*(($qual+1)/($qual+31))/4 || \
    DP > ($DP/2) * (1.7 + 12/($qual+20)) || MQBZ < -(5+DP/20) || RPBZ+SCBZ > 9"
		if (IDVval < 2.f) { isFiltered = true;; }
		if (IMFval < 0.1f) { isFiltered = true;; } //skip indels with high IMF
		if (RPBZval > 6.f || SCBZval > 6.f) { isFiltered = true;; } //skip indels with high RPBZ
		if (RPBZval + SCBZval > 9.f) { isFiltered = true;; }
		if (MQBZval < -(5.f + DPval / 20.f)) { isFiltered = true;; } //skip indels with high MQBZ
		if (QUALval < opts->minCallQual) { isFiltered = true; }
		if (!isFiltered) {
			int a = 0;//DEBUG
		}
	}
	else {
		assert(ref.length() == 1); assert(alt.length() == 1);
		isSnp = true;

		//ok this seems to be a ntVariant 
		//bcftools rec filter: MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 32 || SCBZ > 3
		//or:
		//DP>2*$DP || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)
		if (MQBZval < -3.f) { isFiltered = true;; } //skip snps with high MQBZ
		if (abs(RPBZval) > 3.5f) { isFiltered = true;; } //skip snps with high RPBZ
		if (SCBZval > 2.f + DPval / 30.f) { isFiltered = true;; } //skip snps with high SCBZ
		if (QUALval < opts->minCallQual) { isFiltered = true; }

	}
	return !isFiltered;
}

bool VCFmem::evalVCFentry(options* opts, vcfFields* VF, VCFmem* v2) {
	//in this case the current vcf is converted to a "combined" vcf entry..
	isFiltered = false; altFreq = 0.f;
	bool isF2 = false;//for second read to later decide whom to trust..
	if ((alt == "." || alt == "*" || alt == "<*>") && 
			(v2->alt == "." || v2->alt == "*" || v2->alt == "<*>")) { //nothing needed here...
		isIndel = false; isSnp = false; isFiltered = false; return false;
	}

	if (VF->filter(fields) && VF->filter(v2->fields)) {
		isFiltered = true;
	}
	//combine relevant entries
	//DPval += v2->DPval;
	//SPval = (SPval + v2->SPval) / 2;
	//no eval single entries..

	//filter2 based on the INFO field
	//strand bias
	if (SPval > (30.f + DPval / 2.f)) { isFiltered = true; } //minSP
	if (v2->SPval > (30.f + v2->DPval / 2.f)) { isF2 = true; } //minSP
	if (abs(BQBZval) > (3.1f + DPval / 40.f)) { isFiltered = true; } //prev minBQBZ, now following htslib rec
	if (abs(v2->BQBZval) > (3.1f + v2->DPval / 40.f)) { isF2 = true; } //prev minBQBZ, now following htslib rec

	//custom Filter
	if (MQ0Fval > opts->minMQ0F) { isFiltered = true; }
	if (FSval < opts->minFS) { isFiltered = true; }
	if (v2->MQ0Fval > opts->minMQ0F) { isF2 = true; }
	if (v2->FSval < opts->minFS) { isF2 = true; }


	//assert(DP4[2] + DP4[3] + DP4[1] + DP4[0] == (int)DPval); //is actually not always the same..
	//two ways of calculation altFreq
	if (AF1val > 0.f) {
		altFreq = AF1val;
	} else if ((DP4[2] + DP4[3] + DP4[1] + DP4[0]) > 0) {
		altFreq = float(DP4[2] + DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0]);
	}
	assert(altFreq > 0);

	//combine freqs
	float cFreq((AF1val * DPval + v2->AF1val * v2->DPval) / (v2->DPval + DPval));
	float cFreq2(float(DP4[2] + DP4[3]+ v2->DP4[2] + v2->DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0] + v2->DP4[2] + v2->DP4[3] + v2->DP4[1] + v2->DP4[0]));



	if (isIndel) {
		//https://www.htslib.org/workflow/filter.html
		//DV < 2 ||     IMF < 0.02+(($qual+1)/($qual+31))*(($qual+1)/($qual+31))/4 || \
    DP > ($DP/2) * (1.7 + 12/($qual+20)) || MQBZ < -(5+DP/20) || RPBZ+SCBZ > 9"
		if (IMFval > 0.1f) { isFiltered = true;; } //skip indels with high IMF
		if (RPBZval > 6.f || SCBZval > 6.f) { isFiltered = true;; } //skip indels with high RPBZ
		if (RPBZval + SCBZval > 9.f) { isFiltered = true;; }
		if (MQBZval < -(5.f + DPval / 20.f)) { isFiltered = true;; } //skip indels with high MQBZ
		if (QUALval < opts->minCallQual) { isFiltered = true; }
	}
	else {
		assert(ref.length() == 1); assert(alt.length() == 1);
		isSnp = true;

		//ok this seems to be a ntVariant 
		//bcftools rec filter: MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 32 || SCBZ > 3
		//or:
		//DP>2*$DP || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)
		if (MQBZval < -3.f) { isFiltered = true; } //skip snps with high MQBZ
		if (abs(RPBZval) > 3.5f) { isFiltered = true; } //skip snps with high RPBZ
		if (SCBZval > 2.f + DPval / 30.f) { isFiltered = true;; } //skip snps with high SCBZ
		if (QUALval < opts->minCallQual) { isFiltered = true; }

		if (v2->MQBZval < -3.f) { isF2 = true; } //skip snps with high MQBZ
		if (abs(v2->RPBZval) > 3.5f) { isF2 = true; } //skip snps with high RPBZ
		if (v2->SCBZval > 2.f + DPval / 30.f) { isF2 = true;; } //skip snps with high SCBZ
		if (v2->QUALval < opts->minCallQual) { isF2 = true; }

	}
	//decide on final result
	if (isF2 && isFiltered && alt == v2->alt) {//coincidental? prob not, take this ntVariant
		isFiltered = false;
		altFreq = cFreq;
	}else if (!isF2 && isFiltered) {//second stronger (assumes PB is 2nd)
		if (alt != v2->alt) {
			//transfer values
			alt = v2->alt; altFreq = v2->altFreq; conflict = true;
		}else {//trust because it is the same alt allele
			isFiltered = false; altFreq = cFreq;
		}
	} else if (isF2 && !isFiltered) {//1st stronger (could be ill)
		if (alt != v2->alt) {
			// no transfer needed
			//alt = v2->alt; altFreq = v2->altFreq; conflict = true;
		} else {//trust because it is the same alt allele
			isFiltered = false; altFreq = cFreq;
		}
	} else if (!isF2 && !isFiltered) {
		if (alt != v2->alt) {
			conflict = true;//unreliable site, just remove..
		}
	}


	return !isFiltered;
}

VCFmem::VCFmem(const string& line, vcfFields* VF):
	chrom(""), id(""), ref(""), alt(""), filter(""),
	fieldsSet(false), QUALval(0.f),
	fields(0), DP4(4),
	AF1val(-1.f), FSval(0.f), MQ0Fval(0.f), BQBZval(0.f), SPval(0.f), 
	IDVval(10.f), IMFval(10.f),
	isIndel(false), isSnp(false),
	isFiltered(false), isUnsure(false),
	conflict(false)

{
	// Parse a single VCF line
	std::istringstream ss(line);

	ss >> chrom >> posN >> id >> ref >> alt;
	posN--;//counted differently in vcf to C++
	ss >> QUALval;
	string info, format, xtra;
	//first check if the filter is set to "PASS" or not, if not, skip this line
	ss >> filter >> info >> format >> xtra;
	if (info.substr(0, 5) == "INDEL") {
		isIndel = true;// indelCNT++;
	}
	//cerr << "U"; 
	if (!VF->isSet()) { VF-> parseFields(format); }
	//cerr << "I"; 
	this->splitXtra(xtra,VF);
	//cerr << "O"; 
	this->parseINFO(info);
	//cerr << "P"; 

	return;
}


void VCFmem::parseINFO(const string& info) {
	// Parse the DP4 field from the INFO field
	//DP4 = {0, 0, 0, 0};
	DP4[0] = DP4[1] = DP4[2] = DP4[3] = 0;
	MQ0Fval = BQBZval = SPval =  DPval = RPBZval = SCBZval = 0.f;
	IDVval = IMFval = 10.f; //set to high default value, so that if not set, it is not filtered by the IMF and IDV filter
	FSval = 1.f; AF1val = 0.f; AF2val = 0.f;
	int fieldCnt = 0;
	std::string segment; std::stringstream test(info);
	while (getline(test, segment, ';')) {
		fieldCnt++; //cerr << fieldCnt << " ";
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
			FSval = (float) stod(segment.substr(3));
		}else if (segment.substr(0, 4) == "AF1=") {
			AF1val = (float) stod(segment.substr(4));
		}else if (segment.substr(0, 4) == "AF2=") {
			AF2val = (float) stod(segment.substr(4));
		}else if (segment.substr(0, 3) == "DP=") {
			DPval = (float) stod(segment.substr(3));
		} else if (segment.substr(0, 5) == "MQ0F=") {
			MQ0Fval = (float) stod(segment.substr(5));
		} else if (segment.substr(0, 5) == "BQBZ=") {
			BQBZval = (float) stod(segment.substr(5));
		} else if (segment.substr(0, 3) == "SP=") {
			SPval = (float) stod(segment.substr(3));
		}else if (segment.substr(0, 4) == "IDV=") {
			IDVval = stoi(segment.substr(4));
		}
		else if (segment.substr(0, 4) == "IMF=") {
			IMFval = (float) stod(segment.substr(4));
		}else if (segment.substr(0, 5) == "RPBZ=") {
			RPBZval = (float) stod(segment.substr(5));
		}else if (segment.substr(0, 5) == "SCBZ=") {
			SCBZval = (float) stod(segment.substr(5));
		}else if (segment.substr(0, 5) == "MQBZ=") {
			MQBZval = (float) stod(segment.substr(5));
		}
	}
}


bool VCFmem::splitXtra(string& s, vcfFields* VF) {
	if (fields.size() == 0) { fields.resize(VF->numFields); }
	// Split the string into segments using ':' as the delimiter
	string segment; stringstream test(s);
	int i(0);
	while (getline(test, segment, ':')) {
		assert(i < fields.size());
		fields[i] = segment;
		i++;
	}
	if (i != VF->numFields) {
		cerr << "Warning: number of fields (" << i << ") in VCF xtra field does not match expected number (" << VF->getNumFields() << ")." << " :: " << s << endl;
	}
	return true;
}


//********************************************************************
//********************************************************************
void vcfFields::parseFields(const string& s) {
	if (fieldsSet) { return; }

	// Parse the fields in the VCF file
	string segment; stringstream test(s);
	std::vector<std::string> seglist;

	while (getline(test, segment, ':')) {
		seglist.push_back(segment);
	}
	for (size_t i = 0; i < seglist.size(); i++) {
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
	numFields = seglist.size(); 
	fieldsSet = true;
}



//********************************************************************
//********************************************************************

VCFmem* VCFcollection::recoverVCF(const string& cont, int pos) {
	if (cont != curContig) {
		curContig = cont;
		auto ctgSrch= VMems.find(cont);
		if (ctgSrch == VMems.end()) {
			curMems = nullptr;
			return(nullptr);
		}
		curMems = ctgSrch->second;
	} else if (curMems == nullptr) {
		return(nullptr);
	}
	auto posSrch = curMems->find(pos);
	if (posSrch == curMems->end()) {//pos not called in this vcf..
		return(nullptr);
	}
	return(posSrch->second);
}


void VCFcollection::addVCF(VCFmem* vv) {
	const string cont = vv->getChrom();
	int pos = vv->getPos();
	if (cont != curContig) {
		curContig = cont;
		auto ctgSrch = VMems.find(cont);
		if (ctgSrch == VMems.end()) {
			//create sub collection of vcfs
			curMems = new map2vcfm(0);
			//and register to contig name..
			VMems[curContig] = curMems;
		}
		else {
			curMems = ctgSrch->second;
		}
	}
	auto posSrch = curMems->find(pos);
	if (posSrch == curMems->end()) {//pos not called in this vcf..
		(*curMems)[pos] = vv;
	}else if (posSrch->second->getFreq() < vv->getFreq()) {//pos already called in this vcf..
		delete posSrch->second; (*curMems)[pos] = vv;
	}else if (posSrch->second->getFreq() > vv->getFreq()) {//pos already called in this vcf..
		delete  vv;
	}else if (posSrch->second->isINDEL() && !vv->isINDEL()) {//pos already called in this vcf..
		//currently: overwrite indel with snp call..
		//cerr << "Warning: ntVariant call overwriting indel call at " << cont << ":" << pos << endl;
		delete posSrch->second;
		(*curMems)[pos] = vv;
	}else if (vv->isINDEL()) {//remove indel as not needed..
		delete vv;//just drop
	} else {
		//should never be here.
		throw std::runtime_error("vcf entry already exists: " + cont + to_string(pos));
	}

}

VCFcollection::~VCFcollection() {
	for (auto it : VMems) {
		delete it.second;
	}
}

size_t VCFcollection::size() {
	size_t s(0); 
	if (VMems.size() == 0) { return 0; }
	for (auto X : VMems) { 
		s += X.second->size(); 
	} 
	return s; 
}


