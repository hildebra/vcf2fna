#include "options.h"







vector<string> splitByComma(const string& str) {
	vector<string> result;
	stringstream ss(str);
	string item;
	while (getline(ss, item, ',')) {
		result.push_back(item);
	}
	return result;
}
vector<int> splitByCommaI(const string& str) {
	vector<int> result;
	stringstream ss(str);
	string item;
	while (getline(ss, item, ',')) {
		result.push_back(atoi(item.c_str()));
	}
	return result;
}

string combCommas(const vector<string>& vec) {
	string result;
	for (size_t i = 0; i < vec.size(); ++i) {
		result += vec[i];
		if (i != vec.size() - 1) {
			result += ",";
		}
	}
	return result;
}
string combCommas(const vector<int>& vec) {
	string result;
	for (size_t i = 0; i < vec.size(); ++i) {
		result += to_string(vec[i]);
		if (i != vec.size() - 1) {
			result += ",";
		}
	}
	return result;
}




bool isGZfile(std::string fi) {
	if (fi.length()>=3 && fi.substr(fi.length() - 3) == ".gz") {
		return true;
	}
	return false;
}

istream* openGZUZ(const string& inF) {
	istream* in(nullptr);
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in); 
		//cout << "Straming gzip input on the fly\n";
#else
		cout << "gzip in not supported in your vcf2fna build\n"; exit(50);
#endif
	}
	else { 
		in = new ifstream(inF.c_str()); 
	}

	if (!(*in)) { 
		throw std::runtime_error("Could not open file: " + inF); 
	}
	return in;
}

ostream* writeGZUZ(const string& outF) {
	ostream* out;
	if (isGZfile(outF)) {
#ifdef _gzipread
		out = new ogzstream(outF.c_str(), ios::out);
		//cout << "Writing gzip'd matrix " << outF << endl;
#else
		cout << "gzip out not supported in your vcf2fna build\n"; exit(51);
#endif
	}
	else { out = new ofstream(outF); }

	if (!(*out)) { 
		throw std::runtime_error("Could not open file: " + outF); 
	}
	return out;
}



void helpMsg() {
	cout << "basic usage: ./vcf2fasta -ref [ref.fasta[.gz]] -v [vcf[.gz]] -depthF [.bed[.gz]] -gff [.gff] -oCtg [fna]\n";
	cout << "I/O flags:\n";
	cout << "  -ref [file[.gz]] : reference assembly of metagenome\n";
	cout << "  -v [file[.gz]] : vcf output from samtools mpileup or similar\n";
	cout << "  -depthF [.bed[.gz]] : bed formatted depth of mapping to ref (can be created with e.g. samtools)\n";
	cout << "  -gff [file[.gz]] : genes predicted on ref contigs. These will then be translated.\n";
	cout << "  -tmp [dir] : temporary directory to use\n";
	cout << "  -oCtg [file] : output file for consensus call corrected contigs \n";
	cout << "  -oGeneNT [file] : output file for gene nucleotide sequences\n";
	cout << "  -oGeneAA [file] : output file for gene amino acid sequences\n";
	cout << "  -outType [CNA] : (C): print contigs (N): print genes on ctgs (A): print translated AA of genes\n";
	cout << "  -keepEmptyCtgs : keep empty contigs in output\n";
	cout << "  -keepEmptyGenes : keep empty genes in output\n";
	cout << "SNP/SV call filtering:\n";
	cout << "  -minCallDepth [int] : minimum depth for a call to be considered\n";
	cout << "  -minCallQual [int] : minimum quality for a call to be considered\n";
	cout << "  -seqPlatform [platform] : sequencing platform (ill, PB, ONT)\n";
	cout << "  -minCallQualAdaptive[float] : filter SNPs < Depth * minCallQualAdaptive. Default: 0\n";
	cout << "  -depthFilterScale [float] : if DP < mean contig depth *x, filter. Default: 0.25\n";
	cout << "  -indelRange [int] : SNPs within X nucleotide range of an INDEL will be filtered\n";
	cout << "  -skipINDELs : do not report INDELs in the output sequences\n";
	cout << "  -minorAlleleNoMask : do not mask nucleotides with a minor allele (that passes filtering step)\n";
	cout << "others:\n";
	cout << "  -threads [int] : number of threads to use\n";
	cout << "  -h, --help : print this help message\n";
	cout << "  -v, --version : print version information\n";
}

void versionMsg() {
	cout << "vcf2fasta v" << vcf2fnaVERSION << endl;
	cout << "compiled on: " << __DATE__ << " at " << __TIME__ << endl;
	cout << "compiled with C++ v" << __cplusplus << endl;
	exit(0);
}

void options::announce() {
	cout << "------------------------------------ " << std::endl;
	cout << "vcf2fasta v" << vcf2fnaVERSION << endl;
	print_details();
	cout << "------------------------------------ " << std::endl;
}

options::options(int argc, char** argv) :refFasta(""),
outfna(""), outGeneNT(""), outGeneAA(""),
tmp(""),inVCF(0), depthF(0), gffFile(""), outputTypes(""),
seqPlatform(0),
threads(1), indelRange(5),
minDepthPar(0), minCallQual(20),
minFS(0.01), minMQ0F(0.6), minBQBZ(0.05), minSP(20),
addHDTags(true), skipEmptyContigs(true), skipEmptyGenes(true),
debug1(false), reportINDELs(true), 
maskMinorAllele(true), 
minCallQualAdaptive(0.f), depthFilterScale(0.25f)
{


	bool hasErr = false;


	if (argc == 0) {
		cerr << "Not enough options given to clusterMAGs"; exit(23);
	}//

	if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
		helpMsg();
		exit(0);
	}else if (!strcmp(argv[1], "version") || !strcmp(argv[1], "-version") || !strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) {
		versionMsg();
		exit(0);
	}
	for (int i = 1; i < argc; i++)
	{
		
		
		if (!strcmp(argv[i], "-ref")) {
			refFasta = argv[++i];
		} else if (!strcmp(argv[i], "-oCtg")) {
			outfna = argv[++i];
		} else if (!strcmp(argv[i], "-oGeneNT")) {
			outGeneNT = argv[++i];
		} else if (!strcmp(argv[i], "-oGeneAA")) {
			outGeneAA = argv[++i];
		} else if (!strcmp(argv[i], "-tmp")) {
			tmp = argv[++i];
		} else if (!strcmp(argv[i], "-outType")) {
			outputTypes = argv[++i];
		} else if (!strcmp(argv[i], "-gff")) {
			gffFile = argv[++i];
		} else if (!strcmp(argv[i], "-inVCF")) {
			inVCF = splitByComma(argv[++i]);
		}else if (!strcmp(argv[i], "-seqPlatform")) {
			seqPlatform = splitByComma(argv[++i]);
		} else if (!strcmp(argv[i], "-depthF")) {
			depthF = splitByComma(argv[++i]);
		} else if (!strcmp(argv[i], "-t")) {
			threads = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "-minCallDepth")) {
			vector<string> minDepthParT = splitByComma(argv[++i]);
			minDepthPar.resize(minDepthParT.size());
			for (size_t i = 0; i < minDepthParT.size(); i++) {
				minDepthPar[i] = atoi(minDepthParT[i].c_str());
			}
		}
		else if (!strcmp(argv[i], "-minCallQual")) {  // no swap
			minCallQual = atoi(argv[++i]);
		}
		else if (!strcmp(argv[i], "-indelRange")) {  // no swap
			indelRange = atoi(argv[++i]);
			
		} else if (!strcmp(argv[i], "-minFS")) {
			minFS = stof(argv[++i]);
		} else if (!strcmp(argv[i], "-minMQ0F")) {
			minMQ0F = stof(argv[++i]);
		} else if (!strcmp(argv[i], "-minBQBZ")) {
			minBQBZ = stof(argv[++i]);
		} else if (!strcmp(argv[i], "-minSP")) {
			minSP = stof(argv[++i]);
		}else if (!strcmp(argv[i], "-keepEmptyCtgs")) {
			skipEmptyContigs = false;
		}
		else if (!strcmp(argv[i], "-keepEmptyGenes")) {
			skipEmptyGenes = false;
		}
		else if (!strcmp(argv[i], "-skipINDELs")) {
			reportINDELs = false;
		}
		else if (!strcmp(argv[i], "-minorAlleleNoMask")) {
			maskMinorAllele = false;
		}
		else if (!strcmp(argv[i], "-minCallQualAdaptive")) {
			minCallQualAdaptive = stof(argv[++i]);
		}
		else if (!strcmp(argv[i], "-depthFilterScale")) {
			depthFilterScale = stof(argv[++i]);
		}


		

		else if (!strcmp(argv[i], "-debug1")) {
			debug1 = true;
			cerr << "\n\nDEBUG1 active\n\n";
		}
		else {
			cerr << "Unknown option: " << argv[i] << endl;
			hasErr = true;
		}

	}

	if (seqPlatform.size() == 0 && inVCF.size() == 1) {
		seqPlatform.push_back("ill");//assume illumina by default
	}
	if (minDepthPar.size() == 0 && inVCF.size() == 1) {
		minDepthPar.push_back(2);//assume 2
	}


	// sanity checks
	// we need input
	if (true) {

		if (refFasta == "") {//just set some defaults
			cerr << "-ref must be specified\n";
			hasErr = true;
		}
		/*
		if (outfna == "") {//just set some defaults
			cerr << "-outfna must be specified\n";
			hasErr = true;
		}
		*/
		if (inVCF.size() == 0) {//just set some defaults
			cerr << "-inVCF must be specified\n";
			hasErr = true;
		}
		if (gffFile.empty()) {
			cerr << "-gff must be specified\n";
			hasErr = true;
		}
		if (depthF.size() > 2) {
			cerr << "Max 2 depth/vcf files are allowed!\n";
			hasErr = true;
		}
		if (depthF.size() != inVCF.size()) {
			cerr << "Error: -depthF and -inVCF must have the same num of \",\"\n";
			hasErr = true;
		}
		if (depthF.size() != seqPlatform.size()) {
			cerr << "Error: -depthF and -seqPlatform must have the same num of \",\"\n";
			hasErr = true;
		}
		if (depthF.size() != minDepthPar.size()) {
			cerr << "Error: -depthF and -minCallDepth must have the same num of \",\"\n";
			hasErr = true;
		}

		if (outputTypes == "") {
			if (!outfna.empty()) { outputTypes += "C"; }
			if (!outGeneAA.empty()) { outputTypes += "A"; }
			if (!outGeneNT.empty()) { outputTypes += "N"; }
		}

		if (hasErr) {
			cerr << "Error in option parsing.\nUse \"vcf2fasta -h\" to get full help.\n";
			exit(98);
		}
	}

	announce();


}

void options::print_details() {

	// print run mode:
	cout << "ref assembly, gff file: " << refFasta <<", "<< gffFile << std::endl;
	cout << "seqPlatform:      " << combCommas(seqPlatform) << std::endl;
	cout << "input vcf:    " << combCommas(inVCF) << std::endl;
	cout << "depth file:   " << combCommas(depthF) << std::endl;
	cout << "tmp dir:          " << tmp << std::endl;
	cout << "output types: \"" << outputTypes << "\" ::";
	cout << outfna << " " << outGeneNT << " " << outGeneAA << std::endl;
	if (outfna.empty() && outGeneNT.empty() && outGeneAA.empty()){
		cout << "  (no output files specified - theoretical statistics reported)" << std::endl;
	}
	//cout << "threads:          " << threads << std::endl;
	//cout << "minFS:            " << minFS << std::endl;
	//cout << "minMQ0F:          " << minMQ0F << std::endl;
	//cout << "minBQBZ:          " << minBQBZ << std::endl;
	//cout << "minSP:            " << minSP << std::endl;
	cout << "skipEmptyContigs, reportINDELs: " << skipEmptyContigs << ", " << reportINDELs << std::endl;
	cout << "minCallDepth, minCallQual, maskMinorAllele:     " << combCommas(minDepthPar) << ", "<< minCallQual << ", " << maskMinorAllele<< std::endl;
	cout << "depthFilterScale, minCallQualAdaptive:      " << depthFilterScale << ", " << minCallQualAdaptive << std::endl;
	//cout << "threads:        " << threads << std::endl;
	//cout << "mode:           " << mode  << std::endl;

}
