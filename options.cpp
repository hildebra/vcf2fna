#include "options.h"

bool isGZfile(std::string fi) {
	if (fi.substr(fi.length() - 3) == ".gz") {
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
	cout << "additional flags:\n";
	cout << "  -minCallDepth [int] : minimum depth for a call to be considered\n";
	cout << "  -minCallQual [int] : minimum quality for a call to be considered\n";
	cout << "  -threads [int] : number of threads to use\n";
	cout << "  -tmp [dir] : temporary directory to use\n";
	cout << "  -h, --help : print this help message\n";
	cout << "  -v, --version : print version information\n";
	cout << "  -gff [file] : gff file to use for annotation\n";
	cout << "  -depthF [file] : depth file to use for annotation\n";
	cout << "  -oCtg [file] : output file for consensus call corrected contigs \n";
	cout << "  -oGeneNT [file] : output file for gene nucleotide sequences\n";
	cout << "  -oGeneAA [file] : output file for gene amino acid sequences\n";
	cout << "  -outType [CNA] : (C): print contigs (N): print genes on ctgs (A): print translated AA of genes\n";
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
tmp(""),inVCF(""), depthF(""), gffFile(""), outputTypes("NCA"),
threads(1),
minDepthPar(2), minCallQual(20),
minFS(0.01), minMQ0F(0.6), minBQBZ(0.05), minSP(20),
addHDTags(true), skipEmptyContigs(true)
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
		}else if (!strcmp(argv[i], "-oCtg")) {
			outfna = argv[++i];
		}else if (!strcmp(argv[i], "-oGeneNT")) {
			outGeneNT = argv[++i];
		}else if (!strcmp(argv[i], "-oGeneAA")) {
			outGeneAA = argv[++i];
		}else if (!strcmp(argv[i], "-tmp")) {
			tmp = argv[++i];
		}else if (!strcmp(argv[i], "-outType")) {
			outputTypes = argv[++i];
		}else if (!strcmp(argv[i], "-gff")) {
			gffFile = argv[++i];
		}else if (!strcmp(argv[i], "-inVCF")) {
			inVCF = argv[++i];
		}else if (!strcmp(argv[i], "-depthF")) {
			depthF = argv[++i];
		}		else if (!strcmp(argv[i], "-t")) {
			threads = atoi(argv[++i]);
		}
		else if (!strcmp(argv[i], "-minCallDepth")) {
			minDepthPar = atoi(argv[++i]);
		}
		else if (!strcmp(argv[i], "-minCallQual")) {  // no swap
			minCallQual = atoi(argv[++i]);
		}
		else {
			cerr << "Unknown option: " << argv[i] << endl;
			hasErr = true;
		}

	}


	// sanity checks
	// we need input
	if (true) {

		if (refFasta == "") {//just set some defaults
			cerr << "-ref must be specified\n";
			hasErr = true;
		}
		if (outfna == "") {//just set some defaults
			cerr << "-outfna must be specified\n";
			hasErr = true;
		}
		if (inVCF == "") {//just set some defaults
			cerr << "-inVCF must be specified\n";
			hasErr = true;
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
	cout << "ref assembly:   " << refFasta << std::endl;
	cout << "input vcf:      " << inVCF << std::endl;
	cout << "depth file:     " << depthF << std::endl;
	cout << "gff file:       " << gffFile << std::endl;
	cout << "output types:   " << outputTypes << "::";
	cout << outfna << " " << outGeneNT << " " << outGeneAA << std::endl;
	cout << "tmp dir:        " << tmp << std::endl;
	cout << "threads:        " << threads << std::endl;
	cout << "minFS:         " << minFS << std::endl;
	cout << "minMQ0F:       " << minMQ0F << std::endl;
	cout << "minBQBZ:       " << minBQBZ << std::endl;
	cout << "minSP:         " << minSP << std::endl;
	cout << "skipEmptyContigs: " << skipEmptyContigs << std::endl;	
	cout << "minCallDepth:   " << minDepthPar << std::endl;
	cout << "minCallQual:    " << minCallQual << std::endl;
	//cout << "threads:        " << threads << std::endl;
	//cout << "mode:           " << mode  << std::endl;

}
