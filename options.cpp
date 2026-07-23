#include "options.h"
#include <cerrno>
#include <cctype>
#include <climits>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <unordered_map>

namespace {

const char* const recognizedOptions[] = {
	"-ref", "-inVCF", "-depthF", "-gff", "-oCtg", "-oGeneNT", "-oGeneAA",
	"-tmp", "-outType", "-seqPlatform", "-vcfFilterPolicy", "-t",
	"-minCallDepth", "-minAltReads", "-minAltFreq", "-minCallQual",
	"-indelRange", "-maxRefMismatches", "-minMQ0F", "-minBQBZ", "-minSP", "-minCallQualAdaptive",
	"-depthFilterScale", "-maxDepthFilterScale", "-keepEmptyCtgs",
	"-keepEmptyGenes", "-skipINDELs", "-minorAlleleNoMask", "-debug1",
	"-h", "--help", "-v", "--version"
};

bool isRecognizedOption(const char* value) {
	if (value == nullptr) return false;
	for (const char* option : recognizedOptions) {
		if (strcmp(value, option) == 0) return true;
	}
	return false;
}

const char* takeOptionValue(int& index, int argc, char** argv, bool& hasErr) {
	const string option(argv[index]);
	if (index + 1 >= argc || isRecognizedOption(argv[index + 1])) {
		cerr << "Missing value for " << option << "\n";
		hasErr = true;
		return nullptr;
	}
	++index;
	if (argv[index][0] == '\0') {
		cerr << "Empty value for " << option << "\n";
		hasErr = true;
		return nullptr;
	}
	return argv[index];
}

bool parseInteger(const char* text, const string& option, long long minimum,
	long long maximum, long long& result, bool& hasErr) {
	if (text == nullptr) return false;
	errno = 0;
	char* end = nullptr;
	const long long value = std::strtoll(text, &end, 10);
	if (text == end || end == nullptr || *end != '\0' || errno == ERANGE ||
		value < minimum || value > maximum) {
		cerr << "Invalid value for " << option << ": '" << text << "' (expected integer "
			 << minimum << ".." << maximum << ")\n";
		hasErr = true;
		return false;
	}
	result = value;
	return true;
}

bool parseReal(const char* text, const string& option, float minimum,
	float maximum, float& result, bool& hasErr) {
	if (text == nullptr) return false;
	errno = 0;
	char* end = nullptr;
	const float value = std::strtof(text, &end);
	if (text == end || end == nullptr || *end != '\0' || errno == ERANGE ||
		!std::isfinite(value) || value < minimum || value > maximum) {
		cerr << "Invalid value for " << option << ": '" << text << "' (expected number "
			 << minimum << ".." << maximum << ")\n";
		hasErr = true;
		return false;
	}
	result = value;
	return true;
}

bool validateList(const vector<string>& values, const string& option, bool& hasErr) {
	for (const string& value : values) {
		if (value.empty()) {
			cerr << option << " contains an empty comma-separated item\n";
			hasErr = true;
			return false;
		}
	}
	return true;
}

string normalizedPathKey(const string& value) {
	std::error_code error;
	std::filesystem::path normalized = std::filesystem::weakly_canonical(
		std::filesystem::path(value), error);
	if (error) {
		error.clear();
		normalized = std::filesystem::absolute(std::filesystem::path(value), error);
		if (error) normalized = std::filesystem::path(value);
		normalized = normalized.lexically_normal();
	}
	string key = normalized.generic_string();
#if defined(WIN32) || defined(_WIN32) || (defined(__WIN32) && !defined(__CYGWIN__))
	transform(key.begin(), key.end(), key.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
#endif
	return key;
}

} // namespace

void limitedWarning(const string& category, const string& message) {
	static std::unordered_map<string, size_t> warningCounts;
	constexpr size_t displayLimit = 5;
	size_t& count = warningCounts[category];
	++count;
	if (count <= displayLimit) {
		cerr << "Warning: " << message << '\n';
	} else if (count == displayLimit + 1) {
		cerr << "Info: further warnings of type '" << category
			 << "' will not be shown.\n";
	}
}







vector<string> splitByComma(const string& str) {
	vector<string> result;
	size_t start = 0;
	do {
		const size_t comma = str.find(',', start);
		result.push_back(str.substr(start, comma == string::npos ? string::npos : comma - start));
		if (comma == string::npos) break;
		start = comma + 1;
	} while (true);
	return result;
}
vector<int> splitByCommaI(const string& str) {
	vector<int> result;
	for (const string& item : splitByComma(str)) {
		char* end = nullptr;
		errno = 0;
		const long value = std::strtol(item.c_str(), &end, 10);
		if (item.empty() || end == nullptr || *end != '\0' || errno == ERANGE ||
			value < INT_MIN || value > INT_MAX) {
			throw std::invalid_argument("Invalid comma-separated integer: '" + item + "'");
		}
		result.push_back(static_cast<int>(value));
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
	std::unique_ptr<istream> in;
	if (isGZfile(inF)) {
	#ifdef VCF2FNA_HAS_GZIP
		in.reset(new igzstream(inF.c_str(), ios::in));
		//cout << "Straming gzip input on the fly\n";
	#else
		throw std::runtime_error(
			"gzip input is not enabled in this build (automatic on Linux; otherwise build with VCF2FNA_USE_GZSTREAM): " + inF);
#endif
	}
	else { 
		in.reset(new ifstream(inF.c_str()));
	}

	if (!(*in)) { 
		throw std::runtime_error("Could not open file: " + inF); 
	}
	return in.release();
}

ostream* writeGZUZ(const string& outF) {
	std::unique_ptr<ostream> out;
	if (isGZfile(outF)) {
	#ifdef VCF2FNA_HAS_GZIP
		out.reset(new ogzstream(outF.c_str(), ios::out));
		//cout << "Writing gzip'd matrix " << outF << endl;
	#else
		throw std::runtime_error(
			"gzip output is not enabled in this build (automatic on Linux; otherwise build with VCF2FNA_USE_GZSTREAM): " + outF);
#endif
	}
	else { out.reset(new ofstream(outF)); }

	if (!(*out)) { 
		throw std::runtime_error("Could not open file: " + outF); 
	}
	return out.release();
}



void helpMsg() {
	cout << "Usage: vcf2fna -ref REF.fa -inVCF CALLS.vcf -depthF DEPTH.bed -gff GENES.gff [options]\n\n";
	cout << "Required input:\n";
	cout << "  -ref FILE                 reference metagenome assembly\n";
	cout << "  -inVCF FILE[,FILE]        one or two haploid/metagenomic VCF files\n";
	cout << "  -depthF FILE[,FILE]       matching 0-based, half-open bedGraph depth files;\n";
	cout << "                            contig-grouped, sorted, and non-overlapping\n";
	cout << "  -gff FILE                 GFF3 annotation used for gene output/translation;\n";
	cout << "                            functional CDS/ORF/gene features are retained\n\n";
	cout << "Output:\n";
	cout << "  -oCtg FILE                consensus contig FASTA\n";
	cout << "  -oGeneNT FILE             gene nucleotide FASTA\n";
	cout << "  -oGeneAA FILE             translated gene FASTA\n";
	cout << "                            output paths must be distinct from every input/output\n";
	cout << "  -outType CNA              outputs to write: C=contigs, N=gene NT, A=gene AA\n";
	cout << "                            (inferred from supplied output filenames by default)\n";
	cout << "  -keepEmptyCtgs            retain contigs with no resolvable sequence\n";
	cout << "  -keepEmptyGenes           retain genes with no resolvable sequence\n";
	cout << "  .gz input/output is enabled automatically by Linux builds (zlib required).\n\n";
	cout << "Call and coverage filtering:\n";
	cout << "  -vcfFilterPolicy POLICY   technical (default), all, or ignore; GT is ignored\n";
	cout << "                            technical rejects named failures except clearly\n";
	cout << "                            diploid/population-model codes; all rejects every\n";
	cout << "                            non-PASS code; ignore skips the FILTER column\n";
	cout << "  -minCallDepth INT[,INT]   minimum mapped depth per input (default: 5)\n";
	cout << "  -minAltReads INT          reads required for a credible minor allele (default: 4)\n";
	cout << "  -minAltFreq NUM           frequency required for a credible minor allele (default: 0.05)\n";
	cout << "                            indels/complex replacements retain a conservative 0.10 frequency floor\n";
	cout << "  -minCallQual INT          minimum VCF QUAL (default: 20)\n";
	cout << "                            missing QUAL fails whenever this minimum is nonzero\n";
	cout << "  -minCallQualAdaptive NUM  mean-depth QUAL scaling; 0 disables (default: 0)\n";
	cout << "  -depthFilterScale NUM     minimum/mean contig depth scale (default: 0.25)\n";
	cout << "  -maxDepthFilterScale NUM  maximum/mean contig depth scale; 0 disables (default: 0)\n";
	cout << "                            opt in cautiously: metagenomic abundance is highly uneven\n";
	cout << "  -indelRange INT           filter SNPs near an indel (default: 5)\n";
	cout << "  -maxRefMismatches INT     abort after this many REF/FASTA mismatches (default: 10)\n";
	cout << "  -minMQ0F NUM              maximum zero-mapping-quality fraction (default: 0.1)\n";
	cout << "  -minBQBZ NUM              minimum base-quality-bias Z score (default: -3.1)\n";
	cout << "  -minSP NUM                maximum strand-bias score baseline (default: 40)\n";
	cout << "  -skipINDELs               do not apply indels/complex replacements to output sequences\n";
	cout << "  -minorAlleleNoMask        retain reference at passing minor-allele sites\n\n";
	cout << "GFF translation supports NCBI tables 1, 4, 11, and 25; pseudo features are skipped.\n";
	cout << "FASTA/GFF/VCF sequence matching uses the first whitespace-delimited FASTA token.\n\n";
	cout << "Other:\n";
	cout << "  -seqPlatform LIST         informational ill/PB/ONT labels per VCF; no automatic\n";
	cout << "                            platform-specific threshold changes (default: unspecified)\n";
	cout << "  -tmp DIR                  reserved temporary directory setting\n";
	cout << "  -t INT                    reserved; execution is currently single-threaded\n";
	cout << "  -debug1                   print VCF debug information\n";
	cout << "  -h, --help                print this help message\n";
	cout << "  -v, --version             print version information\n";
}

void versionMsg() {
	cout << "vcf2fna v" << vcf2fnaVERSION << endl;
	cout << "compiled on: " << __DATE__ << " at " << __TIME__ << endl;
	cout << "compiled with C++ v" << __cplusplus << endl;
	exit(0);
}

void options::announce() {
	cout << "------------------------------------ " << std::endl;
	cout << "vcf2fna v" << vcf2fnaVERSION << endl;
	print_details();
	cout << "------------------------------------ " << std::endl;
}

options::options(int argc, char** argv) : refFasta(""),
outfna(""), outGeneNT(""), outGeneAA(""),
tmp(""), inVCF(), depthF(), gffFile(""), outputTypes(""),
seqPlatform(), vcfFilterPolicy("technical"),
threads(1), indelRange(5), maxRefMismatches(10),
minDepthPar(), minAltReads(4), minAltFreq(0.05f), minCallQual(20),
minMQ0F(0.1f), minBQBZ(-3.1f), minSP(40.f),
addHDTags(true), skipEmptyContigs(true), skipEmptyGenes(true),
debug1(false), reportINDELs(true), maskMinorAllele(true),
minCallQualAdaptive(0.f), depthFilterScale(0.25f), maxDepthFilterScale(0.f)
{
	bool hasErr = false;
	if (argc <= 1 || argv == nullptr) {
		cerr << "Not enough arguments. Use \"vcf2fna -h\" for help.\n";
		exit(23);
	}

	for (int i = 1; i < argc; ++i) {
		const string option(argv[i]);
		if (option == "-h" || option == "--help") {
			helpMsg();
			exit(0);
		}
		if (option == "-v" || option == "--version") {
			versionMsg();
		}

		auto readString = [&](string& destination) {
			const char* value = takeOptionValue(i, argc, argv, hasErr);
			if (value != nullptr) destination = value;
		};
		auto readStringList = [&](vector<string>& destination) {
			const char* value = takeOptionValue(i, argc, argv, hasErr);
			if (value != nullptr) destination = splitByComma(value);
		};
		auto readInteger = [&](long long minimum, long long maximum, long long& destination) {
			const char* value = takeOptionValue(i, argc, argv, hasErr);
			parseInteger(value, option, minimum, maximum, destination, hasErr);
		};
		auto readReal = [&](float minimum, float maximum, float& destination) {
			const char* value = takeOptionValue(i, argc, argv, hasErr);
			parseReal(value, option, minimum, maximum, destination, hasErr);
		};

		if (option == "-ref") {
			readString(refFasta);
		} else if (option == "-oCtg") {
			readString(outfna);
		} else if (option == "-oGeneNT") {
			readString(outGeneNT);
		} else if (option == "-oGeneAA") {
			readString(outGeneAA);
		} else if (option == "-tmp") {
			readString(tmp);
		} else if (option == "-outType") {
			readString(outputTypes);
		} else if (option == "-gff") {
			readString(gffFile);
		} else if (option == "-inVCF") {
			readStringList(inVCF);
		} else if (option == "-seqPlatform") {
			readStringList(seqPlatform);
		} else if (option == "-depthF") {
			readStringList(depthF);
		} else if (option == "-vcfFilterPolicy") {
			readString(vcfFilterPolicy);
		} else if (option == "-t") {
			long long value = threads;
			readInteger(1, std::numeric_limits<uint>::max(), value);
			threads = static_cast<uint>(value);
		} else if (option == "-minCallDepth") {
			const char* value = takeOptionValue(i, argc, argv, hasErr);
			if (value != nullptr) {
				const vector<string> values = splitByComma(value);
				validateList(values, option, hasErr);
				minDepthPar.clear();
				for (const string& item : values) {
					long long parsed = 0;
					parseInteger(item.c_str(), option, 0, INT_MAX, parsed, hasErr);
					minDepthPar.push_back(static_cast<int>(parsed));
				}
			}
		} else if (option == "-minAltReads") {
			long long value = minAltReads;
			readInteger(1, INT_MAX, value);
			minAltReads = static_cast<int>(value);
		} else if (option == "-minAltFreq") {
			readReal(0.f, 1.f, minAltFreq);
		} else if (option == "-minCallQual") {
			long long value = minCallQual;
			readInteger(0, INT_MAX, value);
			minCallQual = static_cast<int>(value);
		} else if (option == "-indelRange") {
			long long value = indelRange;
			readInteger(0, INT_MAX, value);
			indelRange = static_cast<uint>(value);
		} else if (option == "-maxRefMismatches") {
			long long value = maxRefMismatches;
			readInteger(0, INT_MAX, value);
			maxRefMismatches = static_cast<int>(value);
		} else if (option == "-minMQ0F") {
			readReal(0.f, 1.f, minMQ0F);
		} else if (option == "-minBQBZ") {
			readReal(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), minBQBZ);
		} else if (option == "-minSP") {
			readReal(0.f, std::numeric_limits<float>::max(), minSP);
		} else if (option == "-minCallQualAdaptive") {
			readReal(0.f, std::numeric_limits<float>::max(), minCallQualAdaptive);
		} else if (option == "-depthFilterScale") {
			readReal(0.f, std::numeric_limits<float>::max(), depthFilterScale);
		} else if (option == "-maxDepthFilterScale") {
			readReal(0.f, std::numeric_limits<float>::max(), maxDepthFilterScale);
		} else if (option == "-keepEmptyCtgs") {
			skipEmptyContigs = false;
		} else if (option == "-keepEmptyGenes") {
			skipEmptyGenes = false;
		} else if (option == "-skipINDELs") {
			reportINDELs = false;
		} else if (option == "-minorAlleleNoMask") {
			maskMinorAllele = false;
		} else if (option == "-debug1") {
			debug1 = true;
			cerr << "\nDEBUG1 active\n";
		} else {
			cerr << "Unknown option: " << option << "\n";
			hasErr = true;
		}
	}

	validateList(inVCF, "-inVCF", hasErr);
	validateList(depthF, "-depthF", hasErr);
	if (!seqPlatform.empty()) validateList(seqPlatform, "-seqPlatform", hasErr);

	if (seqPlatform.empty()) seqPlatform.assign(inVCF.size(), "unspecified");
	if (minDepthPar.empty() && (inVCF.size() == 1 || inVCF.size() == 2)) {
		minDepthPar.assign(inVCF.size(), 5);
	}

	if (refFasta.empty()) {
		cerr << "-ref must be specified\n";
		hasErr = true;
	}
	if (inVCF.empty()) {
		cerr << "-inVCF must be specified\n";
		hasErr = true;
	} else if (inVCF.size() > 2) {
		cerr << "-inVCF accepts at most two files\n";
		hasErr = true;
	}
	if (depthF.empty()) {
		cerr << "-depthF must be specified\n";
		hasErr = true;
	}
	if (gffFile.empty()) {
		cerr << "-gff must be specified\n";
		hasErr = true;
	}
	if (depthF.size() != inVCF.size()) {
		cerr << "-depthF and -inVCF must contain the same number of files\n";
		hasErr = true;
	}
	if (seqPlatform.size() != inVCF.size()) {
		cerr << "-seqPlatform must contain one value per VCF\n";
		hasErr = true;
	}
	if (minDepthPar.size() != depthF.size()) {
		cerr << "-minCallDepth must contain one value per depth file\n";
		hasErr = true;
	}
	for (const string& platform : seqPlatform) {
		if (platform != "ill" && platform != "PB" && platform != "ONT" &&
			platform != "unspecified") {
			cerr << "Unsupported -seqPlatform value: '" << platform
				 << "' (expected ill, PB, ONT, or unspecified)\n";
			hasErr = true;
		}
	}
	if (vcfFilterPolicy != "technical" && vcfFilterPolicy != "all" &&
		vcfFilterPolicy != "ignore") {
		cerr << "Unsupported -vcfFilterPolicy value: '" << vcfFilterPolicy
			 << "' (expected technical, all, or ignore)\n";
		hasErr = true;
	}

	if (outputTypes.empty()) {
		if (!outfna.empty()) outputTypes += "C";
		if (!outGeneNT.empty()) outputTypes += "N";
		if (!outGeneAA.empty()) outputTypes += "A";
	} else {
		set<char> seen;
		for (char outputType : outputTypes) {
			if (outputType != 'C' && outputType != 'N' && outputType != 'A') {
				cerr << "Invalid -outType character: '" << outputType << "' (expected C, N, or A)\n";
				hasErr = true;
			} else if (!seen.insert(outputType).second) {
				cerr << "Duplicate -outType character: '" << outputType << "'\n";
				hasErr = true;
			}
		}
	}
	if (outputTypes.find('C') != string::npos && outfna.empty()) {
		cerr << "-outType C requires -oCtg\n";
		hasErr = true;
	}
	if (outputTypes.find('N') != string::npos && outGeneNT.empty()) {
		cerr << "-outType N requires -oGeneNT\n";
		hasErr = true;
	}
	if (outputTypes.find('A') != string::npos && outGeneAA.empty()) {
		cerr << "-outType A requires -oGeneAA\n";
		hasErr = true;
	}

	set<string> inputPaths;
	auto rememberInput = [&inputPaths](const string& path) {
		if (!path.empty()) inputPaths.insert(normalizedPathKey(path));
	};
	rememberInput(refFasta);
	rememberInput(gffFile);
	for (const string& path : inVCF) rememberInput(path);
	for (const string& path : depthF) rememberInput(path);
	map<string, string> outputOwners;
	auto validateOutputPath = [&](const string& path, const string& optionName) {
		if (path.empty()) return;
		const string key = normalizedPathKey(path);
		if (inputPaths.find(key) != inputPaths.end()) {
			cerr << optionName << " must not overwrite an input file: " << path << "\n";
			hasErr = true;
		}
		auto existing = outputOwners.find(key);
		if (existing != outputOwners.end()) {
			cerr << optionName << " and " << existing->second
				 << " must use different output paths: " << path << "\n";
			hasErr = true;
		} else {
			outputOwners[key] = optionName;
		}
	};
	validateOutputPath(outfna, "-oCtg");
	validateOutputPath(outGeneNT, "-oGeneNT");
	validateOutputPath(outGeneAA, "-oGeneAA");

	if (hasErr) {
		cerr << "Error in option parsing. Use \"vcf2fna -h\" for full help.\n";
		exit(98);
	}

	announce();
}

void options::print_details() {

	// print run mode:
	cout << "ref assembly, gff file: " << refFasta <<", "<< gffFile << std::endl;
	cout << "seqPlatform labels (informational): " << combCommas(seqPlatform) << std::endl;
	cout << "input vcf:        " << combCommas(inVCF) << std::endl;
	cout << "depth file:       " << combCommas(depthF) << std::endl;
	cout << "VCF FILTER policy: " << vcfFilterPolicy << std::endl;
	cout << "tmp dir:          " << tmp << std::endl;
	cout << "output types: \"" << outputTypes << "\" ::";
	cout << outfna << " " << outGeneNT << " " << outGeneAA << std::endl;
	if (outfna.empty() && outGeneNT.empty() && outGeneAA.empty()){
		cout << "  (no output files specified - theoretical statistics reported)" << std::endl;
	}
	cout << "threads (reserved; single-threaded): " << threads << std::endl;
	cout << "skipEmptyContigs, reportINDELs: " << skipEmptyContigs << ", " << reportINDELs << std::endl;
	cout << "minCallDepth, minCallQual, maskMinorAllele: " << combCommas(minDepthPar)
		 << ", " << minCallQual << ", " << maskMinorAllele << std::endl;
	cout << "minAltReads, minAltFreq: " << minAltReads << ", " << minAltFreq << std::endl;
	cout << "depthFilterScale, maxDepthFilterScale, minCallQualAdaptive: "
		 << depthFilterScale << ", " << maxDepthFilterScale << ", " << minCallQualAdaptive << std::endl;
	cout << "minMQ0F, minBQBZ, minSP, maxRefMismatches: " << minMQ0F << ", "
		 << minBQBZ << ", " << minSP << ", " << maxRefMismatches << std::endl;

}
