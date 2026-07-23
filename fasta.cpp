#include "fasta.h"
#include "vcf.h"

#include <cctype>
#include <climits>
#include <limits>
#include <set>
#include <unordered_map>


namespace {

const char* STANDARD_AAS =
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

string trimCopy(const string& value) {
	size_t first = value.find_first_not_of(" \t\r\n");
	if (first == string::npos) {
		return "";
	}
	size_t last = value.find_last_not_of(" \t\r\n");
	return value.substr(first, last - first + 1);
}

string canonicalSequenceId(const string& header) {
	string trimmed = trimCopy(header);
	size_t delimiter = trimmed.find_first_of(" \t");
	return delimiter == string::npos ? trimmed : trimmed.substr(0, delimiter);
}

int hexDigit(char value) {
	if (value >= '0' && value <= '9') return value - '0';
	value = static_cast<char>(std::tolower(static_cast<unsigned char>(value)));
	if (value >= 'a' && value <= 'f') return value - 'a' + 10;
	return -1;
}

string percentDecode(const string& value) {
	string decoded;
	decoded.reserve(value.size());
	for (size_t i = 0; i < value.size(); ++i) {
		if (value[i] == '%' && i + 2 < value.size()) {
			int high = hexDigit(value[i + 1]);
			int low = hexDigit(value[i + 2]);
			if (high >= 0 && low >= 0) {
				decoded.push_back(static_cast<char>((high << 4) | low));
				i += 2;
				continue;
			}
		}
		decoded.push_back(value[i]);
	}
	return decoded;
}

bool parseLongLongStrict(const string& text, long long& value) {
	if (text.empty()) return false;
	size_t used = 0;
	try {
		value = stoll(text, &used);
	} catch (const std::exception&) {
		return false;
	}
	return used == text.size();
}

bool parseIntStrict(const string& text, int& value) {
	long long parsed = 0;
	if (!parseLongLongStrict(text, parsed) || parsed < INT_MIN || parsed > INT_MAX) return false;
	value = static_cast<int>(parsed);
	return true;
}

bool isN(char base) {
	return std::toupper(static_cast<unsigned char>(base)) == 'N';
}

int64_t countDeterminedBases(const string& sequence) {
	return static_cast<int64_t>(std::count_if(sequence.begin(), sequence.end(),
		[](char base) { return !isN(base); }));
}

void addSequenceLength(int64_t& total, size_t length) {
	if (length > static_cast<size_t>(std::numeric_limits<int64_t>::max()) ||
		total > std::numeric_limits<int64_t>::max() - static_cast<int64_t>(length)) {
		throw std::overflow_error("Reference assembly length exceeds the supported 64-bit range");
	}
	total += static_cast<int64_t>(length);
}

int codonBaseIndex(char base) {
	switch (std::toupper(static_cast<unsigned char>(base))) {
	case 'T': return 0;
	case 'C': return 1;
	case 'A': return 2;
	case 'G': return 3;
	default: return -1;
	}
}

bool supportedTranslationTable(int table) {
	return table == 1 || table == 4 || table == 11 || table == 25;
}

bool isInitiatorCodon(const string& codon, int table) {
	if (codon.size() != 3) return false;
	string upper = codon;
	std::transform(upper.begin(), upper.end(), upper.begin(),
		[](char base) { return static_cast<char>(std::toupper(static_cast<unsigned char>(base))); });
	if (table == 1) {
		return upper == "ATG" || upper == "TTG" || upper == "CTG";
	}
	if (table == 4) {
		return upper == "TTA" || upper == "TTG" || upper == "CTG" ||
			upper == "ATT" || upper == "ATC" || upper == "ATA" ||
			upper == "ATG" || upper == "GTG";
	}
	if (table == 11) {
		return upper == "TTG" || upper == "CTG" || upper == "ATT" ||
			upper == "ATC" || upper == "ATA" || upper == "ATG" ||
			upper == "GTG";
	}
	if (table == 25) {
		return upper == "TTG" || upper == "ATG" || upper == "GTG";
	}
	return false;
}

char translateCodon(char p1, char p2, char p3, int table) {
	if (!supportedTranslationTable(table)) {
		throw std::runtime_error("Unsupported NCBI translation table: " + to_string(table) +
			" (supported bacterial-relevant tables: 1, 4, 11, 25)");
	}
	int b1 = codonBaseIndex(p1);
	int b2 = codonBaseIndex(p2);
	int b3 = codonBaseIndex(p3);
	if (b1 < 0 || b2 < 0 || b3 < 0) {
		if (p1 == '-' || p2 == '-' || p3 == '-') return '-';
		if (p1 == '.' || p2 == '.' || p3 == '.') return '.';
		return 'X';
	}
	int index = b1 * 16 + b2 * 4 + b3;
	if (index == 14 && table == 4) return 'W'; // TGA
	if (index == 14 && table == 25) return 'G'; // TGA
	return STANDARD_AAS[index];
}

vector<string> splitTabs(const string& line) {
	vector<string> fields;
	size_t start = 0;
	while (true) {
		size_t tab = line.find('\t', start);
		if (tab == string::npos) {
			fields.push_back(line.substr(start));
			break;
		}
		fields.push_back(line.substr(start, tab - start));
		start = tab + 1;
	}
	return fields;
}

unordered_map<string, string> parseGffAttributes(const string& field) {
	unordered_map<string, string> attributes;
	stringstream input(field);
	string item;
	while (getline(input, item, ';')) {
		item = trimCopy(item);
		if (item.empty()) continue;
		size_t delimiter = item.find('=');
		if (delimiter == string::npos) {
			// Tolerate common GTF-style key "value" attributes.
			delimiter = item.find_first_of(" \t");
		}
		if (delimiter == string::npos) continue;
		string key = trimCopy(item.substr(0, delimiter));
		string value = trimCopy(item.substr(delimiter + 1));
		if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
			value = value.substr(1, value.size() - 2);
		}
		attributes[key] = percentDecode(value);
	}
	return attributes;
}

int translationTableFromText(const string& text, int fallback) {
	size_t key = text.find("transl_table=");
	if (key == string::npos) return fallback;
	key += 13;
	size_t end = key;
	while (end < text.size() && std::isdigit(static_cast<unsigned char>(text[end]))) ++end;
	if (end == key) return fallback;
	return stoi(text.substr(key, end - key));
}

string lowerCopy(string value) {
	std::transform(value.begin(), value.end(), value.begin(),
		[](char c) { return static_cast<char>(std::tolower(static_cast<unsigned char>(c))); });
	return value;
}

bool isCodingFeatureType(const string& type) {
	string lower = lowerCopy(type);
	return lower == "cds" || lower == "orf" || lower == "open_reading_frame" ||
		lower == "protein_coding_gene" || lower == "so:0000316" || lower == "so:0000236";
}

bool isCdsFeatureType(const string& type) {
	const string lower = lowerCopy(type);
	return lower == "cds" || lower == "so:0000316";
}

bool isGeneFeatureType(const string& type) {
	string lower = lowerCopy(type);
	return lower == "gene" || lower == "so:0000704";
}

bool isNonCodingGene(const unordered_map<string, string>& attributes) {
	auto pseudo = attributes.find("pseudo");
	if (pseudo != attributes.end()) {
		const string value = lowerCopy(pseudo->second);
		if (value.empty() || value == "true" || value == "yes" || value == "1") return true;
	}
	if (attributes.find("pseudogene") != attributes.end()) return true;
	static const char* keys[] = { "gene_biotype", "gene_type", "biotype" };
	for (const char* key : keys) {
		auto found = attributes.find(key);
		if (found == attributes.end()) continue;
		string value = lowerCopy(found->second);
		if (value.find("rrna") != string::npos || value.find("trna") != string::npos ||
			value.find("ncrna") != string::npos || value.find("noncoding") != string::npos ||
			value.find("non_coding") != string::npos || value.find("pseudo") != string::npos) {
			return true;
		}
		if (value.find("protein_coding") != string::npos) return false;
	}
	return false;
}

} // namespace


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
	geneID(id), geneStart(sta), geneEnd(end), geneLength(end - sta + 1),
	depthReferenceLength(end - sta + 1),
	segments(1, make_pair(sta, end)), geneStrand(true), translationTable(11), phase(0),
	numOnContig(-1) , accumDepth(0)
{
	if (geneStart < 0 || geneEnd < geneStart) {
		throw std::invalid_argument("Invalid gene coordinates: " + to_string(sta) + "-" + to_string(end));
	}
}

gene::gene(gene* GG) : geneID(GG->geneID), geneStart(GG->geneStart), geneEnd(GG->geneEnd), geneLength(GG->geneLength),
depthReferenceLength(GG->depthReferenceLength),
segments(GG->segments), geneStrand(GG->geneStrand), type(GG->type), translationTable(GG->translationTable), phase(GG->phase),
numOnContig(GG->numOnContig), partial(GG->partial), accumDepth(GG->accumDepth)
{
}

void gene::setSegments(const vector<pair<int, int> >& newSegments) {
	if (newSegments.empty()) return;
	segments = newSegments;
	for (const pair<int, int>& segment : segments) {
		if (segment.first < 0 || segment.second < segment.first) {
			throw std::invalid_argument("Invalid GFF segment coordinates: " +
				to_string(segment.first) + "-" + to_string(segment.second));
		}
	}
	std::sort(segments.begin(), segments.end());
	for (size_t i = 1; i < segments.size(); ++i) {
		if (segments[i].first <= segments[i - 1].second) {
			throw std::invalid_argument("Overlapping or duplicate GFF segments: " +
				to_string(segments[i - 1].first) + "-" + to_string(segments[i - 1].second) +
				" and " + to_string(segments[i].first) + "-" + to_string(segments[i].second));
		}
	}
	recalculateGeometry();
	depthReferenceLength = geneLength;
}

void gene::recalculateGeometry() {
	if (segments.empty()) {
		geneStart = geneEnd = geneLength = 0;
		return;
	}
	geneStart = segments.front().first;
	geneEnd = segments.front().second;
	geneLength = 0;
	for (const pair<int, int>& segment : segments) {
		geneStart = std::min(geneStart, segment.first);
		geneEnd = std::max(geneEnd, segment.second);
		geneLength += segment.second - segment.first + 1;
	}
}

bool gene::fivePrimePartial() const {
	if (partial.size() < 2) return false;
	// Prodigal records partialness at the left and right genomic boundaries.
	return geneStrand ? partial[0] == '1' : partial[1] == '1';
}

int gene::positionOffset(int genomicPos) const {
	int offset = 0;
	if (geneStrand) {
		for (const pair<int, int>& segment : segments) {
			if (genomicPos >= segment.first && genomicPos <= segment.second) {
				return offset + genomicPos - segment.first;
			}
			offset += segment.second - segment.first + 1;
		}
	} else {
		for (auto segment = segments.rbegin(); segment != segments.rend(); ++segment) {
			if (genomicPos >= segment->first && genomicPos <= segment->second) {
				return offset + segment->second - genomicPos;
			}
			offset += segment->second - segment->first + 1;
		}
	}
	return -1;
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
	for (const pair<int, int>& segment : segments) {
		if (segment.first < 0 || segment.second < segment.first ||
			static_cast<size_t>(segment.second) >= seq.length()) {
			throw std::runtime_error("Gene " + geneID + " coordinates " +
				to_string(segment.first) + "-" + to_string(segment.second) +
				" are outside sequence length " + to_string(seq.length()));
		}
		ret += seq.substr(static_cast<size_t>(segment.first),
			static_cast<size_t>(segment.second - segment.first + 1));
	}
	if (!geneStrand){
		reverseComplement(ret);
	} 

	return ret;
}
string gene::geneAA(const string& seq) {
	string ret("");
	if (!supportedTranslationTable(translationTable)) {
		throw std::runtime_error("Gene " + geneID + " requests unsupported NCBI translation table " +
			to_string(translationTable) + " (supported: 1, 4, 11, 25)");
	}
	if (phase < 0 || phase > 2) {
		throw std::runtime_error("Invalid GFF phase " + to_string(phase) + " for gene " + geneID);
	}

	for (size_t i = static_cast<size_t>(phase); i + 2 < seq.length(); i += 3) {
		char P1 = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[i])));
		char P2 = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[i + 1])));
		char P3 = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[i + 2])));
		char protein = translateCodon(P1, P2, P3, translationTable);
		if (i == static_cast<size_t>(phase) && phase == 0 && !fivePrimePartial() &&
			isInitiatorCodon(string() + P1 + P2 + P3, translationTable)) {
			protein = 'M';
		}
		ret += protein ;
		if (protein == '*' ) break;
	}

	return ret;
}



//***************************************************************

fasta::fasta(string s, string h) :
	seq(std::move(s)), mutSeq(""), mutSeqDone(false), header(trimCopy(h)),
	sequenceId(canonicalSequenceId(h)), seqUse(seq.length(), false),
	SNPsCnt(0), UnctCnt(0), SNPsPos(0), SNPfreqs(0), conflictCnt(0),
	INDELcnt(0), INDELpos(0), INDELfreq(0), depthAccum(0), depthAccumBySource(),
	geneCol(new geneCollection())
{
	if (sequenceId.empty()) {
		throw std::invalid_argument("FASTA header has no sequence identifier: " + h);
	}
	if (seq.length() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::overflow_error("Reference sequence exceeds the supported coordinate range: " + sequenceId);
	}
}

bool fasta::referenceAlleleMatches(int pos, const string& ref, string& reason) const {
	if (pos < 0) {
		reason = "negative zero-based coordinate " + to_string(pos);
		return false;
	}
	if (ref.empty() || ref == "." || ref == "*") {
		reason = "missing or unsupported REF allele '" + ref + "'";
		return false;
	}
	size_t start = static_cast<size_t>(pos);
	if (start >= seq.length() || ref.length() > seq.length() - start) {
		reason = "REF span " + to_string(pos) + "-" +
			to_string(pos + static_cast<int>(ref.length()) - 1) +
			" is outside sequence length " + to_string(seq.length());
		return false;
	}
	for (size_t i = 0; i < ref.length(); ++i) {
		char observed = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[start + i])));
		char expected = static_cast<char>(std::toupper(static_cast<unsigned char>(ref[i])));
		// An N in either source denotes an unresolved reference base and is treated as a wildcard.
		if (observed != 'N' && expected != 'N' && observed != expected) {
			reason = "reference has '" + string(1, observed) + "' but VCF REF has '" +
				string(1, expected) + "' at one-based position " + to_string(pos + static_cast<int>(i) + 1);
			return false;
		}
	}
	return true;
}

bool fasta::validateVariantReference(VCFmem* vx, VariantStats* stats, const options* opts) const {
	if (vx == nullptr || stats == nullptr || opts == nullptr) {
		throw std::invalid_argument("Null argument supplied while validating a VCF reference allele");
	}
	string mismatchReason;
	if (referenceAlleleMatches(vx->getPos(), vx->getRef(), mismatchReason)) return true;

	++stats->refMismatch;
	limitedWarning("reference/VCF mismatch",
		"skipping variant on " + sequenceId + " at one-based position " +
		to_string(vx->getPos() + 1) + ": " + mismatchReason + " (mismatch " +
		to_string(stats->refMismatch) + "/" + to_string(opts->maxRefMismatches) + ")");
	if (stats->refMismatch > opts->maxRefMismatches) {
		throw std::runtime_error("Reference/VCF mismatch limit exceeded (" +
			to_string(stats->refMismatch) + " > " + to_string(opts->maxRefMismatches) +
			"). Check that FASTA, VCF, and coordinate conventions match.");
	}
	return false;
}

//cF->ntVariant(posN, ref, "N", -1.f); 
void fasta::ntVariant(VCFmem* vx, VariantStats* Vstats, VCFfilterStats* filStats,
	options* opt) {
	if (vx == nullptr || Vstats == nullptr || filStats == nullptr || opt == nullptr) {
		throw std::invalid_argument("Null argument supplied while applying a VCF variant");
	}

	if (vx->isINDEL()) {
		Vstats->indelCNT++;
	} else if (vx->isSNP()) {
		Vstats->snpCNT++;
	} else {
		return;
	}

	if (!validateVariantReference(vx, Vstats, opt)) return;
	if (vx->conflicted()) {
		++conflictCnt;
		++Vstats->conflictCnt;
		if (vx->isINDEL()) {
			++Vstats->unsrINDEL;
			maskUncertainSpan(vx->getPos(), vx->getUncertainRefLength());
		} else {
			++Vstats->unsrSNP;
			SNP(vx->getPos(), vx->getRef(), "N", -1.f);
		}
		return;
	}

	//don't include if position is not deep enough in first place..
	if (!seqUse[static_cast<size_t>(vx->getPos())]) {
		Vstats->SNPlowCov++;
		return; 
	}
	if (vx->isINDEL()) {
		if (vx->majorAllele()) {
			if (vx->filtered()) {
				++Vstats->unsrINDEL;
				maskUncertainSpan(vx->getPos(), vx->getRef().length());
			}
			else if (opt->reportINDELs) {
				INDEL(vx->getPos(), vx->getRef(), vx->getAlt(), vx->getFreq());
				Vstats->indelUsed++;
			}
		} else if (!vx->filtered() && opt->maskMinorAllele) {
			// A supported but non-major indel leaves the haploid consensus
			// ambiguous across its reference span, just as a minor SNP does.
			++Vstats->unsrINDEL;
			maskUncertainSpan(vx->getPos(), vx->getRef().length());
		}
		if (vx->filtered()) {
			Vstats->indelFILT++;
		}
		return; //skip indels for now
	}

	/// ************************* SNP starts here *************************
	//implement SNPs (maybe)
	if (vx->majorAllele()) {
		if (vx->filtered()) {
			Vstats->majorAlleleFilt++;//snpFILT
			Vstats->unsrSNP++;
			filStats->mumaFilt->addMut(vx->getRef(), vx->getAlt());
			SNP(vx->getPos(), vx->getRef(), "N", -1.f);// vx->getFreq());
		} else {
			Vstats->majorAllele++;//SNPused
			filStats->muma->addMut(vx->getRef(), vx->getAlt());
			SNP(vx->getPos(), vx->getRef(), vx->getAlt(), vx->getFreq());
		}
	}
	else { //minor alleles.. for stats only, not applied to sequence
		if (vx->filtered()) {
			filStats->mumaLowFreqFilt->addMut(vx->getRef(), vx->getAlt());
			Vstats->minorAlleleFilt++;
		} else {
			filStats->mumaLowFreq->addMut(vx->getRef(), vx->getAlt());
			Vstats->minorAllele++;
			if (opt->maskMinorAllele) {
				SNP(vx->getPos(), vx->getRef(), "N", -1.f);//
			}
		}
	}
	/*if (vx->filtered()) {
		Vstats->unsrSNP++;
	}*/


}

//structural variant (indel) of some sort
void fasta::INDEL(int pos, string r, string a, float freq)
{
	string reason;
	if (!referenceAlleleMatches(pos, r, reason)) {
		throw std::logic_error("Attempted to register an unvalidated INDEL on " + sequenceId + ": " + reason);
	}
	if (a.empty()) {
		throw std::invalid_argument("Empty ALT allele for INDEL on " + sequenceId + ":" + to_string(pos + 1));
	}
	mutSeqDone = false;
	INDELcnt++;
	INDELpos.push_back(pos);
	INDELfreq.push_back(freq);
	INDELref.push_back(r);
	INDELalt.push_back(a);//for now only catalogue.. needs to be applied later..
}

void fasta::maskUncertainSpan(int pos, size_t refLength) {
	if (pos < 0 || refLength == 0 || static_cast<size_t>(pos) >= seq.length() ||
		refLength > seq.length() - static_cast<size_t>(pos)) {
		throw std::out_of_range("Cannot mask uncertain variant span on " + sequenceId + ":" +
			to_string(pos + 1));
	}
	for (size_t i = 0; i < refLength; ++i) seq[static_cast<size_t>(pos) + i] = 'N';
	++UnctCnt;
	mutSeqDone = false;
}

//single nucleotide polymorphism
void fasta::SNP(int pos, string r, string a, float freq)
{ 
	string reason;
	if (!referenceAlleleMatches(pos, r, reason)) {
		throw std::logic_error("Attempted to register an unvalidated SNP on " + sequenceId + ": " + reason);
	}
	if (r.size() != 1 || a.size() != 1) {
		throw std::invalid_argument("SNP alleles must each contain exactly one base on " +
			sequenceId + ":" + to_string(pos + 1));
	}
	mutSeqDone = false;
	if (freq < 0.f || isN(a[0])) {
		//maskSeq(pos, pos + 1); //replace with N
		UnctCnt++; //return; 
	} else {
		SNPsCnt++;
		SNPsPos.push_back(pos);
		SNPfreqs.push_back(freq);
	}
	seq[static_cast<size_t>(pos)] = a[0];
}


string fasta::getMutatedHeader(bool hdTags, const string& renderedSeq) {
	string ret = ">" + header;
	if (!hdTags) {
		return ret;
	}
	//add tags for MGTK
	ret += " COV=" + to_string(countDeterminedBases(renderedSeq)) + " REPL=" + to_string(SNPsCnt);
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
	vector<int> sumFreqs(11, 0);
	bool firstVariant = true;
	for (it = SNPsPos.begin(), it2 = SNPfreqs.begin();
		it != SNPsPos.end() && it2 != SNPfreqs.end(); ++it, ++it2) {
		int offset = positionOffset(*it);
		if (offset < 0) continue;
		//fequence vector
		if (!firstVariant) freqV += ",";
		freqV += to_string_with_precision(*it2, 2);
		int bin = static_cast<int>((*it2 * 10.f) + 0.5f);
		bin = std::max(0, std::min(10, bin));
		sumFreqs[static_cast<size_t>(bin)]++;
		//position vector
		if (!firstVariant) hd += ",";
		hd += to_string(offset);
		firstVariant = false;
	}

	int sumMidFreqs(sumFreqs[2] + sumFreqs[3] + sumFreqs[4] + sumFreqs[5] + sumFreqs[6] + sumFreqs[7]);
	int sumFixFreqs(sumFreqs[0] + sumFreqs[1] + sumFreqs[8] + sumFreqs[9] + sumFreqs[10]);
	float oCSP(0.f);
	if (sumMidFreqs + sumFixFreqs > 0) {
		oCSP = float(sumMidFreqs) / float(sumMidFreqs + sumFixFreqs);
	}
	float CSP(0.f);
	if (nonNs == -1) {
		nonNs = static_cast<int>(countDeterminedBases(seq));
	}
	if (nonNs>0) {
		CSP = float(sumMidFreqs) / float(nonNs);
	}
	
	hd += "] F = [" + freqV + "]";
	hd+= " oCSP=" + to_string_with_precision(oCSP,3) +" CSP="+to_string_with_precision(CSP, 3); 
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
	if (INDELpos.empty() || !opts->reportINDELs) {
		mutSeqDone = true;
		return mutSeq;
	}
	//INDELs
	assert(INDELpos.size() == INDELalt.size());
	assert(INDELpos.size() == INDELref.size());
	assert(INDELpos.size() == INDELfreq.size());
	geneCol->prepMuts();

	vector<size_t> order(INDELpos.size());
	for (size_t i = 0; i < order.size(); ++i) order[i] = i;
	std::sort(order.begin(), order.end(), [this](size_t lhs, size_t rhs) {
		return INDELpos[lhs] > INDELpos[rhs];
	});
	for (size_t index : order) { // right-to-left keeps earlier reference coordinates stable
		int pos = INDELpos[index];
		const string& ref = INDELref[index];
		const string& alt = INDELalt[index];
		if (alt == "N" || alt == ".") {
			continue; //skip uncertain indels
		}
		if ((!alt.empty() && alt.front() == '<') || alt == "*") {
			limitedWarning("symbolic indel",
				"symbolic indel " + alt + " at position " + to_string(pos) +
				"; skipping this variant");
			continue; //skip symbolic indels
		}
		if (pos < 0 || static_cast<size_t>(pos) > mutSeq.length() ||
			ref.length() > mutSeq.length() - static_cast<size_t>(pos)) {
			throw std::out_of_range("INDEL on " + sequenceId + " no longer fits the consensus at " +
				to_string(pos + 1));
		}
		if (alt.length() > ref.length() &&
			alt.length() - ref.length() >
				static_cast<size_t>(std::numeric_limits<int>::max()) - mutSeq.length()) {
			throw std::overflow_error("Sequence replacement would exceed the supported coordinate range on " +
				sequenceId + ":" + to_string(pos + 1));
		}
		mutSeq.replace(static_cast<size_t>(pos), ref.length(), alt);
		// POS tags describe the emitted consensus coordinates. Keep registered
		// SNP positions aligned when an upstream indel changes sequence length.
		size_t sharedPrefix = 0;
		while (sharedPrefix < ref.size() && sharedPrefix < alt.size() &&
			std::toupper(static_cast<unsigned char>(ref[sharedPrefix])) ==
			std::toupper(static_cast<unsigned char>(alt[sharedPrefix]))) ++sharedPrefix;
		size_t sharedSuffix = 0;
		while (sharedSuffix < ref.size() - sharedPrefix && sharedSuffix < alt.size() - sharedPrefix &&
			std::toupper(static_cast<unsigned char>(ref[ref.size() - 1 - sharedSuffix])) ==
			std::toupper(static_cast<unsigned char>(alt[alt.size() - 1 - sharedSuffix]))) ++sharedSuffix;
		const int editStart = pos + static_cast<int>(sharedPrefix);
		const int refCoreLength = static_cast<int>(ref.size() - sharedPrefix - sharedSuffix);
		const int altCoreLength = static_cast<int>(alt.size() - sharedPrefix - sharedSuffix);
		const int editEnd = editStart + refCoreLength;
		const int delta = altCoreLength - refCoreLength;
		for (int& snpPosition : SNPsPos) {
			if (snpPosition >= editEnd) snpPosition += delta;
		}
		//now apply to every gene coordinate
		geneCol->correctCoords(pos, ref, alt);
	}
	mutSeqDone = true;



	return mutSeq; 
}  //


int fasta::maskSeq(int start, int end,  bool repl) {
	//mask sequence from start to end (exclusive)
	assert(seq.length() == seqUse.size());

	if (start < 0 || end < 0 || static_cast<size_t>(end) > seq.length() || start >= end) {
		throw std::out_of_range("Invalid depth-mask range: " + to_string(start) + ", " + to_string(end) + ", " + to_string(seq.length()));
	}
	int replCnt = 0;
	for (int i = start; i < end; ++i) {		
		if (seqUse[static_cast<size_t>(i)] != repl) {
			seqUse[static_cast<size_t>(i)] = repl;
			++replCnt;
		}
	}
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
	if (sta < 0 || sto < sta) {
		throw std::invalid_argument("Invalid half-open depth interval: " + to_string(sta) + "-" + to_string(sto));
	}
	for (size_t i = 0; i < genes.size(); ++i) {
		for (const pair<int, int>& segment : genes[i]->segments) {
			// Depth intervals are BED/bedGraph-style [sta, sto); GFF segments are inclusive.
			int overlapStart = max(sta, segment.first);
			int overlapEnd = min(sto, segment.second + 1);
			if (overlapEnd > overlapStart) {
				genes[i]->addAccumDepth(static_cast<int64_t>(overlapEnd - overlapStart) *
					static_cast<int64_t>(depth));
			}
		}
	}
}

void fasta::addDepth(int sta, int sto, int depth, size_t source) {
	if (sta < 0 || sto < sta || depth < 0) {
		throw std::invalid_argument("Invalid depth contribution");
	}
	const int64_t contribution = static_cast<int64_t>(sto - sta) * static_cast<int64_t>(depth);
	depthAccum += contribution;
	if (depthAccumBySource.size() <= source) depthAccumBySource.resize(source + 1, 0);
	depthAccumBySource[source] += contribution;
}
void geneCollection::prepMuts() {
	if (mutsPrepared) { return; }
	mutsPrepared = true; 

	genesMut.clear(); genesMut.resize(genes.size(), nullptr);
	for (size_t i = 0; i < genes.size(); i++) {
		genesMut[i] = new gene(*genes[i]); //copy gene
	}
}


void geneCollection::correctCoords(int pos, const string& ref, const string& alt) {
	size_t sharedPrefix = 0;
	while (sharedPrefix < ref.size() && sharedPrefix < alt.size() &&
		std::toupper(static_cast<unsigned char>(ref[sharedPrefix])) ==
		std::toupper(static_cast<unsigned char>(alt[sharedPrefix]))) ++sharedPrefix;
	size_t sharedSuffix = 0;
	while (sharedSuffix < ref.size() - sharedPrefix && sharedSuffix < alt.size() - sharedPrefix &&
		std::toupper(static_cast<unsigned char>(ref[ref.size() - 1 - sharedSuffix])) ==
		std::toupper(static_cast<unsigned char>(alt[alt.size() - 1 - sharedSuffix]))) ++sharedSuffix;
	const int editStart = pos + static_cast<int>(sharedPrefix);
	const int refCoreLength = static_cast<int>(ref.size() - sharedPrefix - sharedSuffix);
	const int altCoreLength = static_cast<int>(alt.size() - sharedPrefix - sharedSuffix);
	const int editEnd = editStart + refCoreLength;
	const int geneDiff = altCoreLength - refCoreLength;
	if (geneDiff == 0) return; // substitutions do not move feature boundaries
	for (size_t i = 0; i < genesMut.size(); ++i) {
		vector<pair<int, int> > updated;
		for (pair<int, int> segment : genesMut[i]->segments) {
			if (refCoreLength == 0) {
				if (segment.second < editStart) {
					// Insertion is after this feature (including at its final anchor base).
				} else if (segment.first >= editStart) {
					segment.first += geneDiff;
					segment.second += geneDiff;
				} else {
					segment.second += geneDiff; // insertion lies inside the feature
				}
			} else if (segment.second < editStart) {
				// Before the replaced core.
			} else if (segment.first >= editEnd) {
				segment.first += geneDiff;
				segment.second += geneDiff;
			} else {
				// Overlap with the changed core. Never pull a feature boundary back to
				// the shared VCF anchor, which is outside the actual changed interval.
				if (segment.first >= editStart) segment.first = editStart;
				if (segment.second >= editEnd) segment.second += geneDiff;
				else segment.second = editStart + altCoreLength - 1;
			}
			if (segment.second >= segment.first) updated.push_back(segment);
		}
		std::sort(updated.begin(), updated.end());
		vector<pair<int, int> > merged;
		for (const pair<int, int>& segment : updated) {
			if (!merged.empty() &&
				static_cast<int64_t>(segment.first) <= static_cast<int64_t>(merged.back().second) + 1) {
				merged.back().second = max(merged.back().second, segment.second);
			} else {
				merged.push_back(segment);
			}
		}
		genesMut[i]->segments.swap(merged);
		genesMut[i]->recalculateGeometry();
	}
}


void geneCollection::writeAllGenes(options* opts, string& NTs, string& AAs, 
			bool doNT, bool doAA, fasta* fa, OutputStats* Ostats) {
	string ctg = fa->getMutatedSeq(opts);
	if (ctg.empty()) {
		return; //skip empty sequences
	}
	if (opts->skipEmptyContigs && countDeterminedBases(ctg) == 0) {
		return;
	}
	//string ret("");
	for (size_t i = 0; i < genesMut.size(); ++i) {
		string geneSeq = genesMut[i]->geneNT(ctg);
		int geneNonNs = static_cast<int>(countDeterminedBases(geneSeq));
		if (opts->skipEmptyGenes && (geneSeq.length() == 0 || 
			geneNonNs == 0)) {
			continue; //skip empty sequences
		}
		string hd = ">" + fa->sequenceId + "_" + to_string(genes[i]->getIdx()+1);
		if (opts->addHDTags) {
			hd+=genesMut[i]->createHDtag(geneSeq,fa, geneNonNs);
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
			const auto geneXs = count(AAseq.begin(), AAseq.end(), 'X');
			const int64_t geneNonXs = static_cast<int64_t>(AAseq.length()) -
				static_cast<int64_t>(geneXs);
			Ostats->totalGeneAAs += geneNonXs;
		}
		Ostats->totalGenes++;
	}
	
}


string fasta::write(options* opts, OutputStats* Ostats){//
	bool hdTags = opts->addHDTags;
	bool skipEmpty = opts->skipEmptyContigs;
	string ret("");
	prepMuts();
	string renderedSeq = getMutatedSeq(opts);
	const int64_t determinedBases = countDeterminedBases(renderedSeq);
	if (skipEmpty && determinedBases == 0){
		return ret; //skip empty sequences
	}

	Ostats->totalContigs++;
	Ostats->totalCtgNTs += determinedBases;

	ret += getMutatedHeader(hdTags, renderedSeq) + "\n";
	ret += renderedSeq + "\n";
	return ret;
}

int64_t fasta::getNonNcount(options* opts) {
	prepMuts();
	return countDeterminedBases(getMutatedSeq(opts));
}

void fasta::addGene(string id, int sta, int end, string strand, string type, int transTab,
	string partial, int gffPhase, const vector<pair<int, int> >& geneSegments) {
	std::unique_ptr<gene> AG(new gene(id,sta,end));
	AG->setStrand(strand);
	AG->setType(type);
	AG->setTranslationTable(transTab);
	AG->setPhase(gffPhase);
	AG->setSegments(geneSegments);
	if (geneCol->size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::overflow_error("Too many genes on contig " + sequenceId);
	}
	AG->setNumOnContig(static_cast<int>(geneCol->size()));
	AG->setPartial(partial);
	geneCol->push_back(AG.release());
}





refAssembly::refAssembly(options* opt):
	NSeqs(0), refFasta(opt->refFasta),  
	fastas(0),
	hd2ID(0), opts(opt)
{
	const std::string& filename = opt->refFasta;
	cout << "Reading reference assembly from: " << filename << endl;
	istream* file = openGZUZ(filename);
	int64_t totalBp(0);
	std::string line; string sequence(""); string header("");
	while (std::getline(*file, line)) {
		if (!line.empty() && line.back() == '\r') line.pop_back();
		if (line.empty() ) {
			continue; // Skip empty lines and headers
		}
		if (line[0] == '>') {					
			// Store the previous sequence before starting a new one
			if (!header.empty()) {
				fastas.push_back(new fasta(sequence, header));
				addSequenceLength(totalBp, sequence.length());
				sequence.clear();
				NSeqs++;
			}
			header = line.substr(1); // Remove '>' character
			continue;
		}
		if (header.empty()) {
			delete file;
			throw std::runtime_error("Reference FASTA contains sequence data before its first header: " + filename);
		}
		sequence += line;
	}
	delete file;
	//add last readin batch:
	if (!header.empty()) {
		fastas.push_back(new fasta(sequence, header));
		addSequenceLength(totalBp, sequence.length());
		sequence.clear();
		NSeqs++;
	}


	// GFF3/VCF seqids match the first whitespace-delimited token in a FASTA header.
	if (fastas.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::overflow_error("Too many reference sequences");
	}
	for (size_t i = 0; i < fastas.size(); ++i) {
		const string id = fastas[i]->getSequenceId();
		auto existing = hd2ID.find(id);
		if (existing != hd2ID.end() && existing->second != static_cast<int>(i)) {
			throw std::runtime_error("Duplicate FASTA sequence identifier: " + id);
		}
		hd2ID[id] = static_cast<int>(i);
		const string fullHeader = fastas[i]->getHeader();
		if (fullHeader != id && hd2ID.find(fullHeader) == hd2ID.end()) {
			hd2ID[fullHeader] = static_cast<int>(i);
		}
	}

	cout << "Number of sequences: " << NSeqs << " containing "<< totalBp<< "bp."<< endl;
}

refAssembly::~refAssembly() {
	for (fasta* sequence : fastas) {
		delete sequence;
	}
}

void refAssembly::setFasta(string id, fasta* f){
	auto it = hd2ID.find(id);
	if (it != hd2ID.end()) {
		fastas[static_cast<size_t>(it->second)] = f;
	} else {
		throw std::runtime_error("Header not found: " + id);
	}
}

fasta* refAssembly::getFasta(string id)
{ 
	vector<string> candidates;
	candidates.push_back(trimCopy(id));
	candidates.push_back(percentDecode(candidates.front()));
	candidates.push_back(canonicalSequenceId(candidates.front()));
	candidates.push_back(canonicalSequenceId(candidates[1]));
	for (const string& candidate : candidates) {
		auto it = hd2ID.find(candidate);
		if (it != hd2ID.end()) return fastas[static_cast<size_t>(it->second)];
	}
	throw std::runtime_error("Sequence identifier not found in reference FASTA: " + id);
}
bool refAssembly::isSequence(string id) {
	vector<string> candidates;
	candidates.push_back(trimCopy(id));
	candidates.push_back(percentDecode(candidates.front()));
	candidates.push_back(canonicalSequenceId(candidates.front()));
	candidates.push_back(canonicalSequenceId(candidates[1]));
	for (const string& candidate : candidates) {
		if (hd2ID.find(candidate) != hd2ID.end()) return true;
	}
	return false;
}

void refAssembly::readGFF() {
	istream* file = openGZUZ(opts->gffFile);
	cout << "Reading GFF file from: " << opts->gffFile << endl;

	struct SegmentRecord {
		int start;
		int end;
		int phase;
	};
	struct FeatureRecord {
		string chrom;
		string id;
		string strand;
		string type;
		int translationTable;
		bool coding;
		bool nonFunctional;
		bool translExcept;
		bool missingRequiredPhase;
		bool partialLeft;
		bool partialRight;
		size_t sourceLine;
		vector<SegmentRecord> segments;
	};

	vector<FeatureRecord> features;
	unordered_map<string, size_t> featureByKey;
	unordered_map<string, int> sequenceTranslationTables;
	unordered_map<string, vector<string> > featureParents;
	set<string> nonFunctionalFeatureIds;
	set<string> codingFeatureRelations;
	string line;
	int commentTranslationTable = 11;
	size_t lineNumber = 0;
	size_t skippedFeatureRows = 0;
	while (std::getline(*file, line)) {
		++lineNumber;
		if (!line.empty() && line.back() == '\r') line.pop_back();
		if (line.empty() ) {
			continue; // Skip empty lines
		}
		if (line[0] == '#') {
			if (trimCopy(line) == "##FASTA") break;
			commentTranslationTable = translationTableFromText(line, commentTranslationTable);
			continue;
		}

		vector<string> fields = splitTabs(line);
		if (fields.size() != 9) {
			delete file;
			throw std::runtime_error("GFF3 line " + to_string(lineNumber) +
				" must contain exactly 9 tab-delimited columns: " + line);
		}
		string chrom = percentDecode(trimCopy(fields[0]));
		string type = trimCopy(fields[2]);
		string strand = trimCopy(fields[6]);
		auto attributes = parseGffAttributes(fields[8]);
		const bool rowNonFunctional = isNonCodingGene(attributes);
		auto rowId = attributes.find("ID");
		if (rowId != attributes.end() && !rowId->second.empty()) {
			const string featureKey = chrom + "\x1f" + rowId->second;
			if (rowNonFunctional) nonFunctionalFeatureIds.insert(featureKey);
			auto parent = attributes.find("Parent");
			if (parent != attributes.end() && !parent->second.empty()) {
				stringstream parents(parent->second);
				string parentId;
				while (getline(parents, parentId, ',')) {
					parentId = trimCopy(parentId);
					if (!parentId.empty()) featureParents[featureKey].push_back(chrom + "\x1f" + parentId);
				}
			}
		}

		long long startOneBased = 0;
		long long endOneBased = 0;
		if (!parseLongLongStrict(fields[3], startOneBased) ||
			!parseLongLongStrict(fields[4], endOneBased)) {
			delete file;
			throw std::runtime_error("Invalid GFF3 coordinates on line " + to_string(lineNumber) + ": " + line);
		}
		if (startOneBased < 1 || endOneBased < startOneBased || endOneBased > INT_MAX) {
			delete file;
			throw std::runtime_error("Out-of-range GFF3 coordinates on line " + to_string(lineNumber) + ": " + line);
		}
		int start = static_cast<int>(startOneBased - 1);
		int end = static_cast<int>(endOneBased - 1);

		int translationTable = commentTranslationTable;
		auto sequenceTable = sequenceTranslationTables.find(chrom);
		if (sequenceTable != sequenceTranslationTables.end()) translationTable = sequenceTable->second;
		auto tableAttribute = attributes.find("transl_table");
		if (tableAttribute != attributes.end()) {
			if (!parseIntStrict(tableAttribute->second, translationTable)) {
				delete file;
				throw std::runtime_error("Invalid transl_table on GFF3 line " + to_string(lineNumber));
			}
		}

		if (lowerCopy(type) == "region") {
			sequenceTranslationTables[chrom] = translationTable;
			continue;
		}

		bool coding = isCodingFeatureType(type);
		bool geneFallback = isGeneFeatureType(type) && !rowNonFunctional;
		if (!coding && !geneFallback) {
			++skippedFeatureRows;
			continue;
		}
		if (strand != "+" && strand != "-") {
			delete file;
			throw std::runtime_error("Coding/gene feature has invalid strand on GFF3 line " +
				to_string(lineNumber) + ": " + strand);
		}
		const bool rowTranslExcept = attributes.find("transl_except") != attributes.end();

		int gffPhase = 0;
		bool missingRequiredPhase = false;
		if (fields[7] != ".") {
			if (fields[7].size() != 1 || fields[7][0] < '0' || fields[7][0] > '2') {
				delete file;
				throw std::runtime_error("Invalid GFF3 phase on line " + to_string(lineNumber) + ": " + fields[7]);
			}
			gffPhase = fields[7][0] - '0';
		} else if (isCdsFeatureType(type)) {
			missingRequiredPhase = true;
		}

		bool partialLeft = false;
		bool partialRight = false;
		bool rangeSpecified = attributes.find("start_range") != attributes.end() ||
			attributes.find("end_range") != attributes.end();
		auto partialAttribute = attributes.find("partial");
		if (partialAttribute != attributes.end()) {
			string value = lowerCopy(partialAttribute->second);
			if (value.size() >= 2 && (value[0] == '0' || value[0] == '1') &&
				(value[1] == '0' || value[1] == '1')) {
				partialLeft = value[0] == '1';
				partialRight = value[1] == '1';
			} else if ((value == "true" || value == "1" || value == "yes") && !rangeSpecified) {
				partialLeft = partialRight = true;
			}
		}
		if (attributes.find("start_range") != attributes.end()) partialLeft = true;
		if (attributes.find("end_range") != attributes.end()) partialRight = true;

		string id;
		static const char* idKeys[] = { "ID", "locus_tag", "protein_id", "Name", "Parent", "transcript_id", "gene_id" };
		for (const char* key : idKeys) {
			auto found = attributes.find(key);
			if (found != attributes.end() && !found->second.empty()) {
				id = found->second;
				break;
			}
		}
		if (id.empty()) id = "feature_at_line_" + to_string(lineNumber);

		string groupId = id;
		if (coding) {
			static const char* groupingKeys[] = { "Parent", "transcript_id", "protein_id", "locus_tag", "ID", "gene_id" };
			for (const char* key : groupingKeys) {
				auto found = attributes.find(key);
				if (found != attributes.end() && !found->second.empty()) {
					groupId = found->second;
					break;
				}
			}
			if (groupId.find(',') != string::npos) {
				delete file;
				throw std::runtime_error("CDS/ORF with multiple Parent values is not supported on GFF3 line " +
					to_string(lineNumber));
			}
			codingFeatureRelations.insert(chrom + "\x1f" + groupId);
			auto parent = attributes.find("Parent");
			if (parent != attributes.end()) codingFeatureRelations.insert(chrom + "\x1f" + parent->second);
			auto ownId = attributes.find("ID");
			if (ownId != attributes.end()) codingFeatureRelations.insert(chrom + "\x1f" + ownId->second);
		}

		string groupKey = chrom + "\x1f" + lowerCopy(type) + "\x1f" + groupId;
		auto existing = featureByKey.find(groupKey);
		if (existing == featureByKey.end()) {
			FeatureRecord feature;
			feature.chrom = chrom;
			feature.id = groupId;
			feature.strand = strand;
			feature.type = type;
			feature.translationTable = translationTable;
			feature.coding = coding;
			feature.nonFunctional = rowNonFunctional;
			feature.translExcept = rowTranslExcept;
			feature.missingRequiredPhase = missingRequiredPhase;
			feature.partialLeft = partialLeft;
			feature.partialRight = partialRight;
			feature.sourceLine = lineNumber;
			feature.segments.push_back({ start, end, gffPhase });
			featureByKey[groupKey] = features.size();
			features.push_back(feature);
		} else {
			FeatureRecord& feature = features[existing->second];
			if (feature.strand != strand || feature.translationTable != translationTable) {
				delete file;
				throw std::runtime_error("Inconsistent multipart GFF3 feature " + id + " on line " +
					to_string(lineNumber));
			}
			feature.partialLeft = feature.partialLeft || partialLeft;
			feature.partialRight = feature.partialRight || partialRight;
			feature.nonFunctional = feature.nonFunctional || rowNonFunctional;
			feature.translExcept = feature.translExcept || rowTranslExcept;
			feature.missingRequiredPhase = feature.missingRequiredPhase || missingRequiredPhase;
			feature.segments.push_back({ start, end, gffPhase });
		}
	}
	delete file;

	bool inheritedPseudoChanged = true;
	while (inheritedPseudoChanged) {
		inheritedPseudoChanged = false;
		for (const auto& relation : featureParents) {
			if (nonFunctionalFeatureIds.find(relation.first) != nonFunctionalFeatureIds.end()) continue;
			for (const string& parent : relation.second) {
				if (nonFunctionalFeatureIds.find(parent) != nonFunctionalFeatureIds.end()) {
					nonFunctionalFeatureIds.insert(relation.first);
					inheritedPseudoChanged = true;
					break;
				}
			}
		}
	}
	auto featureIsNonFunctional = [&nonFunctionalFeatureIds](const FeatureRecord& feature) {
		return feature.nonFunctional ||
			nonFunctionalFeatureIds.find(feature.chrom + "\x1f" + feature.id) !=
				nonFunctionalFeatureIds.end();
	};
	for (const FeatureRecord& feature : features) {
		if (featureIsNonFunctional(feature)) continue;
		if (feature.missingRequiredPhase) {
			throw std::runtime_error("CDS feature is missing its required GFF3 phase on line " +
				to_string(feature.sourceLine));
		}
		if (opts->outputTypes.find('A') != string::npos &&
			!supportedTranslationTable(feature.translationTable)) {
			throw std::runtime_error("Unsupported NCBI translation table " +
				to_string(feature.translationTable) + " on GFF3 line " +
				to_string(feature.sourceLine) + " (supported: 1, 4, 11, 25)");
		}
		if (opts->outputTypes.find('A') != string::npos && feature.translExcept) {
			throw std::runtime_error("GFF3 transl_except is not yet supported for amino-acid output (line " +
				to_string(feature.sourceLine) + "); refusing to emit a potentially incorrect protein");
		}
	}
	codingFeatureRelations.clear();
	for (const FeatureRecord& feature : features) {
		if (feature.coding && !featureIsNonFunctional(feature)) {
			vector<string> relations(1, feature.chrom + "\x1f" + feature.id);
			set<string> visited;
			while (!relations.empty()) {
				const string relation = relations.back();
				relations.pop_back();
				if (!visited.insert(relation).second) continue;
				codingFeatureRelations.insert(relation);
				auto parents = featureParents.find(relation);
				if (parents != featureParents.end()) {
					relations.insert(relations.end(), parents->second.begin(), parents->second.end());
				}
			}
		}
	}

	size_t addedFeatures = 0;
	size_t fallbackGenesSkipped = 0;
	size_t nonFunctionalFeaturesSkipped = 0;
	set<string> codingFeatureGeometry;
	for (const FeatureRecord& feature : features) {
		if (!feature.coding || feature.segments.empty() || featureIsNonFunctional(feature)) continue;
		int first = feature.segments.front().start;
		int last = feature.segments.front().end;
		for (const SegmentRecord& segment : feature.segments) {
			first = min(first, segment.start);
			last = max(last, segment.end);
		}
		codingFeatureGeometry.insert(feature.chrom + "\x1f" + feature.strand + "\x1f" +
			to_string(first) + "\x1f" + to_string(last));
	}
	for (const FeatureRecord& feature : features) {
		if (featureIsNonFunctional(feature)) {
			++nonFunctionalFeaturesSkipped;
			continue;
		}
		if (!feature.coding) {
			int first = feature.segments.front().start;
			int last = feature.segments.front().end;
			for (const SegmentRecord& segment : feature.segments) {
				first = min(first, segment.start);
				last = max(last, segment.end);
			}
			const string geometry = feature.chrom + "\x1f" + feature.strand + "\x1f" +
				to_string(first) + "\x1f" + to_string(last);
			if (codingFeatureRelations.find(feature.chrom + "\x1f" + feature.id) != codingFeatureRelations.end() ||
				codingFeatureGeometry.find(geometry) != codingFeatureGeometry.end()) {
				++fallbackGenesSkipped;
				continue;
			}
		}
		fasta* target = getFasta(feature.chrom);
		vector<pair<int, int> > segments;
		segments.reserve(feature.segments.size());
		vector<const SegmentRecord*> orientedSegments;
		orientedSegments.reserve(feature.segments.size());
		for (const SegmentRecord& segment : feature.segments) {
			if (segment.start < 0 || segment.end >= target->getLength()) {
				throw std::runtime_error("GFF3 feature " + feature.id + " on " + feature.chrom +
					" lies outside reference sequence length " + to_string(target->getLength()));
			}
			segments.push_back(make_pair(segment.start, segment.end));
			orientedSegments.push_back(&segment);
		}
		std::sort(orientedSegments.begin(), orientedSegments.end(), [&feature](
			const SegmentRecord* lhs, const SegmentRecord* rhs) {
			return feature.strand == "+" ? lhs->start < rhs->start : lhs->end > rhs->end;
		});
		const SegmentRecord* firstOrientedSegment = orientedSegments.empty() ? nullptr : orientedSegments.front();
		if (isCdsFeatureType(feature.type) && firstOrientedSegment != nullptr) {
			const int firstLength = firstOrientedSegment->end - firstOrientedSegment->start + 1;
			if (firstOrientedSegment->phase > firstLength) {
				throw std::runtime_error("GFF3 phase exceeds the first CDS segment length for " + feature.id);
			}
			int64_t translatedBases = static_cast<int64_t>(firstLength - firstOrientedSegment->phase);
			for (size_t i = 1; i < orientedSegments.size(); ++i) {
				const int expectedPhase = static_cast<int>((3 - (translatedBases % 3)) % 3);
				if (orientedSegments[i]->phase != expectedPhase) {
					throw std::runtime_error("Inconsistent multipart CDS phase for " + feature.id +
						": expected " + to_string(expectedPhase) + " but found " +
						to_string(orientedSegments[i]->phase));
				}
				translatedBases += static_cast<int64_t>(
					orientedSegments[i]->end - orientedSegments[i]->start + 1);
			}
		}
		string partial;
		partial += feature.partialLeft ? '1' : '0';
		partial += feature.partialRight ? '1' : '0';
		target->addGene(feature.id, segments.front().first, segments.front().second,
			feature.strand, feature.type, feature.translationTable, partial,
			firstOrientedSegment == nullptr ? 0 : firstOrientedSegment->phase, segments);
		++addedFeatures;
	}
	cout << "GFF features retained: " << addedFeatures << "; skipped non-gene/non-coding rows: "
		<< skippedFeatureRows << "; skipped gene parents where coding features were present: "
		<< fallbackGenesSkipped << "; skipped nonfunctional/pseudo features: "
		<< nonFunctionalFeaturesSkipped << endl;
}


void refAssembly::processDepth(const string& inF, int minDepth, size_t source) {
	
	//open connection to file
	std::unique_ptr<istream> in(openGZUZ(inF));

	cout << "Reading depth file from: " << inF << endl;

	string line; string curChrom(""); //int curChromIdx(-1); string curChromSeq("");
	fasta* curFasta(nullptr); int curChromL(0);  // Initialized to 0
	int64_t totalChromL(0);
	for (const fasta* sequence : fastas) {
		addSequenceLength(totalChromL, static_cast<size_t>(sequence->getLength()));
	}
	int64_t cntPosKept = minDepth == 0 ? totalChromL : 0;
	if (minDepth == 0) {
		// A zero minimum explicitly disables coverage masking for this source,
		// including bases omitted from a sparse bedGraph file (implicit depth 0).
		for (fasta* sequence : fastas) {
			if (sequence->getLength() > 0) sequence->unmaskSeq(0, sequence->getLength());
		}
	}
	vector<int64_t> posKeptPerDepth(10 , 0); // thresholds 0 through 9
	posKeptPerDepth[0] = totalChromL; // omitted intervals have implicit depth zero
	set<string> seenDepthContigs;
	int previousStop = 0;


	geneCollection* genes(nullptr);
	
	//start reading depth file line by line
	//example: "A.14.1640M2__C15_L=476= 25      104     2"
	while (std::getline(*in, line)) {
		if (line.empty()) {
			continue; // Skip empty lines
		}
		if (line[0] == '#' || line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) continue;
		std::istringstream iss(line);
		string header;
		int sta, sto, depth;
		if (!(iss >> header >> sta >> sto >> depth)) {
			throw std::runtime_error("Error reading depth file: " + line + "\n" + inF);
		}
		if (sta < 0 || sto <= sta || depth < 0) {
			throw std::runtime_error("Invalid half-open depth interval or negative depth: " + line + "\n" + inF);
		}
		if (curChrom != header) {//switch to new chromosome
			curChrom = header;
			curFasta = getFasta(curChrom);
			const string canonicalChrom = curFasta->getSequenceId();
			if (!seenDepthContigs.insert(canonicalChrom).second) {
				throw std::runtime_error("Depth contig occurs in multiple blocks (file must be sorted): " + header);
			}
			previousStop = 0;
			curChromL = curFasta->getLength();
			genes = curFasta->getGeneCollection();
		}
		
		//some basic checks:
		if (sto > curChromL) {
			throw std::runtime_error("Error: Stop post seq length: " + line + "\n" + inF);
		}
		if (sta < previousStop) {
			throw std::runtime_error("Depth intervals overlap or are unsorted: " + line + "\n" + inF);
		}
		previousStop = sto;
		curFasta->addDepth(sta, sto, depth, source);
		genes->depthInGenes(sta, sto, depth);
		//replace  seq with N where < depthThreshold
		if (depth >= minDepth) {
			curFasta->unmaskSeq(sta, sto);
			if (minDepth > 0) cntPosKept += static_cast<int64_t>(sto - sta);
		}
		//for stats only, check at which depth positions where kept
		for (size_t depthIdx = 1; depthIdx < posKeptPerDepth.size(); ++depthIdx) {
			if (depth >= static_cast<int>(depthIdx)) {
				posKeptPerDepth[depthIdx] += (sto - sta);
			} else {
				break;
			}
		}
	}
	const int64_t cntPosRm = totalChromL - cntPosKept;
	//general stats on depth filter
	cout << "Depth filter (" << minDepth << "): Number of positions kept: " << cntPosKept << " removed: " << cntPosRm << endl;

	//stats on depth filter per depth bin
	cout << "Theoretical depth filter stats per depth bin, increasing from 0 to "
		<< (posKeptPerDepth.size() - 1) << " (of max "<< totalChromL<<"): " ;
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
		this->processDepth(inF, minDepth, i);
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
		int64_t totalBp(0);
		for (size_t i = 0; i < fastas.size(); ++i) {
			totalBp += fastas[i]->getNonNcount(opts);
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
