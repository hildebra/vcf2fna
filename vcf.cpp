#include "vcf.h"

#include <cctype>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>

namespace {

float missingFloat() {
	return std::numeric_limits<float>::quiet_NaN();
}

bool present(float value) {
	return std::isfinite(value);
}

string lowerCopy(string value) {
	transform(value.begin(), value.end(), value.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return value;
}

bool parseFloatValue(const string& text, float& value) {
	if (text.empty() || text == ".") return false;
	char* end = nullptr;
	const double parsed = std::strtod(text.c_str(), &end);
	if (end == text.c_str() || end == nullptr || *end != '\0' ||
		!std::isfinite(parsed) ||
		parsed > static_cast<double>(std::numeric_limits<float>::max()) ||
		parsed < -static_cast<double>(std::numeric_limits<float>::max())) {
		return false;
	}
	// strtof/std::stof may report a valid subnormal value as ERANGE. Parse as
	// double and narrow deliberately so tiny VCF values (for example FS
	// p-values around 1e-41) remain valid; values below float range become zero.
	value = static_cast<float>(parsed);
	return true;
}

bool parseIntValue(const string& text, int& value) {
	if (text.empty() || text == ".") return false;
	size_t used = 0;
	try {
		long parsed = stol(text, &used);
		if (used != text.size() || parsed < std::numeric_limits<int>::min() ||
			parsed > std::numeric_limits<int>::max()) return false;
		value = static_cast<int>(parsed);
	} catch (const std::exception&) {
		return false;
	}
	return true;
}

bool isSequenceAllele(const string& allele) {
	if (allele.empty()) return false;
	for (char rawBase : allele) {
		const unsigned char c = static_cast<unsigned char>(rawBase);
		char upper = static_cast<char>(std::toupper(c));
		if (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T' && upper != 'N') return false;
	}
	return true;
}

} // namespace


VCFReader::VCFReader(options* opt, refAssembly* R): 
	refG(R), opts(opt), header(""), cF(nullptr),
	vcfFile(opts->inVCF),
	curChrom(""), cntAvContigs(0),
	lnCnt(0), //snpCNT(0), indelCNT(0), snpFILT(0),indelFILT(0), 
	//unsrSNP(0), unsrINDEL(0), SNPused(0), indelUsed(0), conflictCnt(0),
	seenCtgs(0),
	snpRepl(0), snpKept(0),
	VCF2ptr(nullptr), Vstats(nullptr), filterStats(nullptr)
	//VF(nullptr)
{


	Vstats.reset(new VariantStats());
	filterStats.reset(new VCFfilterStats());


	// Open the VCF file

	// A second file is combined as independent evidence; platform order is not assumed.
	if (vcfFile.size()==2) {
		std::unique_ptr<VCFcollection> secondCalls(new VCFcollection());
		std::unique_ptr<istream> in2(openGZUZ(vcfFile[1]));
		read_vcf_file(in2.get(), 1, secondCalls.get());
		VCF2ptr = secondCalls.get();
		std::unique_ptr<istream> in1(openGZUZ(vcfFile[0]));
		read_vcf_file(in1.get(),0,nullptr);
		VCF2ptr = nullptr;

	} else if (vcfFile.size()==1){//just one file
		std::unique_ptr<istream> in(openGZUZ(vcfFile[0]));
		read_vcf_file(in.get(),0,nullptr);
	}
	else {
		throw std::runtime_error("Only one or two VCF files are supported: " + combCommas(vcfFile));
	}

}


void VCFReader::read_vcf_file(istream* fp, int VCFnum, VCFcollection* Vcol) {
	std::string line;
	VCFmulti collector;
	string inputChrom;
	int previousPos = -1;
	bool sawHeader = false;
	size_t expectedColumns = 0;

	while (getline(*fp, line)) {
		if (!line.empty() && line.back() == '\r') line.pop_back();
		//std::cout << line << std::endl;;
		if (line.empty()) continue;
		if (line[0] == '#') {
			if (line.size() > 1 && line[1] == '#') {
				// meta line
				if (line.substr(2, 7) == "contig=") {
					size_t idStart = line.find("ID=");
					size_t idEnd = idStart == string::npos ? string::npos : line.find_first_of(",>", idStart + 3);
					string tmpContig = idStart == string::npos ? "" : line.substr(idStart + 3, idEnd - (idStart + 3));
					cntAvContigs++;//
					
					if (!tmpContig.empty() && !refG->isSequence(tmpContig)) {
						limitedWarning("VCF contig absent from reference",
							"contig " + tmpContig + " not found in reference assembly.");
					}
				}
				//meta.push_back(line);
			} else {
				// header line
				if (sawHeader) {
					throw std::runtime_error("VCF contains more than one #CHROM header line");
				}
				header = line;
				vector<string> headerColumns;
				stringstream headerStream(header);
				string headerColumn;
				while (getline(headerStream, headerColumn, '\t')) headerColumns.push_back(headerColumn);
				static const char* fixedColumns[] = {
					"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
				};
				bool validFixedColumns = headerColumns.size() >= 8;
				for (size_t i = 0; validFixedColumns && i < 8; ++i) {
					validFixedColumns = headerColumns[i] == fixedColumns[i];
				}
				if (!validFixedColumns || (headerColumns.size() != 8 && headerColumns.size() != 10) ||
					(headerColumns.size() == 10 &&
						(headerColumns[8] != "FORMAT" || headerColumns[9].empty()))) {
					throw std::runtime_error(
						"VCF header must contain the eight fixed columns and, optionally, FORMAT plus one pooled sample: " +
						header);
				}
				expectedColumns = headerColumns.size();
				sawHeader = true;
			}
		} else {
			if (!sawHeader) {
				throw std::runtime_error("VCF data record encountered before the #CHROM header line");
			}
			const size_t recordColumns = static_cast<size_t>(
				std::count(line.begin(), line.end(), '\t')) + 1;
			if (recordColumns != expectedColumns) {
				throw std::runtime_error("VCF record column count does not match its #CHROM header");
			}
			lnCnt++;
			if (opts->debug1) cerr << "VCF[" << VCFnum + 1 << "] " << line << '\n';
			
			std::unique_ptr<VCFmem> parsed;
			try {
				parsed.reset(new VCFmem(line, VCFnum));
			} catch (const std::exception& e) {
				throw std::runtime_error("Invalid VCF record at data line " + to_string(lnCnt) + ": " + e.what() + "\n" + line);
			}
			VCFmem* vcf = parsed.get();
			if (vcf->getChrom() != inputChrom) {
				inputChrom = vcf->getChrom();
				previousPos = -1;
			} else if (vcf->getPos() < previousPos) {
				throw std::runtime_error("VCF positions are not sorted within contig " + inputChrom);
			}
			previousPos = vcf->getPos();

			if (Vcol != nullptr) {
				Vcol->addVCF(parsed.release());
				continue;
			}
			//take care of loading correct fasta..
			if (vcf->getChrom() != curChrom) {

				collector.evalVCFs(cF, VCF2ptr, Vstats.get(), filterStats.get(), opts);

				curChrom = vcf->getChrom();
				cF = refG->getFasta(curChrom);
				const string canonicalChrom = cF->getSequenceId();
				if (seenCtgs.find(canonicalChrom) != seenCtgs.end()) {
					throw std::runtime_error("Chromosome occurring multiple times. Is the VCF sorted?: " + curChrom);
				} else {
					seenCtgs[canonicalChrom] = 1;
				}
			}

			collector.addVCFmem(parsed.release());

		}
	}
	if (!sawHeader) {
		throw std::runtime_error("VCF is missing the required #CHROM header line");
	}
	if (Vcol == nullptr) {
		// Contig switches flush previous batches; EOF must flush the final batch.
		collector.evalVCFs(cF, VCF2ptr, Vstats.get(), filterStats.get(), opts);
		// Process contigs called only by the second input. Calls from a paired
		// contig are consumed by evalVCFs, leaving only genuinely unseen contigs.
		if (VCF2ptr != nullptr) {
			for (const string& contig : VCF2ptr->contigs()) {
				VCFmulti secondaryOnly;
				for (VCFmem* call : VCF2ptr->takeContig(contig)) secondaryOnly.addVCFmem(call);
				secondaryOnly.evalVCFs(refG->getFasta(contig), nullptr, Vstats.get(), filterStats.get(), opts);
			}
		}
		cout << endl << "Read " << cntAvContigs << " contigs and " << lnCnt << " entries from VCF file." << endl;
		Vstats->printSNPstats();
/*		cout << endl << "Read " << cntAvContigs << " contigs and " << lnCnt << " entries from VCF file." << endl;
		cout << "  - Found " << Vstats->snpCNT << " SNPs and " << Vstats->indelCNT << " INDELS." << endl;
		cout << "  - Filtered " << Vstats->snpFILT << ";" << Vstats->indelFILT << " entries (SNP; INDEL)." << endl;
		cout << "  - Passed " << Vstats->SNPused << ";" << Vstats->indelUsed << " SNPs and INDELS. Conflicts resolved: "<< Vstats->conflictCnt << endl;
		cout << "  - Unsure " << Vstats->unsrSNP << ";" << Vstats->unsrINDEL << " SNPs and INDELS replaced with N." << endl;
		*/
		filterStats->printStats();
	}
	else {
		cout << "Read and stored " << Vcol->size() << " vcf entries\n";
		cntAvContigs = 0; lnCnt = 0; //reset vars
	}
}


/*
minQual(opts->minCallQual), minDep(opts->minDepthPar),
minFS(opts->minFS), minMQ0F(opts->minMQ0F),
minBQBZ(opts->minBQBZ), minSP(opts->minSP),
*/

bool VCFmem::isNoCall() const {
	return alt == "." || alt == "*" || alt == "<*>" || alt == "<NON_REF>";
}

float VCFmem::evidenceDepth() const {
	if (totalAlleleDepth > 0) return static_cast<float>(totalAlleleDepth);
	if (present(DPval) && DPval >= 0.f) return DPval;
	long long dp4Total = 0;
	for (int value : DP4) dp4Total += value;
	return dp4Total > 0 ? static_cast<float>(dp4Total) : missingFloat();
}

bool VCFmem::callerFilterFails(const options* opts, VCFfilterStats* stats) const {
	if (filter.empty() || filter == "PASS" || filter == ".") return false;
	string policy = lowerCopy(opts->vcfFilterPolicy);
	if (policy == "ignore") {
		stats->VCFfilterIgnored++;
		return false;
	}
	if (policy == "all") {
		stats->VCFfilter++;
		return true;
	}

	// FILTER IDs are caller-defined. In the default "technical" mode, named
	// failures are vetoes except for labels that clearly describe a diploid,
	// germline/somatic, or population-genetic model. Those labels are not a
	// sound reason to discard pooled bacterial evidence and are reported as
	// ignored. The stricter "all" policy is available for caller-tuned VCFs.
	bool callerModelOnly = true;
	stringstream ids(filter);
	string id;
	while (getline(ids, id, ';')) {
		string token = lowerCopy(id);
		token.erase(remove_if(token.begin(), token.end(),
			[](char c) { return c == '_' || c == '-' || c == '.'; }), token.end());
		static const set<string> modelOnlyIds = {
			"excesshet", "excessheterozygosity", "inbreeding", "inbreedingcoeff",
			"hardyweinberg", "hwe", "ploidy", "diploid", "diploidmodel",
			"genotype", "genotypemodel", "germline", "germlinefilter",
			"somatic", "somaticfilter", "allelebalance", "homref",
			"heterozygous", "monomorphic"
		};
		if (modelOnlyIds.find(token) == modelOnlyIds.end()) callerModelOnly = false;
	}
	if (callerModelOnly) {
		stats->VCFfilterIgnored++;
		return false;
	}
	stats->VCFfilter++;
	return true;
}

bool VCFmem::evalVCFentry(options* opts,  float mD, int lastIndelPos, int nextIndelPos,
	VCFfilterStats* filterStats) {
	conflict = false;
	conflictSpan = 0;
	isFiltered = callerFilterFails(opts, filterStats);
	if (isNoCall()) {
		isIndel = false;
		isSnp = false;
		return false;
	}
	if (!isIndel && !isSnp) return false;
	if (alt.find('N') != string::npos) {
		isFiltered = true;
		filterStats->ambiguousAllele++;
	}

	const float depth = evidenceDepth();
	const float qualityForFormula = present(QUALval) ? max(QUALval, 1.f) : 20.f;
	bool depthFailed = !present(depth);

	if (altFreq < opts->minAltFreq || (altDepth >= 0 && altDepth < opts->minAltReads)) {
		depthFailed = true;
	}
	if (altDepth < 0) depthFailed = true;
	if (present(depth)) {
		float absoluteMinimum = sourceFile >= 0 && static_cast<size_t>(sourceFile) < opts->minDepthPar.size()
			? static_cast<float>(opts->minDepthPar[static_cast<size_t>(sourceFile)]) : 0.f;
		float relativeMinimum = mD > 0.f ? min(8.f, mD * opts->depthFilterScale) : 0.f;
		if (depth < max(absoluteMinimum, relativeMinimum)) depthFailed = true;
		if (mD > 0.f && opts->maxDepthFilterScale > 0.f && depth > mD * opts->maxDepthFilterScale) {
			depthFailed = true;
		}
	}
	if ((!present(QUALval) && opts->minCallQual > 0) ||
		(present(QUALval) && QUALval < static_cast<float>(opts->minCallQual))) {
		isFiltered = true;
		filterStats->QUAL++;
	}
	if (opts->minCallQualAdaptive > 0.f && mD > 0.f &&
		(!present(QUALval) || QUALval < min(20.f, mD * opts->minCallQualAdaptive))) {
		isFiltered = true;
		filterStats->QUALadaptive++;
	}

	if (isIndel) {
		if (present(IDVval) && IDVval < static_cast<float>(opts->minAltReads)) depthFailed = true;
		if (present(IMFval) && IMFval < max(0.1f, opts->minAltFreq)) depthFailed = true;
		if (present(MQBZval) && present(depth) && MQBZval < -(5.f + depth / 20.f)) {
			isFiltered = true;
			filterStats->MQBZ++;
		}
		if (present(RPBZval) && present(SCBZval) && RPBZval + SCBZval > 9.f) {
			isFiltered = true;
			filterStats->RPBZ++;
		}
	} else if (isSnp) {
		if ((lastIndelPos >= 0 && std::abs(posN - lastIndelPos) <= static_cast<int>(opts->indelRange)) ||
			(nextIndelPos >= 0 && std::abs(posN - nextIndelPos) <= static_cast<int>(opts->indelRange))) {
			isFiltered = true;
			filterStats->indelProx++;
		}
		if (present(SPval) && SPval > opts->minSP + (present(depth) ? depth / 2.f : 0.f)) {
			isFiltered = true;
			filterStats->SP++;
		}
		// htslib guidance uses only the low BQBZ tail; a large positive value is not equivalent.
		if (present(BQBZval) && BQBZval < opts->minBQBZ - (present(depth) ? depth / 40.f : 0.f)) {
			isFiltered = true;
			filterStats->BQBZ++;
		}
		float mqThreshold = -(3.5f + (present(depth) ? 4.f * depth / qualityForFormula : 0.f));
		if (present(MQBZval) && MQBZval < mqThreshold) {
			isFiltered = true;
			filterStats->MQBZ++;
		}
		float rpThreshold = 3.f + (present(depth) ? 3.f * depth / qualityForFormula : 0.f);
		if (present(RPBZval) && std::abs(RPBZval) > rpThreshold) {
			isFiltered = true;
			filterStats->RPBZ++;
		}
		if (present(SCBZval) && SCBZval > 2.5f + (present(depth) ? depth / 30.f : 0.f)) {
			isFiltered = true;
			filterStats->SCBZ++;
		}
	}

	if (present(MQ0Fval) && MQ0Fval > opts->minMQ0F) {
		isFiltered = true;
		filterStats->MQ0Ffilt++;
	}
	if (depthFailed) {
		isFiltered = true;
		filterStats->DP++;
	}
	return !isFiltered;

#if 0
	isFiltered = false; //altFreq = 0.f;
	if (alt == "." || alt == "*" || alt == "<*>") { //nothing needed here...
		isIndel = false; isSnp = false; isFiltered = false; return false;
	}
	//not used
	//if (VF->filter(fields) ) {isFiltered = true;}
	//filter2 based on the INFO field



	//assert(DP4[2] + DP4[3] + DP4[1] + DP4[0] == (int)DPval); //is actually not always the same..
	//two ways of calculation altFreq
	/*if (AF1val > 0.f) {
		altFreq = AF1val;
	}
	else if ((DP4[2] + DP4[3] + DP4[1] + DP4[0]) > 0) {
		altFreq = float(DP4[2] + DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0]);
	}
	assert(altFreq >= 0);
	*/
	
	if (isIndel) {
		//******************************* INDEL
		//https://www.htslib.org/workflow/filter.html
		// Historical filter sketch retained only in this disabled block.
		// DP > adaptive high-depth threshold; reject poor MQ/read-position support.
		if (IDVval < 2.f) { isFiltered = true;;		}
		if (IMFval < 0.1f) { isFiltered = true;; } //skip indels with high IMF
		if (RPBZval > 6.f || SCBZval > 6.f) { isFiltered = true;; } //skip indels with high RPBZ
		if (RPBZval + SCBZval > 9.f) { isFiltered = true;; }
		if (MQBZval < -(5.f + DPval / 20.f)) { isFiltered = true;; } //skip indels with high MQBZ
		if (QUALval < opts->minCallQual) { isFiltered = true; }
		/*if (!isFiltered) {
			int a = 0;//DEBUG
		}*/
	
	
	}	else if (ref.length() == 1 && alt.length() == 1) {
		//******************************* SNP
		//assert(ref.length() == 1); assert(alt.length() == 1);
		isSnp = true;

		//indel proximity filtering
		if ( (posN - opts->indelRange < lastIndelPos && posN + opts->indelRange > lastIndelPos) ||
			(posN - opts->indelRange < nextIndelPos && posN + opts->indelRange > nextIndelPos)){
			isFiltered = true; filterStats->indelProx++;
		}


		//QUAL from SNP caller
		
		if (opts->minCallQualAdaptive>0.f && DPval > 2.f &&
			QUALval < min(20.f,float(mD* opts->minCallQualAdaptive)) ) {
			isFiltered = true; filterStats->QUALadaptive++; 
		} 
		if (DPval > 2.f &&  QUALval < (float)opts->minCallQual) {
			isFiltered = true; filterStats->QUAL++; 
		}
		if (DPval < min(8.f,float(mD*opts->depthFilterScale))) { isFiltered = true; filterStats->DP++; }

		//mean depth based filtering:
/*		if (mD >= 20.f) {
			//'QUAL<30 || DP<8 || MQ<40 || SB>40'
			if ( DPval < 8.f) { isFiltered = true; filterStats->DP++; }

		} else if (mD >= 5.f) {
			//'QUAL<20 || DP<3 || MQ<40 || SB>40'
			//if (opts->minCallQualAdaptive && QUALval < 10.f) { isFiltered = true; filterStats->QUALadaptive++; }
			if (DPval < 3.f) { isFiltered = true; filterStats->DP++; }

		}	else {
			//'QUAL<20 || DP<2 || MQ<50 || SB>30'
			//if (opts->minCallQualAdaptive && QUALval < 5.f) { isFiltered = true; filterStats->QUALadaptive++; }
		}
*/
		//strand bias
		if (SPval > (40.f + DPval / 2.f)) { isFiltered = true; filterStats->SP++;} //minSP
		//SNP base quality
		if (abs(BQBZval) > (3.1f + DPval / 40.f)) { isFiltered = true; filterStats->BQBZ++;} //prev minBQBZ, now following htslib rec

		//ok this seems to be a ntVariant 
		//bcftools rec filter: MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 32 || SCBZ > 3
		//or:
		//DP>2*$DP || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)
		if (MQBZval < -3.f) { isFiltered = true;; filterStats->MQBZ++;} //skip snps with high MQBZ
		//position in read
		if (abs(RPBZval) > 3.5f) { isFiltered = true;; filterStats->RPBZ++;} //skip snps with high RPBZ
		//soft clips
		if (SCBZval > 2.f + DPval / 30.f) { isFiltered = true;; filterStats->SCBZ++;} //skip snps with high SCBZ
	}

	//custom Filter
	if (MQ0Fval > opts->minMQ0F) { isFiltered = true; filterStats->MQ0Ffilt++; }
	if (FSval < opts->minFS) { isFiltered = true; filterStats->FSfilt++; }
	
	
	if (altFreq >= 0.95f) {//can't imagine a condition where this is not the major allele
		isFiltered = false; 
	}

	return !isFiltered;
#endif
}

bool VCFmem::evalVCFentry(options* opts, VCFmem* v2, float mD, float mD2,
	int lastIndelPos, int nextIndelPos, VCFfilterStats* filterStats) {
	bool pass1 = evalVCFentry(opts, mD, lastIndelPos, nextIndelPos, filterStats);
	bool pass2 = v2->evalVCFentry(opts, mD2, lastIndelPos, nextIndelPos, filterStats);
	return reconcileEvaluated(v2, pass1, pass2);
}

bool VCFmem::reconcileEvaluated(VCFmem* v2, bool pass1, bool pass2) {
	bool call1 = !isNoCall();
	bool call2 = !v2->isNoCall();

	if (!call1 && !call2) {
		isFiltered = false;
		isIndel = false;
		isSnp = false;
		return false;
	}
	if (!call1) {
		copyCallFrom(*v2);
		isFiltered = !pass2;
		return pass2;
	}
	if (!call2) {
		isFiltered = !pass1;
		return pass1;
	}

	if (ref == v2->ref && alt == v2->alt) {
		if (pass1 == pass2) {
			float weight1 = evidenceDepth();
			float weight2 = v2->evidenceDepth();
			if (present(weight1) && present(weight2) && weight1 + weight2 > 0.f) {
				altFreq = (altFreq * weight1 + v2->altFreq * weight2) / (weight1 + weight2);
				if (altDepth >= 0 && v2->altDepth >= 0) {
					if (altDepth > std::numeric_limits<long long>::max() - v2->altDepth) {
						throw std::runtime_error("combined alternate depth overflows a 64-bit integer");
					}
					altDepth += v2->altDepth;
				} else {
					altDepth = -1;
				}
				if (totalAlleleDepth > 0 && v2->totalAlleleDepth > 0) {
					if (totalAlleleDepth > std::numeric_limits<long long>::max() - v2->totalAlleleDepth) {
						throw std::runtime_error("combined total depth overflows a 64-bit integer");
					}
					totalAlleleDepth += v2->totalAlleleDepth;
				} else {
					totalAlleleDepth = 0;
				}
			} else if (!present(weight1) && present(weight2)) {
				copyCallFrom(*v2);
			} else if (!present(weight1) && !present(weight2) && v2->altFreq > altFreq) {
				copyCallFrom(*v2);
			}
		} else if (!pass1 && pass2) {
			copyCallFrom(*v2);
		}
		// Agreement is useful evidence, but it never rescues two failed calls.
		isFiltered = !(pass1 || pass2);
		return !isFiltered;
	}

	if (!pass1 && pass2) {
		copyCallFrom(*v2);
		conflict = false;
		isFiltered = false;
		return true;
	}
	if (pass1 && !pass2) {
		conflict = false;
		isFiltered = false;
		return true;
	}
	if (!pass1 && !pass2) {
		const float firstDepth = evidenceDepth();
		const float secondDepth = v2->evidenceDepth();
		const bool preferSecond =
			(present(secondDepth) && (!present(firstDepth) || secondDepth > firstDepth)) ||
			((!present(firstDepth) || !present(secondDepth) || secondDepth == firstDepth) &&
				(v2->altFreq > altFreq || (v2->altFreq == altFreq && v2->alt < alt)));
		if (preferSecond) copyCallFrom(*v2);
		conflict = false;
		conflictSpan = 0;
		isFiltered = true;
		return false;
	}

	// Two passing but discordant alleles are not arbitrarily assigned to a
	// platform. Preserve the stronger call for diagnostics and mark the site
	// conflicted so the consensus layer masks it.
	const size_t combinedSpan = max(ref.size(), v2->ref.size());
	const float firstDepth = evidenceDepth();
	const float secondDepth = v2->evidenceDepth();
	const bool depthsTie = (present(firstDepth) && present(secondDepth) && secondDepth == firstDepth) ||
		(!present(firstDepth) && !present(secondDepth));
	const bool preferSecond = v2->ref.size() > ref.size() ||
		(v2->ref.size() == ref.size() &&
			((present(secondDepth) && (!present(firstDepth) || secondDepth > firstDepth)) ||
			 (depthsTie && (v2->altFreq > altFreq ||
				(v2->altFreq == altFreq && v2->alt < alt)))));
	if (preferSecond) {
		copyCallFrom(*v2);
	}
	conflict = true;
	conflictSpan = combinedSpan;
	isFiltered = false;
	return false;

#if 0
	//in this case the current vcf is converted to a "combined" vcf entry..
	isFiltered = false; //altFreq = 0.f;
	bool isF2 = false;//for second read to later decide whom to trust..
	if ((alt == "." || alt == "*" || alt == "<*>") && 
			(v2->alt == "." || v2->alt == "*" || v2->alt == "<*>")) { //nothing needed here...
		isIndel = false; isSnp = false; isFiltered = false; return false;
	}

	//if (VF->filter(fields) && VF->filter(v2->fields)) {isFiltered = true;}
	//combine relevant entries
	//DPval += v2->DPval;
	//SPval = (SPval + v2->SPval) / 2;
	//no eval single entries..

	//filter2 based on the INFO field
	//strand bias
	
	if (abs(BQBZval) > (3.1f + DPval / 40.f)) { isFiltered = true; } //prev minBQBZ, now following htslib rec
	if (abs(v2->BQBZval) > (3.1f + v2->DPval / 40.f)) { isF2 = true; } //prev minBQBZ, now following htslib rec

	//custom Filter
	if (MQ0Fval > opts->minMQ0F) { isFiltered = true; }
	if (FSval < opts->minFS) { isFiltered = true; }
	if (v2->MQ0Fval > opts->minMQ0F) { isF2 = true; }
	if (v2->FSval < opts->minFS) { isF2 = true; }


	//assert(DP4[2] + DP4[3] + DP4[1] + DP4[0] == (int)DPval); //is actually not always the same..
	//two ways of calculation altFreq
	/*if (AF1val > 0.f) {
		altFreq = AF1val;
	} else if ((DP4[2] + DP4[3] + DP4[1] + DP4[0]) > 0) {
		altFreq = float(DP4[2] + DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0]);
	}
	assert(altFreq >= 0);
	*/

	//combine freqs
    float cFreq(0.f);
	if ((v2->DPval + DPval) > 0.f) {
		cFreq = (AF1val * DPval + v2->AF1val * v2->DPval) / (v2->DPval + DPval);
		}
	float cFreq2(float(DP4[2] + DP4[3]+ v2->DP4[2] + v2->DP4[3]) / float(DP4[2] + DP4[3] + DP4[1] + DP4[0] + v2->DP4[2] + v2->DP4[3] + v2->DP4[1] + v2->DP4[0]));



	if (isIndel) {
		//https://www.htslib.org/workflow/filter.html
		// Historical filter sketch retained only in this disabled block.
		// DP > adaptive high-depth threshold; reject poor MQ/read-position support.
		if (IMFval > 0.1f) { isFiltered = true;; } //skip indels with high IMF
		if (RPBZval > 6.f || SCBZval > 6.f) { isFiltered = true;; } //skip indels with high RPBZ
		if (RPBZval + SCBZval > 9.f) { isFiltered = true;; }
		if (MQBZval < -(5.f + DPval / 20.f)) { isFiltered = true;; } //skip indels with high MQBZ
		if (QUALval < opts->minCallQual) { isFiltered = true; }
	}
	else if (ref.length() == 1 && alt.length() == 1){
		//assert(ref.length() == 1); assert(alt.length() == 1);
		isSnp = true;
		if (SPval > (40.f + DPval / 2.f)) { isFiltered = true; } //minSP
		if (v2->SPval > (40.f + v2->DPval / 2.f)) { isF2 = true; } //minSP

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
		if (v2->SCBZval > 2.f + v2->DPval / 30.f) { isF2 = true;; } //skip snps with high SCBZ
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
		}else {//trust because it iof ths the same alt allele
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
#endif
}

VCFmem::VCFmem(const string& line, int inputSource):
	chrom(""), id(""), ref(""), alt(""), filter(""),
	fieldsSet(false), QUALval(missingFloat()), posN(-1), sourceFile(inputSource),
	altDepth(-1), totalAlleleDepth(0), DP4(4, 0), altFreq(0.f),
	AF1val(missingFloat()), AF2val(missingFloat()), FSval(missingFloat()),
	MQ0Fval(missingFloat()), BQBZval(missingFloat()), SPval(missingFloat()),
	IDVval(missingFloat()), IMFval(missingFloat()), DPval(missingFloat()),
	RPBZval(missingFloat()), SCBZval(missingFloat()), MQBZval(missingFloat()),
	isIndel(false), isSnp(false),
	isFiltered(false), isUnsure(false),
	conflict(false), conflictSpan(0), VF(new vcfFields(inputSource))

{
	vector<string> columns;
	stringstream row(line);
	string column;
	while (getline(row, column, '\t')) columns.push_back(column);
	if (columns.size() < 8) {
		throw std::runtime_error("expected at least eight tab-delimited VCF columns");
	}
	if (columns.size() > 10) {
		throw std::runtime_error("multiple sample columns are not supported; provide one pooled/metagenomic sample per VCF");
	}
	if (columns.size() == 9) {
		throw std::runtime_error("FORMAT is present without its required sample column");
	}
	chrom = columns[0];
	int vcfPos = 0;
	if (chrom.empty() || !parseIntValue(columns[1], vcfPos) || vcfPos < 1) {
		throw std::runtime_error("invalid CHROM or POS");
	}
	posN = vcfPos - 1;
	id = columns[2];
	ref = columns[3];
	filter = columns[6];
	if (!isSequenceAllele(ref)) throw std::runtime_error("REF must contain only A/C/G/T/N bases");
	if (ref.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
		ref.size() - 1 > static_cast<size_t>(std::numeric_limits<int>::max() - posN)) {
		throw std::runtime_error("REF span exceeds the supported coordinate range");
	}
	transform(ref.begin(), ref.end(), ref.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
	if (columns[5] != "." && !parseFloatValue(columns[5], QUALval)) {
		throw std::runtime_error("invalid QUAL value");
	}
	if (present(QUALval) && QUALval < 0.f) throw std::runtime_error("QUAL must be non-negative");
	if (filter.empty()) throw std::runtime_error("FILTER is empty (use PASS or . when appropriate)");
	parseINFO(columns[7]);
	if (columns.size() >= 10) parseFormat(columns[8], columns[9]);
	if (present(DPval) && DPval < 0.f) throw std::runtime_error("DP must be non-negative");
	if (present(IDVval) && IDVval < 0.f) throw std::runtime_error("IDV must be non-negative");
	if (present(IMFval) && (IMFval < 0.f || IMFval > 1.f)) throw std::runtime_error("IMF must be in [0,1]");
	if (present(AF1val) && (AF1val < 0.f || AF1val > 1.f)) throw std::runtime_error("AF must be in [0,1]");
	if (present(AF2val) && (AF2val < 0.f || AF2val > 1.f)) throw std::runtime_error("AF2 must be in [0,1]");
	if (present(MQ0Fval) && (MQ0Fval < 0.f || MQ0Fval > 1.f)) throw std::runtime_error("MQ0F must be in [0,1]");
	if (present(SPval) && SPval < 0.f) throw std::runtime_error("SP must be non-negative");

	vector<string> alts;
	stringstream altStream(columns[4]);
	while (getline(altStream, column, ',')) alts.push_back(column);
	if (alts.empty()) throw std::runtime_error("ALT is empty");
	chooseAllele(alts);
	return;

#if 0

	VF = new vcfFields( 0);
	// Parse a single VCF line
	std::istringstream ss(line);

	ss >> chrom >> posN >> id >> ref >> alt;
	posN--;//counted differently in vcf to C++
	ss >> QUALval;
	string info, format, xtra;
	//first check if the filter is set to "PASS" or not, if not, skip this line
	ss >> filter >> info >> format >> xtra;
	if (info.substr(0, 5) == "INDEL" && (alt.length()!=1 || ref.length()!=1 || alt == "<DEL>" || alt == "<INS>")) {
		isIndel = true;// indelCNT++;
	}
	//cerr << "U"; 
	if (!VF->isSet()) { VF-> parseFields(format); }
	//cerr << "I"; 
	VF->splitXtra(xtra);
	//cerr << "O"; 
    this->parseINFO(info);

	altFreq = 0.f;


	//check if ADfield > 2, means multi allele..
	if (VF->ADfield.size() > 2 ) {
		//cerr << "Warning: multi-allelic site found, currently not supported. Skipping site.\n" << line << endl;
		string segment; stringstream altT(alt);
		vector<string> alts; int maxAltCnt = 0; int maxAltIdx = 0; int j = -1;
		while (getline(altT, segment, ',')) {
			j++;
			alts.push_back(segment);
			if (VF->ADfield[j+1]> maxAltCnt) {
				maxAltCnt = VF->ADfield[j+1];
				maxAltIdx = j;
			}
		}
		alt = alts[maxAltIdx];
		altFreq = float(maxAltCnt)/float((maxAltCnt+VF->ADfield[0]));

		return;

	}


	// compute initial alt frequency so stored VCF entries have a usable value
	int sumDP4 = DP4[0] + DP4[1] + DP4[2] + DP4[3];
	if (sumDP4 > 0) {
		altFreq = float(DP4[2] + DP4[3]) / float(sumDP4);
	}else if (AF1val >= 0.f) {
		altFreq = AF1val;
	} else {
		altFreq = 0.f;
	}

	//cerr << "P"; 

	return;
#endif
}


void VCFmem::parseFormat(const string& format, const string& sample) {
	if (format.empty() || format == "." || sample.empty() || sample == ".") return;
	VF->parseFields(format);
	if (!VF->splitXtra(sample)) {
		throw std::runtime_error("FORMAT/sample field cardinality or numeric value is invalid");
	}
	auto readFloatField = [&](int index, float& target, const char* name) {
		if (index < 0 || static_cast<size_t>(index) >= VF->fields.size()) return;
		const string& value = VF->fields[static_cast<size_t>(index)];
		if (value.empty() || value == ".") return;
		if (!parseFloatValue(value, target)) throw std::runtime_error(string("invalid FORMAT/") + name);
	};
	readFloatField(VF->DP, DPval, "DP");
	readFloatField(VF->SP, SPval, "SP");
	readFloatField(VF->BQBZ, BQBZval, "BQBZ");
	readFloatField(VF->RPBZ, RPBZval, "RPBZ");
	readFloatField(VF->SCBZ, SCBZval, "SCBZ");
	readFloatField(VF->MQBZ, MQBZval, "MQBZ");
	readFloatField(VF->IDV, IDVval, "IDV");
	readFloatField(VF->IMF, IMFval, "IMF");
	readFloatField(VF->MQ0F, MQ0Fval, "MQ0F");

	if (VF->ADfield.empty() && VF->ADF >= 0 && VF->ADR >= 0 &&
		static_cast<size_t>(VF->ADF) < VF->fields.size() && static_cast<size_t>(VF->ADR) < VF->fields.size()) {
		auto parseDepthVector = [](const string& text, vector<int>& values) {
			if (text.empty() || text == ".") return false;
			stringstream stream(text);
			string token;
			while (getline(stream, token, ',')) {
				int value = 0;
				if (!parseIntValue(token, value) || value < 0) return false;
				values.push_back(value);
			}
			return !values.empty();
		};
		vector<int> forward;
		vector<int> reverse;
		const string& forwardText = VF->fields[static_cast<size_t>(VF->ADF)];
		const string& reverseText = VF->fields[static_cast<size_t>(VF->ADR)];
		const bool forwardMissing = forwardText.empty() || forwardText == ".";
		const bool reverseMissing = reverseText.empty() || reverseText == ".";
		if (forwardMissing != reverseMissing) {
			throw std::runtime_error("FORMAT/ADF and FORMAT/ADR must both be present or both be missing");
		}
		if (!forwardMissing &&
			(!parseDepthVector(forwardText, forward) || !parseDepthVector(reverseText, reverse) ||
			 forward.size() != reverse.size())) {
			throw std::runtime_error("FORMAT/ADF and FORMAT/ADR must contain matching non-negative depth vectors");
		}
		if (!forward.empty()) {
			VF->ADfield.resize(forward.size());
			for (size_t i = 0; i < forward.size(); ++i) {
				if (forward[i] > std::numeric_limits<int>::max() - reverse[i]) {
					throw std::runtime_error("FORMAT/ADF+ADR depth overflows an integer");
				}
				VF->ADfield[i] = forward[i] + reverse[i];
			}
		}
	}
}

void VCFmem::chooseAllele(const vector<string>& rawAlts) {
	vector<string> alts = rawAlts;
	for (string& allele : alts) {
		if (isSequenceAllele(allele)) {
			transform(allele.begin(), allele.end(), allele.begin(),
				[](unsigned char c) { return static_cast<char>(std::toupper(c)); });
		}
	}

	size_t selected = 0;
	if (!VF->ADfield.empty()) {
		if (VF->ADfield.size() != alts.size() + 1) {
			throw std::runtime_error("AD must contain one value for REF and every ALT allele");
		}
		for (int value : VF->ADfield) {
			if (value < 0) throw std::runtime_error("AD values must be non-negative");
			if (totalAlleleDepth > std::numeric_limits<long long>::max() - value) {
				throw std::runtime_error("AD depth sum overflows a 64-bit integer");
			}
			totalAlleleDepth += value;
		}
		bool foundSequenceAlt = false;
		for (size_t i = 0; i < alts.size(); ++i) {
			if (!isSequenceAllele(alts[i])) continue;
			if (!foundSequenceAlt || VF->ADfield[i + 1] > VF->ADfield[selected + 1]) {
				selected = i;
				foundSequenceAlt = true;
			}
		}
		if (!foundSequenceAlt) {
			for (size_t i = 1; i < VF->ADfield.size(); ++i) {
				if (VF->ADfield[i] > VF->ADfield[selected + 1]) selected = i - 1;
			}
		}
		altDepth = VF->ADfield[selected + 1];
	} else if (alts.size() > 1) {
		throw std::runtime_error("multi-allelic records require AD or ADF+ADR evidence");
	}

	alt = alts[selected];
	if (alt.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::runtime_error("selected ALT exceeds the supported coordinate range");
	}
	if (totalAlleleDepth > 0 && altDepth >= 0) {
		altFreq = static_cast<float>(altDepth) / static_cast<float>(totalAlleleDepth);
		if (!present(DPval)) DPval = static_cast<float>(totalAlleleDepth);
	} else {
		long long dp4Total = 0;
		for (int value : DP4) dp4Total += value;
		if (dp4Total > 0) {
			altDepth = static_cast<long long>(DP4[2]) + static_cast<long long>(DP4[3]);
			totalAlleleDepth = dp4Total;
			altFreq = static_cast<float>(altDepth) / static_cast<float>(dp4Total);
			if (!present(DPval)) DPval = static_cast<float>(dp4Total);
		} else if (present(IDVval) || present(IMFval)) {
			if (present(IMFval)) altFreq = IMFval;
			if (present(IDVval)) {
				if (static_cast<long double>(IDVval) >
					static_cast<long double>(std::numeric_limits<long long>::max())) {
					throw std::runtime_error("IDV is too large");
				}
				altDepth = std::llround(IDVval);
			}
			if (present(DPval) && DPval > 0.f) {
				if (static_cast<long double>(DPval) >
					static_cast<long double>(std::numeric_limits<long long>::max())) {
					throw std::runtime_error("DP is too large");
				}
				totalAlleleDepth = std::llround(DPval);
				if (!present(IMFval) && altDepth >= 0) {
					altFreq = static_cast<float>(altDepth) / DPval;
				}
				if (!present(IDVval) && present(IMFval)) {
					altDepth = std::llround(IMFval * DPval);
				}
			}
		} else if (present(AF1val)) {
			altFreq = AF1val;
			if (present(DPval) && DPval > 0.f) {
				if (static_cast<long double>(DPval) >
					static_cast<long double>(std::numeric_limits<long long>::max())) {
					throw std::runtime_error("DP is too large");
				}
				totalAlleleDepth = std::llround(DPval);
				altDepth = std::llround(altFreq * DPval);
			}
		}
	}
	if (altFreq < 0.f || altFreq > 1.f) throw std::runtime_error("allele frequency is outside [0,1]");
	if (!present(IDVval) && altDepth >= 0) IDVval = static_cast<float>(altDepth);
	if (!present(IMFval) && altFreq >= 0.f) IMFval = altFreq;

	isIndel = false;
	isSnp = false;
	if (isNoCall()) return;
	if (!isSequenceAllele(alt)) return; // symbolic/breakend alleles are not sequence-resolved here
	if (ref == alt) throw std::runtime_error("selected ALT allele is identical to REF");
	if (ref.size() == 1 && alt.size() == 1) isSnp = true;
	else isIndel = true; // internal replacement path also handles sequence-resolved MNPs
}

void VCFmem::copyCallFrom(const VCFmem& other) {
	ref = other.ref;
	alt = other.alt;
	filter = other.filter;
	QUALval = other.QUALval;
	DP4 = other.DP4;
	altFreq = other.altFreq;
	AF1val = other.AF1val;
	AF2val = other.AF2val;
	FSval = other.FSval;
	MQ0Fval = other.MQ0Fval;
	BQBZval = other.BQBZval;
	SPval = other.SPval;
	IDVval = other.IDVval;
	IMFval = other.IMFval;
	DPval = other.DPval;
	RPBZval = other.RPBZval;
	SCBZval = other.SCBZval;
	MQBZval = other.MQBZval;
	sourceFile = other.sourceFile;
	altDepth = other.altDepth;
	totalAlleleDepth = other.totalAlleleDepth;
	isIndel = other.isIndel;
	isSnp = other.isSnp;
	isUnsure = other.isUnsure;
}


void VCFmem::parseINFO(const string& info) {
	if (info.empty() || info == ".") return;
	string segment;
	stringstream fieldsStream(info);
	while (getline(fieldsStream, segment, ';')) {
		if (segment.empty()) continue;
		size_t equals = segment.find('=');
		string key = equals == string::npos ? segment : segment.substr(0, equals);
		string value = equals == string::npos ? "" : segment.substr(equals + 1);
		if (value == ".") continue;
		if (key == "DP4") {
			vector<int> values;
			stringstream list(value);
			string token;
			while (getline(list, token, ',')) {
				int parsed = 0;
				if (!parseIntValue(token, parsed) || parsed < 0) throw std::runtime_error("invalid INFO/DP4");
				values.push_back(parsed);
			}
			if (values.size() != 4) throw std::runtime_error("INFO/DP4 must contain exactly four counts");
			DP4 = values;
		} else if (key == "AD") {
			vector<int> values;
			stringstream list(value);
			string token;
			while (getline(list, token, ',')) {
				int parsed = 0;
				if (!parseIntValue(token, parsed) || parsed < 0) throw std::runtime_error("invalid INFO/AD");
				values.push_back(parsed);
			}
			if (!values.empty()) VF->ADfield = values;
		} else {
			float parsed = missingFloat();
			string firstValue = value.substr(0, value.find(','));
			if (key == "FS") {
				if (!parseFloatValue(value, FSval)) throw std::runtime_error("invalid INFO/FS");
			} else if (key == "AF1" || key == "AF") {
				if (!parseFloatValue(firstValue, AF1val)) throw std::runtime_error("invalid INFO/AF");
			} else if (key == "AF2") {
				if (!parseFloatValue(firstValue, AF2val)) throw std::runtime_error("invalid INFO/AF2");
			} else if (key == "DP") {
				if (!parseFloatValue(value, DPval)) throw std::runtime_error("invalid INFO/DP");
			} else if (key == "MQ0F") {
				if (!parseFloatValue(value, MQ0Fval)) throw std::runtime_error("invalid INFO/MQ0F");
			} else if (key == "BQBZ") {
				if (!parseFloatValue(value, BQBZval)) throw std::runtime_error("invalid INFO/BQBZ");
			} else if (key == "SP") {
				if (!parseFloatValue(value, SPval)) throw std::runtime_error("invalid INFO/SP");
			} else if (key == "IDV") {
				if (!parseFloatValue(value, IDVval)) throw std::runtime_error("invalid INFO/IDV");
			} else if (key == "IMF") {
				if (!parseFloatValue(value, IMFval)) throw std::runtime_error("invalid INFO/IMF");
			} else if (key == "RPBZ") {
				if (!parseFloatValue(value, RPBZval)) throw std::runtime_error("invalid INFO/RPBZ");
			} else if (key == "SCBZ") {
				if (!parseFloatValue(value, SCBZval)) throw std::runtime_error("invalid INFO/SCBZ");
			} else if (key == "MQBZ") {
				if (!parseFloatValue(value, MQBZval)) throw std::runtime_error("invalid INFO/MQBZ");
			}
			(void)parsed;
		}
	}
	return;

#if 0
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
#endif
}


bool vcfFields::splitXtra(const string& sample) {
	fields.assign(static_cast<size_t>(numFields), ".");
	string segment;
	stringstream input(sample);
	size_t index = 0;
	while (getline(input, segment, ':')) {
		if (index >= fields.size()) return false;
		fields[index++] = segment.empty() ? "." : segment;
	}

	// FORMAT/AD has precedence over INFO/AD, but a missing FORMAT value must
	// not erase a usable INFO/AD field.
	if (AD >= 0 && static_cast<size_t>(AD) < index) {
		const string& ad = fields[static_cast<size_t>(AD)];
		if (!ad.empty() && ad != ".") {
			vector<int> parsedDepths;
			stringstream values(ad);
			string value;
			while (getline(values, value, ',')) {
				int parsed = 0;
				if (!parseIntValue(value, parsed) || parsed < 0) return false;
				parsedDepths.push_back(parsed);
			}
			if (parsedDepths.empty()) return false;
			ADfield.swap(parsedDepths);
		}
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
	for (size_t i = 0; i < seglist.size(); ++i) {
		if (seglist[i].empty()) throw std::runtime_error("empty FORMAT field name");
		const int fieldIndex = static_cast<int>(i);
		if (seglist[i] == "GT") GT = fieldIndex;
		else if (seglist[i] == "MQ0F") MQ0F = fieldIndex;
		else if (seglist[i] == "IDV") IDV = fieldIndex;
		else if (seglist[i] == "IMF") IMF = fieldIndex;
		else if (seglist[i] == "BQBZ") BQBZ = fieldIndex;
		else if (seglist[i] == "RPBZ") RPBZ = fieldIndex;
		else if (seglist[i] == "SCBZ") SCBZ = fieldIndex;
		else if (seglist[i] == "MQBZ") MQBZ = fieldIndex;
		else if (seglist[i] == "DP") DP = fieldIndex;
		else if (seglist[i] == "SP") SP = fieldIndex;
		else if (seglist[i] == "ADF") ADF = fieldIndex;
		else if (seglist[i] == "ADR") ADR = fieldIndex;
		else if (seglist[i] == "AD") AD = fieldIndex;
	}
	if (seglist.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::runtime_error("too many FORMAT fields");
	}
	numFields = static_cast<int>(seglist.size());
	fieldsSet = true;
	fields.resize(seglist.size());
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
	if (posSrch == curMems->end() || posSrch->second.empty()) {//pos not called in this vcf..
		return(nullptr);
	}
	return(posSrch->second.front());
}

vector<VCFmem*> VCFcollection::takeContig(const string& cont) {
	vector<VCFmem*> calls;
	auto found = VMems.find(cont);
	if (found == VMems.end()) return calls;
	map2vcfm* stored = found->second;
	if (stored != nullptr) {
		size_t callCount = 0;
		for (const auto& entry : *stored) callCount += entry.second.size();
		calls.reserve(callCount);
		for (auto& entry : *stored) {
			for (VCFmem* call : entry.second) calls.push_back(call);
			entry.second.clear();
		}
		delete stored;
	}
	VMems.erase(found);
	if (curContig == cont) {
		curContig.clear();
		curMems = nullptr;
	}
	return calls;
}

vector<string> VCFcollection::contigs() const {
	vector<string> result;
	result.reserve(VMems.size());
	for (const auto& entry : VMems) result.push_back(entry.first);
	sort(result.begin(), result.end());
	return result;
}

void VCFmulti::evalVCFs(fasta* cF, VCFcollection* VCF2ptr, 
	VariantStats* Vstats, VCFfilterStats* filterStats,options* opts) {
	if (mems.empty()) return;
	const string curChrom = mems.front()->getChrom();
	const size_t primaryCount = mems.size();
	vector<std::unique_ptr<VCFmem>> owned;
	owned.reserve(primaryCount);
	for (VCFmem* call : mems) owned.emplace_back(call);
	mems.clear();
	if (VCF2ptr != nullptr) {
		for (VCFmem* call : VCF2ptr->takeContig(curChrom)) owned.emplace_back(call);
		if (cF != nullptr && cF->getSequenceId() != curChrom) {
			for (VCFmem* call : VCF2ptr->takeContig(cF->getSequenceId())) owned.emplace_back(call);
		}
	}
	if (cF == nullptr) return;

	map<int, vector<VCFmem*> > primary;
	map<int, vector<VCFmem*> > secondary;
	for (size_t i = 0; i < owned.size(); ++i) {
		VCFmem* call = owned[i].get();
		map<int, vector<VCFmem*> >& destination = i < primaryCount ? primary : secondary;
		destination[call->getPos()].push_back(call);
	}

	vector<int> positions;
	positions.reserve(primary.size() + secondary.size());
	for (const auto& entry : primary) positions.push_back(entry.first);
	for (const auto& entry : secondary) positions.push_back(entry.first);
	sort(positions.begin(), positions.end());
	positions.erase(unique(positions.begin(), positions.end()), positions.end());

	// Check every sequence-resolved candidate before it can influence proximity,
	// duplicate selection, or two-caller reconciliation.
	set<VCFmem*> invalidReferenceCalls;
	auto validateCandidates = [&](const map<int, vector<VCFmem*> >& calls) {
		for (const auto& entry : calls) {
			for (VCFmem* call : entry.second) {
				if (call != nullptr && (call->isSNP() || call->isINDEL()) &&
					!cF->validateVariantReference(call, Vstats, opts)) {
					invalidReferenceCalls.insert(call);
				}
			}
		}
	};
	validateCandidates(primary);
	validateCandidates(secondary);

	// SNP artefacts can occur near either edge of a deletion. Treat every REF
	// base covered by an indel as an indel position, not just its anchor.
	vector<pair<int, int> > indelPositions; // covered base, anchor position
	auto recordIndelSpan = [&indelPositions, &invalidReferenceCalls](VCFmem* call) {
		if (call == nullptr || !call->isINDEL() || invalidReferenceCalls.count(call) != 0) return;
		const int span = max(1, static_cast<int>(call->getRef().size()));
		for (int offset = 0; offset < span; ++offset) {
			indelPositions.push_back(make_pair(call->getPos() + offset, call->getPos()));
		}
	};
	for (const auto& entry : primary) for (VCFmem* call : entry.second) recordIndelSpan(call);
	for (const auto& entry : secondary) for (VCFmem* call : entry.second) recordIndelSpan(call);
	sort(indelPositions.begin(), indelPositions.end());
	indelPositions.erase(unique(indelPositions.begin(), indelPositions.end()), indelPositions.end());
	vector<VCFmem*> evaluatedResults;
	evaluatedResults.reserve(positions.size());

	for (int position : positions) {
		int previousIndel = -1;
		int nextIndel = -1;
		auto atOrAfter = lower_bound(indelPositions.begin(), indelPositions.end(),
			make_pair(position, std::numeric_limits<int>::min()));
		for (auto right = atOrAfter; right != indelPositions.end(); ++right) {
			if (right->second != position) {
				nextIndel = right->first;
				break;
			}
		}
		auto left = atOrAfter;
		while (left != indelPositions.begin()) {
			--left;
			if (left->second != position) {
				previousIndel = left->first;
				break;
			}
		}

		auto selectCall = [&](const vector<VCFmem*>& candidates) {
			VCFmem* best = nullptr;
			bool bestPass = false;
			for (VCFmem* candidate : candidates) {
				if (candidate == nullptr || (!candidate->isSNP() && !candidate->isINDEL())) continue;
				if (invalidReferenceCalls.count(candidate) != 0) continue;
				const float mean = cF->getAvgDepth(static_cast<size_t>(candidate->getSourceFile()));
				const bool pass = candidate->evalVCFentry(opts, mean,
					previousIndel, nextIndel, filterStats);
				if (best == nullptr || (pass && !bestPass) ||
					(pass == bestPass && candidate->getFreq() > best->getFreq()) ||
					(pass == bestPass && candidate->getFreq() == best->getFreq() &&
						best->isINDEL() && candidate->isSNP())) {
					best = candidate;
					bestPass = pass;
				}
			}
			return make_pair(best, bestPass);
		};

		pair<VCFmem*, bool> firstChoice(nullptr, false);
		pair<VCFmem*, bool> secondChoice(nullptr, false);
		auto firstIt = primary.find(position);
		if (firstIt != primary.end()) firstChoice = selectCall(firstIt->second);
		auto secondIt = secondary.find(position);
		if (secondIt != secondary.end()) secondChoice = selectCall(secondIt->second);

		VCFmem* result = firstChoice.first != nullptr ? firstChoice.first : secondChoice.first;
		if (firstChoice.first != nullptr && secondChoice.first != nullptr) {
			firstChoice.first->reconcileEvaluated(secondChoice.first,
				firstChoice.second, secondChoice.second);
		}
		if (result != nullptr) evaluatedResults.push_back(result);
	}

	// Replacements at different anchors may nevertheless overlap. Applying both
	// is not a defined haplotype operation, so mask the complete reference union
	// rather than guessing an order-dependent consensus.
	vector<bool> suppressed(evaluatedResults.size(), false);
	size_t activeRepresentative = evaluatedResults.size();
	int64_t activeStart = -1;
	int64_t activeEnd = -1;
	auto affectsConsensus = [opts, cF](VCFmem* call) {
		if (call == nullptr) return false;
		if (!cF->isDepthResolved(call->getPos())) return false;
		if (call->conflicted()) return true;
		if (call->isINDEL()) {
			if (call->majorAllele()) return call->filtered() || opts->reportINDELs;
			return !call->filtered() && opts->maskMinorAllele;
		}
		if (call->isSNP()) {
			if (call->majorAllele()) return true;
			return !call->filtered() && opts->maskMinorAllele;
		}
		return false;
	};
	for (size_t i = 0; i < evaluatedResults.size(); ++i) {
		VCFmem* call = evaluatedResults[i];
		if (!affectsConsensus(call)) continue;
		const int64_t start = call->getPos();
		const int64_t end = start + static_cast<int64_t>(call->getUncertainRefLength());
		if (activeRepresentative == evaluatedResults.size() || start >= activeEnd) {
			activeRepresentative = i;
			activeStart = start;
			activeEnd = end;
			continue;
		}

		VCFmem* representative = evaluatedResults[activeRepresentative];
		activeEnd = max(activeEnd, end);
		representative->conflict = true;
		representative->isFiltered = false;
		representative->conflictSpan = static_cast<size_t>(activeEnd - activeStart);
		suppressed[i] = true;
		limitedWarning("overlapping sequence replacements",
			"overlapping sequence replacements on " + representative->getChrom() +
			" near one-based position " + to_string(activeStart + 1) +
			"; masking their combined reference span");
	}

	for (size_t i = 0; i < evaluatedResults.size(); ++i) {
		if (suppressed[i]) {
			if (evaluatedResults[i]->isINDEL()) ++Vstats->indelCNT;
			else if (evaluatedResults[i]->isSNP()) ++Vstats->snpCNT;
			// The representative owns the single conflict/uncertainty region.
			continue;
		}
		cF->ntVariant(evaluatedResults[i], Vstats, filterStats, opts);
	}
}

VCFmulti::~VCFmulti() {
	for (VCFmem* call : mems) delete call;
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
	(*curMems)[pos].push_back(vv);

}

VCFcollection::~VCFcollection() {
    for (auto it : VMems) {
		// delete stored VCFmem pointers inside each map
		map2vcfm* m = it.second;
		if (m) {
			for (auto &p : *m) {
				for (VCFmem* call : p.second) delete call;
			}
			delete m;
		}
	}
}

size_t VCFcollection::size() {
	size_t s(0); 
	if (VMems.size() == 0) { return 0; }
	for (const auto& contig : VMems) {
		for (const auto& position : *contig.second) s += position.second.size();
	} 
	return s; 
}


