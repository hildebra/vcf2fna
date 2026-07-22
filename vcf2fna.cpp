#include "options.h"
#include "fasta.h"
#include "vcf.h"
#include "include/Benchmark.h"

void startMsg() {
	cout << "VCF2Fasta v"<< vcf2fnaVERSION <<" - A tool to convert VCF files to FASTA format\n";
}



int main(int argc, char* argv[])
{
	try {
		if (argc < 2) {
			cerr << "Not enough arguments. Use \"vcf2fna -h\" for getting started.\n";
			return 3;
		}

		Benchmark bench("Time vcf2fasta ");
		bench.start();
		options opts(argc, argv);
		refAssembly refFA(&opts);
		refFA.readGFF();
		refFA.readDepth();
		bench.now_total_time();

		VCFReader vcf(&opts, &refFA);
		bench.now_total_time();
		refFA.writeOutputs();

		bench.stop();
		bench.printResults();
		return 0;
	} catch (const std::exception& error) {
		cerr << "vcf2fna: " << error.what() << '\n';
		return 1;
	}
}
