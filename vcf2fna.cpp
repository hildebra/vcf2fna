#include "options.h"
#include "fasta.h"
#include "vcf.h"
#include "include/Benchmark.h"

void startMsg() {
	cout << "VCF2Fasta v"<< vcf2fnaVERSION <<" - A tool to convert VCF files to FASTA format\n";
}



int main(int argc, char* argv[])
{

if (argc < 2) { cerr << "Not enough arguments. Use \"vcf2fna -h\" for getting started.\n"; exit(3); }

ini_AA();

//clock_t tStart = clock();
Benchmark* bench = new Benchmark("Time vcf2fasta ");
bench->start();

//1 read_vcf_file in the command line arguments
options* opts = new options(argc, argv);


//2 read_vcf_file reference genome
refAssembly* refFA = new refAssembly(opts);
refFA->readDepth();
refFA->readGFF();
bench->now_total_time();


//3 go though the VCF file, reconstruct fasta sequences
VCFReader* vcf = new VCFReader(opts, refFA);
bench->now_total_time();

//writes requested contigs, genes (NT) and genes (AA)
refFA->writeOutputs();

//cleanup
delete vcf;
delete refFA;
delete opts;


bench->stop();
bench->printResults();
delete bench;

}