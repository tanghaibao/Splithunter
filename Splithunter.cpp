#include <getopt.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"

using namespace SeqLib;
using namespace std;


static const char *USAGE_MESSAGE =
"Program: Splithunter \n"
"Contact: Haibao Tang \n"
"Usage: Splithunter [options]\n\n"
"Commands:\n"
"  --verbose,   -v        Set verbose output\n"
"  --bam,       -b <file> Input a BAM\n"
"  --reference, -r <file> Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to <htang@humanlongevity.com>\n\n";

namespace opt {

  static bool verbose = false;
  static char mode = 's';
  static std::string input;
  static std::string reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static std::string fasta; // input is a fasta
  static std::string target; // input target sequence
}

static const char* shortopts = "hvb:r:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "bam",                     no_argument, NULL, 'b' },
  { "reference",               required_argument, NULL, 'r' },
  { NULL, 0, NULL, 0 }
};


int run(int argc, char *argv[]) {
    RefGenome ref;
    cout << argv[1] << endl;
    ref.LoadIndex(argv[1]);

    // get sequence at given locus
    string seq = ref.QueryRegion("chr1", 1000000,1001000);
    cout << seq.substr(0, 10) << endl;

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    UnalignedSequenceVector usv = {{"chr_reg1", seq}};
    bwa.ConstructIndex(usv);

    // align an example string with BWA-MEM
    std::string querySeq = "CAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAG";
    BamRecordVector results;
    // hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
    bwa.AlignSequence(querySeq, "my_seq", results, false, 0.9, 10);

    // print results to stdout
    for (auto& i : results)
        std::cout << i << std::endl;

    return 0;
}

// parse the command line options
void parseOptions(int argc, char** argv, const char* msg) {

  bool die = false;
  bool help = false;

  // get the first argument as input
  if (argc > 1)
    opt::input = std::string(argv[1]);

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'f': opt::mode = 'f'; break;
    case 'F': arg >> opt::fasta; break;
    case 'b': opt::mode = 'b'; break;
    case 'C': opt::mode = 'C'; break;
    case 'T': arg >> opt::target; break;
    case 'G': arg >> opt::reference; break;
    default: die= true;
    }
  }

  if (die || help || (opt::input.empty() && opt::fasta.empty())) {
      std::cerr << "\n" << msg;
      if (die)
	exit(EXIT_FAILURE);
      else
	exit(EXIT_SUCCESS);
    }
}


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        cerr << USAGE_MESSAGE;
        return 0;
    } else {
        run(argc, argv);
    }

    return EXIT_SUCCESS;
}



