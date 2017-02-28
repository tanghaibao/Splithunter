#include <getopt.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"

using namespace SeqLib;
using namespace std;


static const char *USAGE_MESSAGE =
"Program: Splithunter \n"
"Contact: Haibao Tang \n"
"Usage: Splithunter bamfile [options]\n\n"
"Commands:\n"
"  --verbose,   -v        Set verbose output\n"
"  --reference, -r <file> Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to <htang@humanlongevity.com>\n\n";

namespace opt {
  static bool verbose = false;
  static std::string bam;
  static std::string reference = "/mnt/ref/hg38.upper.fa";
}

static const char* shortopts = "hvb:r:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "reference",               required_argument, NULL, 'r' },
  { NULL, 0, NULL, 0 }
};


// Where work is done
int run() {
    cout << opt::bam << endl;
    cout << opt::reference << endl;

    RefGenome ref;
    ref.LoadIndex(opt::reference);

    // get sequence at given locus
    string seq = ref.QueryRegion("chr1", 1000000,1001000);
    cout << seq << endl;

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    UnalignedSequenceVector usv = {{"chr_reg1", seq}};
    bwa.ConstructIndex(usv);

    // align an example string with BWA-MEM
    std::string querySeq = "CGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGC";
    BamRecordVector results;
    // hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
    bwa.AlignSequence(querySeq, "my_seq", results, false, 0.9, 10);

    // print results to stdout
    for (auto& i : results)
        std::cout << i << std::endl;

    return 0;
}


// Parse the command line options
int main(int argc, char** argv) {
    if (argc <= 1) {
        cerr << USAGE_MESSAGE;
        return 0;
    }

    // Get the first argument as input
    if (argc > 1)
        opt::bam = string(argv[1]);

    bool die = false;
    bool help = false;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'v': opt::verbose = true; break;
            case 'r': arg >> opt::reference; break;
            default: die = true;
        }
    }

    run();

    if (die || help || opt::reference.empty() || opt::bam.empty()) {
        std::cerr << "\n" << USAGE_MESSAGE;
        if (die)
	        exit(EXIT_FAILURE);
        else
	        exit(EXIT_SUCCESS);
    }
}
