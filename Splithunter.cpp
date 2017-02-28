#include <getopt.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"

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
    static string bam;
    static string reference = "/mnt/ref/hg38.upper.fa";
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
    string seq = ref.QueryRegion("chr14", 22386000, 22388000);
    cout << seq << endl;

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    UnalignedSequenceVector usv = {{"chr_reg1", seq}};
    bwa.ConstructIndex(usv);

    BamReader br;
    br.Open(opt::bam);
    BamRecord r;
    bool hardclip = false;
    float secondary_cutoff = .9;
    int secondary_cap = 10;

    int counts = 0;
    while (br.GetNextRecord(r)) {
        BamRecordVector results;
        bwa.AlignSequence(r.Sequence(), r.Qname(),
                          results, hardclip, secondary_cutoff, secondary_cap);
        // print results to stdout
        for (auto& i : results)
        {
            cout << i;
            counts++;
        }
    }
    cout << "Number of alignments:" << counts << endl;

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
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'v': opt::verbose = true; break;
            case 'r': arg >> opt::reference; break;
            default: die = true;
        }
    }

    run();

    if (die || help || opt::reference.empty() || opt::bam.empty()) {
        cerr << "\n" << USAGE_MESSAGE;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }
}
