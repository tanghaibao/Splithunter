#include <getopt.h>
#include <json/json.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/GenomicRegion.h"
#include "bedFile/bedFile.h"

using namespace SeqLib;
using namespace std;


static const char *USAGE_MESSAGE =
"Program: BuildDB\n"
"Contact: Haibao Tang\n"
"Usage: BuildDB bedfile [options]\n\n"
"Commands:\n"
"  --verbose,   -v          Set verbose output\n"
"  --reference, -r <file>   Reference genome if using BWA-MEM realignment\n"
"  --outdir,    -s <string> Output directory\n"
"\nReport bugs to <htang@humanlongevity.com>\n\n";

namespace opt {
    static bool verbose = false;
    static string bed = "";
    static string outdir = "data";
    static string reference = "/mnt/ref/hg38.upper.fa";
}

static const char* shortopts = "hvb:o:r:";
static const struct option longopts[] = {
    { "help",       no_argument,       NULL, 'h' },
    { "verbose",    no_argument,       NULL, 'v' },
    { "bed",        required_argument, NULL, 'b' },
    { "outdir",     required_argument, NULL, 'o' },
    { "reference",  required_argument, NULL, 'r' },
    { NULL, 0, NULL, 0 }
};

// Where work is done
int run(BED& bedEntry, RefGenome& ref) {
    string tchr = bedEntry.chrom;
    int32_t tpos1 = bedEntry.start;
    int32_t tpos2 = bedEntry.end;
    string name = bedEntry.name;
    ostringstream ss;
    ss << tchr << ":" << tpos1 << "-" << tpos2;
    const string tchrFull = ss.str();

    cerr << "[    Target ] " << name << " (" << tchrFull << ")" << endl;

    // get sequence at given locus
    string seq = ref.QueryRegion(tchr, tpos1, tpos2);
    //cout << seq << endl;

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    UnalignedSequenceVector usv = {{name, seq}};
    bwa.ConstructIndex(usv);

    // Write to disk
    bwa.WriteIndex(opt::outdir + "/" + name);

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
        opt::bed = string(argv[1]);

    bool die = false;
    bool help = false;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'v': opt::verbose = true; break;
            case 'o': arg >> opt::outdir; break;
            case 'r': arg >> opt::reference; break;
            default: die = true;
        }
    }

    if (die || help || opt::reference.empty() || opt::bed.empty()) {
        cerr << "\n" << USAGE_MESSAGE;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }

    RefGenome ref;
    ref.LoadIndex(opt::reference);

    cerr << "[ BED input ] " << opt::bed << endl;
    cerr << "[ Reference ] " << opt::reference << endl;

    // Parse BEDFILE
    BED bedEntry;
    BedFile bed(opt::bed);
    bed.Open();
    while(bed.GetNextBed(bedEntry)) {
        run(bedEntry, ref);
    }

    return EXIT_SUCCESS;
}
