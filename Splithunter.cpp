#include <getopt.h>
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegion.h"

using namespace SeqLib;
using namespace std;


static const char *USAGE_MESSAGE =
"Program: Splithunter\n"
"Contact: Haibao Tang\n"
"Usage: Splithunter bamfile [options]\n\n"
"Commands:\n"
"  --verbose,   -v        Set verbose output\n"
"  --reference, -r <file> Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to <htang@humanlongevity.com>\n\n";
static const char *DEBUG = "[ DEBUG ] ";

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

    RefGenome ref;
    ref.LoadIndex(opt::reference);

    const string tchr = "chr14";
    const int32_t tpos1 = 22386000;
    const int32_t tpos2 = 22477000;
    const string name = "TRA";
    ostringstream ss;
    ss << tchr << ":" << tpos1 << "-" << tpos2;
    const string tchrFull = ss.str();

    cerr << "[ BAM input ] " << opt::bam << endl;
    cerr << "[ Reference ] " << opt::reference << endl;
    cerr << "[    Target ] " << name << " (" << tchrFull << ")" << endl;

    // get sequence at given locus
    string seq = ref.QueryRegion(tchr, tpos1, tpos2);
    //cout << seq << endl;

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    UnalignedSequenceVector usv = {{name, seq}};
    bwa.ConstructIndex(usv);

    BamReader br;
    br.Open(opt::bam);
    GenomicRegion gr(tchrFull, br.Header());
    br.SetRegion(gr);

    BamRecord r;
    // Global settings
    const bool hardclip = false;
    const float secondary_cutoff = .9;
    const int secondary_cap = 0;
    const int PAD = 30;       // threshold for a significant match
    const int INDEL = 10000;  // threshold for left-right distance

    int totalSR = 0, validSR = 0;
    string leftPart, rightPart;
    unordered_map<string, GenomicRegionVector> cache;

    while (br.GetNextRecord(r)) {
        if (r.DuplicateFlag()) continue;

        totalSR++;
        string readName = r.Qname();
        int32_t readScore = r.NumAlignedBases();

        // Store the pairs in a cache
        if (r.PairedFlag() && (readScore >= r.Length() - PAD)) {
            auto got = cache.find(readName);
            if (got == cache.end()) {
                GenomicRegionVector grv;
                grv.push_back(r.AsGenomicRegion());
                cache.insert({ readName, grv });
            } else {
                (got->second).push_back(r.AsGenomicRegion());
            }
        }

        if (r.NumClip() < PAD || r.NumClip() > r.Length() - PAD) continue;

        BamRecordVector results;
        bwa.AlignSequence(r.Sequence(), readName,
                          results, hardclip, secondary_cutoff, secondary_cap);

        if (results.empty()) continue;

        BamRecord i = results[0];
        int32_t readLength = i.Length();
        if (i.NumClip() < PAD || i.NumClip() > readLength - PAD) continue;

        // Bipartite alignment: 0-index coordinate at this point
        int32_t queryStart = i.AlignmentPosition();
        int32_t queryEnd   = i.AlignmentEndPosition();
        if (queryStart > PAD) {
            leftPart = i.Sequence().substr(0, queryStart);
            rightPart = i.Sequence().substr(queryStart, readLength);
        } else if (queryEnd < readLength - PAD) {
            leftPart = i.Sequence().substr(0, queryEnd);
            rightPart = i.Sequence().substr(queryEnd, readLength);
        } else continue;

        BamRecordVector resultsL, resultsR;
        bwa.AlignSequence(leftPart, readName + "L",
                          resultsL, hardclip, secondary_cutoff, secondary_cap);

        bwa.AlignSequence(rightPart, readName + "R",
                          resultsR, hardclip, secondary_cutoff, secondary_cap);

        int32_t leftScore = 0, rightScore = 0;
        GenomicRegion leftAlign, rightAlign;
        BamRecord leftRec, rightRec;
        if (resultsL.empty() || resultsR.empty()) continue;

        leftRec = resultsL[0];
        rightRec = resultsR[0];

        // Condition 1: Each part is significant
        leftAlign = leftRec.AsGenomicRegion();
        leftRec.GetIntTag("AS", leftScore);
        rightAlign = rightRec.AsGenomicRegion();
        rightRec.GetIntTag("AS", rightScore);
        if (leftScore < PAD || rightScore < PAD) continue;

        // Condition 2: Total score
        int32_t totalScore = leftScore + rightScore;
        if (totalScore < readLength - PAD / 2) continue;

        // Condition 3: Distinct region
        int32_t dist = leftAlign.DistanceBetweenStarts(rightAlign);
        if (dist < INDEL) continue;

        // Verified alignment
        validSR++;
        if (opt::verbose) {
            cout << DEBUG << i;
            cout << DEBUG << "start = " << queryStart << " end = " << queryEnd
                          << " len = " << readLength << endl;
            cout << DEBUG << leftRec;
            cout << DEBUG << rightRec;
            cout << DEBUG << "SR Score    " << leftScore << " + "
                                            << rightScore << " = " << totalScore << endl;
            cout << DEBUG << "SR Distance " << leftAlign << " - "
                                            << rightAlign << " = " << dist << endl;
        }
    }

    int32_t totalSP = 0, validSP = 0;
    // Scan through the pairs
    for (auto i : cache) {
        // Condition 1: Properly paired
        if (i.second.size() != 2) continue;
        totalSP++;

        GenomicRegion leftAlign, rightAlign;
        leftAlign = i.second[0];
        rightAlign = i.second[1];
        int32_t dist = leftAlign.DistanceBetweenStarts(rightAlign);

        // Condition 2: Distinct region
        if (dist < INDEL) continue;
        validSP++;

        if (opt::verbose) {
            cout << DEBUG << i.first << endl;
            cout << DEBUG << "SP Distance " << leftAlign << " - "
                                            << rightAlign << " = " << dist << endl;
        }
    }

    cerr << "[  Total SR ] " << totalSR << endl;
    cerr << "[  Valid SR ] " << validSR << endl;
    cerr << "[  SR ratio ] " << validSR * 1e6 / totalSR << " ppm" << endl;
    cerr << "[  Total SP ] " << totalSP << endl;
    cerr << "[  Valid SP ] " << validSP << endl;
    cerr << "[  SP ratio ] " << validSP * 1e6 / totalSP << " ppm" << endl;
    cerr << endl;

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
    /*
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);
    cout << displayRuntime(start) << endl;
    */

    if (die || help || opt::reference.empty() || opt::bam.empty()) {
        cerr << "\n" << USAGE_MESSAGE;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }
}
