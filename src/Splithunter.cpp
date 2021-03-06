#include <getopt.h>
#include <json/json.h>
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegion.h"
#include "bedFile/bedFile.h"

using namespace SeqLib;
using namespace std;


static const char *USAGE_MESSAGE =
"Program: Splithunter\n"
"Contact: Haibao Tang\n"
"Usage: Splithunter bamfile [options]\n\n"
"Commands:\n"
"  --help,      -h          Display help\n"
"  --verbose,   -v          Set verbose output\n"
"  --samplekey, -s <string> SampleKey, output will be written `samplekey.json`\n"
"  --bedfile,   -b <string> BedFile that contains the loci to go through\n"
"  --locus,     -l <string> Specific locus to compute, one of TRA/TRB/TRG/IGH/IGK/IGL\n"
"\nReport bugs to <htang@humanlongevity.com>\n\n";
static const char *DEBUG = "[ DEBUG ] ";

// Global settings
static const bool hardclip = false;
static const float secondary_cutoff = .9;
static const int secondary_cap = 0;
static const int PAD = 30;       // threshold for a significant match
static const int INDEL = 10000;  // threshold for left-right distance
static const int MINENT = 50;    // threshold for sequence complexity

namespace opt {
    static bool verbose = false;
    static string bam;
    static string bed = "";
    static string locus = "";
    static string path = "";
    static string samplekey = "";
}

static const char* shortopts = "hvb:l:s:";
static const struct option longopts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "verbose",     no_argument,       NULL, 'v' },
    { "bedfile",     required_argument, NULL, 'b' },
    { "locus",       required_argument, NULL, 'l' },
    { "samplekey",   required_argument, NULL, 's' },
    { NULL, 0, NULL, 0 }
};


// https://academic.oup.com/bioinformatics/article/27/6/863/236283/Quality-control-and-preprocessing-of-metagenomic
// Schmieder and Edwards. Quality control and preprocessing of metagenomic datasets. (2011) Bioinformatics
static double entropy(const string& seq) {
    unordered_map<string, int> counts;
    int l = seq.length() - 2;
    if (l <= 0) return 0;

    int k = (l < 64) ? l: 64;
    for (int i = 0; i < l; i++) {
        string trinuc = seq.substr(i, 3);
        auto j = counts.find(trinuc);
        if (j == counts.end()) {
            counts.insert({ trinuc, 1 });
        } else {
            j->second++;
        }
    }

    double res = 0.0;
    for (auto i: counts) {
        double f = (double) i.second / l;
        res += f * log(f) / log(k);
    }
    return res * -100;
}

// Where work is done
int run(BED& bedEntry, Json::Value& root) {
    string tchr = bedEntry.chrom;
    int32_t tpos1 = bedEntry.start;
    int32_t tpos2 = bedEntry.end;
    string name = bedEntry.name;
    ostringstream ss;
    ss << tchr << ":" << tpos1 << "-" << tpos2;
    const string tchrFull = ss.str();

    cerr << "[    Target ] " << name << " (" << tchrFull << ")" << endl;

    // Make an in-memory BWA-MEM index of region
    BWAWrapper bwa;
    bwa.LoadIndex(opt::path + "/data/" + name);

    BamReader br;
    br.Open(opt::bam);
    GenomicRegion gr(tchrFull, br.Header());
    br.SetRegion(gr);

    BamRecord r;
    int totalSR = 0, validSR = 0;
    string leftPart, rightPart;
    // Store the pairs in a cache
    unordered_map<string, BamRecordVector> cache;

    // Details of the alignments
    stringstream SR_details, SP_details;

    while (br.GetNextRecord(r)) {
        if (r.DuplicateFlag()) continue;

        totalSR++;
        string readName = r.Qname();
        int32_t readScore = r.NumAlignedBases();

        // SP Condition 1: Each part is significant
        if (r.PairedFlag() && (readScore >= r.Length() - PAD)) {
            auto got = cache.find(readName);
            if (got == cache.end()) {
                BamRecordVector brv;
                brv.push_back(r);
                cache.insert({ readName, brv });
            } else {
                (got->second).push_back(r);
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

        // SR Condition 1: Each part is significant
        leftAlign = leftRec.AsGenomicRegion();
        leftRec.GetIntTag("AS", leftScore);
        rightAlign = rightRec.AsGenomicRegion();
        rightRec.GetIntTag("AS", rightScore);
        if (leftScore < PAD || rightScore < PAD) continue;

        // SR Condition 2: Total score
        int32_t totalScore = leftScore + rightScore;
        if (totalScore < readLength - PAD / 2) continue;

        // SR Condition 3: Distinct region
        int32_t dist = leftAlign.DistanceBetweenStarts(rightAlign);
        if (dist < INDEL) continue;

        // SR Condition 4: Complexity filter
        double leftEnt = entropy(leftPart);
        double rightEnt = entropy(rightPart);
        if ((leftEnt < MINENT) || (rightEnt < MINENT)) continue;

        // Verified alignment
        validSR++;
        SR_details << leftAlign << "|" << rightAlign << ";";
        if (opt::verbose) {
            cerr << DEBUG << i;
            cerr << DEBUG << "start = " << queryStart << " end = " << queryEnd
                          << " len = " << readLength << endl;
            cerr << DEBUG << leftRec;
            cerr << DEBUG << rightRec;
            cerr << DEBUG << "SR Score    " << leftScore << " + "
                                            << rightScore << " = " << totalScore << endl;
            cerr << DEBUG << "SR Distance " << leftAlign << " - "
                                            << rightAlign << " = " << dist << endl;
            cerr << DEBUG << "Entropy     " << leftEnt << " " << rightEnt << endl;
        }
    }

    int32_t totalSP = 0, validSP = 0;
    // Scan through the pairs
    for (auto i : cache) {
        if (i.second.size() != 2) continue;
        totalSP++;

        // SP Condition 2: Distinct region
        GenomicRegion leftAlign, rightAlign;
        leftAlign = i.second[0].AsGenomicRegion();
        rightAlign = i.second[1].AsGenomicRegion();
        int32_t dist = leftAlign.DistanceBetweenStarts(rightAlign);
        if (dist < INDEL) continue;

        // SP Condition 3: Sequence complexity
        double leftEnt = entropy(i.second[0].Sequence());
        double rightEnt = entropy(i.second[1].Sequence());
        if ((leftEnt < MINENT) || (rightEnt < MINENT)) continue;

        validSP++;
        SP_details << leftAlign << "|" << rightAlign << ";";
        if (opt::verbose) {
            cerr << DEBUG << i.first << endl;
            cerr << DEBUG << "SP Distance " << leftAlign << " - "
                                            << rightAlign << " = " << dist << endl;
        }
    }

    double SR_PPM = totalSR ? validSR * 1e6 / totalSR : 0;
    double SP_PPM = totalSP ? validSP * 1e6 / totalSP : 0;
    cerr << "[  SR ratio ] " << validSR << " / " << totalSR << " = "
                             << SR_PPM << " ppm" << endl;
    cerr << "[  SP ratio ] " << validSP << " / " << totalSP << " = "
                             << SP_PPM << " ppm" << endl;

    root[name + ".SR-SIGNAL"] = validSR;
    root[name + ".SR-TOTAL"] = totalSR;
    root[name + ".SR-PPM"] = SR_PPM;
    root[name + ".SR-DETAILS"] = SR_details.str();
    root[name + ".SP-SIGNAL"] = validSP;
    root[name + ".SP-TOTAL"] = totalSP;
    root[name + ".SP-PPM"] = SP_PPM;
    root[name + ".SP-DETAILS"] = SP_details.str();

    if (opt::samplekey != "") {
        ofstream ofs(opt::samplekey + ".json", ofstream::out);
        ofs << root << endl;
        ofs.close();
    }

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

    string executable(argv[0]);
    opt::path = executable.substr(0, executable.find_last_of("\\/"));

    bool die = false;
    bool help = false;

    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'h': help = true; break;
            case 'v': opt::verbose = true; break;
            case 'b': arg >> opt::bed; break;
            case 'l': arg >> opt::locus; break;
            case 's': arg >> opt::samplekey; break;
            default: die = true;
        }
    }

    if (die || help || opt::bam.empty()) {
        cerr << "\n" << USAGE_MESSAGE;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }

    if (opt::bed.empty()) {
        opt::bed = opt::path + "/data/TR_IG.bed";
    }

    // JSON result object
    Json::Value root;
    root["bam"] = opt::bam;
    root["SampleKey"] = opt::samplekey;

    cerr << "[ BED input ] " << opt::bed << endl;
    cerr << "[ BAM input ] " << opt::bam << endl;

    // Parse BEDFILE
    BED bedEntry;
    BedFile bed(opt::bed);
    bed.Open();
    while(bed.GetNextBed(bedEntry)) {
        if ((opt::locus != "") && (bedEntry.name != opt::locus))
            continue;
        run(bedEntry, root);
    }

    cerr << root << endl;

    return EXIT_SUCCESS;
}
