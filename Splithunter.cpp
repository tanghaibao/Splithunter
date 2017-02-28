#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"

using namespace SeqLib;
using namespace std;


int main(int argc, char* argv[]) {
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
