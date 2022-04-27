#include "fmindex.h"

using namespace std;

vector<pair<string, string>> getPairedReads(const string& base) {
    ifstream ifs(base + ".reads.fasta");
    if (!ifs)
        throw runtime_error("Problem reading: " + base + ".reads.fasta");

    string line;
    bool first = true;
    vector<pair<string, string>> pairedReads;

    while (getline(ifs, line)) {
        if (line[0] == '>')
            continue;
        if (first)
            pairedReads.emplace_back(line, "");
        else
            pairedReads.back().second = line;
        first = !first;
    }

    return pairedReads;
}

int main(int argc, char* argv[]) {
    string base = "testset/CP001363";
    FMIndex index = FMIndex(base, 32, true);
    string text = index.getText();
    auto reads = getPairedReads(base);
}
