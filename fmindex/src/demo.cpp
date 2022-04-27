#include "bandmatrix.h"
#include "bidirectionalfmindex.h"
#include "searchscheme.h"
#include <chrono>

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

    string base = "../testset/CP001363";
    BiFMIndex bifmindex = BiFMIndex(base, 32, true);
    string text = bifmindex.getText();

    const auto reads = getPairedReads(base + "errors");

    length_t i = 0;

    auto start = chrono::high_resolution_clock::now();

    for (const auto& read : reads) {

        if (i % 16 == 0) {
            cout << "Progress: " << i << "/" << reads.size() * 2 << "\r";
            cout.flush();
        }
        bifmindex.naiveApproxMatch(read.first, 4);
        bifmindex.naiveApproxMatch(read.second, 4);

        i += 2;
    }

    /*   SearchScheme ss(bifmindex, "../search_schemes/pigeon/", 4);
      for (const auto& read : reads) {

          if (i % 16 == 0) {
              cout << "Progress: " << i << "/" << reads.size() * 2 << "\r";
              cout.flush();
          }
          ss.matchApprox(read.first);
          ss.matchApprox(read.second);

          i += 2;
      } */

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() * 2 << "/" << reads.size() * 2 << "\n";
    cout << "Total duration: " << fixed << elapsed.count() << "s\n";
}
