#include "fmindex.h"

#include "bandmatrix.h"
#include <fstream>

using namespace std;

// ============================================================================
// Debugging functions
// ============================================================================
ostream& operator<<(ostream& os, const Range& r) {
    os << "Range(" << r.begin << ", " << r.end << ")";
    return os;
}

ostream& operator<<(ostream& os, const TextOcc& t) {
    os << "TextOcc(" << t.range << ", " << t.distance << ")";
    return os;
}

ostream& operator<<(ostream& os, const FMPos& p) {
    os << "FMPos(" << p.range << ", " << p.depth << ")";
    return os;
}

ostream& operator<<(ostream& os, const FMOcc& o) {
    os << "FMOcc(" << o.pos << ", " << o.distance << ")";
    return os;
}
// ============================================================================
// I/O functions
// ============================================================================
bool readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        return false;

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)&buf[0], buf.size());

    return true;
}

bool readArray(const string& filename, vector<length_t>& array) {
    ifstream ifs(filename, ios::binary);
    if (!ifs)
        return false;

    ifs.seekg(0, ios::end);
    array.resize(ifs.tellg() / sizeof(length_t));
    ifs.seekg(0, ios::beg);
    ifs.read((char*)&array[0], array.size() * sizeof(length_t));

    return true;
}

void readSATextMode(const string& filename, vector<length_t>& sa,
                    size_t saSizeHint) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    sa.reserve(saSizeHint);
    length_t el;
    while (ifs >> el)
        sa.push_back(el);
}

void readSA(const string& filename, vector<length_t>& sa, size_t saSizeHint) {
    ifstream ifs(filename, ios::binary);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    size_t numElements = ifs.tellg() / sizeof(length_t);

    if (numElements == saSizeHint) { // file is likely binary
        sa.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, ios::beg);
        ifs.read((char*)sa.data(), sa.size() * sizeof(length_t));
    } else { // try to read SA in text mode
        readSATextMode(filename, sa, saSizeHint);
    }
}

// ============================================================================
// FM Index Construction: week 1
// ============================================================================

void FMIndex::createBWTFromSA(const vector<length_t>& sa) {
    // 5-10 lines of code
    bwt.resize(sa.size());
    for (size_t i=0; i<sa.size(); i++){
        bwt[i] = sa[i] > 0 ? text[sa[i]-1] : '$';
    }
}

void FMIndex::createCounts() {
    // 3 - 7 lines of code
    for (size_t i=0; i<sigma.size(); i++) {
      for (length_t j=0; j<textLength; j++) {
        if(sigma.c2i(text[j]) < (int)i) counts[i]++;
      }
    }
}

void printInfo(const string& info, bool verbose) {
    if (verbose) {
        cout << info << "...";
        cout.flush();
    }
}

void printDone(bool verbose) {
    if (verbose)
        cout << "done" << endl;
}

void FMIndex::read(const string& base, bool verbose) {
    // read the text
    printInfo("Reading: " + base + ".txt", verbose);

    if (!readText(base + ".txt", text))
        throw runtime_error("Problem reading: " + base + ".txt");

    textLength =
        (text[text.size() - 1] == '\n') ? text.size() - 1 : text.size();

    printDone(verbose);

    // read the suffix array
    printInfo("Reading: " + base + ".sa", verbose);

    vector<length_t> SA;
    readSA(base + ".sa", SA, textLength);
    printDone(verbose);

    printInfo("Creating sparse suffix array", verbose);
    sparseSA.createSparseSA(SA);
    printDone(verbose);
    // create BWT from SA
    printInfo("Creating BWT", verbose);
    createBWTFromSA(SA);
    printDone(verbose);

    dollarPos = distance(bwt.begin(), find(bwt.begin(), bwt.end(), '$'));

    // create ALPHABET from text
    sigma = Alphabet<ALPHABET>(text);

    // set counts to zero for each character
    for (length_t i = 0; i < ALPHABET; i++)
        counts[i] = 0;

    printInfo("Create Counts", verbose);
    createCounts();
    printDone(verbose);

    // create the occTable
    // Initialize the bitvectors with size BWT.size()
    printInfo("Create Occ table", verbose);
    occTable.resize(ALPHABET - 1);

    for (auto& bv : occTable)
        bv = Bitvec(textLength + 1);

    for (size_t i = 0; i < bwt.size(); i++) {
        char c = bwt[i];
        if (c != '$') {
            occTable[sigma.c2i(c) - 1][i] = true;
        }
    }

    // index the bitvectors
    for (auto& bv : occTable)
        bv.index();

    printDone(verbose);
    printInfo("FMIndex construction successful", verbose);
    if (verbose)
        cout << endl;
}

// ============================================================================
// FMIndex functionality:  week 1
// ============================================================================

length_t FMIndex::occ(const length_t& charIdx, const length_t& index) const {
    // 2 - 4 lines of code
    if(charIdx == 0) return index<=dollarPos ? 0 : 1;
    length_t occ = 0;
    occ = occTable[charIdx-1].rank(index);
    return occ;
}

length_t FMIndex::findLF(length_t k) const {
    // 1 - 2 lines of code
    return counts[sigma.c2i(bwt[k])] + occ(sigma.c2i(bwt[k]), k);
}

length_t FMIndex::findSA(length_t k) const {
    // 4 - 6 lines of code
    if(sparseSA.hasStored(k)) return sparseSA[k];
    else{
      int i = 0;
      while (true) {
          k = findLF(k);
          i++;
          if(sparseSA.hasStored(k)){
              return ((sparseSA[k]+i) < 0) || ((sparseSA[k]+i) >= text.size()) ? 0 : sparseSA[k]+i;
          }
      }
    }
}

bool FMIndex::addCharLeft(length_t charIdx, const Range& originalRange,
                          Range& newRange) const {
    // 2 - 4 lines of code
    newRange = Range(counts[charIdx] + occ(charIdx, originalRange.getBegin()), counts[charIdx] + occ(charIdx, originalRange.getEnd()));
    return !newRange.empty();
}

// ============================================================================
// FMIndex Integration: week 1
// ============================================================================

#include <iostream>
using namespace std;
#include <string>

vector<length_t> FMIndex::matchExact(const string& str) const {
    // 8 - 12 lines of code
    vector<length_t> result;
    bool matchLeft = false;
    Range range = Range(0, text.size());
    for (int i = (str.size()-1); i >= 0; i--) {
      matchLeft = addCharLeft(sigma.c2i(str[i]), range, range);
      if(!matchLeft) return result;
    }
    for (length_t i = range.getBegin(); i < range.getEnd(); i++) result.push_back(findSA(i));
    return result;
}

#include <cmath>

tuple<length_t, length_t, bool>
FMIndex::bestPairedMatch(const pair<string, string>& reads,
                         const length_t& meanInsSize) const {
    // 15 - 25 lines of code
    vector<length_t> matchesRead1 = matchExact(reads.first), matchesRead2 = matchExact(reads.second), matchesRead1Rev = matchExact(revCompl(reads.first)), matchesRead2Rev = matchExact(revCompl(reads.second));
    int currentBestInsert = text.size();
    length_t currentBestRead1, currentBestRead2;
    bool is2Rev = false;

    for (size_t i = 0; i < matchesRead1.size(); i++) {
      for (size_t j = 0; j < matchesRead2Rev.size(); j++) {
        if(currentBestInsert >= abs(abs((int)(matchesRead2Rev[j] + reads.second.size() - matchesRead1[i]))-(int)meanInsSize)){
          currentBestInsert = abs(abs((int)(matchesRead2Rev[j] + reads.second.size() - matchesRead1[i]))-(int)meanInsSize);
          currentBestRead1 = matchesRead1[i];
          currentBestRead2 = matchesRead2Rev[j];
          is2Rev = 1;
        }
      }
    }
    for (size_t i = 0; i < matchesRead2.size(); i++) {
      for (size_t j = 0; j < matchesRead1Rev.size(); j++) {
        if(currentBestInsert >= abs(abs((int)(matchesRead1Rev[j] + reads.first.size() - matchesRead2[i]))-(int)meanInsSize)){
          currentBestInsert = abs(abs((int)(matchesRead1Rev[j] + reads.first.size() - matchesRead2[i]))-(int)meanInsSize);
          currentBestRead2 = matchesRead2[i];
          currentBestRead1 = matchesRead1Rev[j];
          is2Rev = 0;
        }
      }
    }

    return make_tuple(currentBestRead1,currentBestRead2,is2Rev);
}

// ============================================================================
// FMIndex functionality:  week 2
// ============================================================================
void FMIndex::extendFMPos(const Range& range, const length_t& depth,
                          std::vector<FMPosExt>& stack) const {
    // 4 lines of code
    Range r=range;
    for (length_t i=1; i<sigma.size(); i++){
        if(addCharLeft(i, range, r)) stack.push_back(FMPosExt(sigma.i2c(i), r, depth+1));
    }
}

void FMIndex::convertFMOccToTextOcc(const FMOcc& fmocc,
                                    std::vector<TextOcc>& textOcc) const {
    // 3 - 4 lines of code
    for(length_t i=fmocc.getRange().getBegin(); i<fmocc.getRange().getEnd(); i++ ){
        textOcc.push_back(TextOcc(Range(findSA(i), findSA(i)+fmocc.getDepth()), fmocc.getDistance()));
    }
}

// ============================================================================
// FMIndex Integration:  week 2
// ============================================================================

vector<TextOcc> FMIndex::naiveApproxMatch(const string& pattern,
                                          length_t k) const {

    // create the occurrences vector (to which an FMOcc is pushed)
    vector<FMOcc> occ;

    // create the stack and reserve space
    vector<FMPosExt> stack;
    stack.reserve((pattern.size() + k + 1) * (sigma.size() - 1));

    // Create the matrix for this pattern and edit distance value
    BandedMatrix matrix(pattern.size(), k, 0);

    // Create the first 4 entries in the stack, corresponding to the "A",
    // "C", "G" and "T" strings with depth 1
    extendFMPos(Range(0, text.size()), 0, stack);

    // Create a substring from pattern with backward direction (1 line)
    Substring p(pattern, BACKWARD);

    while (!stack.empty()) {

        // Get the final element from the stack and pop it back (= remove from
        // the stack)
        // Uncomment these lines
        /*  FMPosExt currentPos = stack.back();
         stack.pop_back(); */

        // 7 - 15 lines of code
        throw std::runtime_error(
            "naiveApproxMatch has not been implemented yet");
    }

    return filterRedundantMatches(occ, k);
}

std::vector<TextOcc> FMIndex::filterRedundantMatches(std::vector<FMOcc>& fmocc,
                                                     const length_t& k) const {

    // A) remove doubles in the fmoccurrences
    std::sort(fmocc.begin(), fmocc.end());
    fmocc.erase(std::unique(fmocc.begin(), fmocc.end()), fmocc.end());

    // B) convert fmoccurrences to occurrences in text
    std::vector<TextOcc> textocc;
    for (const auto& f : fmocc) {
        convertFMOccToTextOcc(f, textocc);
    }

    // C) erase doubles from textocc
    std::sort(textocc.begin(), textocc.end());
    textocc.erase(std::unique(textocc.begin(), textocc.end()), textocc.end());

    // D) remove occurrences that just have extra deletions/insertions
    if (textocc.empty()) {
        return textocc;
    }
    std::vector<TextOcc> r;
    r.emplace_back(textocc.front());

    length_t maxDiff = 2 * k;

    for (const auto& o : textocc) {
        length_t prevBegin = r.back().begin();
        length_t prevLength = r.back().getRange().width();
        length_t prevED = r.back().getDistance();

        auto diff = o.begin() - prevBegin;
        if (diff == 0) {
            continue;
        }
        if (diff <= maxDiff) {
            // check if this later occurrence is better than the
            // previous one
            if (o.getDistance() > prevED) {
                continue;
            }
            if (o.getDistance() == prevED &&
                o.getRange().width() >= prevLength) {
                continue;
            }

            // prev was worse so pop_back
            r.pop_back();
        }

        r.emplace_back(o);
    }

    // D) return the occurrences
    return r;
}
