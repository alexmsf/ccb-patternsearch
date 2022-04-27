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
    // SLIDE 32 FMINDEX
    bwt.resize(sa.size());
    for (size_t i=0; i<sa.size(); i++){
        bwt[i] = sa[i] > 0 ? text[sa[i]-1] : '$';
    }
}

void FMIndex::createCounts() {
    // 3 - 7 lines of code
    int index = 0;
    for (size_t i=0; i<textLength; i++) {
        index = sigma.c2i(text[i]);
        for (int j=0; j<index; j++) counts[i]++;
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
    throw runtime_error("occ is not implemented yet!");
}

length_t FMIndex::findLF(length_t k) const {
    // 1 - 2 lines of code
    throw runtime_error("findLF is not implemented yet!");
}

length_t FMIndex::findSA(length_t k) const {
    // 4 - 6 lines of code
    throw runtime_error("findSA is not implemented yet!");
}

bool FMIndex::addCharLeft(length_t charIdx, const Range& originalRange,
                          Range& newRange) const {
    // 2 - 4 lines of code
    throw std::runtime_error("addCharLeft is not implemented yet!");
}

// ============================================================================
// FMIndex Integration: week 1
// ============================================================================

vector<length_t> FMIndex::matchExact(const string& str) const {
    // 8 - 12 lines of code
    throw runtime_error("matchExact is not implemented yet!");
}

tuple<length_t, length_t, bool>
FMIndex::bestPairedMatch(const pair<string, string>& reads,
                         const length_t& insSize) const {
    // 15 - 25 lines of code
    throw runtime_error("bestPairedMatch is not implemented yet!");
}

// ============================================================================
// FMIndex functionality:  week 2
// ============================================================================
void FMIndex::extendFMPos(const Range& range, const length_t& depth,
                          std::vector<FMPosExt>& stack) const {
    // 4 lines of code
    throw std::runtime_error("extendFMPos has not been implemented yet!");
}

void FMIndex::convertFMOccToTextOcc(const FMOcc& fmocc,
                                    std::vector<TextOcc>& textOcc) const {
    // 3 - 4 lines of code
    throw std::runtime_error(
        "ConvertFMOccToTextOcc has not been implemented yet!");
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

    // TODO create the matrix for this pattern and edit distance value (1 line)

    // TODO Create the first 4 entries in the stack, corresponding to the "A",
    // "C", "G" and "T" strings with depth 1

    // create a substring from pattern with backward direction (1 line)
    Substring p(pattern, BACKWARD);

    while (!stack.empty()) {
        // 10 - 15 lines of code
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
