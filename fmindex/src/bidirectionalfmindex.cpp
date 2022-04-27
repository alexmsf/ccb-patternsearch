#include "bidirectionalfmindex.h"

using namespace std;
#include "bandmatrix.h"

ostream& operator<<(ostream& os, const RangePair& r) {
    os << "RangePair(" << r.getBackwardRange() << ", " << r.getForwardRange()
       << ")";
    return os;
}

void BiFMIndex::read(const string& base, bool verbose) {
    // step 1 create the rev bwt
    vector<length_t> revSA;
    readSA(base + ".rev.sa", revSA, textLength);
    string revBWT;
    createRevBWTFromRevSA(revSA, revBWT);

    // step 2 create the cumulative bit vectors
    backwardOccTable = CumulativeBitvectors<ALPHABET>(sigma, bwt);
    forwardOccTable = CumulativeBitvectors<ALPHABET>(sigma, revBWT);
}

void BiFMIndex::createRevBWTFromRevSA(const vector<length_t>& revSA,
                                      string& revBWT) {
    // 5-10 lines of code
    std::cout << "Warning createRevBWTFromRevSA has not been implmented yet!"
              << std::endl;
}

bool BiFMIndex::addCharRight(length_t charIdx, const RangePair& originalRanges,
                             RangePair& newRanges) const {
    // 8 - 12 lines of code
    throw std::runtime_error("addCharRight has not been implemented yet!");
}

bool BiFMIndex::addCharLeft(length_t charIdx, const RangePair& originalRanges,
                            RangePair& newRanges) const {

    // 8 - 12 lines of code
    throw std::runtime_error("addCharLeft has not been implemented yet!");
}

void BiFMIndex::extendFMPos(const RangePair& ranges, const length_t& depth,
                            vector<BiFMPosExt>& stack) {
    // 4 lines of code
    throw std::runtime_error("extendFMPos has not been implemented yet!");
}

RangePair BiFMIndex::matchExactBidirectionally(const Substring& str,
                                               RangePair ranges) const {

    // Match the string in one direction, but keep track of the ranges in both
    // directions. Return the ranges corresponding to the match. If no match can
    // be found the returned ranges should be empty

    assert(dir == str.getDirection());

    // 5 - 10 lines of code
    throw runtime_error(
        "matchExactBidirectionally has not been implemented yet");
}

void BiFMIndex::recApproxMatch(const Search& s, const BiFMOcc& startOcc,
                               vector<FMOcc>& occ,
                               const vector<Substring>& parts, const int& idx) {

    // TODO create the matrix for the current part, with the correct width and
    // intialization value (4 - 5 lines)

    // Create the stack and reserve space
    vector<BiFMPosExt> stack; // stack with positions to visit
    stack.reserve((parts.back().end()) * ALPHABET);

    // TODO set the direction (1 -2 lines)

    // TODO add the children of the start occurrence to the stack, make sure
    // they have depth 1 (row in this matrix)

    // Branch and bound algorithm
    while (!stack.empty()) {
        // 15 - 25 lines of code
        throw runtime_error("recApproxMatch has not been implemented yet");
    }
}
