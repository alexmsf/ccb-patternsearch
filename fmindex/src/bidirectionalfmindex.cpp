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
    originalOccTable = CumulativeBitvectors<ALPHABET>(sigma, bwt);
    reverseOccTable = CumulativeBitvectors<ALPHABET>(sigma, revBWT);
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

    // assert that the direction of str and the BiFMIndex match, leave this line
    // in! If your program blocks on this line it means that somewhere in your
    // code you forgot to correctly set the direction of the index
    assert(dir == str.getDirection());

    // 5 - 10 lines of code
    throw runtime_error(
        "matchExactBidirectionally has not been implemented yet");
}

void BiFMIndex::recApproxMatch(const Search& s, const BiFMOcc& startOcc,
                               vector<FMOcc>& occ,
                               const vector<Substring>& parts, const int& idx) {

    // Fill in the TODO's at lines 70 and 77. You can use the provided code
    // in function FMIndex::naiveApproxMatch in file fmindex.cpp as inspiration
    // Afterwards fill in the while loop at line 88

    // TODO create the matrix for the current part, with the correct width and
    // intialization value (4 - 5 lines)

    // Create the stack and reserve space
    vector<BiFMPosExt> stack; // stack with positions to visit
    stack.reserve((parts.back().end()) * ALPHABET);

    // TODO set the direction (1 - 3 lines)
    // use function setDirection()

    // Add the children of the start occurrence to the stack, make sure
    // they have depth/row 1
    extendFMPos(startOcc.getRanges(), 0, stack);

    // Branch and bound algorithm
    while (!stack.empty()) {

        // Get the final element from the stack and pop it back (= remove from
        // the stack)
        // Uncomment these lines
        /* BiFMPosExt currentPos = stack.back();
        stack.pop_back(); */
        // 10 - 25 lines of code
        throw runtime_error("recApproxMatch has not been implemented yet");
    }
}
