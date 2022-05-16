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
    revBWT.resize(revSA.size());
    revBWT[0] = text[0];
    for (size_t i=1; i<revSA.size(); i++){
        revBWT[i] = revSA[i] > 0 ? text[text.size()-revSA[i]] : '$';
    }
}

bool BiFMIndex::addCharRight(length_t charIdx, const RangePair& originalRanges,
                             RangePair& newRanges) const {
    // 8 - 12 lines of code
    length_t i_ = originalRanges.getForwardRange().getBegin();
    length_t j_ = originalRanges.getForwardRange().getEnd();
    size_t k_ = counts[charIdx] + reverseOccTable.occ(charIdx, i_);
    size_t l_ = counts[charIdx] + reverseOccTable.occ(charIdx, j_);

    length_t i = originalRanges.getBackwardRange().getBegin();
    //length_t j = originalRanges.getBackwardRange().getEnd();
    size_t x = reverseOccTable.cumulocc(charIdx, j_) - reverseOccTable.cumulocc(charIdx, i_);
    size_t y = l_ - k_, k = i + x, l = i + x + y;
    Range newRangeForward = Range(k_, l_);
    Range newRangeBackward = Range(k, l);

    newRanges = RangePair(newRangeBackward, newRangeForward);
    return !newRanges.empty();
}

bool BiFMIndex::addCharLeft(length_t charIdx, const RangePair& originalRanges,
                            RangePair& newRanges) const {

    // 8 - 12 lines of code
    length_t i_ = originalRanges.getBackwardRange().getBegin();
    length_t j_ = originalRanges.getBackwardRange().getEnd();
    size_t k_ = counts[charIdx] + originalOccTable.occ(charIdx, i_);
    size_t l_ = counts[charIdx] + originalOccTable.occ(charIdx, j_);

    length_t i = originalRanges.getForwardRange().getBegin();
    //length_t j = originalRanges.getForwardRange().getEnd();
    size_t x = originalOccTable.cumulocc(charIdx, j_) - originalOccTable.cumulocc(charIdx, i_);
    size_t y = l_ - k_, k = i + x, l = i + x + y;
    Range newRangeBackward = Range(k_, l_);
    Range newRangeForward = Range(k, l);

    newRanges = RangePair(newRangeBackward, newRangeForward);
    return !newRanges.empty();
}

void BiFMIndex::extendFMPos(const RangePair& ranges, const length_t& depth,
                            vector<BiFMPosExt>& stack) {
    // 4 lines of code
    RangePair r=ranges;
    for (length_t i=1; i<sigma.size(); i++){
        if(dir == FORWARD) {
            if(addCharRight(i, ranges, r)) stack.push_back(BiFMPosExt(sigma.i2c(i), r, depth+1));
        }
        else {
            if(addCharLeft(i, ranges, r)) stack.push_back(BiFMPosExt(sigma.i2c(i), r, depth+1));   
        }
    }
}

RangePair BiFMIndex::matchExactBidirectionally(const Substring& str,
                                               RangePair ranges) const {

    // assert that the direction of str and the BiFMIndex match, leave this line
    // in! If your program blocks on this line it means that somewhere in your
    // code you forgot to correctly set the direction of the index
    assert(dir == str.getDirection());

    // 5 - 10 lines of code
    /*for (length_t i = 0; i < str.size(); i++) {
        bool add = dir == FORWARD ? addCharRight(sigma.c2i(str[i]), ranges, ranges) : addCharLeft(sigma.c2i(str[i]), ranges, ranges); 
        if(!add) return RangePair();
    }
    //for (length_t i = range.getBegin(); i < range.getEnd(); i++) result.push_back(findSA(i));
    return ranges;*/
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
    // (length_t patternsize, int W, int startValue)
    
    //BandedMatrix matrix(parts.size(), k, 0);


    // Create the stack and reserve space
    vector<BiFMPosExt> stack; // stack with positions to visit
    stack.reserve((parts.back().end()) * ALPHABET);

    // TODO set the direction (1 - 3 lines)
    // use function setDirection()
    for(size_t i=0; i<=parts.size(); i++) {
        //parts[i].setDirection(s.getDirection(i));
    }
    

    // Add the children of the start occurrence to the stack, make sure
    // they have depth/row 1
    extendFMPos(startOcc.getRanges(), 0, stack);

    // Branch and bound algorithm
    while (!stack.empty()) {

        // Get the final element from the stack and pop it back (= remove from
        // the stack)
        // Uncomment these lines
        //BiFMPosExt currentPos = stack.back();
        //stack.pop_back();
        // 10 - 25 lines of code
        //if(currentPos.)
        //length_t minimalEditDist = matrix.updateMatrixRow(p, currentPos.getRow(), currentPos.getCharacter());
        //if(minimalEditDist<k) extendFMPos(currentPos.getRange(), currentPos.getDepth(), stack);
        //if(matrix.inFinalColumn(currentPos.getRow())) {
        //    length_t matrixValueFinalCol = matrix.getValueInFinalColumn(currentPos.getRow());
        //    if(matrixValueFinalCol<=k) occ.push_back(FMOcc(currentPos, matrixValueFinalCol));
        //}  
        throw runtime_error("recApproxMatch has not been implemented yet");   
    }
}
