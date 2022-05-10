#ifndef BIDIRECTIONALFMINDEX_H
#define BIDIRECTIONALFMINDEX_H

#include "cumulativebitvec.h"
#include "fmindex.h"

// ============================================================================
// CLASS RANGEPAIR: PROVIDED STEP 3
// ============================================================================
class RangePair {
  private:
    Range backwardRange;
    Range forwardRange;

  public:
    /**
     * Default constructor, creates two empty ranges
     */
    RangePair() : backwardRange(Range()), forwardRange(Range()) {
    }
    RangePair(Range backwardRange, Range forwardRange)
        : backwardRange(backwardRange), forwardRange(forwardRange) {
    }
    RangePair(const length_t backwardBegin, const length_t backwardEnd,
              const length_t forwardBegin, const length_t forwardEnd)
        : backwardRange(backwardBegin, backwardEnd),
          forwardRange(forwardBegin, forwardEnd) {
    }

    const Range& getForwardRange() const {
        return forwardRange;
    }
    const Range& getBackwardRange() const {
        return backwardRange;
    }

    /**
     * @returns true if the ranges are empty, false otherwise
     */
    bool empty() const {
        return forwardRange.empty();
    }

    length_t width() const {
        return forwardRange.width();
    }

    /**
     * Operator overloading
     * @returns true if this is equal to rhs
     */
    bool operator==(const RangePair& o) const {
        // only the first range matters as the two ranges imply each other
        return o.getForwardRange() == forwardRange;
    }

    friend std::ostream& operator<<(std::ostream& os, const RangePair& r);
};

class BiFMPos : public FMPos {
  protected:
    Range forwardRange;

  public:
    BiFMPos() : FMPos(), forwardRange() {
    }

    BiFMPos(const RangePair ranges, length_t depth)
        : FMPos(ranges.getBackwardRange(), depth),
          forwardRange(ranges.getForwardRange()) {
    }

    bool operator==(const BiFMPos& rhs) const {
        return getRanges() == rhs.getRanges() && depth == rhs.getDepth();
    }

    const RangePair getRanges() const {
        return RangePair(range, forwardRange);
    }
};

class BiFMPosExt : public FMPosExt {
  private:
    Range forwardRange;

  public:
    /**
     * Create a node of the search tree
     * @param character the character of this node
     * @param range the range over the suffix array
     * @param row the row of this node in thea lignment matrix = depth of
     * this node
     */
    BiFMPosExt(char character, RangePair ranges, length_t row)
        : FMPosExt(character, ranges.getBackwardRange(), row),
          forwardRange(ranges.getForwardRange()) {
    }

    /**
     * Default constructor, this Node will have empty range
     */
    BiFMPosExt() : FMPosExt(), forwardRange() {
    }

    /**
     * Gets the range of this node
     * @returns the range of this node
     */
    const RangePair getRanges() const {
        return RangePair(range, forwardRange);
    }
};

class BiFMOcc : public FMOcc {
  private:
    Range forwardRange;

  public:
    BiFMOcc() : FMOcc(), forwardRange(){};
    /**
     * Make an approximate match in the suffix array
     * @param range the range of this approximate match (range over SA)
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param depth the depth (=length) of this approximate match
     */
    BiFMOcc(RangePair ranges, length_t distance, length_t depth)
        : FMOcc(ranges.getBackwardRange(), distance, depth),
          forwardRange(ranges.getForwardRange()) {
    }
    /**
     * Make an approximate match in the suffix array
     * @param pos, the position in the FMIndex of this approximate match
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    BiFMOcc(BiFMPos pos, length_t distance)
        : FMOcc(pos, distance),
          forwardRange(pos.getRanges().getForwardRange()) {
    }

    const RangePair getRanges() const {
        return RangePair(pos.getRange(), forwardRange);
    }
};

// ============================================================================
// CLASS SEARCH: PROVIDED STEP 3
// ============================================================================
class Search {
  private:
    std::vector<length_t> lowerBounds; // the vector with lower bounds
    std::vector<length_t> upperBounds; // the vector with upper bounds
    std::vector<length_t> order;       // the vector with the order of the parts

    std::vector<Direction>
        directions; // the directions of each phase (calculated from the order)

    Search(std::vector<length_t>& order, std::vector<length_t>& lowerBounds,
           std::vector<length_t>& upperBounds,
           std::vector<Direction>& directions)
        : lowerBounds(lowerBounds), upperBounds(upperBounds), order(order),
          directions(directions) {
    }

    void sanityCheck() const {
        if (!connectivitySatisfied())
            throw std::runtime_error("Search does not satisfy connectivity "
                                     "property!");
        if (!zeroBased())
            throw std::runtime_error("Search is not zero based!");
        if (!noDecreasingInBounds())
            throw std::runtime_error("Illegal bounds in search!");
    }

  public:
    /**
     * Static function to construct a search. The directions and switches of
     * the search are calculated
     * @param order, the order of the search
     * @param lowerBounds, the lower bounds of the search
     * @param upperBounds, the upper bounds of the search
     */
    static Search makeSearch(std::vector<length_t> order,
                             std::vector<length_t> lowerBounds,
                             std::vector<length_t> upperBounds) {
        // check correctness of sizes
        if (order.size() != lowerBounds.size() ||
            order.size() != upperBounds.size()) {
            throw std::runtime_error("Could not create search, the sizes of "
                                     "all vectors are not equal");
        }

        // compute the directions
        std::vector<Direction> directions;
        directions.reserve(order.size());
        directions.push_back((order[1] > order[0]) ? FORWARD : BACKWARD);

        for (length_t i = 1; i < order.size(); i++) {
            Direction d = (order[i] > order[i - 1]) ? FORWARD : BACKWARD;
            directions.push_back(d);
        }

        Search s = Search(order, lowerBounds, upperBounds, directions);
        s.sanityCheck();
        return s;
    }

    /**
     * Sets the directions of the parts to the directions of the search
     * @param parts, the parts to set the direction of
     */
    void setDirectionsInParts(std::vector<Substring>& parts) const {
        // set the directions for the parts
        for (length_t i = 0; i < order.size(); i++) {
            parts[order[i]].setDirection(directions[i]);
        }
    }

    /**
     * @returns the lower bound for the idx'th part
     */
    length_t getLowerBound(length_t idx) const {
        assert(idx < lowerBounds.size());
        return lowerBounds[idx];
    }

    /**
     * @returns the upper bound for the idx'th part
     */
    length_t getUpperBound(length_t idx) const {
        assert(idx < upperBounds.size());
        return upperBounds[idx];
    }
    /**
     * @returns  the idx'th part
     */
    length_t getPart(length_t idx) const {
        assert(idx < order.size());
        return order[idx];
    }

    /**
     * @returns the direction for the idxith part
     */
    Direction getDirection(length_t idx) const {
        assert(idx < directions.size());
        return directions[idx];
    }

    /**
     * Get the number of parts in this search
     * @return the number of parts
     */
    length_t getNumParts() const {
        return order.size();
    }

    /**
     * Checks if the idxth part is the final part of the search
     * @returns true if idx is the final part of the search, false otherwise
     */
    bool isEnd(length_t idx) const {
        return idx == order.size() - 1;
    }

    /**
     * Checks if the connectivity property is satisfied
     * @returns true if the property is satisfied
     */
    bool connectivitySatisfied() const {
        length_t highestSeen = order[0];
        length_t lowestSeen = order[0];
        for (length_t i = 1; i < order.size(); i++) {
            if (order[i] == highestSeen + 1) {
                highestSeen++;
            } else if (order[i] == lowestSeen - 1) {
                lowestSeen--;
            } else {
                return false;
            }
        }
        return true;
    }
    /**
     * Checks if the upper and lower bounds are not decreasing
     * @returns true if the bounds are valid
     */
    bool noDecreasingInBounds() const {
        for (length_t i = 1; i < order.size(); i++) {
            if (lowerBounds[i] < lowerBounds[i - 1]) {
                return false;
            }
            if (upperBounds[i] < upperBounds[i - 1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Check if the search is zero based (order must contain zero)
     * @returns  true if the search is zero based else false
     */
    bool zeroBased() const {
        return *std::min_element(order.begin(), order.end()) == 0;
    }
};

class BiFMIndex : public FMIndex {
  private:
    CumulativeBitvectors<ALPHABET> reverseOccTable;
    CumulativeBitvectors<ALPHABET> originalOccTable;

    // search direction variables
    Direction dir;

    void read(const std::string& base, bool verbose);

  public:
    BiFMIndex(const std::string& base, int sa_sparse = 1, bool verbose = true)
        : FMIndex(base, sa_sparse, verbose) {
        read(base, verbose);
    }

    /**
     * Create the BWT of the reverse text from the SA of the reverse text and
     * the text
     * @param sa the (dense) suffix array of the reversed text
     * @param revBWT [output] the bwt of the reversed text
     */
    void createRevBWTFromRevSA(const std::vector<length_t>& revSA,
                               std::string& revBWT);

    // ============================================================================
    // Functionality
    // ============================================================================

    /**
     * Extends the patern with one character  to the right (find ranges of Pc
     * using ranges of P)
     * @param charIdx the position in alphabet of the character
     * that is added in the front
     * @param originalRanges the ranges of pattern P
     * @param newRanges the ranges of Pc  [output]
     * @return true if the newRanges are not empty and false otherwise
     */
    bool addCharRight(length_t charIdx, const RangePair& originalRanges,
                      RangePair& newRanges) const;

    /**
     * Extends the patern with one character  to the left (find ranges of cP
     * using ranges of P)
     * @param charIdx the position in alphabet of the character
     * that is added in the back
     * @param originalRanges the ranges of pattern P
     * @param newRanges ranges cP  [output]
     * @return true if the newRanges are not empty and false otherwise
     */
    bool addCharLeft(length_t charIdx, const RangePair& originalRanges,
                     RangePair& newRanges) const;

    /**
     * Creates all child positions of the position and pushes them on the
     * stack
     * @param range, the range of the position to get the children of
     * @param depth, the depth of the position to get the children of
     * @param stack, the stack to push the children on
     */
    void extendFMPos(const RangePair& ranges, const length_t& depth,
                     std::vector<BiFMPosExt>& stack);

    /**
     * Creates all child positions of the position and pushes them on the stack
     * @param pos, the position to get the children of
     * @param stack, the stack to push the children on
     */
    void extendFMPos(const BiFMPosExt& pos, std::vector<BiFMPosExt>& stack) {
        extendFMPos(pos.getRanges(), pos.getDepth(), stack);
    }

    /**
     * This function matches a string exactly starting form startRange while
     * keeping track of the ranges in both directions
     * @param string the string to match
     * @param ranges, the start ranges to search, if
     * this range is empty the procedure must search in the whole index
     * @returns the pair of ranges that matches string
     */
    RangePair matchExactBidirectionally(const Substring& str,
                                        RangePair ranges) const;

    /**
     * This function matches a string exactly starting from the empty string,
     * while keeping track of the ranges in both directions
     * @param string the string to match
     * @returns the pair of ranges that matches string
     */
    RangePair matchExactBidirectionally(const Substring& str) const {
        return matchExactBidirectionally(
            str, RangePair(0, textLength, 0, textLength));
    }

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    void setDirection(Direction d) {
        dir = d;
    }

    // ============================================================================
    // INTEGRATION
    // ============================================================================
    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using edit distance metric
     * @param s, the search to follow
     * @param startOcc, the approximate occurrence found for all previous parts
     * of the search
     * @param occ, a vector with occurrences of the complete search, if  such
     * an occurrence is found i,t will be added to this datastructure
     * @param parts the parts of the pattern, with correct direction
     * @param idx, the index of the part to match
     */
    void recApproxMatch(const Search& s, const BiFMOcc& startOcc,
                        std::vector<FMOcc>& occ,
                        const std::vector<Substring>& parts, const int& idx);
};
#endif