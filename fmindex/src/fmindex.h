#ifndef FMINDEX_H
#define FMINDEX_H
#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

#include "alphabet.h"
#include "substring.h"
#include "suffixarray.h"

// ============================================================================
// IO helper functions
// ============================================================================

bool readText(const std::string& filename, std::string& buf);
/**
 * Read a binary file and stores content in array
 * @param filename File name
 * @param array Suffix array (contents will be overwritten)
 * @returns True if successful, false otherwise
 */
bool readArray(const std::string& filename, std::vector<length_t>& array);

void readSA(const std::string& filename, std::vector<length_t>& sa,
            size_t saSizeHint);

// ============================================================================
// CLASS RANGE
// ============================================================================

class Range {
  private:
    length_t begin; // beginning of the range
    length_t end;   // end of the range (non-inclusive)

  public:
    /**
     * Constructor
     * @param b, the beginning of the range
     * @param e, the end of the range (non-inclusive)
     */
    Range(length_t b, length_t e) : begin(b), end(e) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    Range() : begin(0), end(0) {
    }

    length_t getBegin() const {
        return begin;
    }
    length_t getEnd() const {
        return end;
    }
    /**
     * Check if this range is empty
     * @returns true if the range is empty, false otherwise
     */
    bool empty() const {
        return end <= begin;
    }

    /**
     * Gets the width of the range (end - begin)
     * @returns the width of this range
     */
    length_t width() const {
        return (empty()) ? 0 : end - begin;
    }

    /**
     * Operator overloading, two ranges are equal if their begin and end field
     * are equal, or if they are both empty
     */
    bool operator==(const Range& o) const {
        return (o.getBegin() == begin && o.getEnd() == end) ||
               (o.empty() && empty());
    }

    /**
     * Operator overloading
     */
    friend std::ostream& operator<<(std::ostream& os, const Range& r);
};

// ============================================================================
// CLASS FMPos
// ============================================================================

/**
 * A position in the bidirectional FM-index.
 */
class FMPos {
  protected:
    Range range;    // the range over the suffix arrays
    length_t depth; // the depth of the prefix of the suffixes of this position

  public:
    /**
     * Default constructor for empty position (= empty ranges and depth of
     * zero)
     */
    FMPos() : range(Range()), depth(0) {
    }

    /**
     * Constructor
     */
    FMPos(const Range range, length_t depth) : range(range), depth(depth) {
    }

    const Range& getRange() const {
        return range;
    }

    const length_t& getDepth() const {
        return depth;
    }

    /**
     * Operator overloading, two FMPos are equal if their ranges and depth
     * are equal
     * @param rhs the FMPos to compare to this
     * @returns true if this is equal to rhs
     */
    bool operator==(const FMPos& rhs) const {
        return getRange() == rhs.getRange() && depth == rhs.getDepth();
    }
    /**
     * @returns true if the ranges are not empty, false otherwise
     */
    bool isValid() const {
        return !range.empty();
    }

    /**
     * Operator overloading
     */
    friend std::ostream& operator<<(std::ostream& os, const FMPos& r);
};

// ============================================================================
// CLASS FMOcc
// ============================================================================

/**
 * An occurrence (match) in the  FM-index
 */
class FMOcc {
  protected:
    FMPos pos;         // The FM position of this occurrence
    length_t distance; // the edit distance

  public:
    FMOcc() : pos(), distance(0) {
    }
    /**
     * Make an approximate match in the suffix array
     * @param range the range of this approximate match (range over SA)
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param depth the depth (=length) of this approximate match
     */
    FMOcc(Range range, length_t distance, length_t depth)
        : pos(range, depth), distance(distance) {
    }
    /**
     * Make an approximate match in the suffix array
     * @param pos, the position in the FMIndex of this approximate match
     * @param distance the (edit or hamming) distance of this approximate
     * match
     */
    FMOcc(FMPos pos, length_t distance) : pos(pos), distance(distance) {
    }

    const Range& getRange() const {
        return pos.getRange();
    }
    const length_t& getDistance() const {
        return distance;
    }

    const length_t& getDepth() const {
        return pos.getDepth();
    }

    const FMPos& getFMPos() const {
        return pos;
    }

    length_t getWidth() const {
        return pos.getRange().width();
    }

    /**
     * @returns true if the position is valid, false otherwise
     */
    bool isValid() const {
        return pos.isValid();
    }
    /**
     * Operator overloading to sort FMOcc
     * First the FMOcc are sorted on the begin of the range over the suffix
     * array of their position Then they are sorted on their distance score.
     * Lastly they are sorted on their depth.
     * @param rhs the FMOcc to compare to this
     * @returns true if this is smaller than rhs
     */
    bool operator<(const FMOcc& rhs) const {
        if (getRange().getBegin() != rhs.getRange().getBegin()) {
            return getRange().getBegin() < rhs.getRange().getBegin();
        }
        if (distance != rhs.getDistance()) {
            // begin is equal, better ed is smarter
            return distance < rhs.getDistance();
        }
        // shorter read is smaller...
        return getDepth() < rhs.getDepth();
    }
    /**
     * Operatoroverloading
     * Two FMocc are equal if their ranges, distance and depth are all equal
     * @param returns true if this is equal to rhs
     */
    bool operator==(const FMOcc& rhs) {
        return pos == rhs.getFMPos() && distance == rhs.getDistance();
    }

    /**
     * Operator overloading
     */
    friend std::ostream& operator<<(std::ostream& os, const FMOcc& r);
};

// ============================================================================
// CLASS FMPosExt
// ============================================================================
/**
 * A position in the FM index extended with information about the character of
 * this position
 */
class FMPosExt : public FMPos {
  private:
    char c; // the character of this position
  public:
    /**
     * Create a node of the search tree
     * @param character the character of this node
     * @param range the range over the suffix array
     * @param row the row of this node in the alignment matrix = depth of
     * this node
     */
    FMPosExt(char character, Range range, length_t row)
        : FMPos(range, row), c(character) {
    }

    /**
     * Default constructor, this Node will have empty range
     */
    FMPosExt() : FMPos(), c(char(0)) {
    }

    /**
     * Gets the range of this node
     * @returns the range of this node
     */
    const Range& getRange() const {
        return range;
    }

    /**
     * Get the character of this node
     * @returns the character of this node
     */
    const char getCharacter() const {
        return c;
    }

    /**
     * Get the row of this node
     * @returns the row of this node
     */
    length_t getRow() const {
        return depth;
    }
};

// ============================================================================
// CLASS TextOccurrence
// ============================================================================

/**
 * An occurrence in the text. Its range corresponds to the begin and end
 * position (non-inclusive) in the text. The distance is the (edit) distance
 * with the original pattern.
 */
class TextOcc {
  private:
    Range range;       // the range in the text
    length_t distance; // the distance to this range (edit or hamming)

  public:
    /**
     * Constructor
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     */
    TextOcc(Range range, length_t distance) : range(range), distance(distance) {
    }

    /**
     * Constructor
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     * @param CIGAR the CIGAR string of the match
     */
    TextOcc(Range range, length_t distance,
            std::vector<std::pair<char, uint>>& CIGAR)
        : range(range), distance(distance) {
    }

    /**
     * Constructor for an invalid text occurrence (empty range)
     */
    TextOcc() : range(0, 0) {
    }

    const Range getRange() const {
        return range;
    }
    const length_t getDistance() const {
        return distance;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance and finally on their length
     */
    bool operator<(const TextOcc& r) {

        if (range.getBegin() != r.getRange().getBegin()) {
            return range.getBegin() < r.getRange().getBegin();
        }
        // begin is equal, better ed is smarter
        if (distance != r.getDistance()) {
            return distance < r.getDistance();
        }
        // shorter read is smaller...
        return range.width() < r.getRange().width();
    }

    bool operator==(const TextOcc& r) const {
        return r.getRange() == range && r.getDistance() == distance;
    }

    bool isValid() const {
        return !range.empty();
    }

    length_t width() const {
        return range.width();
    }

    length_t begin() const {
        return range.getBegin();
    }

    length_t end() const {
        return range.getEnd();
    }

    friend std::ostream& operator<<(std::ostream& os, const TextOcc& r);
};

// ============================================================================
// CLASS FMINDEX: PROVIDED STEP 1/2/3 (ADAPATED FOR EACH VERSION)
// ============================================================================

class FMIndex {
  protected:
    std::string bwt;                       // the bwt of the text
    std::string text;                      // the original text
    length_t textLength;                   // the length of the text
    std::array<length_t, ALPHABET> counts; // the counts array
    SparseSuffixArray sparseSA; // the suffix array of the reference genome
    Alphabet<ALPHABET> sigma;   // the alphabet

    std::vector<Bitvec> occTable; // the occurrence table
    length_t dollarPos;           // the position of the dollar in the BWT

    // ============================================================================
    // FM Index Construction
    // ============================================================================

    /**
     * Helper function for constructor, reads in the text and suffix array and
     * creates a sparse suffix array, the occTable, the BWT and the counts
     */
    void read(const std::string& base, bool verbose);

  public:
    // ============================================================================
    // FM Index Construction
    // ============================================================================
    /**
     * Constructor
     */
    FMIndex(const std::string& base, int sa_sparse = 1, bool verbose = true)
        : sparseSA(sa_sparse) {
        read(base, verbose);
    }

    /**
     * Create the BWT from the SA and the text
     * @param sa the (dense) suffix array
     */
    void createBWTFromSA(const std::vector<length_t>& sa);

    /**
     * Create the counts vector. Where counts[c] is equal to the number of
     * characters smaller than c in the text
     */
    void createCounts();

    // ============================================================================
    // ACCESSING DATA STRUCTURE
    // ============================================================================
    const std::string& getText() const {
        return text;
    }

    const Alphabet<ALPHABET>& getAlphabet() const {
        return sigma;
    }

    /**
     * Get a substring in the text corresponding to an occurrence in the text
     * @param occ the occurrence in the text
     */
    std::string getSubstr(const TextOcc& occ) const {
        return text.substr(occ.getRange().getBegin(), occ.getRange().width());
    }

    const std::string& getBWT() const {
        return bwt;
    }

    const std::array<length_t, ALPHABET>& getCounts() const {
        return counts;
    };

    /**
     * Takes the reverse complement
     * @param s the string to take the reverse complement of
     * @returns the reverse complement of s
     */
    std::string revCompl(const std::string& s) const {
        std::string r;
        r.reserve(s.size());
        for (auto it = s.crbegin(); it != s.crend(); it++)
            r += sigma.i2c(5 - (sigma.c2i(*it)));
        return r;
    }

    // ============================================================================
    // FMIndex functionalily: week 1
    // ============================================================================

    /**
     * Occ function, selects the correct bitvector and applies a rank operation
     * to find the number of occurrences of the character before the index
     * @param charIdx the index of the character in the alphabet to do the occ
     * query on
     * @param index the index for the occ query
     */
    length_t occ(const length_t& charIdx, const length_t& index) const;
    /**
     * Finds the LF mapping of the character at index k in the bwt string
     * @param k the index to find the LF mapping off
     * @returns the row that is the LF mapping of k.
     */
    length_t findLF(length_t k) const;

    /**
     * Finds the entry in the suffix array of this index.
     * @param k the index to find the entry in the SA off
     * @returns the entry in the SA of the index
     */
    length_t findSA(length_t k) const;

    /**
     * Finds the range of cP using the range of P over the SA
     * @param charIdx the position in alphabet of c
     * @param originalRange the range over the SA of pattern P
     * @param newRange the range over the SA of cP   [output]
     * @return true if the newRange is not empty and false otherwise
     */
    bool addCharLeft(length_t charIdx, const Range& originalRange,
                     Range& newRange) const;

    // ============================================================================
    // INTEGRATION WEEK 1
    // ============================================================================

    /**
     * This function matches a string exactly
     * @param str the string to match
     * @returns the sorted start positions of the exact matches of str in the
     * text
     */
    std::vector<length_t> matchExact(const std::string& str) const;

    /**
     * Finds the best paired match of a pair of reads, given the insertion size.
     * @param reads  a pair of reads to be matched, one against the forward
     * strand and one against the backward strand.
     * @param insSize the average distance between the two extreme ends of the
     * reads closest together.
     */
    std::tuple<length_t, length_t, bool>
    bestPairedMatch(const std::pair<std::string, std::string>& reads,
                    const length_t& insSize) const;

    // ============================================================================
    // FUNCTIONALITY WEEK 2
    // ============================================================================

    /**
     * Creates all child positions of the position and pushes them on the
     * stack
     * @param pos, the position to get the children of
     * @param stack, the stack to push the children on
     */
    void extendFMPos(const FMPos& pos, std::vector<FMPosExt>& stack) const {
        extendFMPos(pos.getRange(), pos.getDepth(), stack);
    }

    /**
     * Creates all child positions of the position and pushes them on the stack
     * @param range, the range of the position to get the children of
     * @param depth, the depth of the position to get the children of
     * @param stack, the stack to push the children on
     */
    void extendFMPos(const Range& range, const length_t& depth,
                     std::vector<FMPosExt>& stack) const;

    void convertFMOccToTextOcc(const FMOcc& fmocc,
                               std::vector<TextOcc>& textOcc) const;

    // ============================================================================
    // INTEGRATION WEEK 2
    // ============================================================================

    /**
     * Matches the pattern approximately. All matches are at most a certain
     * edit distance away from the pattern
     * @param pattern the pattern to match
     * @param k the maximum edit distance
     * @returns a vector with occurrences in the text
     */
    std::vector<TextOcc> naiveApproxMatch(const std::string& pattern,
                                          length_t k) const;

    /**
     * Helper function wich filters out redundant matches
     * @param fmocc a vector with occurrences in the fmindex
     * @param k the maximal allowed edit distance
     * @returns a vector with all non-redundant occurrence in the text
     */
    std::vector<TextOcc> filterRedundantMatches(std::vector<FMOcc>& fmocc,
                                                const length_t& k) const;
};

#endif