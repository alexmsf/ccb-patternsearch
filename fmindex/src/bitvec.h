/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef BITVEC_H
#define BITVEC_H

/**
 * The implementation implements the rank9 algorithm as described in
 * S. Vigna, "Broadword Implementation of Rank/Select Queries", WEA 2008
 * It relies on GCC's __builtin_popcountll, so please build this software
 * using the -mpopcnt flag to enable the SSE 4.2 POPCNT instruction.
 */

#include <string.h>

#include <cassert>
#include <cstdint>
#include <fstream>
#include <vector>

// ============================================================================
// BIT REFERENCE CLASS
// ============================================================================

class Bitref {
  private:
    uint64_t& wordRef; // reference to a word in the bitvector
    uint64_t bitmask;  // bitmask of the form (1 << bitIdx)

  public:
    /**
     * Constructor
     * @param wordRef Reference to a word in the bitvector
     * @param bitmask Bitmask of the form (1 << bitIdx)
     */
    Bitref(uint64_t& wordRef, uint64_t bitmask)
        : wordRef(wordRef), bitmask(bitmask) {
    }

    /**
     * Set bitref to particular value
     * @param val Target value
     * @return Bitref reference after modification
     */
    const Bitref& operator=(bool val) {
        if (val)
            wordRef |= bitmask;
        else
            wordRef &= ~bitmask;
        return *this;
    }

    /**
     * Set bitref to another bitref
     * @param br Another bitref
     * @return Bitref reference
     */
    const Bitref& operator=(const Bitref& br) {
        return this->operator=(bool(br));
    }

    /**
     * Bool conversion operator
     */
    operator bool() const {
        return (wordRef & bitmask) != 0;
    }
};

// ============================================================================
// BIT VECTOR CLASS
// ============================================================================

class Bitvec {
  private:
    uint64_t N;                   // size of the bitvector
    std::vector<uint64_t> bv;     // actual bitvector
    std::vector<uint64_t> counts; // interleaved 1st and 2nd level counts

  public:
    /**
     * Get a bit at a certain position
     * @param p Position
     * @return true or false
     */
    bool operator[](uint64_t p) const {
        assert(p < N);
        uint64_t w = p / 64;
        uint64_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at a certain position
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator[](uint64_t p) {
        assert(p < N);
        uint64_t w = p / 64;
        uint64_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        counts = std::vector<uint64_t>((bv.size() + 7) / 4, 0ull);

        uint64_t countL1 = 0, countL2 = 0;
        for (uint64_t w = 0, q = 0; w < bv.size(); w++) {
            if (w % 8 == 0) { // store the L1 counts
                countL1 += countL2;
                counts[q] = countL1;
                countL2 = __builtin_popcountll(bv[w]);
                q += 2;
            } else { // store the L2 counts
                counts[q - 1] |= (countL2 << (((w % 8) - 1) * 9));
                countL2 += __builtin_popcountll(bv[w]);
            }
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos p)
     * @param p Position
     */
    uint64_t rank(uint64_t p) const {
        // 3 - 5 lines of code
        assert(p < N);
        return firstLevelCounts(p/64)+secondLevelCounts(p/64)+popcount(p/64, p);
    }

    /**
     * Get the first level count
     * @param w the word index to get the first level count of
     */
    uint64_t firstLevelCounts(uint64_t w) const {
        return counts[(w / 8) * 2];
    }

    /**
     * Get the second level count
     * @param w the word index to get the second level count of
     */
    uint64_t secondLevelCounts(uint64_t w) const {
        uint64_t q = (w / 8) * 2; // counts index
        int64_t t = (w % 8) - 1;
        return counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;
    }

    /**
     * Count the number of 1 bits in word w[0...b)
     * @param w the word index
     * @param b the bit offset
     */
    uint64_t popcount(uint64_t w, uint64_t b) const {
        return __builtin_popcountll((bv[w] << 1) << (63 - b));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv.data(), bv.size() * sizeof(uint64_t));
        ofs.write((char*)counts.data(), counts.size() * sizeof(uint64_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));

        bv.resize((N + 63) / 64);
        ifs.read((char*)bv.data(), bv.size() * sizeof(uint64_t));

        counts.resize((bv.size() + 7) / 4);
        ifs.read((char*)counts.data(), counts.size() * sizeof(uint64_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    uint64_t size() const {
        return N;
    }

    /**
     * Default constructor, move constructor and move assignment operator
     */
    Bitvec() : N(0){};
    Bitvec(Bitvec&& rhs) = default;
    Bitvec& operator=(Bitvec&& rhs) = default;

    /**
     * Deleted copy constructor and copy assignment operator
     */
    Bitvec(const Bitvec&) = delete;
    Bitvec& operator=(const Bitvec&) = delete;

    /**
     * Constructor
     * @param N Number of bits in the bitvector
     */
    Bitvec(uint64_t N) : N(N), bv((N + 63) / 64, 0ull) {
    }
};

#endif
