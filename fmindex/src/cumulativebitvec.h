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

#ifndef BWTREPR_H
#define BWTREPR_H

#include <array>
#include <cstdlib>
#include <vector>

#include "alphabet.h"
#include "bitvec.h"

// ============================================================================
// CLASS BWT REPRESENTATION (supports occ(c,k) and cumulocc(c,k) in O(1) time)
// ============================================================================

template <size_t S>          // S is the size of the alphabet (including '$')
class CumulativeBitvectors { // e.g. S = 5 for DNA (A,C,G,T + $)

  private:
    // The '$' character (cIdx == 0) is not encoded in the bitvector.
    // Hence, we use only S-1 bitvectors.

    std::array<Bitvec, S - 1> bvs; // bitvector representations of the BWT
    size_t dollarPos;              // position of the dollar sign

  public:
    /**
     * Default constructor
     */
    CumulativeBitvectors() {
    }

    /**
     * Constructor
     * @param sigma Alphabet
     * @param BWT Burrows-Wheeler transformation
     */
    CumulativeBitvectors(const Alphabet<S>& sigma, const std::string& BWT)
        : dollarPos(BWT.size()) {
        // Initialize the bitvectors with size BWT.size()
        for (auto& bv : bvs)
            bv = Bitvec(BWT.size() + 1);

        for (size_t i = 0; i < BWT.size(); i++) {
            if (BWT[i] == '$') {
                dollarPos = i;
                continue;
            }

            // The $-character (cIdx == 0) is not encoded in the bitvector.
            // Hence, use index cIdx-1 in the bitvector.
            for (size_t cIdx = sigma.c2i(BWT[i]); cIdx < S; cIdx++)
                bvs[cIdx - 1][i] = true;
        }

        // index the bitvectors
        for (auto& bv : bvs)
            bv.index();
    }

    /**
     * Get occurrence count of character c in the range BWT[0...j[
     * @param cIdx Character index
     * @param j index
     * @return occ(c, j)
     */
    size_t occ(int cIdx, size_t j) const {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        // 3 - 6 lines of code
        throw std::runtime_error("occ has not been implemented yet!");
    }

    /**
     * Get cumulative occurrence count of characters SMALLER than c
     * in the range BWT[0...j[
     * @param cIdx Character index
     * @param j index
     * @return cumulocc(c, j)
     */
    size_t cumulocc(int cIdx, size_t j) const {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        // 3 - 6 lines of code
        throw std::runtime_error("cumulocc has not been implemented yet!");
    }
};

#endif