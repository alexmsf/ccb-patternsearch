#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include "bitvec.h"
#include <fstream>
#include <iostream> // used for printing
#include <stdint.h>
#include <string>
#include <vector>

typedef uint32_t length_t;

class SparseSuffixArray {
  private:
    length_t sparsenessFactor; // the sparseness factor
    Bitvec
        bitvector; // Only necessary if every ith suffix is stored (comment out
                   // or delete if every ith entry in suffix array is stored)
    std::vector<length_t> sparseSA; // the sparse suffix array

  public:
    /**
     * Returns true if the value at index i is stored in the sparse suffix
     * array. This depends on which way you chose to sparsify the suffix array
     * @param i the index to check
     * @return true if the value is stored
     */
    bool hasStored(const length_t i) const {
        // 1 line of code
        return ((i % sparsenessFactor) == 0);
    }

    /**
     * Get the value at index i in the suffix array. This depends on the way you
     * chose to sparsify the suffix array. Only call this method if hasStored(i)
     * returned true.
     * @param i the index to get the value of
     * @return the value found at index i in the original suffix array.
     */
    length_t operator[](const length_t i) const {
        assert(hasStored(i));
        // 1 line of code
        return sparseSA[i/sparsenessFactor];
    }

    /**
     * Creates the sparse suffix array from the original suffix array.
     * @param sa the original suffix array
     */
    void createSparseSA(const std::vector<length_t>& sa) {
        sparseSA.resize(sa.size());
         for (length_t i = 0; i < sa.size(); i++) {
             // 2 - 3 lines of code
             if((i % sparsenessFactor) == 0){
               sparseSA[i/sparsenessFactor] = sa[i];
             }
         }
    }

    SparseSuffixArray(length_t sparseness) : sparsenessFactor(sparseness) {
    }
};

#endif
