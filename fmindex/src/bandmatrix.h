#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include "substring.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

typedef uint32_t length_t;

class BandedMatrix {
  private:
    std::vector<length_t> matrix;
    length_t W; // The off diagonal width of this bandmatrix
    length_t m; // number of rows
    length_t n; // number of columns

    length_t colPerRow;

    void initializeMatrix(length_t startValue) {

        // initialize the top row and leftmost column
        for (length_t i = 0; i <= W + 1; i++) {
            operator()(0, i) = i + startValue;
            operator()(i, 0) = i + startValue;
        }

        // set max elements at sides
        // first the elements on rows [1, W]
        for (length_t i = 1; i <= W; i++) {
            // right of band
            operator()(i, i + W + 1) = W + 1 + startValue;
        }

        // then the elements on rows [W + 1, x]
        for (length_t i = W + 1; i + W + 1 < n; i++) {
            // right of band
            operator()(i, i + W + 1) = W + 1 + startValue;
            // left of band
            operator()(i, i - (W + 1)) = W + 1 + startValue;
        }

        // finally the elements on the final rows
        for (length_t i = std::max<int>((int)n - (W + 1), W + 1); i < m; i++) {
            // left of band
            operator()(i, i - (W + 1)) = W + 1 + startValue;
        }
    }

  public:
    /**
     * Constructor
     * @param patternsize, the size of the pattern to match, this will
     * initialize the top row
     * @param W, the maximal width
     * @param startValue the value in the origin
     */
    BandedMatrix(length_t patternsize, int W, int startValue)
        : W(W), n(patternsize + 1) {

        m = patternsize + W + 1;
        colPerRow = (2 * W + 1) + 2;
        matrix.resize(m * colPerRow);
        initializeMatrix(startValue);
        return;
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Element at position (i, j)
     */
    length_t operator()(length_t i, int j) const {
        return matrix[i * colPerRow + j - i + W];
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    length_t& operator()(length_t i, int j) {
        return matrix[i * colPerRow + j - i + W];
    }

    /**
     * Index the matrix
     * @param i Row index
     * @param j Column index
     * @return Element at position (i, j)
     */
    length_t at(length_t i, int j) const {
        return operator()(i, j);
    }

    /**
     * Index the matrix
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    length_t& at(length_t i, int j) {
        return operator()(i, j);
    }

    /**
     * Retrieves the first column that needs to be filled in for the row
     * @param row the row to fill in
     * @returns the first column to fill in
     */
    const int getFirstColumn(int row) const {
        // leftmost cell of band
        return std::max<int>(1, row - W);
    }
    /**
     * Retrieves the last column that needs to be filled in for the row
     * @param row the row to fill in
     * @returns the last column to fill in
     */
    const int getLastColumn(int row) const {
        // rightmost cell of band
        return std::min(n - 1, W + row);
    }

    /**
     * Get the number of rows of the matrix
     */
    length_t getNumberOfRows() const {
        return m;
    }

    /**
     * Update the matrix by calculating the element at position row, column.
     * @param notMatch whether the character was not a match
     * @param row the row of the element to update
     * @param column the column of the element to update
     * @returns the new value at row, column
     */
    length_t updateMatrixCell(bool notMatch, unsigned int row,
                              unsigned int column) {
        // 6 lines of code
        
        if(column==0) at(row, column) = row;
        else if(row==0) at(row, column) = column;
        else {
            int s = notMatch ? 1 : 0;
            at(row, column) = min(min(at(row-1, column-1) + s, at(row, column-1) + 1), at(row-1, column) + 1);
        }
        return at(row, column);
    }

    /**
     * Update the matrix by calculating the elements at row within the band
     * @param pattern the horizontal sequence (with direction correctly set)
     * @param row the row to update
     * @param c the character associated with this row
     * @returns minimal value found at this row
     */
    length_t updateMatrixRow(const Substring& pattern, length_t row, char c) {
        // 5 - 7 lines of code
        
        vector<length_t> rowValues;
        for(int j=0; j<=n; j++){
            bool match = pattern[j-1] == c ? true : false;
            rowValues.push_back(updateMatrixCell(match, row, j));
        }
        return distance(rowValues.begin(), min_element(rowValues.begin(), rowValues.end()));
    }

    /**
     * Indicates whether the row is within the band of the matrix in the
     * final column of the matrix
     * @param row the row to check
     * @returns true if the row is within the band of the matrix in the
     * final column
     */
    bool inFinalColumn(length_t row) const {
        return (length_t)getLastColumn(row) == n - 1;
    }

    /**
     * Finds the value in matrix(row, finalCol). Warning: only call this
     * function if the row is in the band for the final column
     * @param row
     * @returns matrix(row, final column)
     */
    length_t getValueInFinalColumn(length_t row) const {
        assert(inFinalColumn(row));
        return operator()(row, n - 1);
    }

    void printMatrix(length_t maxRow = 500) const {
        length_t mRow = std::min<length_t>(maxRow + 1, m);

        for (length_t i = 0; i < mRow; i++) {
            length_t firstCol = getFirstColumn(i), lastCol = getLastColumn(i);
            std::string rowNumber = ((i < 10) ? "0" : "") + std::to_string(i);
            std::string row = "row " + rowNumber + ": ";

            length_t startCol = 0;

            if (firstCol == 1 && i <= W) {
                row += std::to_string(operator()(i, 0)) + " ";
                startCol = 1;
            }

            for (length_t j = startCol; j < firstCol; j++) {
                row += "  ";
            }

            for (length_t j = firstCol; j <= lastCol; j++) {
                int number = operator()(i, j);
                row += std::to_string(number) + " ";
            }
            for (length_t j = lastCol + 1; j < n; j++) {
                row += "  ";
            }
            std::cout << row << std::endl;
        }
        std::cout << "----------------------------------------------\n";
    }
};

#endif