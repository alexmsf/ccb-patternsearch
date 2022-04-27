/***************************************************************************
 *   Copyright (C) 2020 Jan Fostier (jan.fostier@ugent.be)                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "matrix.h"

using namespace std;

/**
 * Write program usage information to the standard output
 */
void printUsage()
{
        cout << "Usage: affineWRONG input.fasta\n\n";
        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>" << endl;
}

/**
 * Read sequences from a FASTA file
 * @param filename FASTA file filename
 * @param sequences Vector of sequences (output)
 */
void readSequences(const string& filename, vector<string>& sequences)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        string line;
        while (ifs) {
                getline(ifs, line);
                if (line.empty())
                        continue;
                if (line.front() == '>') {
                        sequences.push_back(string());
                        continue;
                }

                sequences.back().append(line);
        }
}

/**
 * Perform INCORRECT global alignment of two sequences using affine gap penalities
 * @param X sequence one
 * @param Y sequence two
 */
void alignAffineWRONG(const string& X, const string& Y)
{
        const int Go = -3;
        const int Ge = -1;
        const int M = 1;
        const int I = -1;

        size_t m = X.length();
        size_t n = Y.length();

        // Initialize an (m+1) x (n+1) matrix T to zero
        // This matrix encodes for each cell (i,j) how it was computed
        // using two bits: 01b (top), 10b (left)
        // For example: value 11b means the cell was computed from
        // either its top and left neighbor
        Matrix T(m+1, n+1);
        for (size_t i = 0; i <= m; i++)
                for (size_t j = 0; j <= n; j++)
                        T(i, j) = 0;

        // initialize an (m+1) x (n+1) matrix S
        Matrix S(m+1, n+1);
        S(0, 0) = 0;

        // initialize first column
        for (size_t i = 1; i <= m; i++) {
                S(i, 0) = Go + i * Ge;
                T(i, 0) |= 0x1;
	}

        // initialize first row
        for (size_t j = 1; j <= n; j++) {
                S(0, j) = Go + j * Ge;
                T(0, j) |= 0x2;
        }

        // fill in the rest of the matrix
        for (size_t i = 1; i <= m; i++) {
                for (size_t j = 1; j <= n; j++) {
                        int diag = S(i-1, j-1) + (X[i-1] == Y[j-1] ? M : I);
                        int gapX = S(i, j-1) + Ge + (T(i, j-1) & 0x1 ? 0 : Go);
                        int gapY = S(i-1, j) + Ge + (T(i-1, j) & 0x2 ? 0 : Go);
                        S(i, j) = max(max(diag, gapX), gapY);
                        if (S(i, j) == gapX) T(i, j) |= 0x1;
                        if (S(i, j) == gapY) T(i, j) |= 0x2;
                }
        }

        for (size_t i = 0; i <= m; i++) {
                for (size_t j = 0; j <= n; j++)
                        cout << S(i, j) << "\t";
		cout << endl;
	}

        // create an alignment
        string alX, alY, mid;

        int i = (int)X.size();
        int j = (int)Y.size();

        while (i > 0 || j > 0) {
                if ((i > 0) && (S(i, j) == S(i-1, j) + Ge + (T(i-1, j) & 0x2 ? 0 : Go))) {
                        alX.push_back(X[i-1]);
                        alY.push_back('-');
                        mid.push_back(' ');
                        i--;
                } else if ((j > 0) && (S(i, j) == S(i, j-1) + Ge + (T(i, j-1) & 0x1 ? 0 : Go))) {
                        alX.push_back('-');
                        alY.push_back(Y[j-1]);
                        mid.push_back(' ');
                        j--;
                } else {
                        alX.push_back(X[i-1]);
                        alY.push_back(Y[j-1]);
                        char c = (X[i-1] == Y[j-1]) ? '|' : '*';
                        mid.push_back(c);
                        i--;
                        j--;
                }
        }

        reverse(alX.begin(), alX.end());
        reverse(alY.begin(), alY.end());
        reverse(mid.begin(), mid.end());

        for (size_t i = 0; i < alX.size(); i += 80) {
                cout << alX.substr(i, 80) << "\n"
                     << mid.substr(i, 80) << "\n"
                     << alY.substr(i, 80) << "\n\n";
        }
        cout << "Alignment score: " << S(m, n) << endl;
}

int main(int argc, char** argv)
{
        if (argc != 2) {
                printUsage();
                return EXIT_FAILURE;
        }

        vector<string> sequences;
        readSequences(argv[1], sequences);

        if (sequences.size() != 2) {
                cerr << "Input FASTA file should contain only two sequences\n";
                return EXIT_FAILURE;
        }

        alignAffineWRONG(sequences[0], sequences[1]);

        return EXIT_SUCCESS;
}
