#ifndef SEARCHSCHEME_H
#define SEARCHSCHEME_H

#include "bidirectionalfmindex.h"
#include <sstream>

#define Pattern std::vector<int>

class SearchScheme {
  private:
    BiFMIndex& index; // reference to the index of the text that is searched
    std::string name; // name of the search scheme

    length_t maxED;
    std::vector<Search> searches;

    static Search makeSearchFromLine(const std::string& line) {

        std::stringstream ss(line);

        std::vector<std::string> tokens;
        std::string token;
        while (ss >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() != 3) {
            throw std::runtime_error("A search should have 3 vectors: order, "
                                     "lower bound and upper bound!");
        }

        // read each vector
        std::vector<std::vector<length_t>> vectors;

        for (const auto& token : tokens) {
            if (token.size() < 2)
                throw std::runtime_error(
                    token + " is not a valid vector representation");
            // remove the end brackets
            std::string removeBrackets = token.substr(1, token.size() - 2);
            std::stringstream st(removeBrackets);

            std::vector<length_t> vector;
            std::string t;
            while (getline(st, t, ',')) {
                vector.emplace_back(stoull(t));
            }
            vectors.emplace_back(vector);
        }

        return Search::makeSearch(vectors[0], vectors[1], vectors[2]);
    }

    void doSearch(std::vector<FMOcc>& occ, const Search& s,
                  const std::vector<RangePair> exactMatchRanges,
                  std::vector<Substring>& parts) const {

        // get the first part of the search (already matched)
        length_t first = s.getPart(0);
        RangePair ranges = exactMatchRanges[first];

        if (ranges.width()) {
            // exact match exists
            length_t exactLength = parts[first].size();

            // prepare the parts of the pattern for this search
            s.setDirectionsInParts(parts);

            // Do exact matching on the next parts of the search, as long as the
            // upper bound stays 0

            length_t idxInSearch = 1;

            while (s.getUpperBound(idxInSearch) == 0) {
                // extend the exact match
                // A) get the part to match exactly
                const auto& part = parts[s.getPart(idxInSearch)];
                // B) set direction in index
                index.setDirection(part.getDirection());
                // C) try to match part exactly
                ranges = index.matchExactBidirectionally(part, ranges);
                if (ranges.empty()) {
                    // search failed
                    return;
                }
                // D) update exact length and idx in search
                exactLength += part.size();
                idxInSearch++;
            }

            // Create a start occurrence corresponding to the exact match
            BiFMOcc startOcc = BiFMOcc(ranges, 0, exactLength);
            // Start the approximate matching phase
            index.recApproxMatch(s, startOcc, occ, parts, idxInSearch);
        }
    }

  public:
    SearchScheme(BiFMIndex& index, const std::string& folder,
                 const length_t maxED)
        : index(index), maxED(maxED) {
        // read the search scheme
        // get the name of the file
        std::string line;
        {
            std::ifstream ifs(folder + "name.txt");
            if (!ifs) {
                throw std::runtime_error(
                    "Problem reading: " + folder +
                    "name.txt\nDid you provide a directory to "
                    "a search scheme without a name file?");
            }
            getline(ifs, line);
            name = line;
            ifs.close();
        }

        std::string searchFile =
            folder + std::to_string(maxED) + "/searches.txt";

        std::ifstream stream_searches(searchFile);
        if (!stream_searches) {
            throw std::runtime_error("Prolbem reading: " + searchFile);
        }

        // read the searches line by line
        while (getline(stream_searches, line)) {
            try {
                searches.push_back(makeSearchFromLine(line));
            } catch (const std::runtime_error& e) {
                throw std::runtime_error(
                    "Something went wrong with processing line: " + line +
                    "\nin file: " + searchFile + e.what());
            }
        }

        stream_searches.close();
    }

    std::vector<TextOcc> matchApprox(const std::string& p) const {
        if (maxED == 0) {
            const auto& pos = index.matchExact(p);
            std::vector<TextOcc> r;
            r.reserve(pos.size());
            for (const auto& po : pos) {
                r.emplace_back(Range(po, po + p.size()), 0);
            }
            return r;
        }

        // create the parts of the pattern
        std::vector<Substring> parts;
        unsigned int numParts = searches[0].getNumParts();

        if (numParts * maxED >= p.size()) {
            // splitting up was not viable -> just search the entire pattern
            std::cerr << "Warning: Naive approx matching was used as "
                         "entered pattern is too short "
                      << p.size() << std::endl;

            return index.naiveApproxMatch(p, maxED);
        }

        // partition the read uniformly
        float fraction = ((float)p.size()) / numParts;
        for (unsigned int i = 0; i < numParts; i++) {
            parts.emplace_back(p, i * fraction, (i + 1) * fraction);
        }
        // set the end of the final part correct to be end of p
        parts.back().setEnd(p.size());

        // calculate the ranges corresponding to the exact match for each part
        std::vector<RangePair> exactMatchRanges;

        // the direction of the initial matching can be forward or backward, but
        // the direction in the parts is by default forward so do this in the
        // forward direction
        index.setDirection(FORWARD);

        for (const auto& part : parts) {
            exactMatchRanges.emplace_back(
                index.matchExactBidirectionally(part));
        }

        std::vector<FMOcc> occ; // the vector with all FM occurrences
        // do each search
        for (const auto& s : searches) {
            doSearch(occ, s, exactMatchRanges, parts);
        }

        return index.filterRedundantMatches(occ, maxED);
    }
};

#endif