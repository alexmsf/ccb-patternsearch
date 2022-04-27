#include "bandmatrix.h"
#include "fmindex.h"
#include "gtest/gtest.h"

using namespace std;

class Week2Test : public ::testing::Test {
  protected:
    static string base;
    static FMIndex fmindex;
    static string text;
};

string Week2Test::base = "../../testset/CP001363";
FMIndex Week2Test::fmindex = FMIndex(base, 32, false);
string Week2Test::text = fmindex.getText();

class FunctionalityTest : public Week2Test {};
class IntegrationTest : public Week2Test {};

TEST(BandedMatrixTest, UpdateMatrixCellTest) {

    BandedMatrix m(20, 4, 0);

    vector<bool> bools = {
        false, false, false, false, false, false, true,  false, false, false,
        false, false, false, false, true,  false, false, false, false, false,
        false, true,  false, false, false, false, false, false, true,  false,
        true,  false, true,  false, false, false, false, false, false, false,
        true,  true,  true,  false, true,  false, true,  false, true,  true,
        true,  false, true,  false, true,  false, false, false, false, false,
        true,  true,  false, true,  false, true,  true,  true,  true,  true,
        false, true,  false, false, false, false, false, true,  false, true,
        false, true,  false, false, false, false, false, false, false, false,
        false, true,  false, false, false, false, true,  false, false, false,
        false, false, true,  false, true,  false, false, true,  false, false,
        false, true,  true,  false, false, false, false, true,  true,  false,
        false, true,  false, true,  true,  false, true,  false, true,  false,
        false, true,  false, false, true,  false, false, false, false, false,
        false, false, true,  false, true,  false, false, true,  true,  true,
        false, false, true,  false, false, true,  false, true,  false, false,
        false, false, true,  true,  false, false, false, true,  false, false};

    for (length_t r = 1, counter = 0; r < m.getNumberOfRows(); r++)
        for (int c = m.getFirstColumn(r); c <= m.getLastColumn(r); c++)
            m.updateMatrixCell(bools[counter++], r, c);

    vector<pair<int, int>> indexes = {
        {2, 3},   {2, 5},   {2, 6},   {3, 1},   {3, 7},   {4, 1},   {4, 4},
        {4, 5},   {4, 8},   {6, 8},   {6, 10},  {7, 4},   {7, 7},   {8, 4},
        {8, 8},   {9, 5},   {9, 6},   {9, 7},   {9, 10},  {9, 13},  {10, 9},
        {10, 14}, {11, 7},  {11, 10}, {11, 11}, {11, 12}, {11, 13}, {11, 15},
        {12, 9},  {13, 9},  {13, 13}, {13, 14}, {14, 12}, {14, 15}, {14, 17},
        {15, 14}, {15, 15}, {16, 18}, {16, 19}, {16, 20}, {17, 14}, {17, 20},
        {18, 16}, {18, 17}, {18, 19}, {18, 20}, {19, 15}, {19, 19}, {20, 16},
        {22, 18}, {22, 19}, {23, 20}};
    vector<int> values = {1, 3, 4, 2, 4, 3, 2, 2, 4, 4, 4, 3, 2, 4, 2, 4, 3, 2,
                          4, 6, 2, 6, 4, 2, 3, 4, 5, 6, 3, 4, 3, 4, 3, 4, 6, 2,
                          3, 5, 6, 6, 4, 6, 4, 3, 4, 5, 5, 3, 5, 5, 4, 4};

    for (unsigned int i = 0; i < values.size(); i++) {
        int r = indexes[i].first, c = indexes[i].second,
            correctValue = values[i];
        EXPECT_EQ(m(r, c), correctValue);
    }
}

TEST(BandedMatrixTest, UpdateMatrixRowTest) {

    length_t k = 4;
    BandedMatrix m(20, k, 0);
    string s = "ACGTACGTAAGGCAGAT";

    Substring pattern(s, BACKWARD);

    EXPECT_EQ(m.updateMatrixRow(pattern, 1, 'A'), 1);
    EXPECT_EQ(m.updateMatrixRow(pattern, 1, 'C'), 1);
    EXPECT_EQ(m.updateMatrixRow(pattern, 1, 'G'), 1);
    EXPECT_EQ(m.updateMatrixRow(pattern, 1, 'T'), 0);

    for (length_t c = 1; c < 5; c++) {
        EXPECT_EQ(m(1, c), c - 1);
    }

    EXPECT_EQ(m.updateMatrixRow(pattern, 2, 'A'), 0);
    EXPECT_EQ(m.updateMatrixRow(pattern, 3, 'G'), 0);

    m.updateMatrixRow(pattern, 4, 'C');

    vector<length_t> values = {3, 2, 1, 1, 1, 2, 3, 4};

    for (int col = 1; col <= m.getLastColumn(4); col++) {
        EXPECT_EQ(values[col - 1], m(4, col));
    }

    EXPECT_EQ(m.updateMatrixRow(pattern, 25, 'A') > k, true);
}

TEST_F(FunctionalityTest, ExtendTest) {

    Range s(0, text.size());

    vector<FMPosExt> stack;

    fmindex.extendFMPos(s, 0, stack);
    EXPECT_EQ(stack.size(), 4);
    fmindex.extendFMPos(Range(1819937, 1822541), 5, stack);
    EXPECT_EQ(stack.size(), 8);
    fmindex.extendFMPos(Range(2523519, 2523522), 9, stack);
    EXPECT_EQ(stack.size(), 10);

    vector<Range> correctRanges1 = {Range(1, 1164824), Range(1164824, 2435826),
                                    Range(2435826, 3706926),
                                    Range(3706926, 4870266)};
    vector<Range> correctRanges2 = {
        Range(474733, 475238), Range(1629564, 1630139), Range(2930559, 2931562),
        Range(4092657, 4093178)};
    vector<Range> correctRanges3 = {Range(612854, 612855),
                                    Range(1801770, 1801772)};

    for (int i = 0; i < 4; i++) {
        const auto& p = stack[i];
        EXPECT_EQ(p.getRange(), correctRanges1[i]);
        EXPECT_EQ(p.getCharacter(), fmindex.getAlphabet().i2c(i + 1));
        EXPECT_EQ(p.getDepth(), 1);
    }
    for (int i = 4; i < 8; i++) {
        const auto& p = stack[i];
        EXPECT_EQ(p.getRange(), correctRanges2[i - 4]);
        EXPECT_EQ(p.getCharacter(), fmindex.getAlphabet().i2c((i - 4) + 1));
        EXPECT_EQ(p.getDepth(), 6);
    }
    for (int i = 8; i < 10; i++) {
        const auto& p = stack[i];
        EXPECT_EQ(p.getRange(), correctRanges3[i - 8]);
        EXPECT_EQ(p.getCharacter(), fmindex.getAlphabet().i2c((i - 8) + 1));
        EXPECT_EQ(p.getDepth(), 10);
    }
}

TEST_F(FunctionalityTest, ConvertTest) {
    FMOcc test = FMOcc(FMPos(Range(2273168, 2273173), 49), 3);
    vector<TextOcc> expected = {TextOcc(Range(294646, 294695), 3),
                                TextOcc(Range(4369381, 4369430), 3),
                                TextOcc(Range(4118332, 4118381), 3),
                                TextOcc(Range(4412763, 4412812), 3),
                                TextOcc(Range(4214416, 4214465), 3)};

    vector<TextOcc> result;
    fmindex.convertFMOccToTextOcc(test, result);

    EXPECT_EQ(expected, result);
}

TEST_F(IntegrationTest, NaiveMatchesTest) {
    string test = "GCGATTATCTCTGTCGGCGACGGTAT";

    auto returned = fmindex.naiveApproxMatch(test, 4);

    vector<TextOcc> expected = {TextOcc(Range(1524, 1553), 3)};

    EXPECT_EQ(expected, returned);

    test = "GCTGAAGTAGGTCCCAAGGGTATGGCTGTTGCCATTAAAGTGGTACGC";
    returned = fmindex.naiveApproxMatch(test, 4);
    expected = {TextOcc(Range(294645, 294695), 2),
                TextOcc(Range(4118331, 4118381), 2),
                TextOcc(Range(4214415, 4214465), 2),
                TextOcc(Range(4369380, 4369430), 2),
                TextOcc(Range(4412762, 4412812), 2)};

    EXPECT_EQ(expected, returned);
}
