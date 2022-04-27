#include "fmindex.h"
#include "gtest/gtest.h"

using namespace std;

class Week1Test : public ::testing::Test {
  protected:
    static string base;
    static FMIndex fmindex;
    static string text;
};

string Week1Test::base = "../../testset/CP001363";
FMIndex Week1Test::fmindex = FMIndex(base, 32, false);
string Week1Test::text = fmindex.getText();

class ConstructionTest : public Week1Test {};
class FunctionalityTest : public Week1Test {};
class IntegrationTest : public Week1Test {};

TEST_F(ConstructionTest, BWTConstructionTEST) {

    string correctBWT;

    readText(base + ".bwt", correctBWT);

    const string& bwt = fmindex.getBWT();

    ASSERT_EQ(bwt.size(), correctBWT.size());

    ASSERT_EQ(bwt, correctBWT);
}

TEST_F(ConstructionTest, CountsConstructionTEST) {

    length_t dollarCounts = 0, aCounts = 1, cCounts = 1164824,
             gCounts = 2435826, tCounts = 3706926;

    const auto& counts = fmindex.getCounts();

    EXPECT_EQ(counts[0], dollarCounts);
    EXPECT_EQ(counts[1], aCounts);
    EXPECT_EQ(counts[2], cCounts);
    EXPECT_EQ(counts[3], gCounts);
    EXPECT_EQ(counts[4], tCounts);
}

TEST(BitvecTest, SmallSize) {
    // small size: only third level counts (popcounts) are used
    size_t bvSize = 29;
    Bitvec bv(bvSize);

    EXPECT_EQ(bv.size(), bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv[i], i % 3 == 0);

    bv.index();

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv.rank(i), (i + 2) / 3);
}

TEST(BitvecTest, MediumSize) {
    // medium size: second and third level counts are used
    size_t bvSize = 389;
    Bitvec bv(bvSize);

    EXPECT_EQ(bv.size(), bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv[i], i % 3 == 0);

    bv.index();

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv.rank(i), (i + 2) / 3);
}

TEST(BitvecTest, LargeSize) {
    // large size: level 1, 2 and 3 counts are used
    size_t bvSize = 71234;
    Bitvec bv(bvSize);

    EXPECT_EQ(bv.size(), bvSize);

    for (size_t i = 0; i < bvSize; i += 3)
        bv[i] = true;

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv[i], i % 3 == 0);

    bv.index();

    for (size_t i = 0; i < bvSize; i++)
        EXPECT_EQ(bv.rank(i), (i + 2) / 3);
}

TEST_F(FunctionalityTest, occTest) {

    length_t dollarPos = 626743;

    EXPECT_EQ(fmindex.occ(0, dollarPos - 1), 0);
    EXPECT_EQ(fmindex.occ(0, dollarPos), 0);
    EXPECT_EQ(fmindex.occ(0, dollarPos + 1), 1);
    EXPECT_EQ(fmindex.occ(0, fmindex.getText().size()), 1);

    EXPECT_EQ(fmindex.occ(2, dollarPos - 1), 141847);
    EXPECT_EQ(fmindex.occ(2, text.size() / 4), 324622);
    EXPECT_EQ(fmindex.occ(2, text.size() / 2), 614863);
    EXPECT_EQ(fmindex.occ(2, fmindex.getText().size()), 1271002);

    EXPECT_EQ(fmindex.occ(1, dollarPos - 1), 200906);
    EXPECT_EQ(fmindex.occ(1, text.size() / 4), 349204);
    EXPECT_EQ(fmindex.occ(1, text.size() / 2), 594505);
    EXPECT_EQ(fmindex.occ(1, fmindex.getText().size()), 1164823);

    EXPECT_EQ(fmindex.occ(3, dollarPos - 1), 150946);
    EXPECT_EQ(fmindex.occ(3, text.size() / 4), 296797);
    EXPECT_EQ(fmindex.occ(3, text.size() / 2), 714164);
    EXPECT_EQ(fmindex.occ(3, fmindex.getText().size()), 1271100);

    EXPECT_EQ(fmindex.occ(4, dollarPos - 1), 133043);
    EXPECT_EQ(fmindex.occ(4, text.size() / 4), 246942);
    EXPECT_EQ(fmindex.occ(4, text.size() / 2), 511600);
    EXPECT_EQ(fmindex.occ(4, fmindex.getText().size()), 1163340);
}

TEST_F(FunctionalityTest, findLFTest) {

    vector<length_t> values = {1,       2435828, 3706935, 2435847, 1164846,
                               1164853, 3706962, 3706969, 2435873, 2435881,
                               22,      1164888, 30,      3706999, 38,
                               41,      48,      1164920, 3707032, 3707039};

    for (length_t i = 0; i < 400; i += 20) {
        auto v = values[i / 20];
        EXPECT_EQ(fmindex.findLF(i), v);
    }
}

TEST_F(FunctionalityTest, findSATest) {

    vector<length_t> values = {4870265, 3487107, 305082,  25114,   3611763,
                               1183065, 3214899, 932651,  1152372, 3183828,
                               1971076, 902938,  4854846, 4047837, 3366531,
                               3027742, 1992765, 5923,    4536251, 782818};

    for (length_t i = 0; i < 400; i += 20) {
        auto v = values[i / 20];
        EXPECT_EQ(fmindex.findSA(i), v);
    }
}

TEST_F(FunctionalityTest, AddCharLeftTest) {

    Range startRange(0, text.size());

    EXPECT_EQ(fmindex.addCharLeft(1, startRange, startRange), true);
    EXPECT_EQ(startRange, Range(1, 1164824));

    EXPECT_EQ(fmindex.addCharLeft(3, startRange, startRange), true);
    EXPECT_EQ(startRange, Range(2435826, 2714511));

    EXPECT_EQ(fmindex.addCharLeft(1, startRange, startRange), true);
    EXPECT_EQ(startRange, Range(594619, 650984));

    EXPECT_EQ(fmindex.addCharLeft(3, startRange, startRange), true);
    EXPECT_EQ(startRange, Range(2580221, 2592119));

    EXPECT_EQ(fmindex.addCharLeft(2, startRange, startRange), true);
    EXPECT_EQ(startRange, Range(1819937, 1822541));

    EXPECT_EQ(fmindex.addCharLeft(1, Range(3355845, 3355846), startRange),
              false);
    EXPECT_EQ(startRange, Range(0, 0));
}

TEST_F(IntegrationTest, matchExactTest) {

    length_t size = 5;
    vector<string> substrings = {
        "AAAAA", "AAAAG", "AAGAT", "ACAGC", "ACGAG", "ACTAC", "AGAGA", "AGAGT",
        "AGATA", "AGCAG", "AGGAG", "AGGAT", "ATAGA", "ATCAT", "ATGGT", "ATTTA",
        "CACAG", "CACCA", "CAGAC", "CAGAG", "CAGCA", "CAGGT", "CCACA", "CCATC",
        "CCCCA", "CCGCG", "CGACA", "CGCAG", "CGCGC", "GAACG", "GAAGT", "GACGA",
        "GAGAC", "GAGAG", "GAGGA", "GAGTA", "GATAC", "GATAG", "GATAT", "GCATA",
        "GGATT", "GGCTT", "GGGAA", "GGGGC", "GGGGG", "GTAGA", "TAAGA", "TACGC",
        "TACTA", "TAGAT", "TATAG", "TCTAG", "TGAGT", "TTATA", "TTTAA", "TTTAG",
        "TTTTT"};

    vector<length_t> expectedNumbers = {
        13047, 7169, 4575, 6136, 1869, 2033, 3128, 2091, 4284,  7393,
        2412,  3828, 2614, 6239, 5019, 5327, 3244, 6955, 5631,  2886,
        11146, 6253, 3553, 6687, 3226, 8540, 5336, 9202, 14779, 5099,
        3616,  5168, 2070, 2282, 2597, 2253, 4500, 3984, 6511,  4562,
        4788,  5457, 4417, 3293, 1952, 2291, 2149, 6654, 1316,  1534,
        2195,  134,  2745, 3459, 6273, 2826, 12952};

    for (length_t i = 0; i < substrings.size(); i++) {
        const auto& sub = substrings[i];
        const auto& pos = fmindex.matchExact(sub);
        EXPECT_EQ(pos.size(), expectedNumbers[i]);
        for (const auto& p : pos) {
            EXPECT_EQ(text.substr(p, size), sub);
        }
    }

    // non-present string
    auto r = fmindex.matchExact("ACCGGATCGTGTGAAGAGGGGAACGTTC");
    EXPECT_EQ(r.size(), 0);
}

vector<pair<string, string>> getPairedReads(const string& base) {
    ifstream ifs(base + ".reads.fasta");
    if (!ifs)
        throw runtime_error("Problem reading: " + base + ".reads.fasta");

    string line;
    bool first = true;
    vector<pair<string, string>> pairedReads;

    while (getline(ifs, line)) {
        if (line[0] == '>')
            continue;
        if (first)
            pairedReads.emplace_back(line, "");
        else
            pairedReads.back().second = line;
        first = !first;
    }

    return pairedReads;
}

TEST_F(IntegrationTest, bestPairedTest) {
    vector<tuple<length_t, length_t, bool>> expected = {
        make_tuple(2884786, 2885379, 1), make_tuple(1020415, 1021187, 1),
        make_tuple(4519136, 4519852, 1), make_tuple(1735711, 1736235, 1),
        make_tuple(4681736, 4681107, 0), make_tuple(4116269, 4115640, 0),
        make_tuple(2717073, 2717728, 1), make_tuple(4192027, 4191289, 0),
        make_tuple(4322813, 4323527, 1), make_tuple(1211085, 1210492, 0),
        make_tuple(4759529, 4758937, 0), make_tuple(779546, 778912, 0),
        make_tuple(1815927, 1816818, 1), make_tuple(307280, 307923, 1),
        make_tuple(4800709, 4801480, 1), make_tuple(82931, 83604, 1),
        make_tuple(41436, 42201, 1),     make_tuple(1112359, 1113178, 1),
        make_tuple(1363082, 1362630, 0), make_tuple(1730220, 1729457, 0)};

    const auto pairedReads = getPairedReads(base);

    for (length_t i = 0; i < expected.size(); i++) {

        EXPECT_EQ(fmindex.bestPairedMatch(pairedReads[i], 800), expected[i]);
    }
}