
#include "bidirectionalfmindex.h"
#include "searchscheme.h"
#include "gtest/gtest.h"

using namespace std;

TEST(CumulativeBitvec, OccTest) {
    // alphabet size of string (including '$') below == 20
    string BWT("Hello,Iamastringwith$alargeralphabetsize");
    Alphabet<32> sigma(BWT);
    CumulativeBitvectors<32> test(sigma, BWT);

    // test the occ(c,k) routine
    vector<size_t> expVal(sigma.size(), 0);
    for (size_t i = 0; i < BWT.size(); i++) {
        for (size_t cIdx = 0; cIdx < sigma.size(); cIdx++)
            EXPECT_EQ(test.occ(cIdx, i), expVal[cIdx]);
        expVal[sigma.c2i(BWT[i])]++;
    }
}

TEST(CumulativeBitvec, CumulOccTest) {
    // alphabet size of string (including '$') below == 20
    string BWT("Hello,Iamastringwith$alargeralphabetsize");
    Alphabet<32> sigma(BWT);
    CumulativeBitvectors<32> test(sigma, BWT);

    // test the cumOcc(c,k) routine
    vector<size_t> expVal = vector<size_t>(sigma.size(), 0);
    for (size_t i = 0; i < BWT.size(); i++) {
        for (size_t cIdx = 0; cIdx < sigma.size(); cIdx++)
            EXPECT_EQ(test.cumulocc(cIdx, i), expVal[cIdx]);

        for (int c = sigma.c2i(BWT[i]) + 1; c < (int)sigma.size(); c++)
            expVal[c]++;
    }
}

class Week3Test : public ::testing::Test {
  protected:
    static string base;
    static BiFMIndex bifmindex;
    static string text;
};

string Week3Test::base = "../../testset/CP001363";
BiFMIndex Week3Test::bifmindex = BiFMIndex(base, 32, false);
string Week3Test::text = bifmindex.getText();

class ConstructionTest : public Week3Test {};
class FunctionalityTest : public Week3Test {};
class IntegrationTest : public Week3Test {
  protected:
    static int maxED;
    static SearchScheme ss;
};

int IntegrationTest::maxED = 4;
SearchScheme IntegrationTest::ss =
    SearchScheme(bifmindex, "../../search_schemes/kuch_k+1/", maxED);

TEST_F(ConstructionTest, CreateRevBWTTest) {

    string substr =
        "AGGCAACGCGAAGCACTCCTACAATCTCGAAGGTCGTCCGTCGGCGAAAAGAATCGCCAGAACTTCGACG"
        "CGCTTTAATTTCCCAAACTCAGCATAGCTACCGGAAAAGATCGAACGGACAATAGAGCCTACTGTTCTTA"
        "TTATAAGCAGGAGAATGACGGACCATGGGATGAATGTAAGCAGTTTGTAAATATCCCCAGATCGCTTACA"
        "GGTTAGTTTTAGGAAAGGTCCACACGAATTAGAAAATCGACATAGTATTAATAGGTTTTAGAAATGAGTC"
        "CGAAAAGCAACCATTTAATGAGCGGTCTGTGTTGTTCCTACAAAGATTCCGCATTCAAGGGAGAGAAGTG"
        "GATACGGTGTAATTAACAGAGTCCCAAATTCAAGGGAAAAGTGTTTTAAACGACATCTTTACATAAAAAG"
        "CCTAGTCATGCAATATAAGCTACAAAACAAAAAGTCAATATAAATGATCATGAATAGATAAAAATTTATA"
        "AAAGTGTTTTGTACAAAAAAGCCTCCAATCAGAAATTACACTTAATATTTT";

    string revBWT;

    vector<length_t> SA;
    readSA(base + ".rev.sa", SA, text.size());

    bifmindex.createRevBWTFromRevSA(SA, revBWT);

    EXPECT_EQ(revBWT.substr(10000, 541), substr);
    EXPECT_EQ(revBWT[0], bifmindex.getText()[0]);
}

TEST_F(FunctionalityTest, AddCharLeftTest) {
    auto textLength = text.size();
    RangePair startRanges(0, textLength, 0, textLength);

    EXPECT_EQ(bifmindex.addCharLeft(1, startRanges, startRanges), true);
    EXPECT_EQ(startRanges, RangePair(Range(1, 1164824), Range(1, 1164824)));

    EXPECT_EQ(bifmindex.addCharLeft(3, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(2435826, 2714511), Range(652880, 931565)));

    EXPECT_EQ(bifmindex.addCharLeft(1, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(594619, 650984), Range(652880, 709245)));

    EXPECT_EQ(bifmindex.addCharLeft(3, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(2580221, 2592119), Range(691062, 702960)));

    EXPECT_EQ(bifmindex.addCharLeft(2, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(1819937, 1822541), Range(694190, 696794)));

    // Try extension that should fail
    EXPECT_EQ(bifmindex.addCharLeft(
                  1,
                  RangePair(Range(3355844, 3355848), Range(4692519, 4692523)),
                  startRanges),
              false);

    // check if range has been made empty after extension failed
    EXPECT_EQ(startRanges, RangePair(0, 0, 0, 0));
}

TEST_F(FunctionalityTest, AddCharRightTest) {
    auto textLength = text.size();
    RangePair startRanges(0, textLength, 0, textLength);

    EXPECT_EQ(bifmindex.addCharRight(1, startRanges, startRanges), true);

    EXPECT_EQ(startRanges, RangePair(Range(1, 1164824), Range(1, 1164824)));

    EXPECT_EQ(bifmindex.addCharRight(3, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(594619, 841220), Range(2435826, 2682427)));

    EXPECT_EQ(bifmindex.addCharRight(1, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(594619, 650984), Range(652880, 709245)));

    EXPECT_EQ(bifmindex.addCharRight(3, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(623814, 634621), Range(2605731, 2616538)));

    EXPECT_EQ(bifmindex.addCharRight(2, startRanges, startRanges), true);
    EXPECT_EQ(startRanges,
              RangePair(Range(626942, 630039), Range(1782121, 1785218)));

    // Try extension that should fail
    EXPECT_EQ(bifmindex.addCharLeft(
                  1,
                  RangePair(Range(3355844, 3355848), Range(4692519, 4692523)),
                  startRanges),
              false);

    // check if range has been made empty after extension failed
    EXPECT_EQ(startRanges, RangePair(0, 0, 0, 0));
}

TEST_F(FunctionalityTest, ExtendTest) {

    RangePair s(0, text.size(), 0, text.size());

    vector<BiFMPosExt> stack;

    bifmindex.setDirection(FORWARD);
    bifmindex.extendFMPos(s, 0, stack);
    EXPECT_EQ(stack.size(), 4);

    vector<RangePair> correctRanges1 = {
        RangePair(Range(1, 1164824), Range(1, 1164824)),
        RangePair(Range(1164824, 2435826), Range(1164824, 2435826)),
        RangePair(Range(2435826, 3706926), Range(2435826, 3706926)),
        RangePair(Range(3706926, 4870266), Range(3706926, 4870266))};

    for (int i = 0; i < 4; i++) {
        const auto& p = stack[i];
        EXPECT_EQ(p.getRanges(), correctRanges1[i]);
        EXPECT_EQ(p.getCharacter(), bifmindex.getAlphabet().i2c(i + 1));
        EXPECT_EQ(p.getDepth(), 1);
    }

    bifmindex.setDirection(BACKWARD);
    bifmindex.extendFMPos(
        RangePair(Range(1819937, 1822541), Range(694190, 696794)), 5, stack);
    EXPECT_EQ(stack.size(), 8);

    vector<RangePair> correctRanges2 = {
        RangePair(Range(474733, 475238), Range(694190, 694695)),
        RangePair(Range(1629564, 1630139), Range(694695, 695270)),
        RangePair(Range(2930559, 2931562), Range(695270, 696273)),
        RangePair(Range(4092657, 4093178), Range(696273, 696794))};

    for (int i = 4; i < 8; i++) {
        const auto& p = stack[i];
        EXPECT_EQ(p.getRanges(), correctRanges2[i - 4]);
        EXPECT_EQ(p.getCharacter(), bifmindex.getAlphabet().i2c((i - 4) + 1));
        EXPECT_EQ(p.getDepth(), 6);
    }

    bifmindex.setDirection(FORWARD);
    bifmindex.extendFMPos(
        RangePair(Range(2523519, 2523522), Range(694387, 694390)), 9, stack);
    EXPECT_EQ(stack.size(), 10);

    vector<RangePair> correctRanges3 = {
        RangePair(Range(2523519, 2523521), Range(196070, 196072)),
        RangePair(Range(2523521, 2523522), Range(1311256, 1311257))};

    for (int i = 8; i < 10; i++) {
        const auto& p = stack[i];
        EXPECT_EQ(p.getRanges(), correctRanges3[i - 8]);
        EXPECT_EQ(p.getCharacter(), bifmindex.getAlphabet().i2c((i - 8) + 1));
        EXPECT_EQ(p.getDepth(), 10);
    }
}

TEST_F(FunctionalityTest, matchBidirectionallyTest) {

    vector<string> strings = {
        "AAAAA", "AAAAG", "AAGAT", "ACAGC", "ACGAG", "ACTAC", "AGAGA", "AGAGT",
        "AGATA", "AGCAG", "AGGAG", "AGGAT", "ATAGA", "ATCAT", "ATGGT", "ATTTA",
        "CACAG", "CACCA", "CAGAC", "CAGAG", "CAGCA", "CAGGT", "CCACA", "CCATC",
        "CCCCA", "CCGCG", "CGACA", "CGCAG", "CGCGC", "GAACG", "GAAGT", "GACGA",
        "GAGAC", "GAGAG", "GAGGA", "GAGTA", "GATAC", "GATAG", "GATAT", "GCATA",
        "GGATT", "GGCTT", "GGGAA", "GGGGC", "GGGGG", "GTAGA", "TAAGA", "TACGC",
        "TACTA", "TAGAT", "TATAG", "TCTAG", "TGAGT", "TTATA", "TTTAA", "TTTAG",
        "TTTTT"};

    vector<Substring> substrings;
    for (const auto& s : strings)
        substrings.emplace_back(s);

    vector<RangePair> expectedOutput = {
        RangePair(Range(2, 13049), Range(2, 13049)),
        RangePair(Range(21693, 28862), Range(2435827, 2442996)),
        RangePair(Range(203731, 208306), Range(3869148, 3873723)),
        RangePair(Range(366020, 372156), Range(1742694, 1748830)),
        RangePair(Range(474733, 476602), Range(2616538, 2618407)),
        RangePair(Range(551028, 553061), Range(1379951, 1381984)),
        RangePair(Range(623814, 626942), Range(691062, 694190)),
        RangePair(Range(632530, 634621), Range(4303768, 4305859)),
        RangePair(Range(634621, 638905), Range(970116, 974400)),
        RangePair(Range(660118, 667511), Range(2551990, 2559383)),
        RangePair(Range(745344, 747756), Range(2627129, 2629541)),
        RangePair(Range(747756, 751584), Range(3912481, 3916309)),
        RangePair(Range(881295, 883909), Range(702960, 705574)),
        RangePair(Range(933817, 940056), Range(3847339, 3853578)),
        RangePair(Range(1062039, 1067058), Range(4452982, 4458001)),
        RangePair(Range(1138736, 1144063), Range(1140263, 1145590)),
        RangePair(Range(1240566, 1243810), Range(2505281, 2508525)),
        RangePair(Range(1246815, 1253770), Range(399322, 406277)),
        RangePair(Range(1300307, 1305938), Range(1305055, 1310686)),
        RangePair(Range(1305938, 1308824), Range(2609700, 2612586)),
        RangePair(Range(1315453, 1326599), Range(485175, 496321)),
        RangePair(Range(1376007, 1382260), Range(4404872, 4411125)),
        RangePair(Range(1492625, 1496178), Range(359563, 363116)),
        RangePair(Range(1549938, 1556625), Range(2182191, 2188878)),
        RangePair(Range(1581006, 1584232), Range(418733, 421959)),
        RangePair(Range(1657296, 1665836), Range(2910276, 2918816)),
        RangePair(Range(1799648, 1804984), Range(371604, 376940)),
        RangePair(Range(1868428, 1877630), Range(2559383, 2568585)),
        RangePair(Range(1941214, 1955993), Range(1887673, 1902452)),
        RangePair(Range(2471908, 2477007), Range(2696718, 2701817)),
        RangePair(Range(2496530, 2500146), Range(4279729, 4283345)),
        RangePair(Range(2547066, 2552234), Range(717068, 722236)),
        RangePair(Range(2584142, 2586212), Range(1310686, 1312756)),
        RangePair(Range(2586212, 2588494), Range(2612586, 2614868)),
        RangePair(Range(2606641, 2609238), Range(795941, 798538)),
        RangePair(Range(2616853, 2619106), Range(1042313, 1044566)),
        RangePair(Range(2634109, 2638609), Range(1373076, 1377576)),
        RangePair(Range(2638609, 2642593), Range(2657285, 2661269)),
        RangePair(Range(2642593, 2649104), Range(3969100, 3975611)),
        RangePair(Range(2782057, 2786619), Range(961873, 966435)),
        RangePair(Range(3208555, 3213343), Range(4578778, 4583566)),
        RangePair(Range(3318855, 3324312), Range(4651018, 4656475)),
        RangePair(Range(3324312, 3328729), Range(231747, 236164)),
        RangePair(Range(3358344, 3361637), Range(2027736, 2031029)),
        RangePair(Range(3361637, 3363589), Range(3291616, 3293568)),
        RangePair(Range(3482524, 3484815), Range(705699, 707990)),
        RangePair(Range(3748822, 3750971), Range(667294, 669443)),
        RangePair(Range(3811603, 3818257), Range(1833703, 1840357)),
        RangePair(Range(3825420, 3826736), Range(1009829, 1011145)),
        RangePair(Range(3841148, 3842682), Range(3883977, 3885511)),
        RangePair(Range(3873445, 3875640), Range(2661269, 2663464)),
        RangePair(Range(4166510, 4166644), Range(2664191, 2664325)),
        RangePair(Range(4272090, 4274835), Range(4309908, 4312653)),
        RangePair(Range(4578277, 4581736), Range(1001114, 1004573)),
        RangePair(Range(4756363, 4762636), Range(331934, 338207)),
        RangePair(Range(4769484, 4772310), Range(2679601, 2682427)),
        RangePair(Range(4857314, 4870266), Range(4857314, 4870266))};

    // match string in one direction
    for (length_t i = 0; i < substrings.size(); i++) {
        auto& sub = substrings[i];
        Direction dir = (i % 2) ? FORWARD : BACKWARD;
        bifmindex.setDirection(dir);
        sub.setDirection(dir);

        EXPECT_EQ(bifmindex.matchExactBidirectionally(sub), expectedOutput[i]);
    }
    cout << endl;

    // Extend the match in the other direction with a new string
    vector<string> strings2 = {
        "AAAAG", "AAGAT", "ACAGC", "ACGAG", "ACTAC", "AGAGA", "AGAGT", "AGATA",
        "AGCAG", "AGGAG", "AGGAT", "ATAGA", "ATCAT", "ATGGT", "ATTTA", "CACAG",
        "CACCA", "CAGAC", "CAGAG", "CAGCA", "CAGGT", "CCACA", "CCATC", "CCCCA",
        "CCGCG", "CGACA", "CGCAG", "CGCGC", "GAACG", "GAAGT", "GACGA", "GAGAC",
        "GAGAG", "GAGGA", "GAGTA", "GATAC", "GATAG", "GATAT", "GCATA", "GGATT",
        "GGCTT", "GGGAA", "GGGGC", "GGGGG", "GTAGA", "TAAGA", "TACGC", "TACTA",
        "TAGAT", "TATAG", "TCTAG", "TGAGT", "TTATA", "TTTAA", "TTTAG", "TTTTT",
        "AAAAA"};

    vector<Substring> substrings2;
    for (const auto& s : strings2)
        substrings2.emplace_back(s);

    vector<RangePair> expectedOutput2 = {
        RangePair(Range(10, 16), Range(2435827, 2435833)),
        RangePair(Range(203768, 203775), Range(2441839, 2441846)),
        RangePair(Range(204231, 204237), Range(1747760, 1747766)),
        RangePair(Range(474917, 474917), Range(1746329, 1746329)),
        RangePair(Range(474968, 474970), Range(1381066, 1381068)),
        RangePair(Range(624331, 624332), Range(1380205, 1380206)),
        RangePair(Range(624397, 624400), Range(4304177, 4304180)),
        RangePair(Range(635627, 635630), Range(4304403, 4304406)),
        RangePair(Range(635650, 635653), Range(2553186, 2553189)),
        RangePair(Range(745806, 745810), Range(2557488, 2557492)),
        RangePair(Range(745860, 745861), Range(3915362, 3915363)),
        RangePair(Range(881780, 881781), Range(3912912, 3912913)),
        RangePair(Range(881866, 881867), Range(3848311, 3848312)),
        RangePair(Range(1062718, 1062719), Range(3853075, 3853076)),
        RangePair(Range(1062934, 1062944), Range(1144947, 1144957)),
        RangePair(Range(1241198, 1241200), Range(1142733, 1142735)),
        RangePair(Range(1241331, 1241340), Range(402686, 402695)),
        RangePair(Range(1301499, 1301515), Range(400777, 400793)),
        RangePair(Range(1301603, 1301608), Range(2610376, 2610381)),
        RangePair(Range(1318173, 1318174), Range(2609942, 2609943)),
        RangePair(Range(1318221, 1318223), Range(4405399, 4405401)),
        RangePair(Range(1493682, 1493686), Range(4405313, 4405317)),
        RangePair(Range(1493783, 1493790), Range(2182612, 2182619)),
        RangePair(Range(1581732, 1581740), Range(2182727, 2182735)),
        RangePair(Range(1581812, 1581817), Range(2910935, 2910940)),
        RangePair(Range(1801654, 1801666), Range(2910847, 2910859)),
        RangePair(Range(1801796, 1801807), Range(2559874, 2559885)),
        RangePair(Range(1946808, 1946856), Range(2562663, 2562711)),
        RangePair(Range(1948989, 1949008), Range(2698521, 2698540)),
        RangePair(Range(2498128, 2498132), Range(2701203, 2701207)),
        RangePair(Range(2498217, 2498224), Range(721478, 721485)),
        RangePair(Range(2585133, 2585137), Range(718159, 718163)),
        RangePair(Range(2585153, 2585153), Range(2613266, 2613266)),
        RangePair(Range(2607925, 2607926), Range(2613049, 2613050)),
        RangePair(Range(2607956, 2607957), Range(1042760, 1042761)),
        RangePair(Range(2636581, 2636582), Range(1042953, 1042954)),
        RangePair(Range(2636595, 2636595), Range(2658470, 2658470)),
        RangePair(Range(2645862, 2645862), Range(2660495, 2660495)),
        RangePair(Range(2646015, 2646018), Range(965633, 965636)),
        RangePair(Range(3210965, 3210970), Range(966239, 966244)),
        RangePair(Range(3211347, 3211358), Range(4656164, 4656175)),
        RangePair(Range(3327410, 3327413), Range(4651205, 4651208)),
        RangePair(Range(3327440, 3327444), Range(2027841, 2027845)),
        RangePair(Range(3363152, 3363152), Range(2029757, 2029757)),
        RangePair(Range(3363170, 3363176), Range(707238, 707244)),
        RangePair(Range(3750290, 3750290), Range(705954, 705954)),
        RangePair(Range(3750404, 3750405), Range(1834747, 1834748)),
        RangePair(Range(3826302, 3826303), Range(1835362, 1835363)),
        RangePair(Range(3826309, 3826312), Range(3884505, 3884508)),
        RangePair(Range(3875229, 3875229), Range(3884704, 3884704)),
        RangePair(Range(3875367, 3875367), Range(2664271, 2664271)),
        RangePair(Range(2623953, 2623953), Range(2664308, 2664308)),
        RangePair(Range(4274534, 4274540), Range(1003957, 1003963)),
        RangePair(Range(4762232, 4762237), Range(1001300, 1001305)),
        RangePair(Range(4762492, 4762493), Range(2679711, 2679712)),
        RangePair(Range(4870117, 4870123), Range(2682421, 2682427)),
        RangePair(Range(4857314, 4857332), Range(13031, 13049))};

    // match extra string in other direction
    for (length_t i = 0; i < substrings2.size(); i++) {
        auto& sub = substrings2[i];
        Direction dir = ((i + 1) % 2) ? FORWARD : BACKWARD;
        bifmindex.setDirection(dir);
        sub.setDirection(dir);

        EXPECT_EQ(bifmindex.matchExactBidirectionally(sub, expectedOutput[i]),
                  expectedOutput2[i]);
    }
}

TEST_F(IntegrationTest, SearchSchemeApproxTest) {

    vector<string> tests = {
        "GCGATTATCTCTGTCGGCGACGGTAT",
        "ACAGAATATAAGTCGCAGACCCATTATACAAAAGGTACGCAGTCACACC",
        "ATAAAGAAAAAGCTTCTCTTCTGGCATGGAGAAAGTATCCGGGTACAGGTA",
        "CCGTGGTGGGGTTCCCGAGCGGCTAAAGGGAGCAGACTGTAAATCTGCCG",
        "GCTGAAGTAGGTCCCAAGGGTATGGCTGTTGCCATTAAAGTGGTACGC",
        "TGCCTCTCCGTCACCGCGTTCTTGCTGAGTACATAGGTGGAAGGTGAGTA",
        "GGCAACTTACTTCCACCATCGCCTCAAAACGCCGCAGCGCCTCCGGGCGG",
        "AGAAGGACCTGCGGATGCGAAGGCGCAGGCAAGAGCGTAATTATCAGAA",
        "CTTATCATTTTTATTTAAGTTTAAATATTTTGATAAATGGTTTTTATTTACT",
        "TCCCTAGCTGGTCTGAGCGGATGACCAGCCACACTGGAACTGAGACACGG",
        "CAGGTTCGAATCTTCATATTGCAGATGCAAAAAAGCGCCTTTAGGCGC",
        "CTTCTCCCTGCGCCTGACGCCTGACGTCCGCTTCCAGTTCAGCCAGCTCCG",
        "GGCCGCCGCCTCGGCAGGGTGAAGCTGATAAGGCCGAAGTAGTTCAGCG",
        "TTCACATCCGACTTGACAGACCGCCTGCGATGCGCTTTACGCCCAGTAAT"};

    for (const auto& test : tests) {
        EXPECT_EQ(bifmindex.naiveApproxMatch(test, maxED),
                  ss.matchApprox(test));
    }
}
