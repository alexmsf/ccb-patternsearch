# WEEK 1

## Introduction
In week 1, you will build a paired-end read mapper, using the FM index.
You are provided with a bacterial genome (CP001363, the search text) and a set of simulated paired-end reads (the search patterns).
In this assignment (week 1), the reads have no errors, and thus we can simply use exact pattern matching functionality.

You are provided with an existing implementation of the FM index in which you will need to implement certain key routines.

## Assignments

### Construction of the FM index

#### Creating the BWT
The first assignment is to create the Burrows-Wheeler transform (BWT) of the search text, given the suffix array (SA) and the search text itself.
Fill in the function:
```C++
void FMIndex::createBWTFromSA(const vector<length_t>& sa)
 ```
in the file `fmindex.cpp`
This function needs to update the `bwt` attribute of the `FMIndex` object. Note that the `bwt` attribute is an empty (i.e., not yet allocated) `std::string`.
You will need to allocate memory for the `bwt` attribute (e.g., using its `resize` member function) and set its characters to the correct values.
You can use the `text` attribute (again a `std::string`) of the `FMIndex` object which stores the search text (the genome) and the `sa` parameter.

---
**PROPERTY**

<img src="https://render.githubusercontent.com/render/math?math=BWT[i] = \begin{cases} \$ \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, \mathrm{if}\, SA[i] = 0 \\ T[SA[i] - 1] \,\,\, \mathrm{otherwise}\end{cases}">

---

Currently, this function throws a runtime error to warn users that the function has not yet been implemented. Remember to comment or delete this line, otherwise the program will exit and the test will fail.

---
**TESTING**

This function is tested in the ConstructionTest.BWTConstructionTEST of week 1. To test your solution rebuild your project (run `make` in the `build` folder) and run the unit tests:
```
cd build/unittest
./week1 
```

This will run all 11 tests for week 1, of which at least 10 will fail. To only run the test for this exercise run:

```
./week1 --gtest_filter=*ConstructionTest.BWTConstructionTEST 
```

You will see some warnings regarding not implemented functions. These will be resolved in later assignments. If your solution was correct this test should pass.

---

#### Creating the counts array

Next, we will create the counts array `counts`. The counts array counts for each character `c` in `sigma` (the alphabet) how many characters in the BWT (or text) are lexicographically smaller than `c`.

You need to complete the function:
```C++
void FMIndex::createCounts()
```
in the file `fmindex.cpp`.
You can assume that `counts[i] = 0` for `0 <= i < sigma.size()` prior to the call to this function. You can also assume that the `text` and `bwt` now contain the correct strings. To convert a character `c` to its character index in the alphabet you can use `sigma.c2i(c)`.

---
**NOTE**

The sentinel character ($) is the lexicographically smallest character corresponding to character index 0.

---

---
**TESTING**

Build your project and execute the following commands to test your solution to this assignment.

```
cd build/unittest
./week1 --gtest_filter=*ConstructionTest.CountsConstructionTEST 
```

---


### Rank operation on bit vectors

In this asignment, the occurrences table is stored as 4 bit vectors, one for each character ACGT. Each bit vector `bv` is associated with one character `c`, where `bv[i] = 1` if and only if `BWT[i] = c`. These bit vectors implement the rank9 algorithm to find the number of times character `c` occurs before index `p`. In other words, `bv.rank(p)` returns the number of 1-bits before index `p`.

You will need to implement the rank function on a bit vector.
Fill in function:
```C++
  size_t rank(size_t p) const
```
in the file `bitvec.h`.

---
**HINT**
These bit vectors make use of the rank9 algorithm. Remember that to find the rank for an index `p`, you need to add the first and second level counts of the word `w` in which `p` is located, as well as the number of 1-bits in word `w` before `p`
The `Bitvector` class already implements the functions:
```C++
  size_t firstLevelCounts(size_t w) const 
  size_t secondLevelCounts(size_t w) const
  size_t popcount(size_t w, size_t b) const 
```
where `w` is the wordindex of index `p` and `b` is the bit corresponding to `p`
Remember that a word is 64 bits long.

---


---
**TESTING**

Build your project and execute the following commands to test your solution to this assignment.

```
cd build/unittest
./week1 --gtest_filter=*BitvecTest*
```

---


### Sparse suffix array
To save memory, the suffix array can be made sparse. To achieve this, we only store a certain fraction of the elements of the suffix array. This is controlled by a sparseness factor `f`. Two options exist:

1. Store every `f`th index of the suffix array, i.e., only the values `SA[i]` for which `i % f = 0`.

2. Store every `f`th suffix of the suffix array, i.e., only the values `SA[i]` for which `SA[i] % f = 0`.

#### 1. Creating the sparse suffix array
You need to create the sparse suffix array. Modify the class `SparseSuffixArray` in `suffixarray.h`. You can choose which way you want to sparsify the suffix array (option 1 or 2). 
Option 1 is easier to implement, but has a `O(n)` worst case time complexity for the `findSA` function (see further).
Option 2 guarantees `O(1)` time-complexity for the `findSA` function, but requires the use of an additional bit vector.

First, complete the function:
```C++
void createSparseSA(const std::vector<length_t>& sa)
```
in the file `suffixArray.h` to create an empty sparse sa vector.
It is your responsibility to store the correct values in the `sparseSA` attribute.
Memory for the `sparseSA` attribute has not yet been allocated.
Again, you could for example use the `sparseSA.resize(...)` routine for this.
Keep in mind that the size of the suffix array is not necessarily some multiple of f.
The sparseness factor f is stored in the `sparsenessFactor` attribute.
If you have chosen option 2 you also need to set `bitvector[i]` to `true` if the element at index `i` in the original suffix array is stored.
If you have chosen option 1, the bit vector is not needed.

Next, fill in functions:
```C++
bool hasStored(const length_t i) const
length_t operator[](const length_t i) const
```

The `hasStored(i)` function should return `true` if the sparse suffix array stores the element that was stored at index `i` in the **original** suffix array.
If you have chosen option 1 then the element is stored if `i % f = 0`.
If you have chosen option 2 then you can query `bitvector[i]` to find out if this element is stored.

The second function is the array index operator `[i]`.
This function will only be called if `hasStored(i)` returns true.
If you have chosen option 1, you can use `i/f` to find the corresponding index in `sparseSA`.
If you have chosen option 2, you can use `bitvector.rank(i)` to find the corresponding index in `sparseSA`.

The correctness of this implementation will be tested later.

### Functionality of the FM index

Before starting this phase make sure that all tests using the following commands pass:
```
cd build/unittest
./week1 --gtest_filter=*ConstructionTest*
```

 
#### LF mapping

Next, we will be implementing the Last-First (LF) mapping.
It maps the character at index `i` in the L column (i.e., the BWT) to the corresponding character in the F column at index `j`.
Said otherwise, it maps index `i` to index `j` if and only if `SA[i] = SA[j] + 1`.
To implement the LF mapping, we need the BWT and the occurrence table, since `LF(i) = counts[sigma.c2i(bwt[i])] + occ(sigma.c2i(bwt[i]), i)`.
Note that `sigma.c2i(c)` converts character `c` ($, A, C, G, T) to its character index (0, 1, 2, 3, 4).
It is the inverse function of `sigma.i2c`.

The `occ` function gives, for each character `c` in the text (ACGT$) and each index `i`, the number of times `c` occurs in the half-open interval `bwt[0, i)`.
To implement this functionality, we make use of 4 bit vectors, one for each of the characters in ACGT.
The vector `occTable` stores these 4 bit vectors, i.e., `occTable[0]` is the bit vector for character 'A', `occTable[1]` is the bit vector for character 'C', etc.
In general, `occTable[charIdx - 1][i]` indicates whether `bwt[i] == sigma.i2c(charIdx)`.
For example for 'A' with charIdx 1, `occTable[0][i] = 1 <=> bwt[i] = 'A'`.
For memory efficieny, we do not use a fifth bit vector for the $-character, as it only has a single 1-bit.
Rather, we simply store the offset of this 1-bit in the `dollarPos` attribute.

To find the number of occurrences of e.g. 'A' before index `i`, we simply count the number of 1-bits before index `i` in the appropriate bitvector (**HINT:** remember the rank function).
If the occ value is queried for $-character (charIdx = 0), you need to asses if the index is greater than the position of this sentinel character in the BWT (stored in `dollarPos`).

1. Implement the `occ` function in the file `fmindex.cpp`.
```C++
length_t FMIndex::occ(const length_t& charIdx, const length_t& index) const 
```

2. You will also need to implement the function `findLF` in `fmindex.cpp`
```C++
length_t FMIndex::findLF(length_t k) const 
```

**HINT:** note that `FMIndex::occ` takes a character **index** as a first argument, not a character itself.

---
**TESTING**

Build your project and execute the following commands to test your solution to this assignment.

```
cd build/unittest
./week1 --gtest_filter=*FunctionalityTest.occTest
./week1 --gtest_filter=*FunctionalityTest.findLFTest
```

#### SA access

Next, we need to look up the element at index `k` in the suffix array.
For this we use the `sparseSA` attribute, which uses your implementation of a sparse suffix array.

You need to check if the element at index `k` in the **original** suffix array is stored in the sparse suffix array (using `hasStored(k)`).
If this is the case, you simply use the `sparseSA[k]` to retrieve the element.
Otherwise, you need to iteratively apply the LF mapping (using `findLF`), until an index is found that is stored.
Once you found such index, say `l`, for which `hasStored(l)` is true, after applying the LF relation `j` times, you can return `sparseSA[l] + j`.
**HINT:** make sure that the returned value is within the range [0, text.size()).

Complete the function:
```C++
length_t FMIndex::findSA(length_t k) const
```
in the file `fmindex.cpp`.

---
**TESTING**

Build your project and execute the following commands to test your solution to this assignment.

```
cd build/unittest
./week1 --gtest_filter=*FunctionalityTest.findSATest
```
---

#### Extending a partial match with a character to the left

Given a pattern `P` and its range `R = [b, e)` over the suffix array that contains the start positions of all suffixes of the text that are prefixed by `P`.
The range `R' = [b', e')` corresponding to the range over the suffix array containing the start positions of all suffixes of the text that are prefixed by `cP` can be computed  from `R`, the `counts` array and `occ` function in the following way: `b' = counts[i] + occ(i, b)` and `e' = counts[i] + occ(i, e)`, where `i` is the character index of `c`.

Implement function `addCharLeft` in the file `fmindex.cpp`:
```C++
bool FMIndex::addCharLeft(length_t charIdx, 
                          const Range& originalRange,
                          Range& newRange) const
```
This function takes the original range over the suffix array of a pattern `P` and calculates the new range of pattern `cP` where `charIdx` is the index of `c`. It returns a bool that indicates whether the new range is valid or not (i.e., the new range is not empty).
At the end of the function the value of parameter `newRange` should be set to the range of pattern `cP`.


---
**HINT**:
To get the begin and end of a Range `range` use
```C++
range.getBegin();
range.getEnd();
```

Remember that the end is **exclusive**.

To set the values `[b, e)` of the Range `newRange` use
```C++
newRange = Range(b, e);
```

---
**TESTING**

Build your project and execute the following commands to test your solution to this assignment.

```
cd build/unittest
./week1 --gtest_filter=*FunctionalityTest.AddCharLeftTest
```

---

### Integrating everything in a read mapper

Before starting this phase make sure that all tests using the following commands pass:
```
cd build/unittest
./week1 --gtest_filter=*FunctionalityTest*
```

#### Exact Matching
Complete the function `matchExact` in the file `fmindex.cpp`.

```C++
vector<length_t> FMIndex::matchExact(const string& str) const
```

This function takes a string as argument and finds all the exact matches of the string in the text.
The function returns a vector with all the start positions of the exact matches in the text.
Start with the range of the empty string i.e., `[0, text.size())`.
By iterating in reverse order over the string (right-to-left) and updating the range (using `addCharLeft`) for each character, find the range over the suffix array which contains all suffixes prefixed by `str`.
To find the actual positions in the text, do a look-up on the suffix array (using `findSA`) for each index in the range.
Remember that you can find the character index of a character `c` by using `sigma.c2i(c)`.

---
**TESTING**
```
cd build/unittest
./week1 --gtest_filter=*IntegrationTest.matchExactTest
```
---



#### Paired-end matching

We now have all the functionality to efficiently solve the bioinformatics read mapping problem.
The goal is to align reads (search patterns) against a reference genome (search text).

DNA is double stranded.
The reference genome that is provided contains the sequence of just one strand (arbitrarily chosen, but let us call it the **forward** strand).
The other strand can be obtained by taking the reverse complement of the forward strand.
The reverse complement of a DNA sequence is obtained by reversing the sequence, followed by taking the complement of each character, i.e. replacing A->T, T->A, C->G and G->C.
For example, the reverse complement of ACCTGAG is CTCAGGT.

Because a read can be sequenced from either strand, in principle, you should try and align each read against both strands.
However, in practice, it is more memory-efficient to construct the FM index for just the forward strand and align twice: once the original read, once the reverse complement of the read.

Read alignment can be ambiguous.
There can be multiple candidate alignment positions in case of repeated subsequences in the reference genome.
There can be multiple candidate alignment positions on the same strand, or even across both strands.
To mitigate this problem, and to find a unique alignment position, reads are provided in pairs where the first read align to one strand (you do not know a priori which one!) and the second read aligns to the **other** strand.
Additionally, the distance between the extreme ends of the reads (called the insert size) is known.
This is illustrated below:

![](../images/paired-end-reads.png?raw=true)

The insert size typically follows a normal distribution for which the mean value is known.
Therefore, given a pair of reads, you should find the alignment positions where both reads map to different strands, and where the distance between the extreme ends of the reads best match the provided mean insert size.

Update the function `bestPairedMatch` in the file `fmindex.cpp`
```C++
tuple<length_t, length_t, bool>
FMIndex::bestPairedMatch(const pair<string, string>& reads,
                         const length_t& meanInsSize) const
```

This function takes a pair of reads and the mean insert size and finds the best paired match.
The function returns a tuple of three elements: the first contains the start position of the match of the first read (of the pair), the second element contains the start position of the match of the second read and the third element is a bool which is true if the first read is matched along the forward strand.
Note that start positions are always reported for the forward strand.
Choose the combination of start positions that is the most likely given the mean insert size. 

---
**EXAMPLE**

Consider a pair of reads and a mean insert size of `900`:

"CCCTGCCCAGGCAACAGCCCCAGCCCCCGAGCTGACTGTTGTCGCAACTA" (read 1 of the pair) and 

"TGAGTTCTGCATTAAAGCGGACACGTCGAAAATACCGTTGCCAACATTTT" (read 2 of the pair).

The length of both reads is 50.

First, we will assume that read 1 maps to the forward strand and read 2 maps to the reverse strand.

Read 1 matches exactly to position `2825189` in the reference genome. 
The reverse complement of read 2 matches exactly to position `2826121` in the reference genome. 
Hence, for this assumption only one combination is possible. 
The insertion size for this combination equals the between the end position of the match of the reverse complement of read 2 (`=(2826121 + 50)`) and the start position of the match of read 1 (`2825189`), which equals `982` .

Next, assume that read 1 maps to the reverse strand.
The reverse complement of read 1 matches eactly to start positions `1062049` and `1290352`. Read 2 matches exactly to position `1061117`. 
Hence for this assumption two combinations are possible. 
The insert size equals the difference between the end position of the match of the reverse complement of read 1 (`1062049 + 49` or `1290352 + 49`)) and the start position of the match of read 2 (`1061117`).
These possible insert sizes are `981` and `229284`.

The insert size that is the closest to our mean insert size (`900`) is `981`. 
Thus, we map the reverse complement of read 1 to position `1062049` and read 2 to position `1061117`. The tuple that the function should return is `{1062049, 1061117, 0}`. The final value is 0 because we mapped the reverse complement of read 1 to the forward strand. 

---
**HINT**:
Use the pre-implemented function 
```C++
 std::string revCompl(const std::string& s) const 
```

to take the reverse complement of a string. 

---
**TESTING**
```
cd build/unittest
./week1 --gtest_filter=*IntegrationTest.bestPairedTest
```
---


After completion of this assignment all unit tests shoud pass:
```
cd build/unittest
./week1 
```

## Dependency Graph

Take a look at the dependency graph of this week's assignment:

![](../images/dependency_graph.png?raw=true)

We hope you can see how the different functions depend on a correct implementation of an earlier function.
