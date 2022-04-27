# WEEK 3

## Background information
In week 2, you focused on a naive backtracking approach for approximate alignment. What you might have noticed is that this approach spends most of its time on branches in the search tree that will turn out to be unsuccesful anyway. Hence, this is not an efficient method. If you run the commands:

```
cd build
./demo
```

a program will start that aligns 9816 reads of length 50 with `k = 4` using the naive approximate matching method to the CP001363 reference genome. This demo takes a considerable amount of time. If you then know that the human reference genome is more than 600 times as large as the CP001363 reference genome and that 9816 is only a very small fraction of the number of reads from a typical sequencing experiment, it becomes clear that this naive approach is unfeasible in real-life scenarios.

Hence, other techniques need to be considered. One of these are the so-called search schemes. These search schemes exploit the fact that the search tree is much denser near the root than deeper down into the tree. This is easy to verify as positions near the root of the search tree correspond to short sequences, which  have a  high probability of occurring at least once in the reference genome. Deeper down into the tree the positions correspond to longer sequences withs subsequent smaller probability of occurring. Put more simply, near the root each position has 4 child positions, while deeper down this is no longer the case. The core idea of search schemes is to get deeper down into the tree in a fast manner (hence without exploring all branches near the root), and only when we are deep enough into the tree we start exploring all possible child positions. 

Search schemes need a bidirectional index, which allows you to switch direction at any point. In this week's assignment, you will extend the FM Index to a bidirectional FM Index, which supports the use of search schemes.


## Provided implementation
You are provided with an existing implementation of the FM-Index (using your solution to week 1/2), an imlementation for search schemes, a bidirectional FM-Index and some helper classes. The new relevant helper classes for this week are the classes `RangePair`, `BiFMPos`, `BiFMOcc`, `BiFMPosExt` and `Search` provided in file `src/bidirectionalfmindex.h` and class `CumulativeBitVectors` in file `cumulativebitvectors.h`.
A `RangePair` object contains two ranges `backwardRange` and `forwardRange`. The backward range corresponds to the range used in traditional backwards matching (= range over the BWT of `T`) and the forward range corresponds to the range over the BWT of `T'`.
A `BiFMPos` object is similar to an `FMPos` object, the only difference is that now an extra forward range is stored. Likewise, a `BiFMPosExt` object is an `FMPosExt` object extended with a forward range and a `BiFMOcc` is an `FMOcc` extended with a forward range. 
A `CumulativeBitVectors` object stores bitvectors for each of the characters in the alphabet (excluding the sentinel) and the position of the sentinel in the BWT (or BWT of reversed text).

## Assignment
In this project you will build a bidirectional FM Index that supports the use of search schemes. 



### Construction of the bidirectional FM-Index
The bidirectional FM Index was first introduced by Lam et al. in 2009. It is realized by also storing the BWT of the reversed text `T'`, in addition a second occurrence table occ' (of the BWT of `T'`) is used. By keeping track of both the range `[i, j)` over the BWT of `T` and `[i', j')` over the BWT of `T'` in a synchronized manner, one can extend a pattern `P` to both `cP` as well as `Pc`.

#### Creating the BWT of the reverse text
To construct the bidirectional FM-Index we need to find the BWT of the reverse text. You are provided with the (dense) suffix array of the reverse text and need to use this and the original text to build the BWT over the reverse text. This exercise is very similar to creating the BWT of the original text from week 1, but with a small twist.

Fill in function:
```C++
void BiFMIndex::createRevBWTFromRevSA(const vector<length_t>& revSA,
                                      string& revBWT)
```

in file `bidirectionalfmindex.cpp`.

You can assume that the revBWT parameter corresponds to the empty string prior to a call to this function. 
At the end of this function the revBWT parameter needs to store the BWT of the reverse text.

This function is tested using the ConstructionTest.CreateRevBWTTest.

```
cd build/unittest
./week2 --gtest_filter=*ConstructionTest.CreateRevBWTTest
```


### PrevOcc tables

For efficient bidirectional extending (see later), the bidirectional FM Index uses a so called previous occ (prevocc) table instead of an occ table. 

Where the occ table stores the number of occurrences (up untill that point) in the BWT of character `c`, the prevocc table stores the number of occurrences of characters `b <= c` (so-called previous characters) up untill that point in the BWT of the text. Similarly, the prevocc' table stores the number of occurrences of characters `b <= c` (so-called previous characters) up untill that point in the BWT of the **reversed** text. These prevocc tables need to support both `occ` and `prevocc` queries.


Formally:

1. Prevocc query: The number of occurrences of characters smaller than `c` up untill index `i` is
```
(c == $) ? (0 : prevocctable(c - 1, i))
```

2. Occ query: The number of occurrences of `c`  up untill index `i` is
```
prevocctable(c, i) - (c == $) ? (0 : prevocctable(c - 1, i))
```

These tables can be stored using bitvectors `bvs`.
To save memory, we opt not to store the bitvector for character $, as this character occurs only once in the BWT. Instead we store the position `dollarPos` of this character. 
The prevocc query with character index `ci` and index `k` then becomes:
1. 0 if `ci == 0`
2. `(k <= dollarPos)` if `ci == 1`
3. `bvs[ci-2].rank(k) + (k <= dollarPos)` if `ci > 1`

The occ query with character index `ci` and index `k` then becomes:
1. `(k <= dollarPos)` if `ci == 0`
2. `bvs[ci-1].rank(k) ` if `ci == 1`
3. `bvs[ci-1].rank(k) - bvs[ci-2].rank(k)` if `ci > 1`

Using this information, complete functions:

```C++
 size_t occ(int cIdx, size_t k) const 
 size_t prevOcc(int cIdx, size_t k) const
```
in file `cumulativebitvec.h`. You can use the `bvs` member of this datastructure.

These functions are tested in the CumulativeBitvec.OccTest en CumulativeBitvec.PrevOccTest.

```
cd build/unittest
./week2 --gtest_filter=*CumulativeBitvec*
```

### Functionality



#### Extending with a single character

Assume you want to extend pattern `P` to `Pc` (= right extension), given ranges `[i, j)` (= the backward range) and `[i', j')` (= the forward range). Computing the new forward range `[k', l')` over the BWT/SA of `T'` is the same as doing a backwards extension on `P'` using the forward range and the `occ` query on the BWT of `T'`.

The non-trivial issue is how to compute the new backward range `[k, l)`. We know that this new range is a subrange of `[i, j)` as suffixes with prefix `Pc` are also suffixes with prefix `P`, and range `[i, j)` corresponds to suffixes with prefix `P`. 
Denote this new subrange as `[k, l) = [i + x, i + x + y)`. To find `x`, one needs to consider that `[i, j)` contains all the indices of suffixes prefixed by `P` in **lexicographical order**, hence range `[i, i + x)` corresponds to suffixes prefixed by `Pb` with `b < c`. The number of suffixes of `T` prefixed by `Pb` equals the number of suffixes of `T'` prefixed by `bP'`
Hence:

<img src="https://render.githubusercontent.com/render/math?math=x = (\mathrm{prevocc'}(b, j') - \mathrm{prevocc'}(b, i'))">

Put differently, `x` equals the number of occurrences of characters smaller than `c` in range `[i', j')` over the BWT of `T'`.
The final unknown is `y`, the width of this new range. The width of the range over the BWT of `T` corresponding to pattern `cP` equals the width of the range over the BWT of `T'` of pattern `Pc`, hence <img src="https://render.githubusercontent.com/render/math?math=y = l' - k' )">

Extending in the backwards direction is analogous.

Fill in functions: 
```C++
bool BiFMIndex::addCharRight(length_t charIdx, const RangePair& originalRanges,
                             RangePair& newRanges) const
bool BiFMIndex::addCharLeft(length_t charIdx, const RangePair& originalRanges,
                             RangePair& newRanges) const
```
in file `bidirectionalfmindex.h`.

These function takes the original ranges over the suffix arrays of `T` and `T'` of a pattern `P` and calculates the new ranges of pattern `Pc` (or `cP`) where charIdx is the index of `c`. It returns a bool that indicates whether the new ranges are valid or not (i.e., the new ranges are not empty).


These functions are tested in FunctionalityTest.AddCharRightTest and FunctionalityTest.AddCharLeftTest
```
cd build/unittest
./week2 --gtest_filter=*ConstructionTest.CreateRevBWTTest
```
---
**HINT**

To get the forward and backward range of a `RangePair` object `r` use:
```C++
r.getForwardRange();
r.getBackwardRange();
```

---

---
**HINT**

You can perform `occ` and `prevocc` queries on the `backwardOccTable` and `forwardOccTable`. For adding a char to the left you need the `backwardOccTable`, as this table stores the occurrences of characters in the BWT of `T`, similary the `forwardOccTable` stores the occurrences of characters in the BWT of `T'`.

---


These functions are tested in FunctionalityTest.AddCharRightTest and FunctionalityTest.AddCharLeftTest
```
cd build/unittest
./week2 --gtest_filter=*FunctionalityTest.AddChar*
```


#### Extending a position

Update the following function in file `bidirectionalfmindex.cpp`:
```C++
void BiFMIndex::extendFMPos(const RangePair& ranges, const length_t& depth,
                            vector<BiFMPosExt>& stack)
```

This function corresponds to the `extendFMPos` function of week 2, only now bidirectional positions are used. The implementations will be very similar. 

You can assume that direction is correctly set and that using `(this->*addChar)` is a call to either `addCharLeft` or `addCharRight`, depending on the direction.


This function is tested in the FunctionalityTest.ExtendTest.

```
cd build/unittest
./week2 --gtest_filter=*FunctionalityTest.ExtendTest
```

#### Exact Matching


Update the following function in file `bidirectionalfmindex.cpp`:
```C++
RangePair BiFMIndex::matchExactBidirectionally(const Substring& str,
                                               RangePair ranges) const
```

This function extends ranges, corresponding to a previously matched pattern, with the given substring, in the current direction of the index (which equals the direction of the substring). If the given ranges are empty, this means that you have to start with ranges corresponding to the empty string.

You can assume that the direction is correctly set and that using `(this->*addChar)` is a call to either `addCharLeft` or `addCharRight`, depending on the direction.

If no exact match can be found empty ranges should be returned.

This function is tested in the FunctionalityTest.matchBidirectionallyTest.
```
cd build/unittest
./week2 --gtest_filter=*FunctionalityTest.matchBidirectionallyTest
```

### Integration



Finally, you will provide the support for search schemes. 

#### Search schemes: Formal definitions

1. **Search**:
A search `S` is a tuple `(pi, L, U)` of strings for a partitioning with `p` parts. `pi`is a permutation string of length `p` over `0, 1, ... p - 1` and defines the order in which the parts of the pattern are searched. `L[i]` and `U[i]` indicate the lower and upper bounds for the cumulative number of errors in all parts up to and including part `pi[i]` for `0 <= i < p`.

2. **Error pattern**:
An error pattern is an integer string `A` of length `p`, where `A[i]` equals the number of errors in in part `pi[i]`. The weight `w` of `A` is defined as <img src="https://render.githubusercontent.com/render/math?math=w = \sum_{i = 0}^{p - 1} A[i] )">.

3. **Error pattern coverage**:
An error pattern `A` is covered by a search `S = (pi, L, U)` if 
<img src="https://render.githubusercontent.com/render/math?math=L[i] \leq (\sum_{j = 0}^{i} \leq A[\pi[j]]) \leq U[i] )">, with `0 <= i < p`

4. **Search scheme for `k` errors**:
A search scheme `SS` is a set of searches for which each possible error pattern `A` with weigh `w <= k` is covered by at least one search.


5. **Connectivity property**:
A search `S = (pi, L, U)` is connected if <img src="https://render.githubusercontent.com/render/math?math=\pi[i] \in \min_{j < i} \pi[j] - 1, \max_{j < i} \pi[j] + 1 )">, with `0 <= i < p`. Only searches that satisfy the connectivity property are valid.

#### Search schemes: Examples

The easiest search scheme to explain is based on the pigeon hole principle. Imagine you want to find every (approximate) occurrence of a read `R` with at most three errors and that you have partitioned `R` into four parts. For each approximate occurrence `O` of `R` with at most three errors one of the four parts must appear **exactly** in `O`.
Finding the exact match of part `P_0` takes exactly `|P_0|` extensions, regardless of the denseness of the region. Hence, we can start by finding the exact match of `P_0`, which takes us deeper into the tree, and only after this part is fully matched we start the approximate matching procedure for `P_1`, `P_2` and `P_3`.
Similarly, we can start another exact search for `P_1` and afterwards extend it with a (forward) approximate match of `P_2` and `P_3`, followed by a (backwards) approximate match of `P_0`.


The pigeon hole search scheme for 3 errors is:
```
(0123, 0000, 0333) 
(1230, 0000, 0333)
(2310, 0000, 0333)
(3210, 0000, 0333)
```

This search scheme is very symmetric, however it is also pretty redundant as many error patterns are covered by more than one search. A more optimized search scheme for 3 errors is:

```
(0123, 0000, 0133)
(1023, 0011, 0133)
(2310, 0000, 0133)
(3210, 0011, 0133)
```

In this project, for all searches `U[0]` will be 0, meaning that each search will start with an exact match of one of the parts. 
You are provided with an implementation of search schemes (searchscheme.h), that can handle any search scheme for which `U[0] == 0` for all searches.

The search scheme handles each search in the same manner. Step 1 is performing an exact bidirectional match, using your implemented function from the previous exercise. After the exact matching stage an approximate stage will start using the following function (which you need to update) in `bidirectionalfmindex.cpp`:
```C++
void BiFMIndex::recApproxMatch(const Search& s, const BiFMOcc& startOcc,
                               vector<FMOcc>& occ,
                               const vector<Substring>& parts, const int& idx)
```

This recursive function goes through each part of the search, for which `U[idx] > 0`. The `startOcc` parameter corresponds to the bidirectional occurrence of the already matched pattern. The `occ` vector corresponds to matches of the **fully** matched pattern. During the execution of this fuction new occurrences can be added to the vector.
The `parts` vector corresponds to the different parts of the pattern. You can assume that the direction of each part corresponds to the correct direction this part has in this search.
The `idx` corresponds to the current index in the search. At this stage `parts[s.getPart(idx)]` will be matched approximately. The current lower and upper bounds are `s.getLowerBound(idx)` and `s.getUpperBound(idx)`. The current direction is `s.getDirecion(idx)`.


This function performs an approximate match, akin to the naive algorithm of week 2, on this part.
For this reason a banded matrix is used. This matrix slightly differs from the one used in week 2, as it has to be initialized with `startOcc.getDistance()` to account for the errors made in previous stages. The width of the band is then `U - startOcc.getDistance()`. Also, you will again need a stack, now with `BiFMPosExt` objects.

Again, you will iterate over the stack and update the matrix for each iteration. This procedure mirrors the one from the exercise in week 2. The difference, however, is when we have reached the end of the part (when we filled in the last column of the matrix). If the value in this column is between the lower and upper bound (**inclusive!**), we have found a valid **partial** occurrence. If this part is the final part of the search we can add this partial occurrence to the occurences vector. Else iff this is not the final part, we must continue the search using this partial occurrence as the start occurrence for the part with the next index.

---

**HINT**

The depth of your new found partial occurrence equals the depth of the start occurrence plus the row in the matrix.
You can know if the current stage corresponds to the end part of the search by applying:

```C++
s.isEnd(idx)
```
which returns true if this the end part. 

---

---

**HINT**
At the start of this function, you must set the direction of the index! For this you can use
```C++
setDirection(dir)
```
where `dir` is the correct direction. 

This also means that after each call to `recApproxMatch`, the status of the direction member could be changed. Remember to reset the direction to the correct direction after a call to `recApproxMatch`.

---

This function is tested using the IntegrationTest.SearchSchemeApproxTest

```
cd build/unittest
./week2 --gtest_filter=*IntegrationTest.SearchSchemeApproxTest
```

---

**NOTE**

This test will take some time as it checks whether the solution using the search schemes equals the solution using the naive method and the naive method is slow.


---


## Demo

Go to file `demo.cpp`, comment out lines 42 - 52 and uncomment lines 54 - 65. Build your project and run `./demo` again from the build folder. You will notice the increase in speed using the search scheme as opposed to the demo using the naive approximate method. This demo makes use of the pigeon hole search scheme. To go even faster, you can use a more optimized search scheme by changing line 54 to:

```C++
 SearchScheme ss(bifmindex, "../search_schemes/kuch_k+1/", 4);
```

---

**NOTE**

Make sure that you have built your project in Release mode to get the best performance. Remember to rebuild the project after every change you have made.

---
