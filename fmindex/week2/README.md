# WEEK 2

## Background information
In week 2, you will extend the functionality of the FM index to also support **approximate** pattern matching.
Given a search pattern P, we are interested in finding all occurrences in the search text within a predefined maximum edit distance of P.
Indeed, because of sequencing errors and because of natural variation, there may be small differences between a read itself and the corresponding location in the reference genome.
In this stage you will use a naive backtracing algorithm to find **all** approximate matches under this bound of a single read.


### Provided implementation
You are provided with an existing implementation of the FM index (using your solution to week 1) and some helper classes. 
The relevant helper classes for this week are the classes `Range`, `FMPos`, `FMOcc`, `FMPosExt` and `TextOcc` provided in file `fmindex.h` and class `Substring` in the file `substring.h`.
An `FMPos` object corresponds to a range over the suffix array and the depth (= number of characters) of the already matched substring. 
An `FMPosExt` object is an extension of this, which also stores the character that was last added to the matched substring. 
An `FMOcc` is an occurrence in the FM index, which is a match to the read, with a range over the suffix array, depth of the match and edit distance to the read. 
A `TextOcc` object stores a range in the text of the match and an edit distance value of the match.
A `Substring` object is a wrapper around a string with a start and end index (exclusive) in that string and a direction. 
This direction can either be forward or backward.
A `Substring` can be indexed using the `[i]` operator with `0 <= i < m` and with `m` the length of the substring.
It will convert index `i` to the appropriate index of the original string, taking into account the position and orientation of the substring.


## Assignments

### Banded matrix

The edit distance is a metric for how **dissimilar** two strings `X` and `Y` are. 
This metric represents the minimum number of edit operations (substitutions, insertions or deletions) required to transform string `X` into string `Y`.
To find the best (= most similar) approximate match using the edit distance metric we need to minimize this score.

Let `M(i, j)` represent the edit distance between the first `i` characters of `X` and the first `j` characters of `Y`.
To calculate `M(i, j)` we make use of following recursion formula:
```
                --
                |  M(i - 1, j - 1) + S
  M(i, j) = min |  M(i, j - 1) + 1
                |  M(i - 1, j) + 1
                --
```

with `S = 1` if character `X[i-1] != Y[j-1]` and `S = 0`, otherwise.

Remember that we initialize `M` as:

```
  M(i, 0) = i
  M(0, j) = j
```

If the maximal allowed edit distance `k` is known beforehand (this is the case in this project), it is not useful to calculate cells for which `|i-j|` is greater than `k`.
Hence, only a band around the diagonal with width `2k + 1` needs to be calculated, as all other cells will have a value of at least `k + 1`.

You are provided with an implementation that initializes the matrix and the band. 


#### 1. Update a single cell
The `BandedMatrix` class stores the values in a vector in a row-major order. 
For each row, at most `(2k + 1) + 2` values are needed (`2k + 1` cells with potentially valid values and 2 extra cells that flank the band with a value of `k + 1`). 
The value at `M(i, j)` can be found and set using the custom `operator()(i, j)` or by using `at(i, j)`.

You need to update the following function in the file `bandmatrix.h`
```C++
length_t updateMatrixCell(bool notMatch, unsigned int row,
                          unsigned int column) 
```
which updates `M(row, column)` based on the `notMatch` argument, which indicates if there is a character **mismatch**.
The function both updates the cell and returns the new value in this cell. 
You can assume that all required values have already been filled in correctly. 


This function is tested in the BandedMatrixTest.UpdateMatrixCellTest.
Run this test by using the following commands.

```
cd build/unittest
./wee2 --gtest_filter=*BandedMatrixTest.UpdateMatrixCellTest
```
---

**HINT**

You can use the function 
```C++
void printMatrix() const
```
on a `BandedMatrix` object to print out the matrix for debugging purposes.

---

#### 2. Update a row
In the same file, you will need to complete the function:
```C++
length_t updateMatrixRow(const Substring& pattern, length_t row, char c)
```

This function updates row `row` of the matrix.
At row `row`, the character to consider is `c`.
The horizonal pattern is provided as a `Substring` object.
You can simply access the character that corresponds to column `j` as `pattern[j - 1]`.

This function returns the **minimal** value found at this row. 
You can assume that the element `at(getFirstColumn(row) - 1, row)` is accessible, and if it is not in the band it has a value greater than the maximal value allowed in the band. 
You can also assume that row `row - 1` is filled in correctly.

Some code is already present that handles the case where `row` exceeds or is equal to the number of rows in the matrix. 


---
**HINT**

Use your solution to the assignment above. 
You can find the first and last column of the band that need to be filled in for the current row by applying the following functions.
```C++
 getFirstColumn(row);
 getLastColumn(row);
```

---

This function is tested in the BandedMatrixTest.UpdateMatrixRowTest.
Run this test by using the following commands.

```
cd build/unittest
./week2 --gtest_filter=*BandedMatrixTest.UpdateMatrixRowTest
```


### Functionality of the FM index
 
#### 1. Extending a position
Update the following function in file `fmindex.cpp`:
```C++
void FMIndex::extendFMPos(const Range& range, const length_t& depth,
                          std::vector<FMPosExt>& stack) const
```

This function finds all children positions of the position with this range and depth and adds them to the stack in order of the alphabet. 
The stack stores `FMPosExt` objects, which have a character, range and depth. 
The depth of child position equals the depth of the parent position plus 1. 
The tests expect the children positions (if valid) to be added in lexicographical order.



---

**NOTE**

Adding a position `f` to the stack can be done by using the `push_back` method:
```C++
stack.push_back(f);
```

After applying this function `f` is now the final element of the stack.

You can create an `FMPosExt` object with character `c`, range `r` and depth `d` using:
```C+++
FMPosExt(c, r, d);
```

---

**HINT**

Remember the `addCharLeft` functionality of week 1, which can find the range of pattern `cS` given the range of pattern `S`. 
If the new range is empty, this function returns false, which means that the child position does not exist.

You can iterate over the alphabet (without the sentinel) by using the following for loop:
```C++
 for (length_t i = 1; i < sigma.size(); i++) {}
```
This gives you the character index of each character in the alphabet. 
To convert to a character use `sigma.i2c(i)`.

---

This function is tested in the FunctionalityTest.ExtendTest.

```
cd build/unittest
./wee2 --gtest_filter=*FunctionalityTest.ExtendTest
```


#### 2. Converting an occurrence in the FM index to an occurrence in the text
In week 1 an occurrence in the text could be represented with the startposition of the occurrence, as the occurrence was an exact match to the read. 
Now, the end position of the occurrence and edit distance to the match need to be stored as well. 
A text occurrence thus stores a range (begin and end position in the text) and an edit distance value. 
Similarly, an occurrence in the FM index corresponds to a range, a depth (length of the match) and an edit distance.
In this assignment you will convert an `FMOcc` object to a vector of `TextOcc` objects.

Update the following function in `fmindex.cpp`.

```C++
void FMIndex::convertFMOccToTextOcc(const FMOcc& fmocc,
                                    std::vector<TextOcc>& textOcc) const 
```

For each SA index in the range of `fmocc` (`fmocc.getRange()`) you can find the corresponding start position in the text using the `findSA` method of week 1. 
You can find the depth of the position/length of the match by using `fmocc.getDepth()`, similarly, you can find the edit distance by applying `fmocc.getDistance()`. 
Creating a `TextOcc` object with a given start position `s`, depth `d` and edit distance `e` is achieved by `TextOcc(Range(s, s + d), e)`.


This function is tested in the FunctionalityTest.ConvertTest.

```
cd build/unittest
./wee2 --gtest_filter=*FunctionalityTest.ConvertTest
```




### Integrating everything in an approximate read mapper

#### Naive backtracking
To efficiently find the best approximate match (under a maximal edit distance `k`) of a short read in a large text, we make use of a backtracking algorithm and dynamic programming. 
In this algorithm we visit every position in the search tree in a depth-first order. 
Each position corresponds to a range `r` over the suffix array and a depth `d`. 
Such that that for each start position `s` found in the suffix array in `r`, `text.substr(findSA(s), d)` is the substring `S` of the text corresponding to this position. 
In other words, the range of the position corresponds to the start positions of string `S` in the text (and hence the width of the range corresponds to the number of times `S` occurs in the text) and the depth of the position corresponds to the length of `S`.

We stop exploring a branch of the search tree if the minimal edit distance between `S` and the pattern `P` that is matched is greater than `k`.
To continue along this branch, we try to extend the position which each character in the alphabet (excluding $). 
Thus, at least 1 and at most 4 new positions need to be visited, whose substrings correspond to `AS`, `CS`, `GS` and `TS` (**NOTE** Since the FM index only allows for prepending a character we are matching in the backwards direction!). 
To keep track of the positions we still need to visit we use a stack. 

To start the process we start with `S = "" ` and extend to strings `A`, `C`, `G` and `T` with their corresponding positions.
These positions are pushed to the stack in lexicographical order.
Then the final element is taken from the stack (in the first case this will be the position corresponding to string `T`). 
We fill in the matrix for this position, since we are using a depth-first approach this boils down to updating a single row using the last added character (the first character of the corresponding substring as matching is done in the backwards direction).
If the minimal edit distance value on this updated row does not exceed `k`, we push the (existing) children of this position to the stack in lexicographical order.
This process continues until there are no more positions on the stack.  


To integrate everything you will make a naive backtracking function:
Update the following function in file `fmindex.cpp`:
```C++
vector<TextOcc> FMIndex::naiveApproxMatch(const string& pattern,
                                          length_t k) const
```



In this function you will need a vector `occ` with `FMOcc` objects that stores every found occurrence in the FM index, a stack `stack` with the `FMPosExt` objects that still need to be checked, a banded matrix `matrix` and a `Substring` `p` with backwards direction created from `pattern`. 
Code to create these objects is already present.


As long as the stack is not empty, you will iterate over the stack, pop-back the final `FMPosExt`object `f`  of the stack and update the matrix for this final object `f`. For this you can use the `updateMatrixRow` method you implemented earlier in this assignment. As arguments you can use the `Substring` object `p`, `f.getCharacter()` and `f.getRow()`.

After updating the matrix, you need to check if the minimal edit distance on this row (the return value of the `updateMatrixRow` method) does not exceed the maximal value `k`. 
If this is the case you can immediately continue to the next position on the stack. 
Else, you need to add the child positions of position `f` to the stack (using the `extendFMPos` method).

If updating the matrxi for this row filled in the last column of the matrix (and hence an alignment to the entire read was found) you need to check if the value in the last column is smaller than or equal to `k`, if this is the case you can add the current position and the found value in the final column to the vector with found `FMOcc` objects. 


Finally, after the stack is empty, you can use the pre-implemented function

```C++
std::vector<TextOcc>
FMIndex::filterRedundantMatches(std::vector<FMOcc>& fmocc,
                                const length_t& k) const
```
to convert the occurrences in the FM index to occurrences in the text. 
Additionally, this function removes equal occurrences and occurrences that are equal to another occurrence except from leading and/or trailing gaps. 


---

**HINT**

The `BandedMatrix` class provides the following public helper methods that can come in handy:

```C++
length_t getNumberOfRows() const
bool inFinalColumn(length_t row) const
length_t getValueInFinalColumn(length_t row) const
```

---

**NOTE**

Creating an `FMOcc` object with a `FMPos` position `p` and an edit distance value `v` is achieved using:

```C++
FMOcc(p, v);
```

---



This function is tested in the IntegrationTest.NaiveMatchesTest.

```
cd build/unittest
./week2 --gtest_filter=*IntegrationTest.NaiveMatchesTest
```


## Dependency Graph

Take a look at the updated dependency graph of the assignment:

![](../images/dependency_graph_week_2.png?raw=true)

The functions in grey are functions that you implemented in week 1, those in white are new functions from this week's assignment.

We hope you can see how the different functions depend on a correct implementation of an earlier function.


