# WEEK 2

## Background information
In week 2, you will extend the functionality of the FM Index to support approximate matching under a predefined maximal edit distance. Approximate matching is necessary since sequencing can introduce errors in the reads and because of natural variation. 
In this stage you will use a naive backtracing algorithm to find **all** approximate matches under this bound of a single read. 

To achieve this you will make use of a banded matrix and dynamic programming. 
 


### Provided implementation
You are provided with an existing implementation of the FM Index (using your solution to week 1) and some helper classes. The relevant helper classes for this week are the classes `Range`, `FMPos`, `FMOcc`, `FMPosExt` and `TextOcc` provided in file `fmindex.h` and class `Substring` in file `substring.h`.
An `FMPos` object corresponds to a range over the suffix array and the depth of the already matched substring. An `FMPosExt` object is an extension of this, which also stores the character that was last added to the matched substring. An `FMOcc` is an occurrence in the FM Index, which is a match to the read, with a range over the suffix array, depth of the match and edit distance to the read. A `TextOcc` object stores a range in the text of the match and an edit distance value of the match.
A `Substring` object is a wrapper around a string with a start and end index (exclusive) in that string and a direction. This direction can either be forward or backward. If the direction is backward the `[i]` operator finds the `i`th character starting from the end index.




## Assignments

### Banded matrix

The edit distance is a metric for how **dissimilar** two strings are. This metric sums the the minimum number of operations (substitutions, insertions or deletions) required to transform one string into the other. To find the best approximate match using the edit distance you want to **minimize** this score.


We make use of a banded matrix `M`. We set our pattern (or read) `P` as the horizontal sequence. The vertical sequence corresponds to a subsring `S` of our text `T`.
Assume, that we have filled in the matrix for substring `S_1` for which `|S_1| = d - 1`. 
This means that we have filled in rows `[0, d - 1]`.
From row `d - 1` we can calcuate the values for substring `(S_1)c`.
We need to keep the previous rows in storage for backtracking purposes.
As we are matching in the backwards (as the FM Index allows only for backwards matching) direction, the character of `S` (a substring of length `d`) at row `d` corresponds to `S[0] = c_s` and the character of `P` (in reverse) at column `j` corresponds to `P[|P| - j] = c_p`

To calculate the value in cell `M(d, j)` we make use of a recursion formula:
```
                --
                |  M(d - 1, j - 1) + (c_p != c_s)
  M(d, j) = min |  M(d, j - 1) + 1
                |  M(d - 1, j) + 1
                --
```

If the maximal allowed edit distance `k` is known beforehand (which is the case in this project), it is not useful to calculate cells for which the minimal value is greater than `k`. Hence, only a band around the diagonal with width `2k + 1` needs to be calculated, as all other cells will have a value of at least `k + 1`.

You are provided with an implementation that initializes the matrix and the band. 


#### 1. Update a single cell
The `BandedMatrix` class stores the values in a vector in a row-major order. For each row at most `(2k + 1) + 2` are needed (`2k + 1` cells with potentially valid values and 2 extra cells that flank the band with a value of `k + 1`). The value at `M(i, j)` can be found and set using the custom `operator()(i, j)` or by using `at(i, j)`.

You will need to update the following function in file `bandmatrix.h`
```C++
length_t updateMatrixCell(bool notMatch, unsigned int row,
                              unsigned int column) 
```
which updates `M(row, column)` based on the `notMatch` argument, which indicates if `c_p != c_s` is true. 
The function both updates the cell and returns the new value in this cell. You can assume the rows `<row` have been filled in correctly. 


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
Additionally, in the same file, you will need to update the function:
```C++
length_t updateMatrixRow(const Substring& pattern, length_t row, char c)
```

This function updates the row of the matrix with character `c = c_s`. It makes use of the `Substring` class, this is a custom class that allows for setting the direction of this string. If the direction is set to backwards `pattern[j]` will find the character at index `|P| - j - 1`. You can assume that the direction of the pattern is set correctly. Hence finding the character that correspond to column `j` is done by using `pattern[j]`

This function returns the minimal value found at this row. If the row exceeds or equals the number of rows in the matrix (and hence this is an empty row), the function should return a value greater than the maximal allowed value in the band (`at(0,0) + W`).
You can assume that the element `at(getFirstColumn(row) - 1, row)` is accessible, and if it is not in the band it has a value greater than the maximal value allowed in the band. You can also assume that the band for row `row - 1` is filled in correctly.

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


### Functionality of the FM Index
 
#### 1. Extending a position
Update the following function if file `fmindex.cpp`:
```C++
void FMIndex::extendFMPos(const Range& range, const length_t& depth,
                          std::vector<FMPosExt>& stack) const
```

This function finds all children positions of the position with this range and depth and adds them to the stack in order of the alphabet. The stack stores `FMPosExt` objects, which have a character, range and depth. The depth of child position equals the depth of the parent position plus 1. The tests expect the children positions (if valid) to be added in lexicographical order.

This function is tested in the FunctionalityTest.ExtendTest.

```
cd build/unittest
./wee2 --gtest_filter=*FunctionalityTest.ExtendTest
```

---
**HINT**

Remember the `addCharLeft` functionality of week 1, which can find the range of pattern `cS` given the range of pattern `S`. If the new range is empty, this function returns false, which means that the child position does not exist.

You can iterate over the alphabet (without the sentinel) by using the following for loop:
```C++
 for (length_t i = 1; i < sigma.size(); i++) {}
```
This gives you the characterindex of each character in the alphabet. To convert to a character use `sigma.i2c(i)`.

---


#### 2. Converting an occurrence in the FM Index to an occurrence in the text
In week 1 an occurrence in the text could be represented with the startposition of the occurrence, as the occurrence was an exact match to the read. Now, the end position of the occurrence and edit distance to the match need to be stored as well. A text occurrence thus stores a range (begin and end position in the text) and an edit distance value. Similarly, an occurrence in the FM Index corresponds to a range, a depth (length of the match) and an edit distance.
In this assignment you will convert an `FMOcc` object to a vector of `TextOcc` objects.

Update the following function in `fmindex.cpp`.

```C++
void FMIndex::convertFMOccToTextOcc(const FMOcc& fmocc,
                                    std::vector<TextOcc>& textOcc) const 
```

For each SA index in the range of `fmocc` (`fmocc.getRange()`) you can find the corresponding start position in the text using the `findSA` method of week 1. You can find the depth of the position/length of the match by using `fmocc.getDepth()`, similarly, you can find the edit distance by applying `fmocc.getDistance()`. Creating a `TextOcc` object with a given start position `s`, depth `d` and edit distance `e` is achieved by `TextOcc(Range(s, s + d), e)`


This function is tested in the FunctionalityTest.ConvertTest.

```
cd build/unittest
./wee2 --gtest_filter=*FunctionalityTest.ConvertTest
```




### Integrating everything in an approximate read mapper

#### Naive backtracking
To efficiently find the best approximate match (under a maximal edit distance `k`) of a short read in a large text, we make use of a backtracking algorithm and dynamic programming. In this algorithm we visit every position in the search tree in a depth-first order. Each position corresponds to a range `r` over the suffix array and a depth `d`. Such that that for each start position `s` found in the suffix array in `r`, `text.substr(findSA(s), d)` is the substring `S` of the text corresponding to this position. 
We stop exploring a branch of the search tree if the minimal edit distance between `S` and the pattern `P` that is matched is greater than `k` or if the depth of our pattern in the text equals the depth of the read we are trying to match plus `k`.
To continue along this branch, we try to extend the position which each character in the alphabet (excluding $). Thus, at least 1 and at most 4 new positions need to be visited, whose substrings correspond to `AS`, `CS`, `GS` and `TS`. 
To keep track of the positions we still need to visit we keep track of a stack of these positions. 

To start the process we start with `S = "" ` and extend to strings `A`, `C`, `G` and `T` with their corresponding positions. Then the position corresponding to string `T` is analyzed and its children are pushed to the stack and so on. This process continues if there are no more positions on the stack.  

---

**NOTE**

Since the FM Index only allows for prepending a character we are matching in the backwards direction!

---


To integrate everything you will make a naive backtracking function:
Update the following function in file `fmindex.cpp`:
```C++
vector<TextOcc> FMIndex::naiveApproxMatch(const string& pattern,
                                          length_t k) const
```



In this function you will need a vector with `FMOcc` objects that stores every found occurrence in the FM Index, a stack with the `FMPosExt` objects that still need to be checked, a banded matrix and a `Substring` with backwards direction created from `pattern`. 
As long as the stack is not empty, you will iterate over the stack, pop-back the final `FMPosExt` `f` object of the stack and update the matrix for this final object `f`. For this you can use the `updateMatrixRow` method you implemented earlier in this assignment. As arguments you can use the `Substring` object you created from pattern, `f.getCharacter()` and `f.getRow()`.
After updating the matrix, you need to check if the minimal edit distance on this row (the return value of the `updateMatrixRow` method) does not exceed the maximal value `k`. If this is the case you can immediately continue to the next position on the stack. 
Else, you need to add the child positions of position `f` to the stack (using the `extendFMPos` method).
If this row filled in the last column of the matrix (and hence an alignment to the entire read was found) you need to check if the value in the last column is smaller than or equal to `k`, if this is the case you can add the current position and the found value in the final column to the vector with found `FMOcc` objects. 


Finally, after the stack is empty, you can use the pre-implemented function

```C++
std::vector<TextOcc>
FMIndex::filterRedundantMatches(std::vector<FMOcc>& fmocc,
                                const length_t& k) const
```
to convert the occurrences in the FM Index to occurrences in the text. Additionally this function removes equal occurrences and occurrences that are equal to another occurrence except from leading and/or trailing gaps. 


---
**HINT**



The `BandedMatrix` class provides the following helper methods that can come in handy

```C++
length_t getNumberOfRows() const
bool inFinalColumn(length_t row) const
length_t getValueInFinalColumn(length_t row) const
```

---


This function is tested in the IntegrationTest.NaiveMatchesTest.

```
cd build/unittest
./wee2 --gtest_filter=*IntegrationTest.NaiveMatchesTest
```


