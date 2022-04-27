# CCB project on pattern matching
Welcome to the CCB project on pattern matching.
This project consists out of three parts (week 1, week 2 and week 3).
During the first week you will build a paired-end read mapper based on the FM index. At this stage, the reads are considered to be error-free and exact pattern matching is used.
In week 2, you will extend this functionality to also support approximate pattern matching, i.e., allowing for a certain number of errors between the search pattern and occurrences.
In week 3, the concept of search schemes is used to accelerate approximate pattern matching.

We will be using a bacterial genome (CP001363) as a search text on the one hand and simulated reads, both error-free reads and reads with errors on the other hand. You can find all required data in the folder `testset/`.

## Getting started

To start the project, clone this repository and build the project. We recommend a Linux system with CMake 3.0 and the GNU GCC compiler. However, since we are using standard C++, you should be able to build the source code on other systems (Windows, MacOS, ...) as well.
```
mkdir build
cd build
cmake ..
make
```

All code is in the `src/` folder.
The folders `week1`, `week2` and `week3` each contain a README.md file with information and instructions for that week.

Throughout the source code, we sometimes make use of the macro ALPHABET which denotes the size of the alphabet and is set to a value of 5 for this project. You can assume that the search text (the genome) and the search patterns (the reads) only contain characters $ACGT.
Additionally, you will see the data-type `length_t`, this represents an unsigned 32-bit integer.

For each assignment, the number of lines of code that is expected is indicated in source code. This is only a rough indication and depending on your coding style, your implementation may require fewer or more lines of code.



## Testing your solution

All required data is provided in the `testset/' folder.  
To test your solution for week X execute the following commands:

```
cd build
cd unittest
./weekX
```
---
**HINT**

The Google Test framework allows you to set filters, such that not each test in the suite is run. More information can be found [here](https://github.com/google/googletest/blob/main/docs/advanced.md#running-a-subset-of-the-tests).

---

## Debugging

You can use the main executable for debugging purposes. Some code to get you started is already present. To run the main executable:

```
cd build
./fmindex
```

## Vectors and strings in STL C++

This project makes use of standard C++ vectors and strings.
Vectors in C++ are standard containers that provide dynamic-array like functionality, with fixed time random access to all elements. Vectors are defined as a template class available in the standard library. A c++ string `std::string` can be seen as a vector of `char` elements.

### Declaring a vector

To declare an empty vector you can use:

 ```C++
 // Define a vector that will contain integers
 std::vector<int> myIntsContainer;
 // define a vecror that wil contain doubles
 std::vector<double> myDoublesContainer;
 ```

 ### Initialising

 To initialise a vector with elements you can use two options:
 ```C++
 // Method 1
std::vector<int> myVec {1,3,5};

// Method 2 (equivalent to method 1)
std::vector<int> myAltVec = {2,4,6};
 ```

 ### Resizing a vector
 You can allocate memory by resizing the vector. 
 ```C++
 std::vector<int> myVec;
 myVec.resize(20)
 ```

 The vector `myVec` now contains 20 times 0 (the default constructed value of an `int`). You can update the value at index `i` (`0 <= i < 20`) with `myVec[i] = newValue`.

 ### Adding elements
 To add an element to the end of the vector you can use:
 ```C++
 std::vector<int> myVec {1,3,5};
 myVec.push_back(4);
 ```
 After this operation `myVec` will now contain elements 1, 3, 5 and 4 (in that order).

 ### Accessing elements
 To access the element at index `i` in the vector `myVec` (provided that `0 <= i < myVec.size()`) use:
 ```C++
  std::vector<int> myVec {1,3,5};
  int i = 2;
  int element = myVec[i];
 ```
 `element` now has value 5.

 Similarly, you can use the `[]` operator to overwrite elements:

 ```C++
  std::vector<int> myVec {1,3,5};
  myVec[1] = 20;
 ```

 The vector `myVec` now contains values 1, 20 and 5.

 ### Iterating over a vector
 To iterate over a vector from front to back you have several options. The following examples will print out the values of the vector on separate lines:
 ```C++
 std::vector<int> myVec {1,3,5};

 // Option 1
 for(int i : myVec){
    std:cout << i << std::endl;
 }

 // Option 2
 for(auto it = myVec.begin(); it != myVec.end(); it++){
     std::cout << *it << std::endl;
 }
 // Option 3
 for(size_t i = 0; i < myVec.size(); i++){
     std::cout << myVec[i] << std::endl;
 }
 ```

 To iterate over a vector in reverse order you can use:

 ```C++
 for(auto it = myVec.rbegin(); it != myVec.rend(); it++){
     std::cout << *it << std::endl;
 }

 for(size_t i = myVec.size(); i-- > 0;){
     std::cout << myVec[i] << std::endl;
 }
 ```

 This will print the elements of `myVec` line by line in reverse order.
 

 
