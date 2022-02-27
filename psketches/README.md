# Parallel sketches: KLL and REQ

Based on the sequential C++ implementation of KLL and REQ algorithms for quantile computation by the Apache Software Foundation - [Apache DataSketches](https://datasketches.apache.org) ([github repo](https://github.com/apache/datasketches-cpp)). The code from the Apache Software Foundation is in the `DataSketches` folder.

To build:

cmake -S . -B build/Release -DCMAKE_BUILD_TYPE=Release
cmake --build build/Release

Executable found in dir 'build/Release'
