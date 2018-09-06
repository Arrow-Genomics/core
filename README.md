# core

#### Preparing
1. First build and install the [Arrow C++](https://github.com/apache/arrow/tree/master/cpp) on the host machine.

2. Put together "seqan" and "core" repos in a folder, i.e, "lib"<br />

3. Create a folder in "lib": /core-build/release<br />

4. Add these packages and libraries in core/CMakeLists.txt file after SeqAn package.
```
find_package( Threads )
link_libraries( ${CMAKE_THREAD_LIBS_INIT} )
find_library(LIB_ARROW arrow)
```

#### Building:
```
cmake ../../core    -DCMAKE_PREFIX_PATH="$HOME/lib/seqan/util/cmake"    -DSEQAN_INCLUDE_PATH="$/lib/seqan/include"
make
```

#### Running:
```
./core /path/to/file/sam(bam)file.sam(.bam)
```
