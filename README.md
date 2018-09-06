# core

Put together "seqan" and "core" repos in a folder, i.e, "lib"<br />
<br />
Create a folder in "lib": /core-build/release<br />

Add these packages and libraries in core/CMakeLists.txt file after SeqAn package.
```
find_package( Threads )
link_libraries( ${CMAKE_THREAD_LIBS_INIT} )
find_library(LIB_ARROW arrow)
```

Building:
```
cmake ../../core    -DCMAKE_PREFIX_PATH="$HOME/lib/seqan/util/cmake"    -DSEQAN_INCLUDE_PATH="$/lib/seqan/include"
make
```

Running:
```
./core /path/to/file/sam(bam)file.sam(.bam)
```
