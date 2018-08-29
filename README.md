# core

Put together "seqan" and "core" in a folder, i.e, "lib"<br />
<br />
Create a folder in "lib": /core-build/release<br />

Building
```
-cmake ../../core    -DCMAKE_PREFIX_PATH="$HOME/lib/seqan/util/cmake"    -DSEQAN_INCLUDE_PATH="$/lib/seqan/include"<br />
-make<br />
```

Running
./core /path/to/file/sam(bam)file.sam(.bam)<br />
