# core

Put together "seqan" and "core" in a folder, i.e, "lib"

Create a folder in "lib": /core-build/release
cmake ../../core    -DCMAKE_PREFIX_PATH="$HOME/lib/seqan/util/cmake"    -DSEQAN_INCLUDE_PATH="$/lib/seqan/include"﻿
make
./core /path/to/file/sam(bam)file.sam(.bam)
