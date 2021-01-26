#!/bin/bash

# gtdbtk data formatting

# input file path
# output file path

# double filtered
cut -f 1,2 $1 | grep -v "user_genome" | gsed -E "s/(.+)\t(d__.*);(p__.*);(c__.*);(o__.*);(f__.*);(g__.*);(s__.*)/\1\t\2\t\3\t\4\t\5\t\6\t\7\t\8/" > $2/gtdbtk.tsv