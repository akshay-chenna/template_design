tr ' ' '\n' < clndx30bb.ndx >> merged30bb.ndx
sed '/^$/d' merged30bb.ndx >> mergedline30bb.txt
csplit -k mergedline30bb.txt /Cluster/ '{*}'
#paste xx{01..60} >> clustrows30bb.txt
