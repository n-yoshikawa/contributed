Instructions on use:

```
./obfragment platinum.sdf > /dev/null 2> frequencies.txt
grep INDEX frequencies.txt | awk {'print $3'} | sort | uniq > fragments.txt
