#!/bin/bash -l

# macOS version of compile.sh
# ulong is defined by glibc on Linux but not on macOS, so it is mapped to macOS's u_long here

cd src
src_files=''
for entry in *.cpp
do
  src_files="$src_files"" ""$entry"
done

command="clang++ -O3 -std=c++11 -pthread -w -include sys/types.h -Dulong=u_long -I.$src_files -o ../MC-DC_Simulator"


eval "$command"
