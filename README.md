# Perfect-Hash-Function
A perfect hash function generated at run time

This does something similar to GNU Perf:  https://www.gnu.org/software/gperf/

But it has a key difference. When GNU Perf generates C code (or other language), which you need to do at compiling time,
my Perfect Hash Function (PHF) generates the hash function at run time.

And PHF is written using modern C++ while GNU Perf was written using old C++ (without using STL).

When GNU Perf contains over 10,000 lines of code, this PHF implementation contains about 1,000 lines of code
