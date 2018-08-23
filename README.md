# Perfect-Hash-Function
Dynamic Perfect Hash Function (dphf ) generate a perfect hash function object according to an user provided array.

Dphf does something similar to GNU Perf:  https://www.gnu.org/software/gperf/

Actually majority of the ideas to generate the object is borrowed from GNU Perf

But dphf is better in the following area:

   1.  When GNU Perf generates C code (or other language),
       dphf generates a perfect hash object on the fly
       according to user input.

   2. Dphf is written using modern C++ while GNU Perf was written
      using old C++ (without using STL).

   3. GNU Perf contains over 10,000 lines of code, dphf contains
      less than 1,000 lines of code


Dphf is a template class which can only accept an class that derives from
dphf_hook.

The constructor of dphf accept a pointer to a vector of template class.

dphf_hook contain 2 pure virtual functions:
   1. compare()
   2. get_key()


One of the use case of dphf is parsing "key-value" pairs. For instance, you need to parse a huge amount of machine generated "key-value" text data. And valid keys are known in advance, you can use dphf to significantly improve the "key" look-up performance.


In order to use dphf

1. include "dphf.hpp"
2. define a class derived from dphf_hook
3. populate a vector of your class object (defined in step 2)
4. construct a dphf object using the vector (created in step 3)
5. using the object to find the desired item.


#include "dphf.hpp"

class word_count: public charles_zhang::dphf_hook {
public:

	...
	virtual std::tuple<const char*, int> get_key() const{
        ...
	}

	virtual bool compare(const char* search_key, int search_key_len) const {
       ...
	}

	virtual ~word_count() {}
...
};


void test() {
	std::vector<word_count> words;
    ...
    charles_zhang::dphf<word_count> word_phf(&words);
    
    ...
	// give the key, locate the item using word_phf
    auto item = word_phf.get_item("system", 6);
	...
}
