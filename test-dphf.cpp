/*
 * test-phf.cpp
 *
 *  Created on: Aug 16, 2018
 *      Author: charleszhang
 */


#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "dphf.hpp"

class word_count: public charles_zhang::dphf_hook {
public:

	word_count(const char* name) : name_(name), len_(strlen(name)) {}

	virtual std::tuple<const char*, int> get_key() const{
		return std::make_tuple(name_, len_);
	}

	virtual bool compare(const char* search_key, int search_key_len) const {
		return len_ == search_key_len &&  !memcmp(name_, search_key, len_);
	}

	virtual ~word_count() {}

	void increase_count() { count_++;}

private:
	const char* name_;
	int len_;
	int count_ {0};
};


static void do_test(const char* file_name, charles_zhang::dphf<word_count>& word_phf);

int main(int argc, char* argv[]) {

	if ( argc != 2 ) {
		printf("%s <filename>\n", argv[0]);
		return 0;
	}

	std::vector<word_count> words = {
			"speeds",
			"ARM",
			"Intel",
			"AMD",
			"CPU",
			"RAM",
			"encryption",
			"system",
			"local",
			"The",
			"from",
			"this"
	};

	charles_zhang::dphf<word_count> word_phf(&words);

	do_test(argv[1], word_phf);

	printf("ok\n");
}


static void do_test(const char* file_name,
		charles_zhang::dphf<word_count>& word_phf) {

	std::ifstream infile(file_name);
	std::string line;
	while ( std::getline(infile, line) ) {
		std::stringstream ss(line);
		std::string word;
		while ( std::getline(ss, word, ' ') )  {
			auto item = word_phf.get_item(word.c_str(), word.length());
			if ( item )
				item->increase_count();
		}
	}
}
