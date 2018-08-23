/*
Generate perfect hash function on demand.

Borrowed ideas from GNU gperf ( https://www.gnu.org/software/gperf/ ).

Gperf generates C code (or code in other languages) and your need to compile
the generated code into your applications.

On the other hand, this dynamic perfect hash function (dphf), will generate
a "perfect hash function" object according to the input vector.

see test-dphf.cpp for an example.

Licensed under the MIT License <http://opensource.org/licenses/MIT>.
Copyright (c) 2018 Charles Zhang (charleszhangsse@gmail.com)

Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this software and associated  documentation files (the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell
copies  of  the Software,  and  to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR
IMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include <set>
#include <unordered_set>
#include <vector>
#include <functional>
#include <tuple>
#include <assert.h>
#include <limits>
#include <algorithm>


namespace charles_zhang
{

using IntSet = std::set<int>;
using IntMultiSet = std::multiset<unsigned int>;

static inline unsigned int ToPowerOf2(unsigned int val) {
	int leading = __builtin_clz(val);
	unsigned int ret = (1<< (32-leading));
	if ( ret == (2 * val) )
		ret = ret>>1;
	return ret;
}


struct EqualClass {

	std::vector<unsigned int> keyIndexes;

	/*
	 * The undetermined selected characters for the keywords in this
	 *  equivalence class, as a canonically reordered multiset.
	 */
	IntMultiSet undeterminedChars;
};


struct Step {
	std::vector<unsigned int> changingChars;

	std::vector<unsigned char> undetermined;

	std::vector<EqualClass> partition;

	unsigned int maxAssoVal;
	unsigned int valMask;
};


class dphf_hook {
public:
	virtual std::tuple<const char*, int> get_key() const = 0;
	virtual bool compare(const char* search_key, int search_key_len) const = 0;
	virtual ~dphf_hook() {}
};



// item_type must implement the methods defined in phf_hook
template <class item_type>
class dphf {
public:

	dphf(std::vector<item_type>* items) {
		items_ = items;
		for ( auto item : *items_ ) {
			const char* key;
			int len;
			std::tie(key, len ) = item.get_key();
			if ( maxKeyLen_ < len )
				maxKeyLen_ = len;
		}
		Prepare();
	}

	item_type* get_item(const char* key, int len) const {
		unsigned int hash = hashFun_(this, key, len);
		if ( hash <= maxHashVal_ ) {
			auto item = lookUp_[hash];
			if ( item == nullptr )
				return nullptr;
			else {
				if ( item->compare(key, len ) )
					return item;
				else
					return nullptr;
			}
		}
		else
			return nullptr;
	}

private:
	std::vector<item_type>* items_;

	unsigned int maxKeyLen_ { 0 };


	/*
	 * one vector<int> element for each keyword.
	 * Each vector<int> contains the "character" used for
	 * calculating hash value for the keyword
	 */
	std::vector<std::vector<unsigned int>> selCharsVec_;

	IntSet posSet_; // position used to calculate hash values
	std::vector<unsigned int> positions_;
	bool includeLastCol_{false};

	unsigned int targetDupCount_{0};

	std::vector<int> alphaInc_;

	/*
	 * for each keyword, we have one vector<unsigned int> object.
	 * In each vector<unsigned int> object, we have alphaSize_ of unsigned int.
	 * vector<unsigned int> object is the occurrence of each used characters
	 * that take part in hash calculation
	 */
	std::vector< std::vector<unsigned int> > occurrences_;

	// a vector that contains characters occurs in the keywords
	std::vector<unsigned int> usedChars_;

	// the values associated with each characters.
	std::vector<unsigned int> assoVals_;

	unsigned int maxAssoVal_  { 0 };
	unsigned int maxHashVal_  { 0 };

	std::vector<Step> steps_;
	std::vector<item_type*> lookUp_;
	std::function<unsigned int(const dphf*, const char*, int)> hashFun_;
	//function<unsigned int(const char*, int)> hashFunc_;


	void Prepare() {
		SetPositions();
		SetAlphaInc();
		selCharsVec_ = GetSelCharsVec(posSet_, &alphaInc_);
		SetAssoVals();
		PopulateLookup();
		PopulateHashFunc();

		//clear objects we do not need any more
		occurrences_.clear();
		occurrences_.shrink_to_fit();

		selCharsVec_.clear();
		selCharsVec_.shrink_to_fit();

		steps_.clear();
		steps_.shrink_to_fit();
	}


	/*
	 * If 2 keywords have the same length and differ in just one position,
	 * and it is not the last character this position is mandatory.
	 */
	IntSet GetMandatoryPositions() {
		IntSet mandatory;
		for ( unsigned int i = 0; i < items_->size() - 1; i++ ) {
			for ( unsigned int j = i + 1; j < items_->size(); j++) {
				int pos = GetDiffPosition((*items_)[i], (*items_)[j]);
				if ( pos != -1 )
					mandatory.insert(pos);
			}
		}
		return mandatory;

	}


	void SetPositions() {
		auto mandatory = GetMandatoryPositions();
		posSet_ = mandatory;
		AddPositions();
		RemovePositions(mandatory);
		Replace2with1(mandatory);

		for ( auto pos : posSet_ ) {
			if ( pos == LastColumn )
				includeLastCol_ = true;
			else
				positions_.push_back(pos);
		}
	}


	void SetAssoVals() {
		SetUsedChars();
		SetOccurrences();
		maxAssoVal_ = ToPowerOf2(items_->size());
		maxHashVal_ = ( maxAssoVal_ - 1 ) * posSet_.size();
		CreateSteps();
		PopulateAssoVals();
		if ( HasCollision() )
			//TODO: throw a proper exception
			throw "incorrect internal logic";
		PopulateMaxHashVal();
	}


	void PopulateLookup() {
		lookUp_.assign(maxHashVal_ + 1, nullptr);
		int index = 0;
		for (const auto& selChars : selCharsVec_ ) {
			unsigned int hash = CalHashVal(selChars);
			lookUp_[hash] = &((*items_)[index]);
			index++;
		}
	}


	/*
	 * If the name in key1 is different from the name in key2 by only one
	 * position and the position is not the last one, then it returns the position,
	 * otherwise, it returns -1;
	 */
	int GetDiffPosition(const item_type& item1, const item_type& item2 ) {
		const char *key1, *key2;
		int        len1,  len2;

		std::tie(key1, len1) = item1.get_key();
		std::tie(key2, len2) = item1.get_key();

		if ( len1 != len2 )
			return -1;

		int ret = -1;
		for (int i = 0; i < len1 - 1 ; i++ ) {
			if ( key1[i] != key2[i] ) {
				if ( ret == -1 )
					ret = i;
				else {
					ret = -1;
					break;
				}
			}
		}

		if ( key1[len1-1] == key1[len2-1] )
			return ret;
		else
			return -1;
	}


	static bool HasKey(const IntSet& positions, unsigned int pos) {
		if ( positions.find(pos) == positions.end() )
			return false;
		else
			return true;
	}


	/*
	 * Add positions as long as this decreases the duplicates count.
	 */
	void AddPositions() {
		targetDupCount_ = CountDuplicate(posSet_);

		while ( 1 ) {
			IntSet best;
			unsigned int bestDupCount = std::numeric_limits<unsigned int>::max();

			for ( int i = maxKeyLen_ - 1; i >= LastColumn; i-- ) {
				if ( HasKey(posSet_, i)  )
					continue;

				IntSet tryPositions = posSet_;
				tryPositions.insert(i);
				unsigned int tryDupCount = CountDuplicate(tryPositions);
				/*
				 * We prefer 'tryPositions' to 'bestPositions'
				 * if it produces less duplicates,
				 * or if it produces the same number of duplicates but with
				 * a more efficient hash function.
				 */
				if (tryDupCount < bestDupCount
						|| ( tryDupCount == bestDupCount && i >= 0) ) {
					best.swap(tryPositions);
					bestDupCount = tryDupCount;
				}
			}

			// Stop adding positions when it gives no improvement.
			if ( bestDupCount >= targetDupCount_ )
				break;

			posSet_.swap(best);
			targetDupCount_ = bestDupCount;
		}
	}


	/*
	 * Remove positions, as long as this doesn't increase
	 * the duplicates count.
	 */
	void RemovePositions(const IntSet& mandatoryPositions) {
		while ( 1 ) {
			IntSet best;
			unsigned int bestDupCount = std::numeric_limits<unsigned int>::max();

			for ( auto iter = posSet_.rbegin(); iter != posSet_.rend(); ++iter) {
				if ( HasKey(mandatoryPositions, *iter) )
					continue;
				IntSet tryPositions = posSet_;
				tryPositions.erase(*iter);
				unsigned int tryDupCount = CountDuplicate(tryPositions);

				/*
				 * We prefer 'try' to 'best' if it produces the same number
				 * of duplicates but with a more efficient hash function.
				 */
				if ( tryDupCount < bestDupCount ||
						( tryDupCount == bestDupCount && *iter == -1) ) {
					best.swap(tryPositions);
					bestDupCount = tryDupCount;
				}
			}

			// Stop removing positions when it gives no improvement.
			if ( bestDupCount > targetDupCount_)
				break;

			posSet_.swap(best);
			targetDupCount_ = bestDupCount;
		}
	}


	/*
	 * Replace two positions by one, as long as this doesn't increase the
	 duplicates count.
	 */
	void Replace2with1(const IntSet& mandatoryPositions) {
		while ( 1 ) {
			IntSet best;
			unsigned int bestDupCount = std::numeric_limits<unsigned int>::max();

			for ( auto iter = posSet_.rbegin(); iter != posSet_.rend(); ++iter) {
				if ( HasKey(mandatoryPositions, *iter) )
					continue;

				for ( auto iter2 = posSet_.rbegin(); iter2 != posSet_.rend(); ++iter2) {
					if ( HasKey(mandatoryPositions, *iter) || iter == iter2 )
						continue;

					for (int k = maxKeyLen_ - 1; k >= 0; k--) {
						if ( HasKey(posSet_, k) )
							continue;

						IntSet tryPositions = posSet_;
						tryPositions.erase(*iter);
						tryPositions.erase(*iter2);
						tryPositions.insert(k);
						unsigned int tryDupCount = CountDuplicate(tryPositions);

						/*
						 * We prefer 'try' to 'best' if it produces less
						 * duplicates, or if it produces the same number of
						 * duplicates but with a more efficient hash function.
						 */
						if (  tryDupCount < bestDupCount ||
								( tryDupCount == bestDupCount &&
									( *iter == -1 || *iter2 == -1 || k >= 0) ) ) {
							best.swap(tryPositions);
							bestDupCount = tryDupCount;
						}
					}
				}
			}

			// Stop replacing positions when it gives no improvement.
			if ( bestDupCount > targetDupCount_)
				break;

			posSet_.swap(best);
			targetDupCount_ = bestDupCount;
		}
	}


	int CountDuplicate(const IntSet& positions) {
		if ( positions.size() == 0 )
			return items_->size();

		std::vector< std::vector<unsigned int> > selChars = GetSelCharsVec(positions, nullptr);
		int curDupCount = 0;
		std::set< std::vector<unsigned int> > countSet;

		for ( auto& chars : selChars ) {
			auto ret = countSet.insert(chars);
			if ( !ret.second )
				curDupCount++;
		}
		return curDupCount;
	}


	int CountDuplicateMultiset(const std::vector<int>& inc) {
		if ( posSet_.size() == 0  )
			return items_->size();

		std::vector< std::vector<unsigned int> > selCharsVec = GetSelCharsVec(posSet_, &inc);
		int dupCount = 0;
		std::set< IntMultiSet > countSet;

		for ( auto& selChars : selCharsVec ) {
			IntMultiSet temp;
			for ( unsigned int i = 0; i < selChars.size(); i++ )
				temp.insert(selChars[i]);
			auto ret = countSet.insert(temp);
			if ( !ret.second )
				dupCount++;
		}
		return dupCount;
	}


	void SetAlphaInc() {
		alphaInc_.assign(maxKeyLen_, 0);
		unsigned int curDupCount = CountDuplicateMultiset(alphaInc_);

		if ( targetDupCount_ >= curDupCount )
			return;

		std::vector<int> best;
		do {
			/*
			 * An increment of 1 is not always enough.
			 * Try higher increments also.
			 */
			for ( int inc = 1; ; inc++ ) {
				unsigned int bestDupCount = std::numeric_limits<unsigned int>::max();
				for ( auto pos : posSet_ ) {
					if ( pos != LastColumn ) {
						std::vector<int> tryInc = alphaInc_;
						tryInc[pos] += inc;
						unsigned int tryDupCount = CountDuplicateMultiset(tryInc);
						if ( tryDupCount < bestDupCount ) {
							bestDupCount = tryDupCount;
							best.swap(tryInc);
						}
					}
				}

				// Stop this round when we got an improvement.
				if ( bestDupCount < curDupCount ) {
					alphaInc_.swap(best);
					curDupCount = bestDupCount;
					break;
				}
			}
		} while ( curDupCount > targetDupCount_ );
	}


	void SetUsedChars() {
		unsigned int counts[alphaSize_];
		memset(counts, 0, sizeof(counts));

		for ( auto& selChars: selCharsVec_ )
			for ( auto c : selChars )
				counts[c]++;

		for (unsigned int i = 0; i < alphaSize_; i++ )
			if ( counts[i] )
				usedChars_.push_back(i);
	}


	void SetOccurrences() {
		occurrences_.reserve(items_->size());
		for ( auto& selChars : selCharsVec_ ) {
			std::vector<unsigned int> temp;
			temp.assign(alphaSize_, 0);
			for ( auto c : selChars )
				temp[c]++;
			occurrences_.push_back(temp);
		}
	}


	unsigned int CalCollisions(const std::vector<EqualClass>& ecs, unsigned int c) {
		unsigned int sum = 0;
		unsigned int splitCardinalities[posSet_.size() + 1];
		for ( auto& ec : ecs ) {
			memset(splitCardinalities, 0, sizeof(splitCardinalities));
			for ( auto index : ec.keyIndexes )
				splitCardinalities[occurrences_[index][c]]++;

			sum += ec.keyIndexes.size() * ec.keyIndexes.size();
			for ( unsigned int i = 0; i <= posSet_.size(); i++ )
				sum -= splitCardinalities[i] * splitCardinalities[i];
		}
		return sum;
	}


	std::vector<EqualClass> ComputePartition(
			const std::vector<unsigned char>& undetermined) const {
		std::vector<EqualClass> partition;
		IntMultiSet undeterminedChars;
		unsigned int index = 0;

		for ( auto& selChars : selCharsVec_ ) {
			undeterminedChars.clear();
			for ( auto c : selChars )
				if ( undetermined[c] )
					undeterminedChars.insert(c);

			EqualClass* ec = nullptr;

			for ( auto& item : partition ) {
				if ( item.undeterminedChars == undeterminedChars ) {
					ec = &item;
					break;
				}
			}

			if ( !ec ) {
				EqualClass newEc;
				newEc.undeterminedChars.swap(undeterminedChars);
				newEc.keyIndexes.push_back(index);
				partition.push_back(newEc);
			}
			else
				ec->keyIndexes.push_back(index);
			index++;
		}

		return partition;
	}


	// Check whether adding c to "determined" characters changing
	// the partition
	bool IsPartitionSame(const std::vector<EqualClass>& ecs, unsigned int c) {
		for ( auto& ec: ecs ) {
			unsigned int firstCount = occurrences_[ec.keyIndexes[0]][c];
			for (unsigned int i = 1; i < ec.keyIndexes.size(); i++ )
				if ( occurrences_[ec.keyIndexes[i]][c] != firstCount )
					return false;
		}
		return true;
	}


	void CreateSteps() {
		std::vector<unsigned char> undetermined;
		undetermined.assign(alphaSize_, false);

		while ( 1 ) {
			std::vector<EqualClass> partition = ComputePartition(undetermined);
			unsigned int best = 0;
			unsigned int bestCollisions = std::numeric_limits<unsigned int>::max();

			for ( auto c : usedChars_ ) {
				if ( undetermined[c] )
					continue;
				unsigned int collisions = CalCollisions(partition, c);
				if ( collisions < bestCollisions ) {
					best = c;
					bestCollisions = collisions;
				}
			}

			/*
			 * All used characters are determined.  We are
			 *  at the starting situation and don't need any more step.
			 */
			if ( bestCollisions == std::numeric_limits<unsigned int>::max() )
				break;

			Step step;
			step.undetermined = undetermined;
			step.partition.swap(partition);

			undetermined[best] = true;
			step.changingChars.push_back(best);
			partition = ComputePartition(undetermined);

			/*
			 * If adding new characters do not change the partition, we add them
			 * in the same step
			 */
			for ( auto c : usedChars_) {
				if ( !undetermined[c] && IsPartitionSame(partition, c) ) {
					undetermined[c] = true;
					step.changingChars.push_back(c);
				}
			}
			step.maxAssoVal = maxAssoVal_;
			step.valMask = maxAssoVal_ - 1;
			steps_.push_back(step);
		}
		std::reverse(steps_.begin(), steps_.end());
	}


	void PopulateAssoVals() {
		assoVals_.assign(alphaSize_, 0);
		for (Step& step : steps_ ) {
			if ( step.changingChars.size() > 1 )
				PopulateAssoValuesRandom(step);
			else
				PopulateAssoValuesJump(step);
		}
	}


	void PopulateAssoValuesRandom(const Step& step) {
		srand(time(NULL));
		bool looping = true;
		while ( looping ) {
			for ( auto c : step.changingChars ) {
				assoVals_[c] = ( assoVals_[c] + rand() ) & step.valMask;
				if ( IsAssoValOk(step) ) {
					looping = false;
					break;
				}
			}
		}
	}


	void PopulateAssoValuesJump(Step& step) {
		assert(step.changingChars.size() == 1 );
		unsigned int bound = 0;
		unsigned int c = step.changingChars[0];
		do {
			if ( bound == step.maxAssoVal ) {
				step.maxAssoVal *= 2;
				step.valMask = step.maxAssoVal - 1;
				if ( step.maxAssoVal > maxAssoVal_ )
					maxAssoVal_  = step.maxAssoVal;
			}
			assoVals_[c] = bound & step.valMask;
			bound++;
		} while ( !IsAssoValOk(step) );
	}


	void PopulateMaxHashVal() {
		for (const auto& selChars : selCharsVec_ ) {
			unsigned int hash = CalHashVal(selChars);
			if ( maxHashVal_ < hash )
				maxHashVal_ = hash;
		}
	}


	bool IsAssoValOk(const Step& step) {
		for ( auto& cls : step.partition ) {
			if ( cls.keyIndexes.size() <= 1 )
				continue;
			std::unordered_set<unsigned int> hashSet;
			for ( auto index : cls.keyIndexes ) {
				auto& selChars = selCharsVec_[index];
				unsigned int hashVal = 0;
				for ( auto c : selChars ) {
					if ( !step.undetermined[c] )
						hashVal += assoVals_[c];
				}
				auto res = hashSet.insert(hashVal);
				if ( !res.second )
					return false;
			}
		}
		return true;
	}


	unsigned int CalHashVal(const std::vector<unsigned int>& selChars) {
		unsigned int hashVal = 0;
		for ( auto c : selChars )
			hashVal += assoVals_[c];
		return hashVal;
	}


	bool HasCollision() {
		std::unordered_set<unsigned int> hashSet;
		for (const auto& selChars : selCharsVec_ ) {
			unsigned int hash = CalHashVal(selChars);
			auto res = hashSet.insert(hash);
			if ( !res.second )
				return true;
		}
		return false;
	}


	std::vector<std::vector<unsigned int>> GetSelCharsVec(const IntSet& positions,
			const std::vector<int>* alphaInc) {
		std::vector< std::vector<unsigned int> > selChars;
		selChars.reserve(items_->size());

		for ( auto item : *items_ ) {
			std::vector<unsigned int> chars;
			const char* key;
			int key_len;
			std::tie(key, key_len) = item.get_key();
			for ( auto pos : positions ) {
				if ( pos == LastColumn )
					chars.push_back(key[key_len - 1]);
				else if ( pos < key_len - 1 ) {
					if ( alphaInc )
						chars.push_back(key[pos] + (*alphaInc)[pos]);
					else
						chars.push_back(key[pos]);
				}
			}
			selChars.push_back(chars);
		}
		return selChars;

	}


	void PopulateHashFunc() {
		/*IntSet pos_124 {-1, 2, 4};
		IntSet pos136  { 1, 3, 6};

		if ( positionSet_ == pos_124 )
			hashFun_= &Phf::CalHash_124Column;
		else if ( positionSet_ == pos136 )
			hashFun_= &Phf::CalHash136Column;
		else*/
			hashFun_= &dphf::CalHashGeneric;
	}


	unsigned int CalHashGeneric(const char* key, int len) const {
		unsigned int hash;
		int assoIndex;

		if ( includeLastCol_ ) {
			int pos = len - 1;
			assoIndex = key[pos];
			hash = assoVals_[assoIndex];
		}
		else
			hash = 0;

		for ( auto pos : posSet_ ) {
			if ( pos >= len || pos < 0 )
				continue;
			int assoIndex = key[pos] + alphaInc_[pos];
			hash += assoVals_[assoIndex];
		}
		return hash;
	}


	// when positions is -1, 2, 4, we use this function
	unsigned int CalHash_124Column(const char* key, int len) const {
		if ( len > 4  )
			return assoVals_[key[len-1]] +
					assoVals_[key[4] + alphaInc_[4]] +
					assoVals_[key[2] + alphaInc_[2]];
		else {
			switch ( len ) {
			case 4:
				return assoVals_[key[3]] +
							assoVals_[key[2] + alphaInc_[2]];
			case 3:
				return assoVals_[key[2] ];
			default:
				return assoVals_[key[1] + alphaInc_[1]];
			}
		}

	}


	// when positions is 1, 3, 6, we use this function
	unsigned int CalHash136Column(const char* key, int len) const {
		if ( len > 6  )
			return assoVals_[key[6] + alphaInc_[6]] +
					assoVals_[key[3] + alphaInc_[3]] +
					assoVals_[key[1] + alphaInc_[1]];
		else {
			if ( len > 3 )
				return assoVals_[key[3] + alphaInc_[3]]
								 + assoVals_[key[1] + alphaInc_[1]];
			else
				return assoVals_[key[1] + alphaInc_[1]];
		}
	}

	static const int LastColumn = -1;
	static const unsigned int alphaSize_ = 257;

};

}

