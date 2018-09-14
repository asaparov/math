/**
 * \file features.h
 *
 * This file contains structures for representing sets of discrete features,
 * where each feature_set has `n` features, and each feature is of type
 * `unsigned int`.
 *
 *  <!-- Created on: Dec 1, 2016
 *           Author: asaparov -->
 */

#ifndef FEATURES_H_
#define FEATURES_H_

#include <core/array.h>

using namespace core;

/**
 * This struct represents a sequence of features, with length
 * feature_set::feature_count, each feature having type `unsigned int`.
 */
struct feature_set {
	/**
	 * The native array that stores the feature values.
	 */
	unsigned int* features;

	/**
	 * The number of features and the length of feature_set::features,
	 * feature_set::excluded_counts, and the first dimension of
	 * feature_set::excluded.
	 */
	unsigned int feature_count;

	/**
	 * A two-dimensional array containing *excluded* feature values. That is,
	 * `excluded[i]` contains the set of excluded feature values for the
	 * feature at index `i`. Each `excluded[i]` is a sorted array of
	 * `unsigned int` feature values containing distinct elements. The length
	 * of `excluded[i]` is given by `feature_set::excluded_counts[i]`. For each
	 * `j` such that `feature_set::excluded_counts[j] == 0`, `excluded[j]` is
	 * `NULL`.
	 */
	unsigned int** excluded;

	/**
	 * An array that stores the length of each `feature_set::excluded[i]`.
	 */
	unsigned int* excluded_counts;

	/**
	 * Initializes this feature set with the given `feature_count`.
	 * Each element in feature_set::features is uninitialized, whereas
	 * feature_set::excluded_counts is initialized to all zeros.
	 */
	feature_set(unsigned int feature_count) : feature_count(feature_count) {
		if (!initialize()) exit(EXIT_FAILURE);
	}

	/**
	 * Initializes this feature set by copying from the given feature_set `src`.
	 */
	feature_set(const feature_set& src) : feature_count(src.feature_count) {
		if (!initialize(src)) exit(EXIT_FAILURE);
	}

	~feature_set() { free(); }

	/**
	 * Returns the feature value at the given `index`.
	 */
	inline unsigned int operator [] (unsigned int index) const {
		return features[index];
	}

	/**
	 * Returns the feature value at the given `index`.
	 */
	inline unsigned int& operator [] (unsigned int index) {
		return features[index];
	}

	/**
	 * Sets the feature value at the given `index`.
	 */
	inline void set_feature(unsigned int index, unsigned int value) {
		features[index] = value;
	}

	/**
	 * Returns `true` if the given feature `value` is excluded from the feature
	 * at the given `index`. Otherwise, this function returns `false`.
	 */
	inline bool is_excluded(unsigned int index, unsigned int value) const {
		return index_of(value, excluded[index], excluded_counts[index]) < excluded_counts[index];
	}

	/**
	 * Checks that the excluded array at `feature_set::excluded[index]` has
	 * sufficient `capacity`. If that array is `NULL`, this function will
	 * initialize it (leaving its elements uninitialized).
	 */
	bool ensure_excluded_capacity(unsigned int index, unsigned int capacity) {
		if (excluded_counts[index] == 0) excluded[index] = NULL;
		unsigned int* new_excluded = (unsigned int*) realloc(
			excluded[index], sizeof(unsigned int) * capacity);
		if (new_excluded == NULL) {
			fprintf(stderr, "feature_list.excluded_unsorted ERROR: Out of memory.\n");
			return false;
		}
		excluded[index] = new_excluded;
		return true;
	}

	/**
	 * Initializes the excluded array at `feature_set::excluded[index]` and
	 * copies the contents from `src`. This function assumes
	 * `feature_set::excluded[index]` is previously uninitialized.
	 */
	bool set_excluded(unsigned int index, const unsigned int* src, unsigned int count) {
		if (count == 0) return true;
		excluded[index] = (unsigned int*) malloc(sizeof(unsigned int) * count);
		if (excluded[index] == NULL) {
			fprintf(stderr, "feature_set.set_excluded ERROR: Out of memory.\n");
			return false;
		}
		excluded_counts[index] = count;
		memcpy(excluded[index], src, sizeof(unsigned int) * count);
		return true;
	}

	/**
	 * Excludes the given `item` from the feature at the given `index`, without
	 * sorting the excluded array `feature_set::excluded[index]`.
	 * `feature_set::excluded_counts[index]` is also incremented. This function
	 * assumes the excluded array does not already contain `item` and has
	 * sufficient capacity.
	 */
	void exclude_unsorted(unsigned int index, unsigned int item) {
		excluded[index][excluded_counts[index]] = item;
		excluded_counts[index]++;
	}

	/**
	 * Sorts the excluded array at `feature_set::excluded[index]`. This
	 * function assumes the excluded array is not `NULL`.
	 */
	void sort_excluded(unsigned int index) {
		if (excluded_counts[index] == 0)
			core::free(excluded[index]);
		else insertion_sort(excluded[index], excluded_counts[index]);
	}

	/**
	 * Evaluates the hash function for the given feature_set `set`.
	 */
	static inline unsigned int hash(const feature_set& set) {
		unsigned int hash = default_hash(set.features, set.feature_count);
		for (unsigned int i = 0; i < set.feature_count; i++)
			if (set.excluded_counts[i] > 0)
				hash ^= default_hash(set.excluded[i], set.excluded_counts[i]);
		return hash;
	}

	/**
	 * Moves the feature_set in `src` to `dst`. Note that this function does
	 * not initialize the fields in `dst` and copy the contents from the
	 * corresponding fields in `src` into `dst`. Rather, this function simply
	 * copies the pointers.
	 */
	static inline void move(const feature_set& src, feature_set& dst) {
		dst.features = src.features;
		dst.feature_count = src.feature_count;
		dst.excluded = src.excluded;
		dst.excluded_counts = src.excluded_counts;
	}

	/**
	 * Swaps the underlying pointers in the feature_sets `src` and `dst`.
	 */
	static inline void swap(feature_set& first, feature_set& second) {
		core::swap(first.features, second.features);
		core::swap(first.feature_count, second.feature_count);
		core::swap(first.excluded, second.excluded);
		core::swap(first.excluded_counts, second.excluded_counts);
	}

	/**
	 * Returns whether feature_set::features is `NULL` in the given `set`.
	 */
	static inline bool is_empty(const feature_set& set) {
		return set.features == NULL;
	}

	/**
	 * Sets feature_set::features to `NULL` in the given `set`.
	 */
	static inline void set_empty(feature_set& set) {
		set.features = NULL;
	}

	/**
	 * Sets feature_set::features to `NULL` in every element of the given array `sets`.
	 */
	static inline void set_empty(feature_set* sets, unsigned int length) {
		for (unsigned int i = 0; i < length; i++)
			sets[i].features = NULL;
	}

	/**
	 * Frees the underlying arrays in the given feature_set `set`.
	 */
	static inline void free(feature_set& set) {
		set.free();
	}

private:
	inline bool initialize() {
		features = (unsigned int*) malloc(sizeof(unsigned int) * feature_count);
		if (features == NULL) {
			fprintf(stderr, "feature_set.initialize ERROR: Insufficient memory for feature vector.\n");
			return false;
		}
		excluded = (unsigned int**) malloc(sizeof(unsigned int*) * feature_count);
		if (excluded == NULL) {
			fprintf(stderr, "feature_set.initialize ERROR: Insufficient memory for excluded sets.\n");
			core::free(features);
			return false;
		}
		excluded_counts = (unsigned int*) calloc(sizeof(unsigned int), feature_count);
		if (excluded_counts == NULL) {
			fprintf(stderr, "feature_set.initialize ERROR: Insufficient memory for excluded count array.\n");
			core::free(features); core::free(excluded);
			return false;
		}
		return true;
	}

	inline bool initialize(const feature_set& src) {
		if (!initialize()) return false;
		memcpy(features, src.features, sizeof(unsigned int) * feature_count);
		for (unsigned int i = 0; i < src.feature_count; i++) {
			if (src.excluded_counts[i] > 0) {
				excluded[i] = (unsigned int*) malloc(sizeof(unsigned int) * src.excluded_counts[i]);
				if (excluded[i] == NULL) {
					fprintf(stderr, "feature_set.initialize ERROR: Insufficient memory for excluded set at index %u.\n", i);
					free(); return false;
				}
				memcpy(excluded[i], src.excluded[i], sizeof(unsigned int) * src.excluded_counts[i]);
				excluded_counts[i] = src.excluded_counts[i];
			}
		}
		return true;
	}

	inline void free() {
		core::free(features);
		for (unsigned int i = 0; i < feature_count; i++) {
			if (excluded_counts[i] > 0)
				core::free(excluded[i]);
		}
		core::free(excluded);
		core::free(excluded_counts);
	}

	friend bool init(feature_set&, unsigned int);
	friend bool init(feature_set&, const feature_set&);
};

/**
 * Initializes the given feature_set `set` with the given `feature_count`. Each
 * element in feature_set::features is uninitialized, whereas
 * feature_set::excluded_counts is initialized to all zeros.
 */
inline bool init(feature_set& set, unsigned int feature_count) {
	set.feature_count = feature_count;
	return set.initialize();
}

/**
 * Initializes the given feature_set `set` by copying from the given feature_set `src`.
 */
inline bool init(feature_set& set, const feature_set& src) {
	set.feature_count = src.feature_count;
	return set.initialize(src);
}

/**
 * Returns whether the feature_set in `first` is equivalent to the one in
 * `second`. For equivalence, the feature sets must have the same sequence of
 * features, as well as the same set of excluded features at every index.
 */
inline bool operator == (const feature_set& first, const feature_set& second) {
	if (first.features == NULL || first.feature_count != second.feature_count) return false;
	for (unsigned int i = 0; i < first.feature_count; i++)
		if (first.features[i] != second.features[i]) return false;

	for (unsigned int i = 0; i < first.feature_count; i++) {
		if (first.excluded_counts[i] != second.excluded_counts[i])
			return false;
		if (first.excluded_counts[i] > 0) {
			for (unsigned int j = 0; j < first.excluded_counts[i]; j++)
				if (first.excluded[i][j] != second.excluded[i][j])
					return false;
		}
	}
	return true;
}

/**
 * A sorter for feature sets that enables lexicographical sorting according to
 * the sequence of feature values.
 */
struct feature_set_sorter {
	unsigned int depth;

	/**
	 * Initializes the feature_set_sorter with the given `depth`. This
	 * feature_set_sorter sorts feature set objects that have length `depth - 1`.
	 */
	feature_set_sorter(unsigned int depth) : depth(depth) { }
};

/**
 * Compares the given feature set objects `first` and `second` using the
 * feature_set_sorter `sorter`.
 * \tparam FeatureSet a feature set type that implements the function
 * 		`T get_feature(unsigned int)` where `T` is any type that satisfies
 * 		[LessThanComparable](https://en.cppreference.com/w/cpp/named_req/LessThanComparable).
 */
template<typename FeatureSet>
inline bool less_than(
	const FeatureSet& first,
	const FeatureSet& second,
	const feature_set_sorter& sorter)
{
	for (unsigned int i = 0; i < sorter.depth - 1; i++) {
		if (first.get_feature(i) < second.get_feature(i))
			return true;
		else if (second.get_feature(i) < first.get_feature(i))
			return false;
	}
	return false;
}

/**
 * This structure represents a feature_set weighted by a numerical value, such
 * as a probability.
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
struct weighted_feature_set {
	/**
	 * The underlying feature_set.
	 */
	feature_set features;

	/**
	 * The weight of weighted_feature_set::features.
	 */
	V log_probability;

	/**
	 * Returns the feature value at the given `index`.
	 */
	inline unsigned int get_feature(unsigned int index) const {
		return features[index];
	}

	/**
	 * Sets the feature value at the given `index` to `feature`.
	 */
	inline void set_feature(unsigned int index, unsigned int feature) {
		features[index] = feature;
	}

	/**
	 * Checks that the excluded array at `feature_set::excluded[index]` in
	 * `weighted_feature_set::features` has sufficient capacity. If that array
	 * is `NULL`, this function will initialize it (leaving its elements
	 * uninitialized).
	 */
	inline bool ensure_excluded_capacity(unsigned int index, unsigned int capacity) {
		return features.ensure_excluded_capacity(index, capacity);
	}

	/**
	 * Initializes the excluded array at `feature_set::excluded[index]` in
	 * `weighted_feature_set::features` and copies the contents from `src`.
	 * This function assumes `feature_set::excluded[index]` is previously
	 * uninitialized.
	 */
	inline bool set_excluded(unsigned int index, const unsigned int* src, unsigned int count) {
		return features.set_excluded(index, src, count);
	}

	/**
	 * Excludes the given item from the feature at the given `index`, without
	 * sorting the excluded array `feature_set::excluded[index]` in
	 * `weighted_feature_set::features`. `feature_set::excluded_counts[index]`
	 * is also incremented. This function assumes the excluded array does not
	 * already contain `item` and has sufficient capacity.
	 */
	inline void exclude_unsorted(unsigned int index, unsigned int item) {
		return features.exclude_unsorted(index, item);
	}

	/**
	 * Sorts the excluded array at `feature_set::excluded[index]` in
	 * `weighted_feature_set::features`. This function assumes the excluded
	 * array is not `NULL`.
	 */
	inline void sort_excluded(unsigned int index) {
		return features.sort_excluded(index);
	}

	/**
	 * Returns the natural exponent of the weight weighted_feature_set::log_probability.
	 */
	inline V get_probability() const {
		return exp(log_probability);
	}

	/**
	 * Sets the weight weighted_feature_set::log_probability to the natural
	 * logarithm of the given `probability`.
	 */
	inline void set_probability(V probability) {
		log_probability = log(probability);
	}

	/**
	 * Moves the given weighted_feature_set in `src` into `dst`.
	 */
	static inline void move(const weighted_feature_set<V>& src, weighted_feature_set<V>& dst) {
		dst.log_probability = src.log_probability;
		feature_set::move(src.features, dst.features);
	}

	/**
	 * Swaps the weighted_feature_set structures in `first` and `second`.
	 */
	static inline void swap(weighted_feature_set<V>& first, weighted_feature_set<V>& second) {
		feature_set::swap(first.features, second.features);
		core::swap(first.log_probability, second.log_probability);
	}

	/**
	 * Frees the given weighted_feature_set `set`.
	 */
	static inline void free(weighted_feature_set<V>& set) {
		core::free(set.features);
	}
};

/**
 * Initializes the given weighted_feature_set `set` with the given
 * `feature_count`. Each element in feature_set::features in
 * weighted_feature_set::features is uninitialized, whereas
 * feature_set::excluded_counts is initialized to all zeros. The weight
 * weighted_feature_set::log_probability is uninitialized.
 */
template<typename V>
inline bool init(weighted_feature_set<V>& set, unsigned int feature_count) {
	return init(set.features, feature_count);
}

/**
 * Initializes the given weighted_feature_set `set` with the given
 * feature_set `src` and weight `log_probability`.
 */
template<typename V>
inline bool init(weighted_feature_set<V>& set, const feature_set& src, V log_probability) {
	set.log_probability = log_probability;
	return init(set.features, src);
}

/**
 * Returns whether weighted_feature_set::log_probability of `first` is *greater than* that of `second`.
 */
template<typename V>
inline bool operator < (const weighted_feature_set<V>& first, const weighted_feature_set<V>& second) {
	return first.log_probability > second.log_probability;
}

/**
 * Returns the weight weighted_feature_set::log_probability of `set`.
 */
template<typename V>
inline double log_probability(const weighted_feature_set<V>& set) {
	return set.log_probability;
}

#endif /* FEATURES_H_ */
