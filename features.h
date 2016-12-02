/**
 * features.h
 *
 *  Created on: Dec 1, 2016
 *      Author: asaparov
 */

#ifndef FEATURES_H_
#define FEATURES_H_

#include <core/array.h>

struct feature_set {
	unsigned int* features;
	unsigned int feature_count;

	/* each is a sorted set of distinct ids */
	unsigned int** excluded;
	unsigned int* excluded_counts;

	feature_set(unsigned int feature_count) : feature_count(feature_count) {
		if (!initialize()) exit(EXIT_FAILURE);
	}

	feature_set(const feature_set& src) : feature_count(src.feature_count) {
		if (!initialize(src)) exit(EXIT_FAILURE);
	}

	~feature_set() { free(); }

	inline unsigned int operator [] (unsigned int index) const {
		return features[index];
	}

	inline unsigned int& operator [] (unsigned int index) {
		return features[index];
	}

	inline void set_feature(unsigned int index, unsigned int value) {
		features[index] = value;
	}

	inline bool is_excluded(unsigned int index, unsigned int value) const {
		return index_of(value, excluded[index], excluded_counts[index]) < excluded_counts[index];
	}

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

	/* NOTE: this function assumes the excluded array at the given index is uninitialized */
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

	/* NOTE: this function assumes the excluded array has sufficient capacity */
	void exclude_unsorted(unsigned int index, unsigned int item) {
		excluded[index][excluded_counts[index]] = item;
		excluded_counts[index]++;
	}

	/* NOTE: this assumes the excluded array is initialized */
	void sort_excluded(unsigned int index) {
		if (excluded_counts[index] == 0)
			core::free(excluded[index]);
		else insertion_sort(excluded[index], excluded_counts[index]);
	}

	static inline unsigned int hash(const feature_set& set) {
		unsigned int hash = default_hash(set.features, set.feature_count);
		for (unsigned int i = 0; i < set.feature_count; i++)
			if (set.excluded_counts[i] > 0)
				hash ^= default_hash(set.excluded[i], set.excluded_counts[i]);
		return hash;
	}

	static inline void move(const feature_set& src, feature_set& dst) {
		dst.features = src.features;
		dst.feature_count = src.feature_count;
		dst.excluded = src.excluded;
		dst.excluded_counts = src.excluded_counts;
	}

	static inline void swap(feature_set& first, feature_set& second) {
		core::swap(first.features, second.features);
		core::swap(first.feature_count, second.feature_count);
		core::swap(first.excluded, second.excluded);
		core::swap(first.excluded_counts, second.excluded_counts);
	}

	static inline bool is_empty(const feature_set& set) {
		return set.features == NULL;
	}

	static inline void set_empty(feature_set& set) {
		set.features = NULL;
	}

	static inline void set_empty(feature_set* sets, unsigned int length) {
		for (unsigned int i = 0; i < length; i++)
			sets[i].features = NULL;
	}

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

inline bool init(feature_set& set, unsigned int feature_count) {
	set.feature_count = feature_count;
	return set.initialize();
}

inline bool init(feature_set& set, const feature_set& src) {
	set.feature_count = src.feature_count;
	return set.initialize(src);
}

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

struct feature_set_sorter {
	unsigned int depth;
	feature_set_sorter(unsigned int depth) : depth(depth) { }
};

template<typename FeatureSet>
inline bool less_than(
	const FeatureSet& first,
	const FeatureSet& second,
	const feature_set_sorter& sorter)
{
	for (unsigned int i = 0; i < sorter.depth - 1; i++) {
		if (first.get_feature(i) < second.get_feature(i))
			return true;
		else if (first.get_feature(i) > second.get_feature(i))
			return false;
	}
	return false;
}

template<typename V>
struct weighted_feature_set {
	feature_set features;
	V log_probability;

	inline unsigned int get_feature(unsigned int index) const {
		return features[index];
	}

	inline void set_feature(unsigned int index, unsigned int feature) {
		features[index] = feature;
	}

	inline bool ensure_excluded_capacity(unsigned int index, unsigned int capacity) {
		return features.ensure_excluded_capacity(index, capacity);
	}

	/* NOTE: this function assumes the excluded array at the given index is uninitialized */
	inline bool set_excluded(unsigned int index, const unsigned int* src, unsigned int count) {
		return features.set_excluded(index, src, count);
	}

	/* NOTE: this function assumes the excluded array has sufficient capacity */
	inline void exclude_unsorted(unsigned int index, unsigned int item) {
		return features.exclude_unsorted(index, item);
	}

	inline void sort_excluded(unsigned int index) {
		return features.sort_excluded(index);
	}

	inline V get_probability() const {
		return exp(log_probability);
	}

	inline void set_probability(V probability) {
		log_probability = log(probability);
	}

	static inline void move(const weighted_feature_set<V>& src, weighted_feature_set<V>& dst) {
		dst.log_probability = src.log_probability;
		feature_set::move(src.features, dst.features);
	}

	static inline void swap(weighted_feature_set<V>& first, weighted_feature_set<V>& second) {
		feature_set::swap(first.features, second.features);
		core::swap(first.log_probability, second.log_probability);
	}

	static inline void free(weighted_feature_set<V>& posterior) {
		core::free(posterior.features);
	}
};

template<typename V>
inline bool init(weighted_feature_set<V>& posterior, unsigned int depth) {
	return init(posterior.features, depth);
}

template<typename V>
inline bool init(weighted_feature_set<V>& posterior, const feature_set& src, V log_probability) {
	posterior.log_probability = log_probability;
	return init(posterior.features, src);
}

template<typename V>
inline bool operator < (const weighted_feature_set<V>& first, const weighted_feature_set<V>& second) {
	return first.log_probability > second.log_probability;
}

template<typename V>
inline double log_probability(const weighted_feature_set<V>& set) {
	return set.log_probability;
}

#endif /* FEATURES_H_ */
