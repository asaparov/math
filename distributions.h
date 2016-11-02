/**
 * distributions.h
 *
 *  Created on: Jul 31, 2015
 *      Author: asaparov
 */

#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_

#include <core/random.h>

#include "histogram.h"
#include "sparse_vector.h"

/* TODO: add documentation */
#define CATEGORICAL_MAX_THRESHOLD 16.0
#define CATEGORICAL_MIN_THRESHOLD 1.0e-5

/* forward declarations */
template<typename V> struct dirichlet;

template<typename V>
struct array_categorical {
	V* log_probabilities;
	V* probabilities; /* unnormalized */
	V maximum; /* of log_probabilities */

	array_categorical(unsigned int length) : log_probabilities(NULL), probabilities(NULL)
	{
		if (!initialize(length))
			exit(EXIT_FAILURE);
	}

	~array_categorical() { free(); }

	inline bool resize(unsigned int new_length) {
		return core::resize(log_probabilities, new_length)
			&& core::resize(probabilities, new_length);
	}

	void place(unsigned int i, const V& log_probability) {
		log_probabilities[i] = log_probability;
	}

	void renormalize(unsigned int length) {
		maximum = max(log_probabilities, length);
		for (unsigned int i = 0; i < length; i++)
			probabilities[i] = exp(log_probabilities[i] - maximum);
	}

	static inline void move(const array_categorical<V>& src, array_categorical<V>& dst) {
		dst.log_probabilities = src.log_probabilities;
		dst.probabilities = src.probabilities;
		dst.maximum = src.maximum;
	}

	static inline void swap(array_categorical<V>& first, array_categorical<V>& second) {
		core::swap(first.log_probabilities, second.log_probabilities);
		core::swap(first.probabilities, second.probabilities);
		core::swap(first.maximum, second.maximum);
	}

	/* NOTE: this function assumes that the type V has constant size */
	static inline long unsigned int size_of(const array_categorical<V>& categorical, unsigned int length) {
		return sizeof(V) * (2 * length + 1);
	}

	static inline void free(array_categorical<V>& categorical) {
		categorical.free();
	}

private:
	inline bool initialize(unsigned int length) {
		log_probabilities = (V*) malloc(sizeof(V) * length);
		if (log_probabilities == NULL) {
			fprintf(stderr, "array_categorical ERROR: Out of memory.\n");
			return false;
		}
		probabilities = (V*) malloc(sizeof(V) * length);
		if (probabilities == NULL) {
			fprintf(stderr, "array_categorical ERROR: Out of memory.\n");
			core::free(probabilities);
			return false;
		}
		return true;
	}

	inline void free() {
		core::free(log_probabilities);
		core::free(probabilities);
	}

	template<typename A>
	friend bool init(array_categorical<A>&, unsigned int);
};

template<typename V>
bool init(array_categorical<V>& categorical, unsigned int length) {
	return categorical.initialize(length);
}

template<typename V>
struct symmetric_dirichlet {
	typedef V value_type;

	V pi;
	unsigned int atom_count;

	V log_prob;
	V total;

	symmetric_dirichlet(const symmetric_dirichlet<V>& prior) :
		pi(prior.pi), atom_count(prior.atom_count), log_prob(-log(atom_count)), total(prior.pi * atom_count) { }

	symmetric_dirichlet(unsigned int atom_count, const V& prior) :
		pi(prior), atom_count(atom_count), log_prob(-log(atom_count)), total(prior * atom_count) { }

	inline void ensure_atom_count(unsigned int new_atom_count) {
		if (new_atom_count > atom_count)
			atom_count = new_atom_count;
		total = pi * atom_count;
	}

	inline V sum() const {
		return total;
	}

	inline V max() const {
		return pi;
	}

	inline V get_at_index(unsigned int index) const {
		return pi;
	}

	inline V get_for_atom(unsigned int atom) const {
		return pi;
	}

	template<typename K>
	inline void add(sparse_vector<K, V>& vector) const {
		vector.unspecified = pi;
	}

	inline void print(FILE* out) const {
		fprintf(out, "pi: %lf\n", pi);
		fprintf(out, "atom count: %u\n", atom_count);
	}

	inline const symmetric_dirichlet<V>& get_parameters() const {
		return *this;
	}

	template<typename Destination>
	inline void sample(Destination& dst) const {
		std::gamma_distribution<V> distribution(pi, 1.0);

		V sum = 0.0;
		for (unsigned int i = 0; i < atom_count; i++) {
			V value = distribution(engine);
			dst.set(i, value);
			sum += value;
		}

		for (unsigned int i = 0; i < atom_count; i++)
			dst.set(i, dst.get(i) / sum);
	}

	V log_probability(unsigned int item) const {
		return log_prob;
	}

	template<typename Metric>
	static inline long unsigned int size_of(const symmetric_dirichlet<V>& distribution, const Metric& metric) {
		return core::size_of(distribution.atom_count) + core::size_of(distribution.log_prob)
			 + core::size_of(distribution.pi) + core::size_of(distribution.total);
	}

	static inline void free(symmetric_dirichlet<V>& distribution) { }

private:
	bool initialize(unsigned int length, const V& pi_src) {
		pi = pi_src;
		atom_count = length;
		log_prob = -log(atom_count);
		total = pi * atom_count;
		return true;
	}

	bool initialize(unsigned int length, const V* pi_src) {
		fprintf(stderr, "symmetric_dirichlet_prior.initialize ERROR: "
				"Dense parameterization for the prior is not supported"
				" with a symmetric Dirichlet prior.");
		return false;
	}

	template<typename A>
	friend bool init(
			symmetric_dirichlet<A>& prior,
			const symmetric_dirichlet<A>& src);
};

template<typename V>
inline bool init(symmetric_dirichlet<V>& distribution,
		const symmetric_dirichlet<V>& src)
{
	return distribution.initialize(src.atom_count, src.pi);
}

template<typename V>
inline bool init(symmetric_dirichlet<V>& distribution, const dirichlet<V>& src) {
	fprintf(stderr, "init ERROR: Unsupported initialization of "
			"symmetric_dirichlet_prior with a dirichlet_prior argument.\n");
	exit(EXIT_FAILURE);
}

template<typename V>
inline bool read(symmetric_dirichlet<V>& distribution, FILE* in) {
	if (!read(distribution.pi, in)) return false;
	if (!read(distribution.atom_count, in)) return false;
	distribution.log_prob = -log(distribution.atom_count);
	return true;
}

template<typename V>
inline bool write(const symmetric_dirichlet<V>& distribution, FILE* out) {
	if (!write(distribution.pi, out)) return false;
	return write(distribution.atom_count, out);
}

template<typename V>
struct dirichlet {
	typedef V value_type;

	V* pi;
	V pi_sum;
	unsigned int atom_count;

	dirichlet(const symmetric_dirichlet<V>& prior) {
		if (!initialize(prior.atom_count, prior.prior))
			exit(EXIT_FAILURE);
	}

	dirichlet(const dirichlet<V>& prior) {
		if (!initialize(prior.atom_count, prior.pi))
			exit(EXIT_FAILURE);
	}

	dirichlet(unsigned int atom_count, const V& prior) {
		if (!initialize(atom_count, prior))
			exit(EXIT_FAILURE);
	}

	dirichlet(unsigned int atom_count, const V* prior) {
		if (!initialize(atom_count, prior))
			exit(EXIT_FAILURE);
	}

	~dirichlet() { free(); }

	inline void ensure_atom_count(unsigned int new_atom_count) {
		if (new_atom_count <= atom_count)
			return;
		fprintf(stderr, "dirichlet.ensure_atom_count ERROR: This is not implemented.\n");
	}

	inline V sum() const {
		return pi_sum;
	}

	inline V get_at_index(unsigned int index) const {
		return pi[index];
	}

	inline V get_for_atom(unsigned int atom) const {
		return pi[atom - 1];
	}

	template<typename K>
	void add(sparse_vector<K, V>& vector) const {
		for (unsigned int i = 0; i < vector.length; i++)
			vector.set(i + 1, pi[i]);
	}

	void print(FILE* out) const {
		if (atom_count == 0) {
			fprintf(out, "pi: []\n");
			return;
		}
		fprintf(out, "pi: [%lf", pi[0]);
		for (unsigned int i = 1; i < atom_count; i++)
			fprintf(out, ", %lf", pi[i]);
		fprintf(out, "]\n");
		fprintf(out, "atom count: %u\n", atom_count);
	}

	inline const dirichlet<V>& get_parameters() const {
		return *this;
	}

	template<typename Destination>
	inline void sample(Destination& dst) const {
		V sum = 0.0;
		for (unsigned int i = 0; i < atom_count; i++) {
			std::gamma_distribution<V> distribution(pi[i], 1.0);
			V value = distribution(engine);
			dst.set(i, value);
			sum += value;
		}

		for (unsigned int i = 0; i < atom_count; i++)
			dst.set(i, dst.get(i) / sum);
	}

	V log_probability(unsigned int atom) const {
		return log(pi[atom - 1]) - log(pi_sum);
	}

	/* NOTE: this function assumes that the type V has constant size */
	template<typename Metric>
	static inline long unsigned int size_of(const dirichlet<V>& distribution, const Metric& metric) {
		return core::size_of(distribution.atom_count) + core::size_of(distribution.pi_sum) + sizeof(V) * distribution.atom_count;
	}

	static inline void free(dirichlet<V>& distribution) {
		distribution.free();
	}

private:
	bool initialize(unsigned int length, const V& pi_src) {
		atom_count = length;
		pi = (V*) malloc(sizeof(V) * atom_count);
		if (pi == NULL) {
			fprintf(stderr, "dirichlet_prior.initialize ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < atom_count; i++)
			pi[i] = pi_src;
		pi_sum = pi_src * atom_count;
		return true;
	}

	bool initialize(unsigned int length, const V* pi_src) {
		atom_count = length;
		pi = (V*) malloc(sizeof(V) * atom_count);
		if (pi == NULL) {
			fprintf(stderr, "dirichlet_prior.initialize ERROR: Out of memory.\n");
			return false;
		}
		pi_sum = 0.0;
		for (unsigned int i = 0; i < atom_count; i++) {
			pi[i] = pi_src[i];
			pi_sum += pi_src[i];
		}
		return true;
	}

	inline void free() {
		core::free(pi);
	}

	template<typename K>
	friend bool init(dirichlet<K>&, const symmetric_dirichlet<K>&);

	template<typename K>
	friend bool init(dirichlet<K>&, const dirichlet<K>&);
};

template<typename V>
inline bool init(dirichlet<V>& distribution,
		const symmetric_dirichlet<V>& src)
{
	return distribution.initialize(src.atom_count, src.pi);
}

template<typename V>
inline bool init(dirichlet<V>& distribution, const dirichlet<V>& src)
{
	return distribution.initialize(src.atom_count, src.pi);
}

template<typename V>
inline bool read(dirichlet<V>& distribution, FILE* in) {
	if (!read(distribution.atom_count, in))
		return false;
	distribution.pi = (V*) malloc(sizeof(V) * distribution.atom_count);
	if (distribution.pi == NULL) return false;
	if (!read(distribution.pi, in, distribution.atom_count)) {
		free(distribution.pi);
		return false;
	}
	distribution.pi_sum = 0.0;
	for (unsigned int i = 0; i < distribution.atom_count; i++)
		distribution.pi_sum += distribution.pi[i];
	return true;
}

template<typename V>
inline bool write(const dirichlet<V>& distribution, FILE* out) {
	if (!write(distribution.atom_count, out)) return false;
	return write(distribution.pi, out, distribution.atom_count);
}


/**
 * Some useful type traits for Dirichlet distribution structures.
 */

template<typename T>
struct is_dirichlet : std::false_type { };

template<typename V>
struct is_dirichlet<symmetric_dirichlet<V>> : std::true_type { };

template<typename V>
struct is_dirichlet<dirichlet<V>> : std::true_type { };

/* a categorical distribution represented as an array of probabilities */
template<typename V>
struct dense_categorical
{
	typedef V value_type;

	V* phi;
	unsigned int atom_count;

	dense_categorical(unsigned int atom_count) : atom_count(atom_count) {
		phi = (V*) calloc(atom_count, sizeof(V));
		if (phi == NULL) {
			fprintf(stderr, "dense_categorical ERROR: Insufficient memory for phi.\n");
			exit(EXIT_FAILURE);
		}
	}

	dense_categorical(const dense_categorical& src) : atom_count(src.atom_count) {
		phi = (V*) malloc(atom_count * sizeof(V));
		if (phi == NULL) {
			fprintf(stderr, "dense_categorical ERROR: Insufficient memory for phi.\n");
			exit(EXIT_FAILURE);
		}
		memcpy(phi, src.phi, sizeof(V) * atom_count);
	}

	~dense_categorical() { free(); }

	inline V get(unsigned int index) const {
		return phi[index];
	}

	inline V get_for_atom(unsigned int index) const {
		return phi[index - 1];
	}

	inline V sum() const {
		return 1.0;
	}

	void ensure_atom_count(unsigned int new_atom_count) {
		if (new_atom_count <= atom_count)
			return;
		fprintf(stderr, "dense_categorical.set_atom_count ERROR: This is not implemented.\n");
	}

	inline V probability(unsigned int item) const {
#if !defined(NDEBUG)
		if (item == 0) {
			fprintf(stderr, "dense_categorical.conditional ERROR: Given item is zero.\n");
			return std::numeric_limits<V>::signaling_NaN();
		}
#endif
		return phi[item - 1];
	}

	inline V probability(const array_histogram<unsigned int>& items) const {
		V value = 1.0;
		for (unsigned int i = 0; i < items.counts.size; i++)
			value *= pow(probability(items.counts.keys[i]), items.counts.values[i]);
		return value;
	}

	inline V log_probability(unsigned int item) const {
		return log(phi[item - 1]);
	}

	inline V log_probability(const array_histogram<unsigned int>& items) const {
		V value = 0.0;
		for (const auto& entry : items.counts)
			value = log_probability(entry.key) * entry.value;
		return value;
	}

	template<typename PriorDistribution>
	static inline V probability(const PriorDistribution& prior, unsigned int item) {
		return prior.get_for_atom(item) / prior.sum();
	}

	template<typename PriorDistribution>
	static inline V log_probability(const PriorDistribution& prior, unsigned int item) {
		return prior.log_probability(item);
	}

	template<typename PriorDistribution>
	static inline V log_probability(const PriorDistribution& prior, const array_histogram<unsigned int>& items)
	{
		double log_probability = 0.0;
		for (unsigned int i = 0; i < items.counts.size; i++)
			log_probability += log_rising_factorial(
					prior.get_for_atom(items.counts.keys[i]), items.counts.values[i]);
		return log_probability - log_rising_factorial(prior.sum(), items.total());
	}

	template<typename PriorDistribution>
	static inline V conditional(const PriorDistribution& prior,
			unsigned int item, const array_histogram<unsigned int>& conditioned)
	{
		for (unsigned int i = 0; i < conditioned.counts.size; i++) {
			if (conditioned.counts.keys[i] == item) {
				return (prior.get_for_atom(item) + conditioned.counts.values[i])
						/ (prior.sum() + conditioned.total());
			}
		}

		return (prior.get_for_atom(item)) / (prior.sum() + conditioned.total());
	}

	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			unsigned int item, const array_histogram<unsigned int>& conditioned)
	{
		for (unsigned int i = 0; i < conditioned.counts.size; i++) {
			if (conditioned.counts.keys[i] == item) {
				return log((prior.get_for_atom(item) + conditioned.counts.values[i]))
						- log(prior.sum() + conditioned.total());
			}
		}

		return log(prior.get_for_atom(item)) - log(prior.sum() + conditioned.total());
	}

	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			unsigned int item, const hash_histogram<unsigned int>& conditioned)
	{
		bool contains;
		unsigned int count = conditioned.counts.get(item, contains);

		if (contains) {
			return fasterlog(prior.get_for_atom(item) + count);
		} else {
			return fasterlog(prior.get_for_atom(item));
		}
	}

	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			unsigned int holdout, const array_histogram<unsigned int>& conditioned)
	{
		for (unsigned int i = 0; i < conditioned.counts.size; i++) {
			if (conditioned.counts.keys[i] == holdout) {
				return log(prior.get_for_atom(holdout) + conditioned.counts.values[i] - 1)
						- log(prior.sum() + conditioned.total() - 1);
			}
		}

		return log(prior.get_for_atom(holdout)) - log(prior.sum() + conditioned.total());
	}

	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			unsigned int holdout, const hash_histogram<unsigned int>& conditioned)
	{
		bool contains;
		unsigned int count = conditioned.counts.get(holdout, contains);

		double to_return;
		if (contains) {
			return log(prior.get_for_atom(holdout) + count - 1);
		} else {
			return log(prior.get_for_atom(holdout));
		}
	}

	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			const array_histogram<unsigned int>& items, const array_histogram<unsigned int>& conditioned)
	{
		V log_probability = 0.0;
		unsigned int i = 0, j = 0;
		while (i < items.counts.size && j < conditioned.counts.size) {
			if (items.counts.keys[i] == conditioned.counts.keys[j]) {
				log_probability +=
						log_rising_factorial(prior.get_for_atom(items.counts.keys[i])
						+ conditioned.counts.values[j], items.counts.values[i]);
				i++; j++;
			} else {
				j++;
			}
		}

		while (i < items.counts.size) {
			log_probability += log_rising_factorial(
					prior.get_for_atom(items.counts.keys[i]), items.counts.values[i]);
			i++;
		}

		return log_probability - log_rising_factorial(prior.sum() + conditioned.total(), items.total());
	}

	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			const array_histogram<unsigned int>& items, const hash_histogram<unsigned int>& conditioned)
	{
		V log_probability = 0.0;
		for (unsigned int i = 0; i < items.counts.size; i++) {
			bool contains;
			unsigned int count = conditioned.counts.get(items.counts.keys[i], contains);
			if (contains) {
				log_probability += log_rising_factorial(
						prior.get_for_atom(items.counts.keys[i]) + count, items.counts.values[i]);
			} else {
				log_probability += log_rising_factorial(
						prior.get_for_atom(items.counts.keys[i]), items.counts.values[i]);
			}
		}

		return log_probability;
	}

	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			const array_histogram<unsigned int>& holdout, const array_histogram<unsigned int>& conditioned)
	{
		V log_probability = 0.0;
		unsigned int i = 0, j = 0;
		while (i < holdout.counts.size && j < conditioned.counts.size) {
			if (holdout.counts.keys[i] == conditioned.counts.keys[j]) {
#if !defined(NDEBUG)
				if (conditioned.counts.values[j] < holdout.counts.values[i])
					fprintf(stderr, "dense_categorical.conditional_without WARNING:"
							" 'holdout' is not a subset of 'conditioned'.\n");
#endif
				log_probability += log_rising_factorial(prior.get_for_atom(holdout.counts.keys[i])
						+ conditioned.counts.values[j] - holdout.counts.values[i], holdout.counts.values[i]);
				i++; j++;
			} else {
				j++;
			}
		}

		return log_probability - log_rising_factorial(
				prior.sum() + conditioned.total() - holdout.total(), holdout.total());
	}

	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			const array_histogram<unsigned int>& holdout, const hash_histogram<unsigned int>& conditioned)
	{
		V log_probability = 0.0;
		for (unsigned int i = 0; i < holdout.counts.size; i++) {
			bool contains;
			unsigned int count = conditioned.counts.get(holdout.counts.keys[i], contains);
			if (contains) {
#if !defined(NDEBUG)
				if (count < holdout.counts.values[i])
					fprintf(stderr, "dense_categorical.conditional_without WARNING:"
							" 'holdout' is not a subset of 'conditioned'.\n");
#endif
				log_probability += log_rising_factorial(prior.get_for_atom(holdout.counts.keys[i])
						+ count - holdout.counts.values[i], holdout.counts.values[i]);
			}
		}

		return log_probability;
	}

	static inline void move(const dense_categorical<V>& src, dense_categorical<V>& dst) {
		dst.phi = src.phi;
		dst.atom_count = src.atom_count;
	}

	/* NOTE: this function assumes that the type V has constant size */
	template<typename Metric>
	static inline long unsigned int size_of(const dense_categorical<V>& distribution, const Metric& metric) {
		return sizeof(V) * distribution.atom_count
			 + core::size_of(distribution.atom_count);
	}

	static inline void free(dense_categorical<V>& distribution) {
		distribution.free();
	}

	inline unsigned int get_parameters() const {
		return atom_count;
	}

private:
	inline void free() {
		core::free(phi);
	}
};

template<typename V>
inline bool init(dense_categorical<V>& distribution, unsigned int atom_count) {
	distribution.atom_count = atom_count;
	distribution.phi = (V*) calloc(atom_count, sizeof(V));
	if (distribution.phi == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for phi in the dense_categorical.\n");
		return false;
	}
	return true;
}

template<typename V>
inline bool init(dense_categorical<V>& distribution, const dense_categorical<V>& src) {
	distribution.atom_count = src.atom_count;
	distribution.phi = (V*) malloc(src.atom_count * sizeof(V));
	if (distribution.phi == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for phi in the dense_categorical.\n");
		return false;
	}
	memcpy(distribution.phi, src.phi, sizeof(V) * src.atom_count);
	return true;
}

template<typename V, typename Stream>
bool read(dense_categorical<V>& distribution, Stream& stream) {
	if (!read(distribution.atom_count, stream))
		return false;
	distribution.phi = (V*) malloc(distribution.atom_count * sizeof(V));
	if (distribution.phi == NULL) {
		fprintf(stderr, "read ERROR: Insufficient memory for phi in the dense_categorical.\n");
		return false;
	}
	return read(distribution.phi, stream, distribution.atom_count);
}

template<typename V, typename Stream>
bool write(const dense_categorical<V>& distribution, Stream& stream) {
	return write(distribution.atom_count, stream)
		&& write(distribution.phi, stream, distribution.atom_count);
}


template<typename K, typename V>
struct sparse_categorical
{
	typedef V value_type;

	hash_map<K, pair<V, V>> probabilities;
	unsigned int atom_count;

	V prob;
	V dense_prob;
	V log_prob;

	sparse_categorical(unsigned int atom_count) :
		probabilities(16), atom_count(atom_count),
		prob(1.0 / atom_count), dense_prob(0.0), log_prob(-log(atom_count)) { }
	~sparse_categorical() { free(); }

	bool set(const K& key, const V& probability) {
		if (!probabilities.check_size())
			return false;

		bool contains; unsigned int index;
		pair<V, V>& value = probabilities.get(key, contains, index);
		if (!contains) {
			probabilities.table.keys[index] = key;
			probabilities.table.size++;
			value.key = probability;
			value.value = log(probability);
		} else dense_prob -= value.key;

		dense_prob += probability;
		if (probabilities.table.size >= atom_count) {
			prob = 0.0;
			log_prob = -std::numeric_limits<V>::infinity();
		} else {
			prob = (1.0 - dense_prob) / (atom_count - probabilities.table.size);
			log_prob = log(prob);
		}
		return true;
	}

	inline V probability(const K& observation) const {
		bool contains;
		const pair<V, V>& entry = probabilities.get(observation, contains);
		if (contains) return entry.key;
		else return prob;
	}

	inline V log_probability(const K& observation) const {
		bool contains;
		const pair<V, V>& entry = probabilities.get(observation, contains);
		if (contains) return entry.value;
		else return log_prob;
	}

private:
	inline void free() {
		for (auto entry : probabilities)
			core::free(entry.key);
	}
};


/* the degenerate distribution */
template<typename K>
struct constant
{
	/* returns true if the probability is 1, and false otherwise */
	template<typename PriorDistribution>
	static inline bool conditional(const PriorDistribution& prior,
			const K& item, const array_histogram<K>& conditioned)
	{
		return item == conditioned.counts.keys[0];
	}

	template<typename PriorDistribution>
	static inline bool conditional(const PriorDistribution& prior,
			const array_histogram<K>& items, const array_histogram<K>& conditioned)
	{
		return items.counts.keys[0] == conditioned.counts.keys[0];
	}

	template<typename PriorDistribution>
	static inline double log_conditional(const PriorDistribution& prior,
			const K& item, const array_histogram<K>& conditioned)
	{
		if (item == conditioned.counts.keys[0])
			return 0.0;
		else return -std::numeric_limits<double>::infinity();
	}

	template<typename PriorDistribution>
	static inline double log_conditional(const PriorDistribution& prior,
			const array_histogram<K>& items, const array_histogram<K>& conditioned)
	{
		if (items.counts.keys[0] == conditioned.counts.keys[0])
			return 0.0;
		else return -std::numeric_limits<double>::infinity();
	}

	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type probability(
		const PriorDistribution& prior, const K& item)
	{
		return prior.probability(item);
	}

	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type probability(
		const PriorDistribution& prior, const array_histogram<K>& items)
	{
		return prior.probability(items);
	}

	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type log_probability(
		const PriorDistribution& prior, const K& item)
	{
		return prior.log_probability(item);
	}

	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type log_probability(
		const PriorDistribution& prior, const array_histogram<K>& items)
	{
		return prior.log_probability(items);
	}
};


/* a uniform distribution */
template<typename V>
struct uniform_distribution {
	typedef V value_type;

	V prob;
	V log_prob;

	uniform_distribution(unsigned int count) :
		prob(1.0 / count), log_prob(-log((V) count)) { }

	template<typename K>
	inline V probability(const K& observation) const {
		return prob;
	}

	template<typename K>
	inline V log_probability(const K& observation) const {
		return log_prob;
	}
};


/* distribution over a sequences of independent events */
template<typename ElementDistribution>
struct sequence_distribution
{
	typedef typename ElementDistribution::value_type V;

	V end_probability;
	V log_end_probability;
	V log_not_end_probability;
	ElementDistribution element_distribution;

	sequence_distribution(ElementDistribution& element_distribution, double end_probability) :
		end_probability(end_probability),
		log_end_probability(log(end_probability)),
		log_not_end_probability(log(1.0 - end_probability)),
		element_distribution(element_distribution)
	{ }

	template<typename SequenceType>
	inline double probability(const SequenceType& sequence) const {
		if (sequence.length == 0) return 0.0;
		double product = element_distribution.probability(sequence[0]);
		for (unsigned int i = 1; i < sequence.length; i++)
			product *= element_distribution.probability(sequence[i]) * (1.0 - end_probability);
		return product * end_probability;
	}

	template<typename SequenceType>
	inline double log_probability(const SequenceType& sequence) const {
		if (sequence.length == 0) return -std::numeric_limits<double>::infinity();
		double sum = element_distribution.log_probability(sequence[0]);
		for (unsigned int i = 1; i < sequence.length; i++)
			sum += element_distribution.log_probability(sequence[i]) + log_not_end_probability;
		return sum + log_end_probability;
	}

	static inline void free(sequence_distribution<ElementDistribution>& distribution) {
		core::free(distribution.element_distribution);
	}
};

template<typename ElementDistribution>
inline bool init(sequence_distribution<ElementDistribution>& distribution,
		const sequence_distribution<ElementDistribution>& src)
{
	distribution.end_probability = src.end_probability;
	distribution.log_end_probability = src.log_end_probability;
	distribution.log_not_end_probability = src.log_not_end_probability;
	return init(distribution.element_distribution, src.element_distribution);
}

template<typename ElementDistribution, typename Stream>
bool read(sequence_distribution<ElementDistribution>& distribution, Stream& stream)
{
	if (!read(distribution.end_probability, stream))
		return false;
	distribution.log_end_probability = log(distribution.end_probability);
	distribution.log_not_end_probability = log(1.0 - distribution.end_probability);
	return read(distribution.element_distribution, stream);
}

template<typename ElementDistribution, typename Stream>
bool write(const sequence_distribution<ElementDistribution>& distribution, Stream& stream)
{
	return write(distribution.end_probability, stream)
		&& write(distribution.element_distribution, stream);
}

#endif /* DISTRIBUTIONS_H_ */
