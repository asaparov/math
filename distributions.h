/**
 * \file distributions.h
 *
 *  <!-- Created on: Jul 31, 2015
 *           Author: asaparov -->
 */

#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_

#include <core/random.h>

#include "multiset.h"
#include "sparse_vector.h"
#include "log.h"

/* TODO: add documentation */
#define CATEGORICAL_MAX_THRESHOLD 16.0
#define CATEGORICAL_MIN_THRESHOLD 1.0e-5

/* forward declarations */
#if !defined(DOXYGEN_IGNORE)
template<typename V> struct dirichlet;
#endif

/**
 * Returns the natural logarithm of the
 * [rising factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials):
 * \f[ \log a^{(n)} = \log \left\{ \prod_{i=0}^{n-1} a^i \right\}, \f]
 * where `base` is \f$ a \f$ and `exponent` is \f$ n \f$.
 */
inline double log_rising_factorial(double base, unsigned int exponent) {
	return lgamma(base + exponent) - lgamma(base);
}

/* \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic). */
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

/**
 * A structure representing a symmetric Dirichlet distribution, which is a
 * regular [Dirichlet distribution](https://en.wikipedia.org/wiki/Dirichlet_distribution),
 * where all the elements of the concentration parameter vector *pi* have the
 * same value.
 *
 * The following example constructs a Dirichlet distribution with dimension 3,
 * and concentration parameter [0.1, 0.1, 0.1]. It proceeds to generate two
 * samples from this distribution. The expected output is
 * `[0.004223, 0.995777, 0.000000] [0.964576, 0.035400, 0.000024]`.
 *
 * ```{.cpp}
 * #include <math/distributions.h>
 * using namespace core;
 *
 * template<typename V>
 * struct custom_vector {
 * 	V* elements;
 *
 * 	custom_vector(unsigned int length) {
 * 		elements = (V*) malloc(sizeof(V) * length);
 * 	}
 *
 * 	~custom_vector() { free(elements); }
 *
 * 	inline V get(unsigned int index) const {
 * 		return elements[index];
 * 	}
 *
 * 	void set(unsigned int index, const V& value) {
 * 		elements[index] = value;
 * 	}
 * };
 *
 * int main() {
 * 	set_seed(100);
 * 	symmetric_dirichlet<double> dir(3, 0.1);
 *
 * 	custom_vector<double> output(3);
 * 	for (unsigned int i = 0; i < 2; i++) {
 * 		sample(dir, output);
 * 		print(output.elements, 3, stdout); print(' ', stdout);
 * 	}
 * }
 * ```
 *
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
struct symmetric_dirichlet {
	/**
	 * The type of the probabilities.
	 */
	typedef V value_type;

	/**
	 * The value of each element of the concentration parameter.
	 */
	V pi;

	/**
	 * The number of dimensions of the Dirichlet distribution.
	 */
	unsigned int atom_count;

	/**
	 * The negative logarithm of `symmetric_dirichlet::atom_count`.
	 */
	V log_prob;

	/**
	 * The sum of the elements in the concentration parameter vector. In the
	 * case of a symmetric Dirichlet distribution, this simply
	 * `symmetric_dirichlet::pi * symmetric_dirichlet::atom_count`.
	 */
	V total;

	/**
	 * Constructs a symmetric Dirichlet distribution by copying the fields from
	 * the given symmetric Dirichlet distribution.
	 */
	symmetric_dirichlet(const symmetric_dirichlet<V>& prior) :
		pi(prior.pi), atom_count(prior.atom_count), log_prob(-log(atom_count)), total(prior.pi * atom_count) { }

	/**
	 * Constructs a symmetric Dirichlet distribution with the given dimension
	 * and concentration parameter.
	 */
	symmetric_dirichlet(unsigned int atom_count, const V& prior) :
		pi(prior), atom_count(atom_count), log_prob(-log(atom_count)), total(prior * atom_count) { }

	/**
	 * Checks if symmetric_dirichlet::atom_count is smaller than the given
	 * `new_atom_count`. If so, symmetric_dirichlet::atom_count is set to
	 * `new_atom_count` and symmetric_dirichlet::total is updated accordingly.
	 */
	inline void ensure_atom_count(unsigned int new_atom_count) {
		if (new_atom_count > atom_count)
			atom_count = new_atom_count;
		total = pi * atom_count;
	}

	/**
	 * Returns the sum of the elements in the concentration parameter vector
	 * (symmetric_dirichlet::total).
	 */
	inline V sum() const {
		return total;
	}

	/**
	 * Returns the maximum of the elements in the concentration parameter
	 * vector. In the case of the symmetric Dirichlet distribution, this is
	 * simply symmetric_dirichlet::pi.
	 */
	inline V max() const {
		return pi;
	}

	/**
	 * Returns the element at the given `index` in the concentration parameter
	 * vector. In the case of the symmetric Dirichlet distribution, this
	 * function always returns symmetric_dirichlet::pi.
	 */
	inline V get_at_index(unsigned int index) const {
		return pi;
	}

	/**
	 * Returns the element at the given `atom` in the concentration parameter
	 * vector. The atom is a non-empty element drawn from a categorical
	 * distribution with a Dirichlet prior. In the case of the symmetric
	 * Dirichlet distribution, this function always returns
	 * symmetric_dirichlet::pi.
	 */
	template<typename K>
	inline V get_for_atom(const K& atom) const {
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

	/**
	 * Returns the parameters that may be used to construct this distribution,
	 * using either the constructor or init.
	 */
	inline const symmetric_dirichlet<V>& get_parameters() const {
		return *this;
	}

	/**
	 * Samples from this symmetric Dirichlet distribution and puts the result in `dst`.
	 * \tparam Destination a vector type that implements the public member
	 * 		function `set(unsigned int, const V&)`.
	 */
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

	/**
	 * Returns the log probability of a single observation `atom`, drawn from a
	 * Dirichlet-categorical distribution, where the Dirichlet is represented
	 * by this object. In the case of the symmetric Dirichlet distribution,
	 * this function always returns symmetric_dirichlet::log_prob.
	 */
	V log_probability(unsigned int item) const {
		return log_prob;
	}

	template<typename Metric>
	static inline long unsigned int size_of(const symmetric_dirichlet<V>& distribution, const Metric& metric) {
		return core::size_of(distribution.atom_count) + core::size_of(distribution.log_prob)
			 + core::size_of(distribution.pi) + core::size_of(distribution.total);
	}

	/**
	 * Frees the given symmetric Dirichlet distribution. Since
	 * symmetric_dirichlet does not allocate additional memory, this function
	 * is a no-op.
	 */
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

/**
 * Initializes the given symmetric_dirichlet `distribution` by copying its
 * fields from the given symmetric_dirichlet `src`.
 */
template<typename V>
inline bool init(symmetric_dirichlet<V>& distribution,
		const symmetric_dirichlet<V>& src)
{
	return distribution.initialize(src.atom_count, src.pi);
}

/**
 * Since a symmetric Dirichlet is a special case of a Dirichlet distribution,
 * this function prints an error and exits.
 */
template<typename V>
inline bool init(symmetric_dirichlet<V>& distribution, const dirichlet<V>& src) {
	fprintf(stderr, "init ERROR: Unsupported initialization of "
			"symmetric_dirichlet_prior with a dirichlet_prior argument.\n");
	exit(EXIT_FAILURE);
}

/**
 * Reads a symmetric_dirichlet `distribution` from `in`.
 */
template<typename V>
inline bool read(symmetric_dirichlet<V>& distribution, FILE* in) {
	if (!read(distribution.pi, in)) return false;
	if (!read(distribution.atom_count, in)) return false;
	distribution.log_prob = -log(distribution.atom_count);
	return true;
}

/**
 * Writes the symmetric_dirichlet `distribution` to `out`.
 */
template<typename V>
inline bool write(const symmetric_dirichlet<V>& distribution, FILE* out) {
	if (!write(distribution.pi, out)) return false;
	return write(distribution.atom_count, out);
}

/**
 * Samples from the given symmetric Dirichlet `distribution` and puts the result in `dst`.
 * \tparam Destination a vector type that implements the public member function
 * 		`set(unsigned int, const V&)`.
 */
template<typename V, typename Destination>
inline bool sample(const symmetric_dirichlet<V>& distribution, Destination& dst) {
	distribution.sample(dst);
	return true;
}

/**
 * A structure representing a finite [Dirichlet distribution](https://en.wikipedia.org/wiki/Dirichlet_distribution),
 * which is a generalization of the symmetric Dirichlet distribution (see struct symmetric_dirichlet).
 *
 * The following example constructs a Dirichlet distribution with dimension 3,
 * and concentration parameter [10, 1, 2]. It proceeds to generate two samples
 * from this distribution. The expected output is
 * `[0.491806, 0.023958, 0.484236] [0.804916, 0.019964, 0.175120]`.
 *
 * ```{.cpp}
 * #include <math/distributions.h>
 * using namespace core;
 *
 * template<typename V>
 * struct custom_vector {
 * 	V* elements;
 *
 * 	custom_vector(unsigned int length) {
 * 		elements = (V*) malloc(sizeof(V) * length);
 * 	}
 *
 * 	~custom_vector() { free(elements); }
 *
 * 	inline V get(unsigned int index) const {
 * 		return elements[index];
 * 	}
 *
 * 	void set(unsigned int index, const V& value) {
 * 		elements[index] = value;
 * 	}
 * };
 *
 * int main() {
 * 	set_seed(100);
 * 	double alpha[] = {10.0, 1.0, 2.0};
 * 	dirichlet<double> dir(3, alpha);
 *
 * 	custom_vector<double> output(3);
 * 	for (unsigned int i = 0; i < 2; i++) {
 * 		sample(dir, output);
 * 		print(output.elements, 3, stdout); print(' ', stdout);
 * 	}
 * }
 * ```
 *
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
struct dirichlet {
	/**
	 * The type of the probabilities.
	 */
	typedef V value_type;

	/**
	 * The concentration parameter.
	 */
	V* pi;

	/**
	 * The sum of all elements in the concentration parameter dirichlet::pi.
	 */
	V pi_sum;

	/**
	 * The number of dimensions of this Dirichlet distribution.
	 */
	unsigned int atom_count;

	/**
	 * Constructs a Dirichlet distribution by copying the fields from the given
	 * symmetric Dirichlet distribution.
	 */
	dirichlet(const symmetric_dirichlet<V>& prior) {
		if (!initialize(prior.atom_count, prior.prior))
			exit(EXIT_FAILURE);
	}

	/**
	 * Constructs a Dirichlet distribution by copying the fields from the given
	 * Dirichlet distribution.
	 */
	dirichlet(const dirichlet<V>& prior) {
		if (!initialize(prior.atom_count, prior.pi))
			exit(EXIT_FAILURE);
	}

	/**
	 * Constructs a Dirichlet distribution with the given dimension
	 * `atom_count` and symmetric concentration parameter `prior`.
	 */
	dirichlet(unsigned int atom_count, const V& prior) {
		if (!initialize(atom_count, prior))
			exit(EXIT_FAILURE);
	}

	/**
	 * Constructs a Dirichlet distribution with the given dimension
	 * `atom_count` and concentration parameter vector `prior`.
	 */
	dirichlet(unsigned int atom_count, const V* prior) {
		if (!initialize(atom_count, prior))
			exit(EXIT_FAILURE);
	}

	~dirichlet() { free(); }

	/**
	 * Checks if dirichlet::atom_count is smaller than the given
	 * `new_atom_count`. If so, dirichlet::atom_count is set to
	 * `new_atom_count`.
	 */
	inline void ensure_atom_count(unsigned int new_atom_count) {
		if (new_atom_count <= atom_count)
			return;
		fprintf(stderr, "dirichlet.ensure_atom_count ERROR: This is not implemented.\n");
	}

	/**
	 * Returns the sum of the elements in the concentration parameter vector
	 * (dirichlet::pi_sum).
	 */
	inline V sum() const {
		return pi_sum;
	}

	/**
	 * Returns the element at the given `index` in the concentration parameter
	 * vector. This function does not perform any bounds checking.
	 */
	inline V get_at_index(unsigned int index) const {
		return pi[index];
	}

	/**
	 * Returns the element at the given `atom` in the concentration parameter
	 * vector. The atom is a non-zero unsigned integer drawn from a categorical
	 * distribution with a Dirichlet prior. Thus, the atom `n` corresponds to
	 * the index `n - 1`. This function does not perform any bounds checking.
	 */
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

	/**
	 * Returns the parameters that may be used to construct this distribution,
	 * using either the constructor or init.
	 */
	inline const dirichlet<V>& get_parameters() const {
		return *this;
	}

	/**
	 * Samples from this Dirichlet distribution and puts the result in `dst`.
	 * \tparam Destination a vector type that implements the public member
	 * 		function `set(unsigned int, const V&)`.
	 */
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

	/**
	 * Returns the log probability of a single observation `atom`, drawn from a
	 * Dirichlet-categorical distribution, where the Dirichlet is represented
	 * by this object.
	 */
	V log_probability(unsigned int atom) const {
		return log(pi[atom - 1]) - log(pi_sum);
	}

	/* NOTE: this function assumes that the type V has constant size */
	template<typename Metric>
	static inline long unsigned int size_of(const dirichlet<V>& distribution, const Metric& metric) {
		return core::size_of(distribution.atom_count) + core::size_of(distribution.pi_sum) + sizeof(V) * distribution.atom_count;
	}

	/**
	 * Frees dirichlet::pi in the given Dirichlet distribution.
	 */
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

/**
 * Initializes the given Dirichlet `distribution` by copying its fields from
 * the given symmetric_dirichlet distribution `src`.
 */
template<typename V>
inline bool init(dirichlet<V>& distribution,
		const symmetric_dirichlet<V>& src)
{
	return distribution.initialize(src.atom_count, src.pi);
}

/**
 * Initializes the given Dirichlet `distribution` by copying its fields from
 * the given Dirichlet distribution `src`.
 */
template<typename V>
inline bool init(dirichlet<V>& distribution, const dirichlet<V>& src)
{
	return distribution.initialize(src.atom_count, src.pi);
}

/**
 * Reads a dirichlet `distribution` from `in`.
 */
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

/**
 * Writes the given dirichlet `distribution` to `out`.
 */
template<typename V>
inline bool write(const dirichlet<V>& distribution, FILE* out) {
	if (!write(distribution.atom_count, out)) return false;
	return write(distribution.pi, out, distribution.atom_count);
}

/**
 * Samples from the given Dirichlet `distribution` and puts the result in `dst`.
 * \tparam Destination a vector type that implements the public member function
 * 		`set(unsigned int, const V&)`.
 */
template<typename V, typename Destination>
inline bool sample(const dirichlet<V>& distribution, Destination& dst) {
	distribution.sample(dst);
	return true;
}


/**
 * <!-- Some useful type traits for Dirichlet distribution structures. -->
 */

/**
 * A type trait that is [true_type](http://en.cppreference.com/w/cpp/types/integral_constant)
 * if and only if `T` is either symmetric_dirichlet or dirichlet.
 */
template<typename T>
struct is_dirichlet : std::false_type { };

template<typename V>
struct is_dirichlet<symmetric_dirichlet<V>> : std::true_type { };

template<typename V>
struct is_dirichlet<dirichlet<V>> : std::true_type { };


/**
 * This struct represents a categorical distribution, where the probabilities
 * are represented as a sum of uniform and a non-uniform component. The
 * non-uniform component is stored as a core::hash_map from atoms to
 * probabilities. The dense_categorical struct instead stores all probabilities
 * contiguously as a single native array. dense_categorical should be used if
 * the dimension is small or if the probabilities cannot be easily represented
 * sparsely as a sum of uniform and non-uniform components. Unlike
 * dense_categorical, the observations do not necessarily have type `unsigned
 * int`, and can have generic type `K`.
 *
 * The following is an example where a sparse_categorical distribution is
 * constructed over the domain <code>{'a', 'b', 'c', 'd', 'e'}</code>. The probability of
 * each event is specified and 10 samples are drawn. The expected output is
 * `d, b, d, d, b, c, a, b, e, e,`.
 * ```{.cpp}
 * #include <math/distributions.h>
 * using namespace core;
 *
 * int main() {
 * 	set_seed(100);
 * 	sparse_categorical<char, double> categorical(5);
 * 	categorical.set('a', 0.1);
 * 	categorical.set('b', 0.4);
 * 	categorical.set('c', 0.1);
 * 	categorical.set('d', 0.2);
 * 	categorical.set('e', 0.2);
 *
 * 	for (unsigned int i = 0; i < 10; i++) {
 * 		char c;
 * 		sample(categorical, c);
 * 		printf("%c, ", c);
 * 	}
 * }
 * ```
 *
 * \tparam K the generic type of the observations. `K` must satisfy either:
 * 		1. [is_fundamental](http://en.cppreference.com/w/cpp/types/is_fundamental),
 * 		2. [is_enum](http://en.cppreference.com/w/cpp/types/is_enum),
 * 		3. [is_pointer](http://en.cppreference.com/w/cpp/types/is_pointer),
 * 		4. implements the public static method `unsigned int hash(const T&)`,
 * 			the public static method `void is_empty(const T&)`, implements the
 * 			operators `==`, satisfies [CopyAssignable](https://en.cppreference.com/w/cpp/named_req/CopyAssignable),
 * 			and core::is_moveable. **NOTE:** The first argument to the `==`
 * 			operator may be empty.
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename K, typename V>
struct sparse_categorical
{
	/**
	 * The type of the probabilities.
	 */
	typedef V value_type;

	/**
	 * A hash_map that encodes the non-uniform component of the categorical
	 * distribution. It maps from atoms to pairs, where the first entry in the
	 * pair contains the probability and the second entry contains the log
	 * probability.
	 */
	hash_map<K, pair<V, V>> probabilities;

	/**
	 * The number of dimensions of the categorical distribution.
	 */
	unsigned int atom_count;

	/**
	 * Stores the probability of every atom in the uniform component of the
	 * categorical distribution (i.e. every atom that is not a key in
	 * sparse_categorical::probabilities).
	 */
	V prob;

	/**
	 * Stores the total probability mass in the non-uniform component of the
	 * categorical distribution (i.e. the sum of the probabilities in
	 * sparse_categorical::probabilities).
	 */
	V dense_prob;

	/**
	 * The natural logarithm of sparse_categorical::prob.
	 */
	V log_prob;

	/**
	 * Initializes this categorical distribution with the given dimension
	 * `atom_count`, setting all probabilities to zero.
	 */
	sparse_categorical(unsigned int atom_count) :
		probabilities(16), atom_count(atom_count),
		prob(1.0 / atom_count), dense_prob(0.0), log_prob(-log(atom_count)) { }

	/**
	 * Initializes this categorical distribution by copying its fields from the
	 * given sparse_categorical distribution `src`.
	 */
	sparse_categorical(const sparse_categorical<K, V>& src) : probabilities(src.probabilities.table.capacity),
			atom_count(src.atom_count), prob(src.prob), dense_prob(src.dense_prob), log_prob(src.log_prob)
	{
		if (!initialize(src))
			exit(EXIT_FAILURE);
	}

	~sparse_categorical() { free(); }

	/**
	 * Sets the `probability` of the given observation `key`. This function
	 * will make the key part of the non-uniform component of the categorical
	 * distribution.
	 */
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

	/**
	 * Computes the conditional probability of observing the given `item` drawn
	 * from this constant distribution, *conditioned* on a collection of
	 * observations `conditioned`. This function assumes all the elements in
	 * `conditioned` are identical, and that `item` and `conditioned` have
	 * non-zero probability according to `prior`.
	 * \returns `true` if `item` is equivalent to the first element in `conditioned`.
	 * \returns `false` otherwise.
	 * \tparam PriorDistribution the type of the prior distribution. This function does not use `prior`.
	 */
	template<typename PriorDistribution>
	static inline V conditional(const PriorDistribution& prior,
			const K& item, const array_multiset<K>& conditioned)
	{
		for (unsigned int i = 0; i < conditioned.counts.size; i++) {
			if (conditioned.counts.keys[i] == item) {
				return (prior.get_for_atom(item) + conditioned.counts.values[i])
						/ (prior.sum() + conditioned.total());
			}
		}

		return (prior.get_for_atom(item)) / (prior.sum() + conditioned.total());
	}

	/**
	 * Returns the log probability of observing the given `item`, drawn from a
	 * categorical distribution, which is itself drawn from the given `prior`
	 * distribution, *conditioned* on the set of observations `conditioned`. It
	 * is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			const K& item, const array_multiset<K>& conditioned)
	{
		for (unsigned int i = 0; i < conditioned.counts.size; i++) {
			if (conditioned.counts.keys[i] == item) {
				return log((prior.get_for_atom(item) + conditioned.counts.values[i]))
						- log(prior.sum() + conditioned.total());
			}
		}

		return log(prior.get_for_atom(item)) - log(prior.sum() + conditioned.total());
	}

	/**
	 * Returns the probability of observing the given collection of `items`,
	 * each drawn independently and identically from a categorical
	 * distribution, which is itself drawn from the given `prior` distribution,
	 * *conditioned* on the set of observations `conditioned`. It is assumed
	 * the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			const array_multiset<K>& items, const array_multiset<K>& conditioned)
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

	/**
	 * Returns the probability of the given observation `key`.
	 */
	inline V probability(const K& observation) const {
		bool contains;
		const pair<V, V>& entry = probabilities.get(observation, contains);
		if (contains) return entry.key;
		else return prob;
	}

	/**
	 * Returns the log probability of the given observation `key`.
	 */
	inline V log_probability(const K& observation) const {
		bool contains;
		const pair<V, V>& entry = probabilities.get(observation, contains);
		if (contains) return entry.value;
		else return log_prob;
	}

	/**
	 * Frees the given sparse_categorical `distribution` by releasing the
	 * memory resources associated with sparse_categorical::probabilities,
	 * along with all of its elements.
	 */
	static inline void free(sparse_categorical<K, V>& distribution) {
		distribution.free();
		core::free(distribution.probabilities);
	}

private:
	inline bool initialize(const sparse_categorical<K, V>& src) {
		for (unsigned int i = 0; i < src.probabilities.table.capacity; i++) {
			if (!is_empty(src.probabilities.table.keys[i])) {
				if (!init(probabilities.table.keys[i], src.probabilities.table.keys[i])) {
					set_empty(probabilities.table.keys[i]);
					return false;
				}
				probabilities.values[i] = src.probabilities.values[i];
				probabilities.table.size++;
			}
		}
		return true;
	}

	inline void free() {
		for (auto entry : probabilities)
			core::free(entry.key);
	}

	template<typename A, typename B>
	friend bool init(sparse_categorical<A, B>&, const sparse_categorical<A, B>&);
};

/**
 * Initializes the given sparse_categorical `distribution` by copying its
 * fields from the given sparse_categorical distribution `src`.
 */
template<typename K, typename V>
inline bool init(sparse_categorical<K, V>& distribution, const sparse_categorical<K, V>& src) {
	distribution.atom_count = src.atom_count;
	distribution.prob = src.prob;
	distribution.log_prob = src.log_prob;
	distribution.dense_prob = src.dense_prob;
	if (!hash_map_init(distribution.probabilities, src.probabilities.table.capacity)) {
		fprintf(stderr, "init ERROR: Unable to initialize hash_map in sparse_categorical.\n");
		return false;
	} else if (!distribution.initialize(src)) {
		core::free(distribution.probabilities);
		return false;
	}
	return true;
}

/**
 * Draws a sample from the given sparse_categorical `distribution` and stores
 * it in `output`. It is possible that the observation will be sampled from the
 * uniform component of the sparse categorical distribution, in which case this
 * function will print an error and return `false`.
 * \tparam K satisfies core::is_copyable.
 */
template<typename K, typename V>
inline bool sample(const sparse_categorical<K, V>& distribution, K& output)
{
	V random = sample_uniform<V>();
	if (random < distribution.dense_prob
	 || distribution.probabilities.table.size >= distribution.atom_count)
	{
		V aggregator = 0.0;
		const K* last = NULL;
		for (const auto& entry : distribution.probabilities) {
			last = &entry.key;
			aggregator += entry.value.key;
			if (random < aggregator)
				return copy(entry.key, output);
		}
		return copy(*last, output);
	} else {
		fprintf(stderr, "sample ERROR: We sampled an object "
				"not belonging to the hash_map of known objects.");
		return false;
	}
}

/**
 * This struct represents a categorical distribution, where the probabilities
 * are stored contiguously in a native array of type `V`. The
 * sparse_categorical struct instead stores only the non-uniform probabilities
 * in a core::hash_map. sparse_categorical should be used if the dimension is
 * large and the probabilities can be represented sparsely as a sum of uniform
 * and non-uniform components.
 * **NOTE:** This distribution assumes the observations have type
 * `unsigned int` and take values in `{1, ..., dense_categorical::atom_count}`.
 *
 * In the following example, we define a categorical distribution over the
 * ASCII characters (encoded as positive integers from 1 through 256).
 * ASCII characters (encoded as positive integers from 1 through 256). Positive
 * probability is only assigned to the characters <code>{'a', 'b', 'c', 'd', 'e'}</code>,
 * and we draw 10 samples from the distribution. The expected output of the
 * example is `a, c, b, b, d, e, b, e, b, b, `.
 * ```{.cpp}
 * #include <math/distributions.h>
 * using namespace core;
 *
 * template<typename V>
 * inline void set_probability(
 * 		dense_categorical<V>& categorical,
 * 		unsigned int atom, const V& probability)
 * {
 * 	categorical.phi[atom - 1] = probability;
 * }
 *
 * int main() {
 * 	set_seed(100);
 * 	dense_categorical<double> categorical(256);
 * 	set_probability(categorical, 'a', 0.1);
 * 	set_probability(categorical, 'b', 0.4);
 * 	set_probability(categorical, 'c', 0.1);
 * 	set_probability(categorical, 'd', 0.2);
 * 	set_probability(categorical, 'e', 0.2);
 *
 * 	for (unsigned int i = 0; i < 10; i++) {
 * 		unsigned int output;
 * 		sample(categorical, output);
 * 		printf("%c, ", output);
 * 	}
 * }
 * ```
 *
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
struct dense_categorical
{
	/**
	 * The type of the probabilities.
	 */
	typedef V value_type;

	/**
	 * The native array containing the probabilities of each event.
	 */
	V* phi;

	/**
	 * The number of dimensions of the categorical distribution.
	 */
	unsigned int atom_count;

	/**
	 * Initializes this categorical distribution with the given dimension
	 * `atom_count`, setting all probabilities to zero.
	 */
	dense_categorical(unsigned int atom_count) : atom_count(atom_count) {
		phi = (V*) calloc(atom_count, sizeof(V));
		if (phi == NULL) {
			fprintf(stderr, "dense_categorical ERROR: Insufficient memory for phi.\n");
			exit(EXIT_FAILURE);
		}
	}

	/**
	 * Initializes this categorical distribution by copying the fields from the
	 * given dense_categorical distribution `src`.
	 */
	dense_categorical(const dense_categorical& src) : atom_count(src.atom_count) {
		phi = (V*) malloc(atom_count * sizeof(V));
		if (phi == NULL) {
			fprintf(stderr, "dense_categorical ERROR: Insufficient memory for phi.\n");
			exit(EXIT_FAILURE);
		}
		memcpy(phi, src.phi, sizeof(V) * atom_count);
	}

	~dense_categorical() { free(); }

	/**
	 * Returns the element at the given `index` of the probability vector
	 * dense_categorical::phi. This is equivalent to the probability of
	 * observing `index + 1` drawn from this distribution.
	 */
	inline V get(unsigned int index) const {
		return phi[index];
	}

	/**
	 * Returns the element at the index `atom - 1` of the probability vector
	 * dense_categorical::phi. This is equivalent to the probability of
	 * observing `atom` drawn from this distribution.
	 */
	inline V get_for_atom(unsigned int atom) const {
		return phi[atom - 1];
	}

	/**
	 * It is assumed that the probabilities in dense_categorical::phi sum to 1,
	 * and so this function always returns 1.0.
	 */
	inline V sum() const {
		return 1.0;
	}

	void ensure_atom_count(unsigned int new_atom_count) {
		if (new_atom_count <= atom_count)
			return;
		fprintf(stderr, "dense_categorical.set_atom_count ERROR: This is not implemented.\n");
	}

	/**
	 * Returns the probability of observing `item` drawn from this distribution.
	 */
	inline V probability(unsigned int item) const {
#if !defined(NDEBUG)
		if (item == 0) {
			fprintf(stderr, "dense_categorical.conditional ERROR: Given item is zero.\n");
			return std::numeric_limits<V>::signaling_NaN();
		}
#endif
		return phi[item - 1];
	}

	/**
	 * Returns the joint probability of observing the given collection of
	 * `items`, drawn independently and identically from this distribution.
	 */
	inline V probability(const array_multiset<unsigned int>& items) const {
		V value = 1.0;
		for (unsigned int i = 0; i < items.counts.size; i++)
			value *= pow(probability(items.counts.keys[i]), items.counts.values[i]);
		return value;
	}

	/**
	 * Returns the log probability of observing `item` drawn from this distribution.
	 */
	inline V log_probability(unsigned int item) const {
		return log(phi[item - 1]);
	}

	/**
	 * Returns the joint log probability of observing the given collection of
	 * `items`, drawn independently and identically from this distribution.
	 */
	inline V log_probability(const array_multiset<unsigned int>& items) const {
		V value = 0.0;
		for (const auto& entry : items.counts)
			value = log_probability(entry.key) * entry.value;
		return value;
	}

	/**
	 * Returns the probability of observing the given `item`, drawn from a
	 * categorical distribution, which is itself drawn from the given `prior`
	 * distribution. It is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V probability(const PriorDistribution& prior, unsigned int item) {
		return prior.get_for_atom(item) / prior.sum();
	}

	/**
	 * Returns the log probability of observing the given `item`, drawn from a
	 * categorical distribution, which is itself drawn from the given `prior`
	 * distribution. It is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_probability(const PriorDistribution& prior, unsigned int item) {
		return prior.log_probability(item);
	}

	/**
	 * Returns the probability of observing the given collection of `items`,
	 * each drawn independently and identically from a categorical
	 * distribution, which is itself drawn from the given `prior` distribution.
	 * It is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_probability(const PriorDistribution& prior, const array_multiset<unsigned int>& items)
	{
		V log_probability = 0.0;
		for (unsigned int i = 0; i < items.counts.size; i++)
			log_probability += log_rising_factorial(
					prior.get_for_atom(items.counts.keys[i]), items.counts.values[i]);
		return log_probability - log_rising_factorial(prior.sum(), items.total());
	}

	/**
	 * Returns the probability of observing the given `item`, drawn from a
	 * categorical distribution, which is itself drawn from the given `prior`
	 * distribution, *conditioned* on the set of observations `conditioned`. It
	 * is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V conditional(const PriorDistribution& prior,
			unsigned int item, const array_multiset<unsigned int>& conditioned)
	{
		return sparse_categorical<unsigned int, V>::conditional(prior, item, conditioned);
	}

	/**
	 * Returns the log probability of observing the given `item`, drawn from a
	 * categorical distribution, which is itself drawn from the given `prior`
	 * distribution, *conditioned* on the set of observations `conditioned`. It
	 * is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			unsigned int item, const array_multiset<unsigned int>& conditioned)
	{
		return sparse_categorical<unsigned int, V>::log_conditional(prior, item, conditioned);
	}

	/**
	 * Returns the log probability of observing the given `item`, drawn from a
	 * categorical distribution, which is itself drawn from the given `prior`
	 * distribution, *conditioned* on the set of observations `conditioned`. It
	 * is assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			unsigned int item, const hash_multiset<unsigned int>& conditioned)
	{
		bool contains;
		unsigned int count = conditioned.counts.get(item, contains);

		if (contains) {
			return fasterlog(prior.get_for_atom(item) + count);
		} else {
			return fasterlog(prior.get_for_atom(item));
		}
	}

	/**
	 * Returns the log probability of observing the given item `holdout`, drawn
	 * from a categorical distribution, which is itself drawn from the given
	 * `prior` distribution, *conditioned* on the set of observations
	 * `conditioned` \ { `holdout` } (the set of conditioned observations with
	 * `holdout` being subtracted). It is assumed the given `prior` is a
	 * Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			unsigned int holdout, const array_multiset<unsigned int>& conditioned)
	{
		for (unsigned int i = 0; i < conditioned.counts.size; i++) {
			if (conditioned.counts.keys[i] == holdout) {
				return log(prior.get_for_atom(holdout) + conditioned.counts.values[i] - 1)
						- log(prior.sum() + conditioned.total() - 1);
			}
		}

		return log(prior.get_for_atom(holdout)) - log(prior.sum() + conditioned.total());
	}

	/**
	 * Returns the log probability of observing the given item `holdout`, drawn
	 * from a categorical distribution, which is itself drawn from the given
	 * `prior` distribution, *conditioned* on the set of observations
	 * `conditioned` \ { `holdout` } (the set of conditioned observations with
	 * `holdout` being subtracted). It is assumed the given `prior` is a
	 * Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			unsigned int holdout, const hash_multiset<unsigned int>& conditioned)
	{
		bool contains;
		unsigned int count = conditioned.counts.get(holdout, contains);

		if (contains) {
			return log(prior.get_for_atom(holdout) + count - 1);
		} else {
			return log(prior.get_for_atom(holdout));
		}
	}

	/**
	 * Returns the log probability of observing the given collection of
	 * `items`, each drawn independently and identically from a categorical
	 * distribution, which is itself drawn from the given `prior` distribution,
	 * *conditioned* on the set of observations `conditioned`. It is assumed
	 * the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			const array_multiset<unsigned int>& items, const array_multiset<unsigned int>& conditioned)
	{
		return sparse_categorical<unsigned int, V>::log_conditional(prior, items, conditioned);
	}

	/**
	 * Returns the log probability of observing the given collection of
	 * `items`, each drawn independently and identically from a categorical
	 * distribution, which is itself drawn from the given `prior` distribution,
	 * *conditioned* on the set of observations `conditioned`. It is assumed
	 * the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional(const PriorDistribution& prior,
			const array_multiset<unsigned int>& items, const hash_multiset<unsigned int>& conditioned)
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

	/**
	 * Returns the log probability of observing the given collection of items
	 * `holdout`, each drawn independently and identically from a categorical
	 * distribution, which is itself drawn from the given `prior` distribution,
	 * *conditioned* on the set of observations `conditioned` \ `holdout` (the
	 * set of conditioned observations with `holdout` being subtracted). It is
	 * assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			const array_multiset<unsigned int>& holdout, const array_multiset<unsigned int>& conditioned)
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

	/**
	 * Returns the log probability of observing the given collection of items
	 * `holdout`, each drawn independently and identically from a categorical
	 * distribution, which is itself drawn from the given `prior` distribution,
	 * *conditioned* on the set of observations `conditioned` \ `holdout` (the
	 * set of conditioned observations with `holdout` being subtracted). It is
	 * assumed the given `prior` is a Dirichlet.
	 * \tparam PriorDistribution a distribution type with public member
	 * 		functions `V get_for_atom(unsigned int)` and `V sum()`.
	 */
	template<typename PriorDistribution>
	static inline V log_conditional_without(const PriorDistribution& prior,
			const array_multiset<unsigned int>& holdout, const hash_multiset<unsigned int>& conditioned)
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

	/**
	 * Moves the dense_categorical distribution from `src` into `dst`. This
	 * function simply copies the pointers, and does not initialize a new
	 * probability array.
	 */
	static inline void move(const dense_categorical<V>& src, dense_categorical<V>& dst) {
		dst.phi = src.phi;
		dst.atom_count = src.atom_count;
	}

	/**
	 * Copies the dense_categorical distribution from `src` into `dst`. This
	 * function initializes a new probability array in `dst` and copies into it
	 * the contents from the array in `src`.
	 */
	static inline bool copy(const dense_categorical<V>& src, dense_categorical<V>& dst) {
		dst.phi = (V*) malloc(sizeof(V) * src.atom_count);
		if (dst.phi == NULL) {
			fprintf(stderr, "dense_categorical.copy ERROR: Out of memory.\n");
			return false;
		}
		memcpy(dst.phi, src.phi, src.atom_count * sizeof(V));
		dst.atom_count = src.atom_count;
		return true;
	}

	/* NOTE: this function assumes that the type V has constant size */
	template<typename Metric>
	static inline long unsigned int size_of(const dense_categorical<V>& distribution, const Metric& metric) {
		return sizeof(V) * distribution.atom_count
			 + core::size_of(distribution.atom_count);
	}

	/**
	 * Frees the given categorical distribution, by releasing the resources
	 * associated with the native array dense_categorical::phi.
	 */
	static inline void free(dense_categorical<V>& distribution) {
		distribution.free();
	}

	/**
	 * Returns the parameters that may be used to construct this distribution,
	 * using either the constructor or init. This particular function simply
	 * returns dense_categorical::atom_count.
	 */
	inline unsigned int get_parameters() const {
		return atom_count;
	}

private:
	inline void free() {
		core::free(phi);
	}
};

/**
 * Initializes the given categorical `distribution` with the given dimension
 * `atom_count`, setting all probabilities to zero.
 */
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

/**
 * Initializes the given categorical `distribution`, copying its fields from
 * the given dense_categorical distribution `src`.
 */
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

/**
 * Reads a dense_categorical `distribution` from `stream`.
 * \tparam Stream satisfies is_readable.
 */
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

/**
 * Writes the given dense_categorical `distribution` to `stream`.
 * \tparam Stream satisfies is_writeable.
 */
template<typename V, typename Stream>
bool write(const dense_categorical<V>& distribution, Stream& stream) {
	return write(distribution.atom_count, stream)
		&& write(distribution.phi, stream, distribution.atom_count);
}

/**
 * Draws a sample from the given dense_categorical `distribution` and writes it to `output`.
 */
template<typename V>
inline bool sample(const dense_categorical<V>& distribution, unsigned int& output) {
	output = sample_categorical(distribution.phi, 1.0, distribution.atom_count) + 1;
	return true;
}


/**
 * A struct that represents the [degenerate/constant distribution](https://en.wikipedia.org/wiki/Degenerate_distribution).
 * \tparam K the type of the observation, which must implement the operator `==`.
 */
template<typename K>
struct constant
{
	/**
	 * Computes the conditional probability of observing the given `item` drawn
	 * from this constant distribution, *conditioned* on a collection of
	 * observations `conditioned`. This function assumes all the elements in
	 * `conditioned` are identical, and that `item` and `conditioned` have
	 * non-zero probability according to `prior`.
	 * \returns `true` if `item` is equivalent to the first element in `conditioned`.
	 * \returns `false` otherwise.
	 * \tparam PriorDistribution the type of the prior distribution. This function does not use `prior`.
	 */
	template<typename PriorDistribution>
	static inline bool conditional(const PriorDistribution& prior,
			const K& item, const array_multiset<K>& conditioned)
	{
		return item == conditioned.counts.keys[0];
	}

	/**
	 * Computes the conditional probability of observing the given collection
	 * of `items`, each drawn independently and identically from this constant
	 * distribution, *conditioned* on a collection of observations
	 * `conditioned`. This function assumes all the elements in `items` are
	 * identical, all the elements in `conditioned` are identical, and that
	 * `items` and `conditioned` have non-zero probability according to `prior`.
	 * \returns `true` if the first element in `items` is equivalent to the first element in `conditioned`.
	 * \returns `false` otherwise.
	 * \tparam PriorDistribution the type of the prior distribution. This function does not use `prior`.
	 */
	template<typename PriorDistribution>
	static inline bool conditional(const PriorDistribution& prior,
			const array_multiset<K>& items, const array_multiset<K>& conditioned)
	{
		return items.counts.keys[0] == conditioned.counts.keys[0];
	}

	/**
	 * Computes the conditional log probability of observing the given `item`
	 * drawn from this constant distribution, *conditioned* on a collection of
	 * observations `conditioned`. This function assumes all the elements in
	 * `conditioned` are identical, and that `item` and `conditioned` have
	 * non-zero probability according to `prior`.
	 * \returns `0` if `item` is equivalent to the first element in `conditioned`.
	 * \returns `-inf` otherwise.
	 * \tparam PriorDistribution the type of the prior distribution. This function does not use `prior`.
	 */
	template<typename PriorDistribution>
	static inline double log_conditional(const PriorDistribution& prior,
			const K& item, const array_multiset<K>& conditioned)
	{
		if (item == conditioned.counts.keys[0])
			return 0.0;
		else return -std::numeric_limits<double>::infinity();
	}

	/**
	 * Computes the conditional probability of observing the given collection
	 * of `items`, each drawn independently and identically from this constant
	 * distribution, *conditioned* on a collection of observations
	 * `conditioned`. This function assumes all the elements in `items` are
	 * identical, all the elements in `conditioned` are identical, and that
	 * `items` and `conditioned` have non-zero probability according to `prior`.
	 * \returns `0` if the first element in `items` is equivalent to the first element in `conditioned`.
	 * \returns `-inf` otherwise.
	 * \tparam PriorDistribution the type of the prior distribution. This function does not use `prior`.
	 */
	template<typename PriorDistribution>
	static inline double log_conditional(const PriorDistribution& prior,
			const array_multiset<K>& items, const array_multiset<K>& conditioned)
	{
		if (items.counts.keys[0] == conditioned.counts.keys[0])
			return 0.0;
		else return -std::numeric_limits<double>::infinity();
	}

	/**
	 * Returns the probability of observing the given `item`, drawn from a
	 * constant distribution, which is itself drawn from the given `prior`
	 * distribution.
	 * \tparam PriorDistribution a distribution type that contains the typedef
	 * `value_type` implements the public member function `value_type probability(const K&)`.
	 */
	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type probability(
		const PriorDistribution& prior, const K& item)
	{
		return prior.probability(item);
	}

	/**
	 * Returns the probability of observing the given collection of `items`,
	 * each drawn independently and identically from a constant distribution,
	 * which is itself drawn from the given `prior` distribution.
	 * \tparam PriorDistribution a distribution type that contains the typedef
	 * 		`value_type` implements the public member function `value_type
	 * 		probability(const array_multiset<K>&)`.
	 */
	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type probability(
		const PriorDistribution& prior, const array_multiset<K>& items)
	{
		return prior.probability(items);
	}

	/**
	 * Returns the log probability of observing the given `item`, drawn from a
	 * constant distribution, which is itself drawn from the given `prior`
	 * distribution.
	 * \tparam PriorDistribution a distribution type that contains the typedef
	 * `value_type` implements the public member function `value_type log_probability(const K&)`.
	 */
	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type log_probability(
		const PriorDistribution& prior, const K& item)
	{
		return prior.log_probability(item);
	}

	/**
	 * Returns the log probability of observing the given collection of
	 * `items`, each drawn independently and identically from a constant
	 * distribution, which is itself drawn from the given `prior` distribution.
	 * \tparam PriorDistribution a distribution type that contains the typedef
	 * 		`value_type` implements the public member function `value_type
	 * 		log_probability(const array_multiset<K>&)`.
	 */
	template<typename PriorDistribution>
	static inline typename PriorDistribution::value_type log_probability(
		const PriorDistribution& prior, const array_multiset<K>& items)
	{
		return prior.log_probability(items);
	}

	/**
	 * Samples an observation from a constant distribution, which is itself
	 * drawn from a `prior` distribution, *conditioned* on the fact that the
	 * given set of `observations` was sampled from the constant distribution.
	 * The sample is written to `sample`. This function assumes the
	 * elements in `observations` are identical, and have non-zero probability
	 * according to `prior`.
	 * \tparam K satisfies core::is_copyable.
	 * \tparam PriorDistribution the type of the prior distribution. This function does not use `prior`.
	 */
	template<typename PriorDistribution>
	static inline bool sample(const PriorDistribution& prior,
			const array_multiset<K>& observations, K& sample)
	{
		return copy(observations.counts.keys[0], sample);
	}
};


/**
 * This struct represents a discrete uniform distribution, where probabilities
 * have type `V`.
 */
template<typename V>
struct uniform_distribution {
	/**
	 * The type of the probabilities.
	 */
	typedef V value_type;

	/**
	 * The probability of each event.
	 */
	V prob;

	/**
	 * The log probability of each event.
	 */
	V log_prob;

	/**
	 * Constructs a uniform distribution with the given dimension/number of
	 * distinct events `count`.
	 */
	uniform_distribution(unsigned int count) :
		prob(1.0 / count), log_prob(-log((V) count)) { }

	/**
	 * Returns the probability of the given `observation`, which is
	 * equivalently uniform_distribution::prob.
	 */
	template<typename K>
	inline V probability(const K& observation) const {
		return prob;
	}

	/**
	 * Returns the log probability of the given `observation`, which is
	 * equivalently uniform_distribution::log_prob.
	 */
	template<typename K>
	inline V log_probability(const K& observation) const {
		return log_prob;
	}

	static inline void free(const uniform_distribution<V>& distribution) { }
};

template<typename V>
inline bool init(uniform_distribution<V>& distribution, unsigned int count) {
	distribution.prob = 1.0 / count;
	distribution.log_prob = -log((V) count);
	return true;
}


/**
 * This struct represents a distribution over sequences of events, where each
 * event is drawn independently and identically from `ElementDistribution`,
 * until a special "stop event" is drawn with probability
 * sequence_distribution::end_probability. We assume the first event is not the
 * stop event.
 *
 * The following example constructs a sequence_distribution where
 * the ElementDistribution is a sparse_categorical distribution over the
 * characters <code>{'a', 'b', 'c', 'd', 'e'}</code>. Five samples are drawn from the
 * sequence distribution. The expected output is
 * <code>'ddbc', 'eecdecdd', 'bbbb', 'dd', 'babababc', </code>.
 *
 * ```{.cpp}
 * #include <math/distributions.h>
 * using namespace core;
 *
 * int main() {
 * 	set_seed(100);
 * 	sparse_categorical<char, double> categorical(5);
 * 	categorical.set('a', 0.1);
 * 	categorical.set('b', 0.4);
 * 	categorical.set('c', 0.1);
 * 	categorical.set('d', 0.2);
 * 	categorical.set('e', 0.2);
 *
 * 	sequence_distribution<sparse_categorical<char, double>> seq(categorical, 0.2);
 *
 * 	for (unsigned int i = 0; i < 5; i++) {
 * 		string s;
 * 		sample(seq, s);
 * 		print("'", stdout);
 * 		print(s, stdout);
 * 		print("', ", stdout);
 * 	}
 * }
 * ```
 *
 * \tparam ElementDistribution a distribution type that defines a public
 * 		typedef `value_type` that indicates the type of the probabilities. Some
 * 		operations in this class also require the public member functions
 * 		`value_type probability(const T&)`, `value_type log_probability(const T&)`.
 */
template<typename ElementDistribution>
struct sequence_distribution
{
	/**
	 * The type of the probabilities.
	 */
	typedef typename ElementDistribution::value_type V;

	/**
	 * After the first event in the sequence, this is the probability that no
	 * further events are drawn, and the sequence ends.
	 */
	V end_probability;

	/**
	 * The natural logarithm of sequence_distribution::end_probability.
	 */
	V log_end_probability;

	/**
	 * The natural logarithm of (1 - sequence_distribution::end_probability).
	 */
	V log_not_end_probability;

	/**
	 * The distribution from which each element in the sequence is drawn
	 * independently and identically.
	 */
	ElementDistribution element_distribution;

	/**
	 * Constructs a sequence_distribution with the given
	 * sequence_distribution::element_distribution and
	 * sequence_distribution::end_probability.
	 * \tparam ElementDistribution a distribution type that can be constructed
	 * 		using a single parameter with type ElementDistribution&.
	 */
	sequence_distribution(ElementDistribution& element_distribution, V end_probability) :
		end_probability(end_probability),
		log_end_probability(log(end_probability)),
		log_not_end_probability(log(1.0 - end_probability)),
		element_distribution(element_distribution)
	{ }

	/**
	 * Returns the probability of the given observation `sequence` under this
	 * sequence distribution.
	 * \tparam SequenceType a sequence type that contains the integral field
	 * 		`length` that indicates the number of elements in the sequence, and
	 * 		implements the operator `[]` that returns an element at a
	 * 		particular index of the sequence.
	 */
	template<typename SequenceType>
	inline V probability(const SequenceType& sequence) const {
		if (sequence.length == 0) return 0.0;
		V product = element_distribution.probability(sequence[0]);
		for (unsigned int i = 1; i < sequence.length; i++)
			product *= element_distribution.probability(sequence[i]) * (1.0 - end_probability);
		return product * end_probability;
	}

	/**
	 * Returns the log probability of the given observation `sequence` under
	 * this sequence distribution.
	 * \tparam SequenceType a sequence type that contains the integral field
	 * 		`length` that indicates the number of elements in the sequence, and
	 * 		implements the operator `[]` that returns an element at a
	 * 		particular index of the sequence.
	 */
	template<typename SequenceType>
	inline V log_probability(const SequenceType& sequence) const {
		if (sequence.length == 0) return -std::numeric_limits<V>::infinity();
		V sum = element_distribution.log_probability(sequence[0]);
		for (unsigned int i = 1; i < sequence.length; i++)
			sum += element_distribution.log_probability(sequence[i]) + log_not_end_probability;
		return sum + log_end_probability;
	}

	/**
	 * Copies the sequence_distribution `src` into `dst`.
	 * \tparam ElementDistribution satisfies core::is_copyable.
	 */
	static inline bool copy(
			const sequence_distribution<ElementDistribution>& src,
			sequence_distribution<ElementDistribution>& dst)
	{
		dst.end_probability = src.end_probability;
		dst.log_end_probability = src.log_end_probability;
		dst.log_not_end_probability = src.log_not_end_probability;
		return core::copy(src.element_distribution, dst.element_distribution);
	}

	/**
	 * Frees the given sequence `distribution`.
	 * \tparam ElementDistribution the type of the distribution of each element
	 * 		in the sequence. The function `core::free` is called with an
	 * 		argument of type `ElementDistribution&`.
	 */
	static inline void free(sequence_distribution<ElementDistribution>& distribution) {
		core::free(distribution.element_distribution);
	}
};

/**
 * Initializes the given sequence_distribution `distribution` with the given
 * sequence_distribution `src`.
 * \tparam ElementDistribution a distribution type for which the function
 * 		`bool init(ElementDistribution&, const ElementDistribution)` is implemented.
 */
template<typename ElementDistribution>
inline bool init(sequence_distribution<ElementDistribution>& distribution,
		const sequence_distribution<ElementDistribution>& src)
{
	distribution.end_probability = src.end_probability;
	distribution.log_end_probability = src.log_end_probability;
	distribution.log_not_end_probability = src.log_not_end_probability;
	return init(distribution.element_distribution, src.element_distribution);
}

/**
 * Reads the given sequence_distribution `distribution` from `stream`.
 * \tparam Stream satisfies is_readable.
 */
template<typename ElementDistribution, typename Stream>
bool read(sequence_distribution<ElementDistribution>& distribution, Stream& stream)
{
	if (!read(distribution.end_probability, stream))
		return false;
	distribution.log_end_probability = log(distribution.end_probability);
	distribution.log_not_end_probability = log(1.0 - distribution.end_probability);
	return read(distribution.element_distribution, stream);
}

/**
 * Writes the given sequence_distribution `distribution` to `stream`.
 * \tparam Stream satisfies is_readable.
 */
template<typename ElementDistribution, typename Stream>
bool write(const sequence_distribution<ElementDistribution>& distribution, Stream& stream)
{
	return write(distribution.end_probability, stream)
		&& write(distribution.element_distribution, stream);
}

/**
 * Samples from the given sequence_distribution `distribution` and stores the result in `output`.
 * \tparam SequenceType a sequence type for which the function
 * 		`bool init(SequenceType&, unsigned int)` is implemented, which
 * 		initializes the sequence with the given length. The operator `[]` must
 * 		also be implemented, which returns a reference to the specified index
 * 		in the sequence with type `T`, so that the function
 * 		`bool sample(const ElementDistribution&, T&)` can be called. Finally,
 * 		the function `core::free` may be called on the elements of the
 * 		sequence, if the sampling fails.
 */
template<typename ElementDistribution, typename SequenceType>
bool sample(const sequence_distribution<ElementDistribution>& distribution, SequenceType& output)
{
	/* first sample the length of the sequence */
	unsigned int length = sample_geometric(distribution.end_probability) + 1;
	if (!init(output, length)) return false;
	for (unsigned int i = 0; i < length; i++) {
		if (!sample(distribution.element_distribution, output[i])) {
			for (unsigned int j = 0; j < i; j++)
				free(output[j]);
			return false;
		}
	}
	return true;
}

#endif /* DISTRIBUTIONS_H_ */
