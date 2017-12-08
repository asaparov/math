/**
 * \file log.h
 * 
 * Contains useful functions for performing arithmetic in log space, such as
 * computing the maximum of an array, computing sum-exp, log-sum-exp, and
 * normalization. This file also contains log_cache, which caches the
 * logarithms of positive integers. All functions in this file avoid loss of
 * precision by performing operations in log space when feasible.
 *
 * In the following example, the array
 * `{ -100.0, -100.0, -100.0 + log(2.0), -100.0 + log(4.0) }` of log
 * probabilities is normalized and printed. The expected output is
 * `[0.125000, 0.125000, 0.250000, 0.500000]`.
 *
 * ```{.cpp}
 * #include <math/log.h>
 * #include <core/io.h>
 * using namespace core;
 *
 * int main() {
 * 	double a[] = { -100.0, -100.0, -100.0 + log(2.0), -100.0 + log(4.0) };
 * 	normalize_exp(a, 4);
 * 	print(a, 4, stdout);
 * }
 * ```
 *
 *  <!-- Created on: Aug 1, 2015
 *           Author: asaparov -->
 */

#ifndef LOG_H_
#define LOG_H_

#include <core/array.h>

#include <math.h>
#include <mutex>

/**
 * This value is used for `log(0)` in some operations to avoid floating-point
 * infinities and improve performance.
 */
#define LOG_ZERO -1.0e16

/**
 * Returns `src`. This function is defined so that the function in this file
 * `max`, `sum_exp`, and `normalize` behave correctly when `T` is `float`.
 */
inline float log_probability(float src) { return src; }

/**
 * Returns `src`. This function is defined so that the function in this file
 * `max`, `sum_exp`, and `normalize` behave correctly when `T` is `double`.
 */
inline double log_probability(double src) { return src; }

template<typename T>
struct default_value_type { typedef double type; };

template<>
struct default_value_type<float> { typedef float type; };

/**
 * Computes the maximum log probability of the elements in `src` with the given
 * `length`, only inspecting the elements at every `skip` locations, until
 * `length` elements have been considered. For example, if `skip = 2`, this
 * function will find the maximum of the elements at every *even* index. The
 * computed maximum is stored in `max`.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V>
inline void max(V& max, const T* src, unsigned int length, unsigned int skip)
{
	for (unsigned int i = 0; i < length; i++)
		if (log_probability(src[i * skip]) > max)
			max = log_probability(src[i * skip]);
}

/**
 * Computes the maximum sum of the log probability of the elements in `src` and
 * the elements in `prior` with the given `length`, only inspecting the
 * elements in `src` at every `skip` locations, until `length` elements have
 * been considered. For example, if `skip = 2`, this function will find the
 * maximum of the elements at every *even* index. Note that `skip` does not
 * affect the inspection of elements in `prior`. The computed maximum is stored
 * in `max`.
 * 
 * So more precisely, this function computes:
 * \f[ max_i \{ f(\text{src}[i \cdot \text{skip}]) + g(\text{prior}[i]) \}. \f]
 * where \f$ f(\cdot) \f$ is the function `V log_probability(const T&)` and
 * \f$ g(\cdot) \f$ is the function `V log_probability(const P&)`.
 *
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam P the type of every element in `prior`, for which the function `V log_probability(const P&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename P, typename V>
inline void max(V& max, const T* src, const P* prior, unsigned int length, unsigned int skip)
{
	for (unsigned int i = 0; i < length; i++) {
		V value = log_probability(src[i * skip]) + log_probability(prior[i]);
		if (value > max) max = value;
	}
}

/**
 * Returns the maximum log probability of the elements in `src` with the given
 * `length`, only inspecting the elements at every `skip` locations, until
 * `length` elements have been considered. For example, if `skip = 2`, this
 * function will find the maximum of the elements at every *even* index.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities. By default, if `T` is a floating-point type, `V` is the same type.
 */
template<typename T, typename V = typename default_value_type<T>::type>
inline V max(const T* src, unsigned int length, unsigned int skip)
{
	V maximum = log_probability(src[0]);
	if (length > 1)
		max(maximum, src + skip, length - 1, skip);
	return maximum;
}

/**
 * Returns the maximum sum of the log probability of the elements in `src` and
 * the elements in `prior` with the given `length`, only inspecting the
 * elements in `src` at every `skip` locations, until `length` elements have
 * been considered. For example, if `skip = 2`, this function will find the
 * maximum of the elements at every *even* index. Note that `skip` does not
 * affect the inspection of elements in `prior`. The computed maximum is stored
 * in `max`.
 *
 * So more precisely, this function computes:
 * \f[ max_i \{ f(\text{src}[i \cdot \text{skip}]) + g(\text{prior}[i]) \}. \f]
 * where \f$ f(\cdot) \f$ is the function `V log_probability(const T&)` and
 * \f$ g(\cdot) \f$ is the function `V log_probability(const P&)`.
 *
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam P the type of every element in `prior`, for which the function `V log_probability(const P&)` is defined.
 * \tparam V the type of the log probabilities. By default, if `T` is a floating-point type, `V` is the same type.
 */
template<typename T, typename P, typename V = typename default_value_type<T>::type>
inline V max(const T* src, const P* prior,
		unsigned int length, unsigned int skip)
{
	V maximum = src[0] + prior[0];
	if (length > 1)
		max(maximum, src + skip, prior + 1, length - 1, skip);
	return maximum;
}

/**
 * Returns the maximum log probability of the given native array `src`.
 * \tparam T the type of every element in `src`, for which the function
 * 		`V log_probability(const T&)` is defined, where `V` is the return type
 * 		of this function.
 */
template<typename T>
inline auto max(const T* src, unsigned int length) -> decltype(max(src, 1, length)) {
	return max(src, length, 1);
}

/**
 * Returns the maximum log probability of the given core::array `src`.
 * \tparam T the type of every element in `src`, for which the function
 * 		`V log_probability(const T&)` is defined, where `V` is the return type
 * 		of this function.
 */
template<typename T>
inline auto max(const core::array<T>& src) -> decltype(max(src.data, src.length)) {
	return max(src.data, src.length);
}

/**
 * This function returns
 * \f[ \sum_{i=0}^{\text{length}-1} \exp\{ f(\text{src}[i * \text{skip}]) - \text{shift} \} \f]
 * where \f$ f(\cdot) \f$ is the function `V log_probability(const T&)`.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V>
inline V sumexp(const T* src, unsigned int length, unsigned int skip, const V& shift) {
	V sum = 0.0;
	for (unsigned int i = 0; i < length; i++)
		sum += exp(log_probability(src[i * skip]) - shift);
	return sum;
}

/**
 * This function returns
 * \f[ \sum_{i=0}^{\text{length}-1} \exp\{ f(\text{src}[i * \text{skip}]) + g(\text{prior}[i]) - \text{shift} \} \f]
 * where \f$ f(\cdot) \f$ is the function `V log_probability(const T&)` and
 * \f$ g(\cdot) \f$ is the function `V log_probability(const P&)`.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam P the type of every element in `prior`, for which the function `V log_probability(const P&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename P, typename V>
inline V sumexp(const T* src, const P* prior, unsigned int length, unsigned int skip, const V& shift) {
	V sum = 0.0;
	for (unsigned int i = 0; i < length; i++)
		sum += exp(log_probability(src[i * skip]) + log_probability(prior[i]) - shift);
	return sum;
}

/**
 * Copies the log probabilities from `src` into `dst`, subtracting
 * `normalization` from each element. This function assumes `dst` is
 * initialized and has sufficient capacity. This function only considers
 * elements at every `skip` positions, until `length` elements have been
 * considered. For example, if `skip = 2`, only the elements at even indices
 * are copied. Note that `src` and `dst` can be the same.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V>
inline void normalize(const T* src, V* dst, unsigned int length, unsigned int skip, const V& normalization)
{
	for (unsigned int i = 0; i < length; i++)
		dst[i * skip] = log_probability(src[i * skip]) - normalization;
}

/**
 * Copies the log probabilities from `src` into `dst`, subtracting
 * `normalization` from each element, and taking the natural exponent of the
 * result. This function assumes `dst` is initialized and has sufficient
 * capacity. This function only considers elements at every `skip` positions,
 * until `length` elements have been considered.
 * For example, if `skip = 2`, only the elements at even indices are copied.
 * Note that `src` and `dst` can be the same.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V>
inline void normalize_exp(const T* src, V* dst,
		unsigned int length, unsigned int skip, const V& normalization)
{
	for (unsigned int i = 0; i < length; i++)
		dst[i * skip] = exp(log_probability(src[i * skip]) - normalization);
}

/**
 * This function computes
 * \f[ \exp\{ f(\text{src}[i * \text{skip}]) + g(\text{prior}[i]) - \text{normalization} \} \f]
 * and stores the result in `dst[i]`, where `i` ranges from `0` to
 * `length - 1`, and \f$ f(\cdot) \f$ is the function
 * `V log_probability(const T&)` and \f$ g(\cdot) \f$ is the function
 * `V log_probability(const P&)`. Note that `src` and `dst` can be the same.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam P the type of every element in `prior`, for which the function `V log_probability(const P&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename P, typename V>
inline void normalize_exp(const T* src, const P* prior, V* dst,
		unsigned int length, unsigned int skip, const V& normalization)
{
	for (unsigned int i = 0; i < length; i++)
		dst[i * skip] = exp(log_probability(src[i * skip]) + log_probability(prior[i]) - normalization);
}

/**
 * Returns the natural logarithm of the sum of the natural exponent of the
 * elements in `src`, skipping every `skip` elements, until `length` elements
 * have been considered.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V = typename default_value_type<T>::type>
inline V logsumexp(const T* src, unsigned int length, unsigned int skip = 1) {
	V maximum = max(src, length);
	return log(sumexp(src, length, skip, maximum)) + maximum;
}

/**
 * Returns the natural logarithm of the sum of the natural exponent of the
 * elements in `src`, with the given `prior`, skipping every `skip` elements,
 * until `length` elements have been considered. Note that `skip` does not
 * affect the inspection of the elements in `prior`.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam P the type of every element in `prior`, for which the function `V log_probability(const P&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename P, typename V = typename default_value_type<T>::type>
inline V logsumexp(const T* src, const P* prior, unsigned int length, unsigned int skip = 1) {
	V maximum = max(src, prior, length);
	return log(sumexp(src, prior, length, skip, maximum)) + maximum;
}

/**
 * Computes the *normalized log probabilities* of the elements in `src` and
 * stores them in `dst`, skipping every `skip` elements, until `length`
 * elements have been considered. Note that `src` and `dst` can be the same.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V>
void normalize(const T* src, V* dst, unsigned int length, unsigned int skip = 1) {
	normalize(src, dst, length, skip, logsumexp(src, length, skip));
}

/**
 * Computes the *normalized probabilities* of the elements in `src` and stores
 * them in `dst`, skipping every `skip` elements, until `length` elements have
 * been considered. Note that `src` and `dst` can be the same.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename V>
void normalize_exp(const T* src, V* dst, unsigned int length, unsigned int skip = 1) {
	normalize_exp(src, dst, length, skip, logsumexp(src, length, skip));
}

/**
 * Computes the *normalized probabilities* of the elements in `src` and stores
 * them in `dst`, skipping every `skip` elements, until `length` elements have
 * been considered. Note that `skip` does not affect the inspection of the
 * elements in `prior`. Also note that `src` and `dst` can be the same.
 * \tparam T the type of every element in `src`, for which the function `V log_probability(const T&)` is defined.
 * \tparam P the type of every element in `prior`, for which the function `V log_probability(const P&)` is defined.
 * \tparam V the type of the log probabilities.
 */
template<typename T, typename P, typename V>
void normalize_exp(const T* src, const P* prior, V* dst, unsigned int length, unsigned int skip = 1) {
	normalize_exp(src, prior, dst, length, skip, logsumexp(src, prior, length, skip));
}

/**
 * Normalizes the elements in `src`. Like all functions in this file, loss of
 * precision is avoided by performing operations in log space when feasible.
 */
template<typename T>
inline void normalize_exp(T* x, unsigned int length) {
	normalize_exp(x, x, length, 1);
}

/**
 * Returns the natural logarithm of the sum of the natural exponents of
 * `first` and `second`.
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
inline V logsumexp(const V& first, const V& second) {
	if (first > second)
		return first + log(1 + exp(second - first));
	else return second + log(1 + exp(first - second));
}

/**
 * Returns the natural logarithm of `exp(first) - exp(second)`.
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
inline V logdiffexp(const V& first, const V& second) {
	if (first > second)
		return first + log(1 - exp(second - first));
	else return second + log(-1 + exp(first - second));
}

/**
 * This structure caches natural logarithms of positive integers. An
 * [std::mutex](http://en.cppreference.com/w/cpp/thread/mutex) is used to allow
 * safe simultaneous use of this cache from multiple threads.
 *
 * The following example displays a typical use-case of this cache. The
 * expected output is
 * `log(2) = 0.693147, log(183) = 5.209486, sum from i=1 to 999 of log(i) = 5905.220423`.
 *
 * ```{.cpp}
 * #include <math/log.h>
 * using namespace core;
 *
 * int main() {
 * 	log_cache<double> cache(1000);
 * 	printf("log(2) = %lf, ", cache.get(2));
 * 	printf("log(183) = %lf, ", cache.get(183));
 *
 * 	double sum = 0.0;
 * 	for (unsigned int i = 1; i < 1000; i++)
 * 		sum += cache.get(i);
 * 	printf("sum from i=1 to 999 of log(i) = %lf\n", sum);
 * }
 * ```
 *
 * \tparam V satisfies [is_arithmetic](http://en.cppreference.com/w/cpp/types/is_arithmetic).
 */
template<typename V>
struct log_cache {
	V* values;
	unsigned int size;
	std::mutex lock;

	/**
	 * Constructs the cache with the given `initial_size`. The natural
	 * logarithms up to and includig `log(initial_size - 1)` are precomputed
	 * and stored.
	 */
	log_cache(unsigned int initial_size) : values(NULL), size(0) {
		if (!resize<true>(initial_size))
			exit(EXIT_FAILURE);
	}

	/**
	 * Checks that the natural logarithms up to and including
	 * `log(requested_size - 1)` are computed and stored.
	 */
	inline bool ensure_size(unsigned int requested_size) {
		if (requested_size > size) {
			lock.lock();
			bool result = resize(requested_size);
			lock.unlock();
			return result;
		}
		return true;
	}

	/**
	 * Returns the natural logarithm of `x`. This function assumes that
	 * `log(x)` has already been computed, either by the constructor or by
	 * ensure_size, and no bounds checking is performed.
	 */
	inline V get(unsigned int x) const {
		return values[x];
	}

	/**
	 * Returns a static instance of a log_cache, with type `V`.
	 */
	static inline log_cache<V>& instance() {
		static log_cache<V> cache(16);
		return cache;
	}

private:
	template<bool FirstResize = false>
	inline bool resize(unsigned int new_size) {
		V* new_values = (V*) realloc(values, sizeof(V) * new_size);
		if (new_values == NULL) {
			fprintf(stderr, "log_cache.resize ERROR: Out of memory.\n");
			return false;
		}
		values = new_values;
		if (FirstResize) {
			values[0] = LOG_ZERO;
			size = 1;
		}
		for (unsigned int i = size; i < new_size; i++)
			values[i] = log(i);
		size = new_size;
		return true;
	}
};

#endif /* LOG_H_ */
