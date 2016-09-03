/**
 * log.h - Contains useful functions for performing arithmetic in log space.
 *
 *  Created on: Aug 1, 2015
 *      Author: asaparov
 */

#ifndef LOG_H_
#define LOG_H_

#include <math.h>
#include <mutex>

/* avoid floating-point infinities for performance */
#define LOG_ZERO -1.0e16

inline float log_probability(float src) { return src; }
inline double log_probability(double src) { return src; }

template<typename T>
struct default_value_type { typedef double type; };

template<>
struct default_value_type<float> { typedef float type; };

template<typename T, typename V>
inline void max(V& max, const T* src, unsigned int length, unsigned int skip)
{
	for (unsigned int i = 0; i < length; i++)
		if (log_probability(src[i * skip]) > max)
			max = log_probability(src[i * skip]);
}

template<typename T, typename P, typename V>
inline void max(V& max, const T* src, const P* prior, unsigned int length, unsigned int skip)
{
	for (unsigned int i = 0; i < length; i++) {
		V value = log_probability(src[i * skip]) + log_probability(prior[i]);
		if (value > max) max = value;
	}
}

template<typename T, typename V = typename default_value_type<T>::type>
inline V max(const T* src, unsigned int length, unsigned int skip)
{
	V maximum = log_probability(src[0]);
	if (length > 1)
		max(maximum, src + skip, length - 1, skip);
	return maximum;
}

template<typename T, typename P, typename V = typename default_value_type<T>::type>
inline V max(const T* src, const P* prior,
		unsigned int length, unsigned int skip)
{
	V maximum = src[0] + prior[0];
	if (length > 1)
		max(maximum, src + skip, prior + 1, length - 1, skip);
	return maximum;
}

template<typename T>
inline auto max(const T* src, unsigned int length) -> decltype(max(src, 1, length)) {
	return max(src, length, 1);
}

template<typename T, typename V>
inline V sumexp(const T* src, unsigned int length, unsigned int skip, const V& shift) {
	V sum = 0.0;
	for (unsigned int i = 0; i < length; i++)
		sum += exp(log_probability(src[i * skip]) - shift);
	return sum;
}

template<typename T, typename P, typename V>
inline V sumexp(const T* src, const P* prior, unsigned int length, unsigned int skip, const V& shift) {
	V sum = 0.0;
	for (unsigned int i = 0; i < length; i++)
		sum += exp(log_probability(src[i * skip]) + log_probability(prior[i]) - shift);
	return sum;
}

template<typename T, typename V>
inline void normalize(const T* src, V* dst, unsigned int length, unsigned int skip, const V& normalization)
{
	for (unsigned int i = 0; i < length; i++)
		dst[i * skip] = log_probability(src[i * skip]) - normalization;
}

template<typename T, typename V>
inline void normalize_exp(const T* src, V* dst,
		unsigned int length, unsigned int skip, const V& normalization)
{
	for (unsigned int i = 0; i < length; i++)
		dst[i * skip] = exp(log_probability(src[i * skip]) - normalization);
}

template<typename T, typename P, typename V>
inline void normalize_exp(const T* src, const P* prior, V* dst,
		unsigned int length, unsigned int skip, const V& normalization)
{
	for (unsigned int i = 0; i < length; i++)
		dst[i * skip] = exp(log_probability(src[i * skip]) + log_probability(prior[i]) - normalization);
}

template<typename T, typename V = typename default_value_type<T>::type>
inline V logsumexp(const T* src, unsigned int length, unsigned int skip = 1) {
	V maximum = max(src, length);
	return log(sumexp(src, length, skip, maximum)) + maximum;
}

template<typename T, typename P, typename V = typename default_value_type<T>::type>
inline V logsumexp(const T* src, const P* prior, unsigned int length, unsigned int skip = 1) {
	V maximum = max(src, prior, length);
	return log(sumexp(src, prior, length, skip, maximum)) + maximum;
}

template<typename T, typename V>
void normalize(const T* src, V* dst, unsigned int length, unsigned int skip = 1) {
	normalize(src, dst, length, skip, logsumexp(src, length, skip));
}

template<typename T, typename V>
void normalize_exp(const T* src, V* dst, unsigned int length, unsigned int skip = 1) {
	normalize_exp(src, dst, length, skip, logsumexp(src, length, skip));
}

template<typename T, typename P, typename V>
void normalize_exp(const T* src, const P* prior, V* dst, unsigned int length, unsigned int skip = 1) {
	normalize_exp(src, prior, dst, length, skip, logsumexp(src, prior, length, skip));
}

template<typename T>
inline void normalize_exp(T* x, unsigned int length) {
	normalize_exp(x, x, length, 1);
}

template<typename V>
inline V logsumexp(const V& first, const V& second) {
	if (first > second)
		return first + log(1 + exp(second - first));
	else return second + log(1 + exp(first - second));
}

template<typename V>
inline V logdiffexp(const V& first, const V& second) {
	if (first > second)
		return first + log(1 - exp(second - first));
	else return second + log(-1 + exp(first - second));
}

/* caches logarithms of positive integers */
template<typename V>
struct log_cache {
	V* values;
	unsigned int size;
	std::mutex lock;

	log_cache(unsigned int initial_size) : values(NULL), size(0) {
		if (!resize<true>(initial_size))
			exit(EXIT_FAILURE);
	}

	inline bool ensure_size(unsigned int requested_size) {
		if (requested_size > size) {
			lock.lock();
			bool result = resize(requested_size);
			lock.unlock();
			return result;
		}
		return true;
	}

	inline V get(unsigned int x) const {
		return values[x];
	}

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
