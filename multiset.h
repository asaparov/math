/**
 * \file multiset.h
 *
 * This file contains implementations for array_multiset and hash_multiset,
 * which are structures that count the number of occurrences of each distinct
 * element in a set. array_multiset implements this using a core::array_map,
 * where the elements are the keys, and the frequencies are the values with
 * type `unsigned int`. hash_multiset implements this using a core::hash_map.
 *
 * In the following example, both an array_multiset and a hash_multiset are
 * constructed. The same set of elements are added to both, and the two are
 * printed to stdout. The expected output is
 * `{e:1, i:2, r:0, c:3, q:1}  {c:3, e:1, i:2, q:1, r:0}`.
 *
 * ```{.cpp}
 * #include <math/multiset.h>
 * using namespace core;
 *
 * template<typename Multiset>
 * void add_elements(Multiset& m) {
 * 	m.add('c'); m.add('c');
 * 	m.add('e'); m.add('q');
 * 	m.add('r'); m.add('i');
 * 	m.add('i'); m.add('c');
 * 	m.remove('r');
 * }
 *
 * int main() {
 * 	hash_multiset<char> first(16);
 * 	add_elements(first);
 * 	print(first, stdout); print("  ", stdout);
 *
 * 	array_multiset<char> second(8);
 * 	add_elements(second);
 * 	print(second, stdout);
 * }
 * ```
 *
 *  <!-- Created on: Dec 14, 2015
 *           Author: asaparov -->
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <core/map.h>
#include <core/io.h>

using namespace core;

/**
 * A multiset structure that keeps track of the number of occurrences of
 * distinct element in a set, implemented using a core::array_map, where the
 * elements are the keys, and the values are their frequencies. The entries in
 * the underlying core::array_map are sorted according to the keys.
 *
 * hash_multiset implements the same abstract data type using a core::hash_map
 * and should be used if the number of distinct elements is expected to be
 * large.
 *
 * \tparam T satisfies [LessThanComparable](https://en.cppreference.com/w/cpp/named_req/LessThanComparable)
 * 		and [CopyAssignable](https://en.cppreference.com/w/cpp/named_req/CopyAssignable).
 */
template<typename T>
struct array_multiset {
	/**
	 * The underlying array_map.
	 */
	array_map<T, unsigned int> counts;

	/**
	 * The sum of the values in array_multiset::counts (i.e. the total number
	 * of occurrences of all elements).
	 */
	unsigned int sum;

	/**
	 * Constructs an empty array_multiset with the given `initial_capacity` for
	 * the underlying array_map array_multiset::counts.
	 */
	explicit array_multiset(unsigned int initial_capacity) : counts(initial_capacity), sum(0) { }
	~array_multiset() { free(); }

	/**
	 * Returns the total number of occurrences of all elements (array_multiset::sum).
	 */
	inline unsigned int total() const {
		return sum;
	}

	/**
	 * Adds the given `item` to the multiset with frequency `count`. The given
	 * `start_index` may be provided to improve performance. However, this
	 * function assumes that if `item` already exists in this multiset, its
	 * underlying index in array_multiset::counts is not smaller than
	 * `start_index`.
	 *
	 * This function sorts array_multiset::counts. If the user wishes to add
	 * multiple elements, they should consider using add_unsorted or the
	 * overload of add that takes an array_multiset as its parameter.
	 */
	inline bool add(const T& item, unsigned int count = 1, unsigned int start_index = 0) {
		if (!add_unsorted(item, count, start_index))
			return false;
		if (counts.size > 1)
			sort(counts.keys, counts.values, (unsigned int) counts.size, default_sorter()); /* TODO: test performance of different sorts */
		return true;
	}

	/**
	 * Adds the given `item` to the multiset with frequency `count`. The given
	 * `start_index` may be provided to improve performance. However, this
	 * function assumes that if `item` already exists in this multiset, its
	 * underlying index in array_multiset::counts is not smaller than
	 * `start_index`.
	 *
	 * This function leaves array_multiset::counts unsorted, and so the user
	 * must sort array_multiset::counts after all elements are added.
	 */
	bool add_unsorted(const T& item, unsigned int count = 1, unsigned int start_index = 0) {
		unsigned int index = counts.index_of(item, start_index);
		if (index < counts.size) {
			counts.values[index] += count;
			sum += count;
			return true;
		}
		if (!counts.ensure_capacity((unsigned int) counts.size + 1)) {
			fprintf(stderr, "array_multiset.add WARNING: Unable to expand array_map.\n");
			return false;
		}
		counts.keys[counts.size] = item;
		counts.values[counts.size] = count;
		counts.size++;
		sum += count;
		return true;
	}

	/**
	 * Adds the given multiset `items` to this multiset.
	 */
	bool add(const array_multiset<T>& items) {
		unsigned int i = 0, j = 0;
		unsigned int new_length = (unsigned int) counts.size;
		while (i < counts.size && j < items.counts.size) {
			if (counts.keys[i] == items.counts.keys[j]) {
				counts.values[i] += items.counts.values[j];
				i++; j++;
			} else if (counts.keys[i] < items.counts.keys[j]) {
				i++;
			} else {
				if (!counts.ensure_capacity(new_length + 1)) {
					fprintf(stderr, "array_multiset.add ERROR: Unable to expand underlying map.\n");
					return false;
				}
				counts.keys[new_length] = items.counts.keys[j];
				counts.values[new_length] = items.counts.values[j];
				new_length++; j++;
			}
		}

		/* add the leftover elements from 'items' */
		if (!counts.ensure_capacity(new_length + (unsigned int) items.counts.size - j)) {
			fprintf(stderr, "array_multiset.add ERROR: Unable to expand array_map.\n");
			return false;
		}
		while (j < items.counts.size) {
			counts.keys[new_length] = items.counts.keys[j];
			counts.values[new_length] = items.counts.values[j];
			new_length++; j++;
		}

		counts.size = new_length;
		if (counts.size > 1)
			sort(counts.keys, counts.values, (unsigned int) counts.size, default_sorter());
		sum += items.sum;
		return true;
	}

	/**
	 * Removes the given `item` from the multiset. The given `start_index` may
	 * be provided to improve performance. However, this function assumes that
	 * `item` exists in the multiset with non-zero frequency, and the
	 * underlying index of `item` in array_multiset::counts is not smaller than
	 * `start_index`.
	 */
	unsigned int subtract(const T& item, unsigned int start_index = 0)
	{
		unsigned int index = counts.index_of(item, start_index);
		if (index < counts.size) {
#if !defined(NDEBUG)
			if (counts.values[index] == 0) {
				fprintf(stderr, "array_multiset.subtract WARNING: Attempted "
						"to remove more items from a bin than it contains.\n");
			} else {
				counts.values[index]--;
				sum--;
			}
#else
			counts.values[index]--;
			sum--;
#endif
			return index;
		}

		fprintf(stderr, "array_multiset.subtract WARNING: No such item.\n");
		return (unsigned int) counts.size;
	}

	/**
	 * Removes the given multiset `items` from this multiset. This function
	 * assumes that the given multiset `items` is a subset of this multiset
	 * (i.e. this multiset contains all the keys in `items` with frequencies
	 * at least as large).
	 */
	void subtract(const array_multiset<T>& items)
	{
		unsigned int i = 0, j = 0;
		while (i < counts.size && j < items.counts.size) {
			if (counts.keys[i] == items.counts.keys[j]) {
#if !defined(NDEBUG)
				if (counts.values[i] < items.counts.values[j]) {
					fprintf(stderr, "array_multiset.subtract WARNING: Attempted "
							"to remove more items from a bin than it contains.\n");
					counts.values[i] = 0;
				} else counts.values[i] -= items.counts.values[j];
#else
				counts.values[i] -= items.counts.values[j];
#endif
				i++; j++;
			} else i++;
		}

#if !defined(NDEBUG)
		if (j != items.counts.size)
			fprintf(stderr, "array_multiset.subtract WARNING: Missing bin in multiset.\n");
#endif
		sum -= items.sum;
	}

	/**
	 * Removes all elements from this multiset. This function also frees the
	 * keys in the underlying array_map by calling core::free on each element.
	 */
	inline void clear() {
		for (auto entry : counts)
			core::free(entry.key);
		counts.clear();
		sum = 0;
	}

	/**
	 * This function removes keys from array_multiset::counts whose recorded
	 * frequencies are `0`. If the user is frequently removing items from this
	 * multiset, this function should be called periodically to clean up the
	 * empty entries.
	 */
	void remove_zeros(unsigned int start_index = 0) {
		unsigned int new_length = start_index;
		for (unsigned int i = start_index; i < counts.size; i++) {
			if (counts.values[i] != 0) {
				core::move(counts.keys[i], counts.keys[new_length]);
				counts.values[new_length] = counts.values[i];
				new_length++;
			} else {
				core::free(counts.keys[i]);
			}
		}
		counts.size = new_length;
	}

	/**
	 * Returns true if and only if the underlying array_map
	 * array_multiset::counts is sorted.
	 */
	inline bool is_sorted() const {
		for (unsigned int i = 1; i < counts.size; i++)
			if (counts.keys[i] < counts.keys[i - 1])
				return false;
		return true;
	}

	/**
	 * Returns the hash function evaluated on the given array_multiset `key`.
	 */
	static inline unsigned int hash(const array_multiset<T>& key) {
		return default_hash(key.counts.keys, key.counts.size)
				^ default_hash(key.counts.values, key.counts.size);
	}

	/**
	 * Moves the array_multiset `src` into `dst`. Note that this function
	 * merely copies pointers, and not contents.
	 */
	static inline void move(const array_multiset<T>& src, array_multiset<T>& dst) {
		array_map<T, unsigned int>::move(src.counts, dst.counts);
		dst.sum = src.sum;
	}

	/**
	 * Swaps the contents of the given array_multisets `first` and `second`.
	 */
	static inline void swap(array_multiset<T>& first, array_multiset<T>& second) {
		core::swap(first.counts, second.counts);
		core::swap(first.sum, second.sum);
	}

	/**
	 * Copies the contents of the array_multiset `src` into `dst`.
	 */
	static inline bool copy(const array_multiset<T>& src, array_multiset<T>& dst) {
		if (!init(dst, (unsigned int) src.counts.size)) {
			fprintf(stderr, "array_multiset.copy ERROR: Unable to initialize destination multiset.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.counts.size; i++) {
			dst.counts.keys[i] = src.counts.keys[i];
			dst.counts.values[i] = src.counts.values[i];
		}
		dst.counts.size = src.counts.size;
		dst.sum = src.sum;
		return true;
	}

	template<typename Metric = default_metric>
	static inline long unsigned int size_of(const array_multiset<T>& s, const Metric& metric) {
		return core::size_of(s.counts, metric) + core::size_of(s.sum);
	}

	/**
	 * Frees the given array_multiset `s`. This function also frees the keys in
	 * the underlying array_map by calling core::free on each element.
	 */
	static inline void free(array_multiset<T>& s) {
		s.free();
		core::free(s.counts);
	}

private:
	inline void free() {
		for (auto entry : counts)
			core::free(entry.key);
	}
};

/**
 * Initializes an empty array_multiset `s` with the given `initial_capacity`
 * for the underlying array_map array_multiset::counts.
 */
template<typename T>
inline bool init(array_multiset<T>& s, unsigned int initial_capacity) {
	s.sum = 0;
	return array_map_init(s.counts, initial_capacity);
}

/**
 * Initializes the given array_multiset `s` by copying the contents from the
 * given array_multiset `src`.
 */
template<typename T>
inline bool init(array_multiset<T>& s, const array_multiset<T>& src) {
	if (!array_map_init(s.counts, (unsigned int) src.counts.size))
		return false;
	for (unsigned int i = 0; i < src.counts.size; i++) {
		s.counts.keys[i] = src.counts.keys[i];
		s.counts.values[i] = src.counts.values[i];
	}
	s.counts.size = src.counts.size;
	s.sum = src.sum;
	return true;
}

/**
 * Returns true if and only if the array_multiset `first` is equivalent to `second`.
 */
template<typename T>
inline bool operator == (const array_multiset<T>& first, const array_multiset<T>& second) {
	if (first.sum != second.sum)
		return false;
	unsigned int i = 0, j = 0;
	while (i < first.counts.size && j < second.counts.size) {
		if (first.counts.keys[i] == second.counts.keys[j]) {
			if (first.counts.values[i] != second.counts.values[j])
				return false;
			i++; j++;
		} else if (first.counts.keys[i] < second.counts.keys[j]) {
			if (first.counts.values[i] > 0)
				return false;
			i++;
		} else {
			if (second.counts.values[j] > 0)
				return false;
			j++;
		}
	}

	while (i < first.counts.size) {
		if (first.counts.values[i] > 0)
			return false;
		i++;
	} while (j < second.counts.size) {
		if (second.counts.values[j] > 0)
			return false;
		j++;
	}
	return true;
}

/**
 * Returns false if and only if the array_multiset `first` is equivalent to `second`.
 */
template<typename T>
inline bool operator != (const array_multiset<T>& first, const array_multiset<T>& second) {
	return !(first == second);
}

/**
 * Reads an array_multiset structure `s` from `in`.
 * \param reader a scribe that is passed to `read` for core::array_map.
 */
template<typename T, typename... Reader>
inline bool read(array_multiset<T>& s, FILE* in, Reader&&... reader) {
	if (!read(s.counts, in, std::forward<Reader>(reader)...))
		return false;
	s.sum = 0;
	for (unsigned int i = 0; i < s.counts.size; i++)
		s.sum += s.counts.values[i];
	if (s.counts.size > 1)
		sort(s.counts.keys, s.counts.values, (unsigned int) s.counts.size, default_sorter());
	return true;
}

/**
 * Writes the given array_multiset structure `s` to `out`.
 * \param writer a scribe that is passed to `write` for core::array_map.
 */
template<typename T, typename... Writer>
inline bool write(const array_multiset<T>& s, FILE* out, Writer&&... writer) {
	return write(s.counts, out, std::forward<Writer>(writer)...);
}

/**
 * Prints the given array_multiset structure `s` to `out`.
 * \param printer a scribe that is passed to `print` for core::array_map.
 */
template<typename T, typename... Printer>
inline void print(const array_multiset<T>& s, FILE* out, Printer&&... printer) {
	fputc('{', out);
	if (s.counts.size == 0) {
		fputc('}', out);
		return;
	}
	print(s.counts.keys[0], out, std::forward<Printer>(printer)...);
	fprintf(out, ":%u", s.counts.values[0]);
	for (unsigned int i = 1; i < s.counts.size; i++) {
		fprintf(out, ", ");
		print(s.counts.keys[i], out, std::forward<Printer>(printer)...);
		fprintf(out, ":%u", s.counts.values[i]);
	}
	fputc('}', out);
}


/**
 * A multiset structure that keeps track of the number of occurrences of
 * distinct elements in a set, implemented using a core::hash_map, where the
 * elements are the keys, and the values are their frequencies.
 *
 * array_multiset implements the same abstract data type using a
 * core::array_map and should be used if the number of distinct elements is
 * expected to be small.
 *
 * \tparam T the generic type of the elements. `T` must satisfy either:
 * 		1. [is_fundamental](http://en.cppreference.com/w/cpp/types/is_fundamental),
 * 		2. [is_enum](http://en.cppreference.com/w/cpp/types/is_enum),
 * 		3. [is_pointer](http://en.cppreference.com/w/cpp/types/is_pointer),
 * 		4. implements the public static method `unsigned int hash(const T&)`,
 * 			the public static method `void is_empty(const T&)`, implements the
 * 			operators `==`, satisfies [CopyAssignable](https://en.cppreference.com/w/cpp/named_req/CopyAssignable),
 * 			and core::is_moveable. **NOTE:** The first argument to the `==`
 * 			operator may be empty.
 */
template<typename T>
struct hash_multiset {
	/**
	 * The underlying hash_map.
	 */
	hash_map<T, unsigned int> counts;

	/**
	 * The sum of the values in hash_multiset::counts (i.e. the total number of
	 * occurrences of all elements).
	 */
	unsigned int sum;

	/**
	 * Constructs an empty hash_multiset with the given `initial_capacity` for
	 * the underlying hash_map hash_multiset::counts.
	 */
	explicit hash_multiset(unsigned int initial_capacity) : counts(initial_capacity), sum(0.0) { }
	~hash_multiset() { free(); }

	/**
	 * Returns the total number of occurrences of all elements (hash_multiset::sum).
	 */
	inline unsigned int total() const {
		return sum;
	}

	/**
	 * Adds the given item to the multiset.
	 */
	bool add(const T& item) {
		if (!counts.check_size()) {
			fprintf(stderr, "hash_multiset.add WARNING: Unable to expand hash_map.\n");
			return false;
		}

		bool contains;
		unsigned int index;
		unsigned int& count = counts.get(item, contains, index);
		if (!contains) {
			counts.table.keys[index] = item;
			counts.table.size++;
			count = 1;
		} else {
			count++;
		}
		sum++;
		return true;
	}

	/**
	 * Adds the given multiset of items to this multiset.
	 */
	bool add(const array_multiset<T>& items) {
		if (!counts.check_size(counts.table.size + items.counts.size)) {
			fprintf(stderr, "hash_multiset.add WARNING: Unable to expand hash_map.\n");
			return false;
		}

		bool contains;
		unsigned int index;
		for (unsigned int i = 0; i < items.counts.size; i++) {
			unsigned int& count = counts.get(items.counts.keys[i], contains, index);
			if (!contains) {
				counts.table.keys[index] = items.counts.keys[i];
				counts.table.size++;
				count = items.counts.values[i];
			} else {
				count += items.counts.values[i];
			}
		}
		sum += items.sum;
		return true;
	}

	/**
	 * Removes the given item from the multiset. This function assumes that
	 * item exists in the multiset with non-zero frequency.
	 */
	void subtract(const T& item)
	{
#if !defined(NDEBUG)
		bool contains;
		unsigned int& count = counts.get(item, contains);
		if (!contains) {
			fprintf(stderr, "hash_multiset.subtract WARNING: Attempted "
					"to remove more items from a bin than it contains.\n");
			return;
		}
#else
		unsigned int& count = counts.get(item);
#endif
		count--;
		sum--;
	}

	/**
	 * Removes the given multiset `items` from this multiset. This function
	 * assumes that the given multiset `items` is a subset of this multiset
	 * (i.e. this multiset contains all the keys in items with frequencies at
	 * least as large).
	 */
	void subtract(const array_multiset<T>& items)
	{
		for (unsigned int i = 0; i < items.counts.size; i++) {
			unsigned int& count = counts.get(items.counts.keys[i]);
#if !defined(NDEBUG)
			if (count < items.counts.values[i]) {
				fprintf(stderr, "hash_multiset.subtract WARNING: Attempted "
						"to remove more items from a bin than it contains.\n");
				count = 0;
			} else count -= items.counts.values[i];
#else
			count -= items.counts.values[i];
#endif
		}
		sum -= items.sum;
	}

	/**
	 * Moves the hash_multiset `src` into `dst`. Note that this function merely
	 * copies pointers, and not contents.
	 */
	static inline void move(const hash_multiset<T>& src, hash_multiset<T>& dst) {
		hash_map<T, unsigned int>::move(src.counts, dst.counts);
		dst.sum = src.sum;
	}

	/**
	 * Copies the contents of the array_multiset `src` into `dst`.
	 */
	static inline bool copy(const hash_multiset<T>& src, hash_multiset<T>& dst) {
		dst.sum = src.sum;
		return hash_map<T, unsigned int>::copy(src.counts, dst.counts);
	}

	template<typename Metric>
	static inline long unsigned int size_of(const hash_multiset<T>& s, const Metric& metric) {
		return size_of(s.counts, metric) + size_of(s.sum);
	}

	/**
	 * Frees the given hash_multiset `s`. This function also frees the keys in
	 * the underlying hash_map by calling core::free on each element.
	 */
	static inline void free(hash_multiset<T>& s) {
		free();
		core::free(s.counts);
	}

private:
	inline void free() {
		for (auto entry : counts)
			core::free(entry.key);
	}
};

/**
 * Initializes an empty hash_multiset `s` with the given `initial_capacity` for
 * the underlying hash_map hash_multiset::counts.
 */
template<typename T>
inline bool init(hash_multiset<T>& s, unsigned int initial_capacity) {
	s.sum = 0;
	return hash_map_init(s.counts, initial_capacity);
}

/**
 * Reads a hash_multiset `s` from `in`.
 * \param reader a scribe that is passed to `read` for core::hash_map.
 */
template<typename T, typename Reader>
inline bool read(hash_multiset<T>& s, FILE* in, Reader& reader) {
	if (!read(s.counts, in, reader))
		return false;
	s.sum = 0;
	for (unsigned int i = 0; i < s.counts.table.capacity; i++)
		if (!is_empty(s.counts.table.keys[i]))
			s.sum += s.counts.values[i];
	return true;
}

/**
 * Writes the given hash_multiset `s` to `out`.
 * \param writer a scribe that is passed to `write` for core::hash_map.
 */
template<typename T, typename Writer>
inline bool write(const hash_multiset<T>& s, FILE* out, Writer& writer) {
	return write(s.counts, out, writer);
}

/**
 * Prints the given hash_multiset `s` to `out`.
 * \param printer a scribe that is passed to `print` for core::hash_map.
 */
template<typename T, typename... Printer>
inline void print(const hash_multiset<T>& s, FILE* out, Printer&&... printer) {
	fputc('{', out);
	bool first = true;
	for (unsigned int i = 0; i < s.counts.table.capacity; i++) {
		if (is_empty(s.counts.table.keys[i]))
			continue;
		if (!first)
			fprintf(out, ", ");
		first = false;

		print(s.counts.table.keys[i], out, std::forward(printer)...);
		fprintf(out, ":%u", s.counts.values[i]);
	}
	fputc('}', out);
}

#endif /* HISTOGRAM_H_ */
