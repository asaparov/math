/**
 * histogram.h
 *
 *  Created on: Dec 14, 2015
 *      Author: asaparov
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "core/map.h"
#include <core/io.h>

using namespace core;

template<typename T>
struct array_histogram {
	array_map<T, unsigned int> counts;
	unsigned int sum;

	explicit array_histogram(unsigned int initial_capacity) : counts(initial_capacity), sum(0) { }
	~array_histogram() { free(); }

	inline unsigned int total() const {
		return sum;
	}

	inline bool add(const T& item, unsigned int count = 1, unsigned int start_index = 0) {
		if (!add_unsorted(item, count, start_index))
			return false;
		if (counts.size > 1)
			sort(counts.keys, counts.values, (unsigned int) counts.size, dummy_sorter()); /* TODO: test performance of different sorts */
		return true;
	}

	bool add_unsorted(const T& item, unsigned int count = 1, unsigned int start_index = 0) {
		unsigned int index = counts.index_of(item, start_index);
		if (index < counts.size) {
			counts.values[index] += count;
			sum += count;
			return true;
		}
		if (!counts.ensure_capacity((unsigned int) counts.size + 1)) {
			fprintf(stderr, "array_histogram.add WARNING: Unable to expand array_map.\n");
			return false;
		}
		counts.keys[counts.size] = item;
		counts.values[counts.size] = count;
		counts.size++;
		sum += count;
		return true;
	}

	bool add(const array_histogram<T>& items) {
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
					fprintf(stderr, "array_histogram.add ERROR: Unable to expand underlying map.\n");
					return false;
				}
				counts.keys[new_length] = items.counts.keys[j];
				counts.values[new_length] = items.counts.values[j];
				new_length++; j++;
			}
		}

		/* add the leftover elements from 'items' */
		if (!counts.ensure_capacity(new_length + (unsigned int) items.counts.size - j)) {
			fprintf(stderr, "array_histogram.add ERROR: Unable to expand array_map.\n");
			return false;
		}
		while (j < items.counts.size) {
			counts.keys[new_length] = items.counts.keys[j];
			counts.values[new_length] = items.counts.values[j];
			new_length++; j++;
		}

		counts.size = new_length;
		if (counts.size > 1)
			sort(counts.keys, counts.values, (unsigned int) counts.size, dummy_sorter());
		sum += items.sum;
		return true;
	}

	unsigned int subtract(const T& item, unsigned int start_index = 0)
	{
		unsigned int index = counts.index_of(item, start_index);
		if (index < counts.size) {
#if !defined(NDEBUG)
			if (counts.values[index] == 0) {
				fprintf(stderr, "array_histogram.subtract WARNING: Attempted "
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

		fprintf(stderr, "array_histogram.subtract WARNING: No such item.\n");
		return (unsigned int) counts.size;
	}

	/* this function assumes that 'histogram' contains all the keys in 'items' */
	void subtract(const array_histogram<T>& items)
	{
		unsigned int i = 0, j = 0;
		while (i < counts.size && j < items.counts.size) {
			if (counts.keys[i] == items.counts.keys[j]) {
#if !defined(NDEBUG)
				if (counts.values[i] < items.counts.values[j]) {
					fprintf(stderr, "array_histogram.subtract WARNING: Attempted "
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
			fprintf(stderr, "array_histogram.subtract WARNING: Missing bin in histogram.\n");
#endif
		sum -= items.sum;
	}

	inline void clear() {
		counts.clear();
		sum = 0;
	}

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

	inline bool is_sorted() const {
		for (unsigned int i = 1; i < counts.size; i++)
			if (counts.keys[i] < counts.keys[i - 1])
				return false;
		return true;
	}

	static inline unsigned int hash(const array_histogram<T>& key) {
		return default_hash(key.counts.keys, key.counts.size)
				^ default_hash(key.counts.values, key.counts.size);
	}

	static inline void move(const array_histogram<T>& src, array_histogram<T>& dst) {
		array_map<T, unsigned int>::move(src.counts, dst.counts);
		dst.sum = src.sum;
	}

	static inline bool copy(const array_histogram<T>& src, array_histogram<T>& dst) {
		if (!init(dst, (unsigned int) src.counts.size)) {
			fprintf(stderr, "array_histogram.copy ERROR: Unable to initialize destination histogram.\n");
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

	template<typename Metric = dummy_metric>
	static inline long unsigned int size_of(const array_histogram<T>& h, const Metric& metric) {
		return core::size_of(h.counts, metric) + core::size_of(h.sum);
	}

	static inline void free(array_histogram<T>& h) {
		h.free();
		core::free(h.counts);
	}

private:
	inline void free() {
		for (auto entry : counts)
			core::free(entry.key);
	}
};

template<typename T>
inline bool init(array_histogram<T>& h, unsigned int initial_capacity) {
	h.sum = 0;
	return array_map_init(h.counts, initial_capacity);
}

template<typename T>
inline bool init(array_histogram<T>& h, const array_histogram<T>& src) {
	if (!array_map_init(h.counts, (unsigned int) src.counts.size))
		return false;
	for (unsigned int i = 0; i < src.counts.size; i++) {
		h.counts.keys[i] = src.counts.keys[i];
		h.counts.values[i] = src.counts.values[i];
	}
	h.counts.size = src.counts.size;
	h.sum = src.sum;
	return true;
}

template<typename T>
inline bool operator == (const array_histogram<T>& first, const array_histogram<T>& second) {
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

template<typename T>
inline bool operator != (const array_histogram<T>& first, const array_histogram<T>& second) {
	return !(first == second);
}

template<typename T, typename... Reader>
inline bool read(array_histogram<T>& h, FILE* in, Reader&&... reader) {
	if (!read(h.counts, in, std::forward<Reader>(reader)...))
		return false;
	h.sum = 0;
	for (unsigned int i = 0; i < h.counts.size; i++)
		h.sum += h.counts.values[i];
	return true;
}

template<typename T, typename... Writer>
inline bool write(const array_histogram<T>& h, FILE* out, Writer&&... writer) {
	return write(h.counts, out, std::forward<Writer>(writer)...);
}

template<typename T, typename... Printer>
inline void print(const array_histogram<T>& h, FILE* out, Printer&&... printer) {
	fputc('{', out);
	if (h.counts.size == 0) {
		fputc('}', out);
		return;
	}
	print(h.counts.keys[0], out, std::forward<Printer>(printer)...);
	fprintf(out, " : %u", h.counts.values[0]);
	for (unsigned int i = 1; i < h.counts.size; i++) {
		fprintf(out, ", ");
		print(h.counts.keys[i], out, std::forward<Printer>(printer)...);
		fprintf(out, " : %u", h.counts.values[i]);
	}
	fputc('}', out);
}

template<typename T>
inline void print(const array_histogram<T>& h, FILE* out) {
	dummy_scribe scribe;
	print(h, out, scribe);
}

template<typename T>
struct hash_histogram {
	hash_map<T, unsigned int> counts;
	unsigned int sum;

	explicit hash_histogram(unsigned int initial_capacity) : counts(initial_capacity), sum(0.0) { }

	inline unsigned int total() const {
		return sum;
	}

	bool add(const T& item) {
		if (!counts.check_size()) {
			fprintf(stderr, "hash_histogram.add WARNING: Unable to expand hash_map.\n");
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

	bool add(const array_histogram<T>& items) {
		if (!counts.check_size(counts.table.size + items.counts.size)) {
			fprintf(stderr, "hash_histogram.add WARNING: Unable to expand hash_map.\n");
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

	void subtract(const T& item)
	{
#if !defined(NDEBUG)
		bool contains;
		unsigned int& count = counts.get(item, contains);
		if (!contains) {
			fprintf(stderr, "hash_histogram.subtract WARNING: Attempted "
					"to remove more items from a bin than it contains.\n");
			return;
		}
#else
		unsigned int& count = counts.get(item);
#endif
		count--;
		sum--;
	}

	/* this function assumes that 'histogram' contains all the keys in 'items' */
	void subtract(const array_histogram<T>& items)
	{
		for (unsigned int i = 0; i < items.counts.size; i++) {
			unsigned int& count = counts.get(items.counts.keys[i]);
#if !defined(NDEBUG)
			if (count < items.counts.values[i]) {
				fprintf(stderr, "hash_histogram.subtract WARNING: Attempted "
						"to remove more items from a bin than it contains.\n");
				count = 0;
			} else count -= items.counts.values[i];
#else
			count -= items.counts.values[i];
#endif
		}
		sum -= items.sum;
	}

	static inline void move(const hash_histogram<T>& src, hash_histogram<T>& dst) {
		hash_map<T, unsigned int>::move(src.counts, dst.counts);
		dst.sum = src.sum;
	}

	static inline bool copy(const hash_histogram<T>& src, hash_histogram<T>& dst) {
		dst.sum = src.sum;
		return hash_map<T, unsigned int>::copy(src.counts, dst.counts);
	}

	template<typename Metric>
	static inline long unsigned int size_of(const hash_histogram<T>& h, const Metric& metric) {
		return size_of(h.counts, metric) + size_of(h.sum);
	}

	static inline void free(hash_histogram<T>& h) {
		core::free(h.counts);
	}
};

template<typename T>
inline bool init(hash_histogram<T>& h, unsigned int initial_capacity) {
	h.sum = 0;
	return hash_map_init(h.counts, initial_capacity);
}

template<typename T, typename Reader>
inline bool read(hash_histogram<T>& h, FILE* in, Reader& reader) {
	if (!read(h.counts, in, reader))
		return false;
	h.sum = 0;
	for (unsigned int i = 0; i < h.counts.table.capacity; i++)
		if (!is_empty(h.counts.table.keys[i]))
			h.sum += h.counts.values[i];
	return true;
}

template<typename T, typename Writer>
inline bool write(const hash_histogram<T>& h, FILE* out, Writer& writer) {
	return write(h.counts, out, writer);
}

template<typename T, typename Printer>
inline void print(const hash_histogram<T>& h, FILE* out, Printer& printer) {
	fputc('{', out);
	bool first = true;
	for (unsigned int i = 0; i < h.counts.table.capacity; i++) {
		if (is_empty(h.counts.table.keys[i]))
			continue;
		if (!first)
			fprintf(out, ", ");
		first = false;

		print(h.counts.table.keys[i], out, printer);
		fprintf(out, " : %u", h.counts.values[i]);
	}
	fputc('}', out);
}

template<typename T>
inline void print(const hash_histogram<T>& h, FILE* out) {
	dummy_scribe scribe;
	print(out, h, scribe);
}

#endif /* HISTOGRAM_H_ */
