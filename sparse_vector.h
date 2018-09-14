#ifndef SPARSE_VECTOR_H_
#define SPARSE_VECTOR_H_

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include <core/array.h>
#include <core/map.h>
#include <core/utility.h>

using namespace core;

template<typename K, typename V = double>
struct list
{
	std::vector<std::pair<K, V>> elements;

	list(unsigned int length) {
		elements.reserve(length);
	}

	list(const hash_map<K, V>& src)
	{
		elements.reserve(src.table.size);
		add(src);
	}

	void add(const hash_map<K, V>& src)
	{
		for (unsigned int i = 0; i < src.table.capacity; i++)
			if (!is_empty(src.table.keys[i]))
				elements.push_back(std::make_pair(src.table.keys[i], src.values[i]));
	}

	inline K& get_key(unsigned int index) {
		return elements.data()[index].first;
	}

	inline V& get_value(unsigned int index) {
		return elements.data()[index].second;
	}
};

template<typename K, typename V = double>
struct sparse_vector
{
	unsigned int length;
	hash_map<K, V> elements;
	V unspecified;

	sparse_vector(unsigned int length, const V& unspecified) :
			length(length), elements(128), /* TODO: maybe this shouldn't be fixed? */
			unspecified(unspecified) { }

	inline void set(K index, V value) {
		elements.put(index, value);
	}

	inline void add(K index, V value)
	{
		bool contains;
		V& prev = elements.get(index, contains);
		if (contains)
			prev += value;
		else elements.put(index, unspecified + value);
	}

	inline void add(const sparse_vector<K, V>& vector)
	{
#if !defined(NDEBUG)
		if (vector.length != length)
			fprintf(stderr, "sparse_vector.add WARNING: Vector lengths do not match.\n");
#endif
		unspecified += vector.unspecified;
		for (unsigned int i = 0; i < vector.elements.table.capacity; i++)
			if (!is_empty(vector.elements.table.keys[i]))
				add(vector.elements.table.keys[i], vector.elements.values[i]);
	}

	inline const V get(K index) const
	{
		bool contains;
		V value = elements.get(index, contains);
		if (contains)
			return value;
		else return unspecified;
	}

	inline void clear() {
		elements.clear();
	}

	V sum() const
	{
		V total = 0;
		for (unsigned int i = 0; i < elements.table.capacity; i++)
			if (!is_empty(elements.table.keys[i]))
				total += elements.values[i];

		return total + unspecified * (length - elements.table.size);
	}

	void get_values(K* keys, V* values, unsigned int length) const {
		for (unsigned int i = 0; i < length; i++)
			values[i] = get(keys[i]);
	}

	void set_length(unsigned int new_length) {
		length = new_length;
	}

	void renormalize()
	{
		if (elements.table.size == 0) {
			unspecified = 1.0 / length;
			return;
		}

		V factor = sum();
		unspecified /= factor;

		for (unsigned int i = 0; i < elements.table.capacity; i++)
			if (!is_empty(elements.table.keys[i]))
				elements.values[i] /= factor;
	}

	static void swap(sparse_vector<K, V>& first, sparse_vector<K, V>& second) {
		hash_map<K, V>::swap(first.elements, second.elements);
		core::swap(first.length, second.length);
		core::swap(first.unspecified, second.unspecified);
	}

	static void move(const sparse_vector<K, V>& src, sparse_vector<K, V>& dst) {
		dst.length = src.length;
		dst.unspecified = src.unspecified;
		hash_map<K, V>::move(src.elements, dst.elements);
	}

	static void sparse_vector_free(sparse_vector<K, V>& vector) {
		core::free(vector.elements);
	}
};

template<typename K, typename V>
bool sparse_vector_init(sparse_vector<K, V>& vector,
		unsigned int length, const V& unspecified)
{
	vector.length = length;
	vector.unspecified = unspecified;
	/* TODO: maybe this shouldn't be fixed? */
	return hash_map_init(vector.elements, 32);
}

template<typename K, typename V, typename KeyReader>
bool read(sparse_vector<K, V>& vector, FILE* in, KeyReader& key_reader) {
	bool success = true;
	success &= read(vector.length, in);
	success &= read(vector.unspecified, in);
	success &= read(vector.elements, in, key_reader);
	return success;
}

template<typename K, typename V, typename KeyWriter>
bool write(const sparse_vector<K, V>& vector, FILE* out, KeyWriter& key_writer) {
	bool success = true;
	success &= write(vector.length, out);
	success &= write(vector.unspecified, out);
	success &= write(vector.elements, out, key_writer);
	return success;
}

template<typename K, typename V>
void print(const sparse_vector<K, V>& vector, FILE* out, bool new_line)
{
	if (vector.elements.table.size == 0)
		fprintf(out, "{}");
	else {
		unsigned int i = 0;
		while (is_empty(vector.elements.table.keys[i])) i++;
		fprintf(out, "{");
		print(vector.elements.table.keys[i], out, false);
		fprintf(out, " : ");
		print(vector.elements.values[i], out, false);

		for (i++; i < vector.elements.table.capacity; i++) {
			if (is_empty(vector.elements.table.keys[i]))
				continue;
			fprintf(out, ", ");
			print(vector.elements.table.keys[i], out, false);
			fprintf(out, " : ");
			print(vector.elements.values[i], out, false);
		}
		fprintf(out, "}");
	}
	fprintf(out, " contains %u items.", vector.length);
	if (new_line) fprintf(out, "\n");
}

#endif /* SPARSE_VECTOR_H_ */
