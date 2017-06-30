Documentation and code examples are available at <http://asaparov.org/docs/math>.

This repository contains a handful of data structures and algorithms for common mathematical operations in scientific computing.

### Dependencies

To use the code, simply download the files into a folder named "math". Add this folder to the include path for the compiler, for example by using the `-I` flag.

This library depends on [core](https://github.com/asaparov/core) and the data structures and procedures follow the same paradigm as core. Otherwise, there are no dependencies on external libraries. The code makes use of `C++11` and is regularly tested with `gcc 6` but I have previously compiled it with `gcc 4.8`, `clang 4.0`, and `Microsoft Visual C++ 14.0 (2015)`. The code is intended to be platform-independent, so please create an issue if there are any compilation bugs.

### Overview

This library contains the following files:
 - [distributions.h](http://asaparov.org/docs/math/distributions.h.html) contains structures that represent various probability distributions, such as the Dirichlet (both [symmetric](http://asaparov.org/docs/math/distributions.h.html#struct symmetric_dirichlet) and [general](http://asaparov.org/docs/math/distributions.h.html#struct dirichlet)), categorical (implemented using a [dense array](http://asaparov.org/docs/math/distributions.h.html#struct dense_categorical) or [sparsely using a hash_map](http://asaparov.org/docs/math/distributions.h.html#struct sparse_categorical)), [discrete uniform](http://asaparov.org/docs/math/distributions.h.html#struct uniform_distribution), and [sequence](http://asaparov.org/docs/math/distributions.h.html#struct uniform_distribution) distributions. This file contains functions for computing the probabilities and log probabilities of events according to these distributions, as well as generating samples from them.
 - [features.h](http://asaparov.org/docs/math/features.h.html) contains structures for representing feature vectors, which are commonly used in machine learning.
 - [histogram.h](http://asaparov.org/docs/math/histogram.h.html) contains histogram data structures, which are structures that count the number of occurrences of distinct elements in a set.
 - [log.h](http://asaparov.org/docs/math/log.h.html) contains functions to perform arithmetic in log space, while avoiding loss of floating-point precision.
