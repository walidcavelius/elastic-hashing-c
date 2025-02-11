#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t start;
  size_t size;
  size_t filled_slots;
  double epsilon;
} SubArray;

typedef struct {
  size_t start_index;
  size_t size;
  size_t current;
} Batch;

typedef struct {
  size_t *table;
  size_t n;
  double delta;
  size_t items;
  SubArray *arrays;
  size_t num_arrays;
  Batch current_batch;
  size_t beta;
} ElasticHash;

// implementation of φ(i,j) injection function from paper
static size_t phi(size_t i, size_t j) {
  size_t result = 1; // start with 1 as most significant bit

  // interleave j with 1s
  while (j > 0) {
    result = (result << 2) | ((j & 1) << 1) | 1;
    j >>= 1;
  }

  result = (result << 1);

  // i bits
  result = (result << (size_t)log2(i + 1)) | i;

  return result;
}

ElasticHash *elastic_hash_init(size_t n, double delta) {
  ElasticHash *eh = malloc(sizeof(ElasticHash));
  if (!eh)
    return NULL;

  eh->n = n;
  eh->delta = delta;
  eh->items = 0;
  eh->beta = (size_t)(2 * log2(1 / delta)); // as in the paper

  eh->table = calloc(n, sizeof(size_t));
  if (!eh->table) {
    free(eh);
    return NULL;
  }

  // α = ⌈4 log δ⁻¹ + 10⌉
  eh->num_arrays = (size_t)ceil(4 * log2(1 / delta) + 10);
  eh->arrays = malloc(eh->num_arrays * sizeof(SubArray));
  if (!eh->arrays) {
    free(eh->table);
    free(eh);
    return NULL;
  }

  // geometrically decreasing sizes
  size_t remaining = n;
  size_t start = 0;
  for (size_t i = 0; i < eh->num_arrays; i++) {
    size_t ai = (remaining * 3) / 4; // ai+1 = 3ai/4
    if (ai < eh->beta)
      ai = eh->beta; // at least β slots

    eh->arrays[i].start = start;
    eh->arrays[i].size = ai - (ai % eh->beta);
    eh->arrays[i].filled_slots = 0;
    eh->arrays[i].epsilon = 1.0;

    start += eh->arrays[i].size;
    remaining -= eh->arrays[i].size;

    if (remaining < eh->beta)
      break;
  }

  eh->current_batch.start_index = 0;
  eh->current_batch.size = (size_t)(0.75 * eh->arrays[0].size);
  eh->current_batch.current = 0;

  return eh;
}

// j th probe for key in array i
static size_t probe_sequence(const ElasticHash *eh, size_t key, size_t i,
                             size_t j) {
  size_t hash = phi(i, j);
  hash ^= key;
  return eh->arrays[i].start + (hash % eh->arrays[i].size);
}

// f(ε)
static size_t compute_f(double epsilon, double delta) {
  double log_eps = log2(1.0 / epsilon);
  double log_delta = log2(1.0 / delta);
  return (size_t)(4 * fmin(log_eps * log_eps, log_delta));
}

bool elastic_hash_insert(ElasticHash *eh, size_t key) {
  if (eh->items >= eh->n * (1 - eh->delta))
    return false;

  size_t batch_i = (size_t)log2(eh->current_batch.start_index + 1);

  // three cases from the paper
  for (size_t i = batch_i; i < eh->num_arrays; i++) {
    SubArray *current = &eh->arrays[i];
    current->epsilon = 1.0 - ((double)current->filled_slots / current->size);

    // case 1: if εᵢ > δ/2 and εᵢ₊₁ > 0.25
    if (i < eh->num_arrays - 1 && current->epsilon > eh->delta / 2 &&
        eh->arrays[i + 1].epsilon > 0.25) {

      size_t f_eps = compute_f(current->epsilon, eh->delta);
      for (size_t j = 0; j < f_eps; j++) {
        size_t pos = probe_sequence(eh, key, i, j);
        if (eh->table[pos] == 0) {
          eh->table[pos] = key;
          current->filled_slots++;
          eh->items++;
          return true;
        }
      }

      continue;
    }

    // case 2: if εᵢ ≤ δ/2
    if (current->epsilon <= eh->delta / 2) {
      if (i + 1 < eh->num_arrays) {
        continue;
      }
    }

    // case 3: if εᵢ₊₁ ≤ 0.25
    if (i + 1 < eh->num_arrays && eh->arrays[i + 1].epsilon <= 0.25) {
      for (size_t j = 0;; j++) {
        size_t pos = probe_sequence(eh, key, i, j);
        if (eh->table[pos] == 0) {
          eh->table[pos] = key;
          current->filled_slots++;
          eh->items++;
          return true;
        }
      }
    }
  }

  eh->current_batch.current++;
  if (eh->current_batch.current >= eh->current_batch.size) {
    size_t i = (size_t)log2(eh->current_batch.start_index + 1);
    if (i + 1 < eh->num_arrays) {
      eh->current_batch.start_index += eh->current_batch.size;
      eh->current_batch.size =
          eh->arrays[i + 1].size - (size_t)(0.75 * eh->arrays[i + 1].size);
      eh->current_batch.current = 0;
    }
  }

  return false;
}

ssize_t elastic_hash_search(const ElasticHash *eh, size_t key) {
  for (size_t i = 0; i < eh->num_arrays; i++) {
    const SubArray *current = &eh->arrays[i];
    size_t f_eps = compute_f(current->epsilon, eh->delta);

    for (size_t j = 0; j < f_eps; j++) {
      size_t pos = probe_sequence(eh, key, i, j);
      if (eh->table[pos] == key) {
        return (ssize_t)pos;
      }
      if (eh->table[pos] == 0)
        break;
    }
  }
  return -1;
}

void elastic_hash_destroy(ElasticHash *eh) {
  if (eh) {
    free(eh->table);
    free(eh->arrays);
    free(eh);
  }
}

// exemple use
int main() {
  // initialize with n=1024 and δ=0.1
  ElasticHash *eh = elastic_hash_init(1024, 0.1);

  // insert some keys
  for (size_t i = 1; i <= 900; i++) {
    if (!elastic_hash_insert(eh, i)) {
      printf("Insertion failed at %zu\n", i);
      break;
    }
  }

  // search for a key
  ssize_t pos = elastic_hash_search(eh, 42);
  if (pos >= 0) {
    printf("Found key 42 at position %zd\n", pos);
  }

  elastic_hash_destroy(eh);
  return 0;
}
