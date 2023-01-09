# ALBEDO: Tacit Programming Guide & High-Performance Array Architecture

J is a dialect of APL designed by Kenneth E. Iverson that achieves high throughput by expressing computational graphs as point-free (tacit) compositions of verbs, adverbs, and conjunctions.

This guide explains the tacit programming idioms, rank conjunction semantics, and memory models used throughout ALBEDO.

---

## 1. The Anatomy of Verb Trains: Forks and Hooks

Tacit programs eliminate explicit named arguments, compiling directly into internal bytecode trees that vectorize across SIMD instructions and cache lines.

### The 3-Verb Fork: `(f g h)`
A 3-verb sequence `(f g h)` forms a fork. The outer verbs `f` and `h` process the arguments, and the central verb `g` combines their results.

- **Monadic Fork**:
  $$(f \; g \; h) \; y = (f \; y) \; g \; (h \; y)$$
  
  *Example in ALBEDO (`core/stats.ijs`)*:
  Computing the arithmetic mean vector:
  ```j
  mean_vector =: +/ % #
  ```
  Here:
  - `f` is `+/` (sum items)
  - `g` is `%` (divide)
  - `h` is `#` (tally / length)
  Evaluation: `(+/ y) % (# y)`.

- **Dyadic Fork**:
  $$x \; (f \; g \; h) \; y = (x \; f \; y) \; g \; (x \; h \; y)$$

  *Example in ALBEDO (`core/linalg.ijs`)*:
  Matrix multiplication and projection:
  ```j
  proj =: [ +/ .* %.@]
  ```

### The 2-Verb Hook: `(f g)`
A 2-verb sequence `(f g)` forms a hook:

- **Monadic Hook**:
  $$(f \; g) \; y = y \; f \; (g \; y)$$

  *Example in ALBEDO (`core/stats.ijs`)*:
  Centering data by subtracting its mean:
  ```j
  center =: - mean_vector
  ```
  Evaluation: `y - (mean_vector y)`.

- **Dyadic Hook**:
  $$x \; (f \; g) \; y = x \; f \; (g \; y)$$

---

## 2. The Rank Conjunction `"`: Dimensional Control

In J, every noun has a rank (0 for scalars, 1 for vectors, 2 for matrices, 3+ for tensors). Every verb has an intrinsic rank specifying the sub-arrays it naturally operates on:
- Monadic addition `+` has rank `0 0 0` (scalar operation).
- Matrix inverse `%.` has rank `2 2 2` (matrix operation).

The rank conjunction `v " r` forces verb `v` to operate on sub-arrays of rank `r`:

### Vector Normalization across Matrix Rows
Given matrix `M` of shape `(N, D)`, to compute Euclidean norm for each row:
```j
row_norms =: albedo_linalg_norm2 "1 M
```
Here `"1` instructs `albedo_linalg_norm2` to treat each row (rank 1 sub-array) as an independent argument.

### Column Normalization Agreement
When dividing an $(N \times D)$ matrix by a length-$D$ scale vector:
```j
M_scaled =: M %"1 scale_vec
```
Without `"1`, dyadic `%` would attempt to match the $N$ items of $M$ with the $D$ items of `scale_vec`, throwing a `|length error`. With `%"1`, each length-$D$ row of $M$ is divided element-wise by `scale_vec`.

---

## 3. The Under Conjunction `&.`: Structural Inverses

The Under conjunction `f &. g` applies transformation `g`, executes verb `f` in the transformed domain, and automatically applies the inverse `g^:_1`:

$$(f \; \&. \; g) \; y = g^{-1}(f(g(y)))$$

### Example: Log-Domain Geometric Mean
```j
geom_mean =: (+/ % #) &. ^.
```
Execution trace:
1. Apply `^.` (natural logarithm) to each element.
2. Compute arithmetic mean `+/ % #` in log-space.
3. Automatically apply `^` (exponential, the inverse of `^.`).

---

## 4. In-Place Mutation and Sparse Updating

In J, the amend adverb `}` replaces elements in-place without copying unchanged memory regions:

### Syntax
- `new_val (<indices)} array`

### Examples in ALBEDO
1. **Setting Matrix Diagonal**:
   ```j
   set_diag =: 4 : 0
     diag_vals =. x
     A =. y
     N =. {. $ A
     coords =. (i. N) ,. (i. N)
     diag_vals (<"1 coords)} A
   )
   ```

2. **Zeroing Noise Self-Affinity Diagonal**:
   ```j
   W_zero_diag =: (1.0 - albedo_linalg_eye N) * W
   ```

---

## 5. Tensor Permutations via Dyadic Transpose `|:`

Standard multidimensional programming requires nested for-loops to reorder indices. J executes arbitrary axis permutations in single SIMD strides via `perm |: tensor`.

### Example: Permuting a 3D Hyperspectral Cube
Given cube $\mathcal{X}$ of shape `(B, H, W)`:
```j
BSQ_to_BIL =: 1 0 2 |: Cube   NB. Result: (H, B, W)
BSQ_to_BIP =: 1 2 0 |: Cube   NB. Result: (H, W, B)
BIP_to_BSQ =: 2 0 1 |: Cube   NB. Result: (B, H, W)
```
These primitives enable $O(1)$ stride adjustments that eliminate memory buffer allocation during spectral data ingest and pipeline handoffs.
