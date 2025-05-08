# Ribos Fragments Simulation

**Date**: November 2024

## Inputs and Outputs

### Inputs

#### Raw Data

For each cell line, we are given two tables: `Treated` and `NonTreated`, with the following structure:

| **Length** ($\ell_i$) | **Count** ($c_i$) |
|-----------------------|-------------------|
| $\ell_1$              | $c_1$             |
| $\ell_2$              | $c_2$             |
| $\vdots$              | $\vdots$         |
| $\ell_n$              | $c_n$             |

Where the values $\ell_i$ are distinct integer lengths, and $c_i$ is the corresponding count of fragments of that length.

#### Parameters

- `Total`: Total number of fragments to consider during the simulation.
- `n_reps`: Number of repetitions of the simulation to run while doing the search.
- `threshold_N`: Percentual threshold to finish the search.
- `min_size`: The minimum length in the tables `Treatment` and `NonTreatment`:

  $$
  \texttt{min\_size} = \min(\ell_i)
  $$

### Outputs

- `N`: Number of cuts at which the binary search stopped.
- `Fragments_try`: List of `n_reps` dataframes (tables) from the last round of simulations with structure as above, satisfying $\sum c_i \approx \texttt{Total}$.
- `total_cuts_try`: List of `n_reps` integers, counting how many cuts were introduced in each simulation. Should equal `N` in the updated model.
- `new_lengths_try`: List of `n_reps` lists of fragment lengths after the last simulation round. These lengths are not unique and can repeat.
- `medians_try`: Median of each of the tables in `Fragments_try`. Should be close to `Goal_Median` if the simulation is successful.

## The Algorithms

### Data Preparation

From the raw data, we compute:

1. `Initial_Median`: Median of the `NonTreated` table:

   $$
   \texttt{Initial\_Median} = \ell_k \quad \text{where } k \text{ is the maximal index such that } \frac{\sum_{i=1}^k c_i}{\sum_{i=1}^n c_i} \leq 0.5
   $$

2. `Goal_Median`: Median of the `Treated` table.

3. `Fragments`: A dataframe with the data from `NonTreated` normalized to contain approximately `Total` elements:

   $$
   F(\ell_i) = \left\lfloor \frac{\texttt{NonTreated}(\ell_i)}{\sum_j \texttt{NonTreated}(\ell_j)} \times \texttt{Total} \right\rfloor
   $$

Where $F(\ell_i)$ is the count corresponding to length $\ell_i$ in the normalized table.

### Main Simulation

**Inputs**: $F = \texttt{Fragments}$, $N$ = number of cuts, $\mu = \texttt{min\_size}$.

Steps:

1. Compute weights: $w_i = \ell_i \times F(\ell_i)$.
2. Sample $N$ elements from $\{\ell_1, \ldots, \ell_n\}$ with probabilities:

   $$
   \mathbb{P}(\ell_i) = \frac{w_i}{\sum_{j=1}^n w_j}
   $$

3. Initialize:
   - `new_lengths` = []
   - `total_cuts` = 0

4. For each unique length $\ell$ in the sample:
   - Let $c = E(\ell)$ = number of times $\ell$ appears.
   - Max cuts per fragment: $m = \left\lfloor \frac{\ell}{\mu} \right\rfloor - 1$
   - Number of such fragments: $f = F(\ell)$
   - Max possible cuts: $M = m \times f$
   - Sample $\min(M, c)$ indices from the multiset $\{1, 1, \ldots, f\}$ (each repeated $m$ times), **without** replacement.
   - For each unique index $i$:
     - Cut using `cut_given_length` and update:
       - `new_lengths`
       - `total_cuts`
       - Table $F$ with:

         $$
         F(\ell) \leftarrow F(\ell) - 1; \quad F(l_j) \leftarrow F(l_j) + 1
         $$

**Returns**: $F$, `total_cuts`, `new_lengths`

### `cut_given_length`

**Input**: $\ell$ = fragment length, $\mu$ = min size, $m$ = number of cuts

Steps:

1. Compute:

   $$
   \texttt{extra\_space} = \ell - \mu(m+1)
   $$

2. Sample $m$ integers $x_1, \dots, x_m$ (with replacement) from $\{0, \dots, \texttt{extra\_space}\}$.
3. Order $x_i$ so that $x_1 \leq \cdots \leq x_m$, and define $x_0 = 0$, $x_{m+1} = \texttt{extra\_space}$.
4. Compute lengths:

   $$
   l_i = \mu + (x_i - x_{i-1})
   $$

**Returns**: $\{l_1, \dots, l_{m+1}\}$

### Binary Search for $N$

**Input**: `Goal_Median`, $F$ = `Fragments`, `Total`, `n_reps`, `threshold_N`

Steps:

1. Initialize:
   - $N_{\min} = 0$, $N_{\max} = \texttt{Total}$, $N = \texttt{Total} / 5$

2. Loop until convergence:

   - Run the simulation `n_reps` times
   - Collect:
     - `medians_try`, `Fragments_try`, `total_cuts_try`, `new_lengths_try`
   - Compute:

     $$
     \texttt{Average\_Median} = \text{mean}(\texttt{medians\_try})
     $$

   - Decision:
     - If $|\texttt{Average\_Median} - \texttt{Goal\_Median}| \leq 1$: break
     - If $\texttt{Average\_Median} < \texttt{Goal\_Median}$:
       - $N_{\max} = N$
       - $N_{\text{new}} = \frac{N_{\min} + N_{\max}}{2}$
     - If $\texttt{Average\_Median} > \texttt{Goal\_Median}$:
       - $N_{\min} = N$
       - $N_{\text{new}} = \frac{N_{\min} + N_{\max}}{2}$
     - If $\frac{|N_{\text{new}} - N|}{\texttt{Total}} \leq \texttt{threshold\_N}$: break

**Return**: `N`, `medians_try`, `Fragments_try`, `total_cuts_try`, `new_lengths_try`
