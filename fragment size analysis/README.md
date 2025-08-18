# Ribos Fragments Simulation

# Using the code

1. The Jupyter Notebook `Fragment_Simulation_All.ipynb` runs the simulation for all the gels and celltypes. It takes as arguments:


  -`All_main_files`: A list with the locations of the files. In this case all files are in the folder `fragments_data`. The data should be in .cvs format and follow the naming convention: `gel_name`+`case`.csv, where `case` is either `treated` or `nontreated`. 
  
  -`All_N_initials`: A list of initial values for `N` (number of breaks to simulate) to accelerate the binary search. The first time we ran this algorithm we chose $N=2\times 10^6$ for all cases, but the uploaded file starts the search at the values of $N$ obtained during our simulations. 
  
  -`data_labels`: A list of the column names in the .csv files. It has the name of the fragment size column (in this case `0`) and the names of the celltypes. Similarly `labels` includes just the names of the celltypes

2. Running this notebook produces the files stored in the folder `output`, for each celltype and gel, it stores the average simulated distribution (named as `Mean_Distribution_`+(celltype)+`_Total_`+(size of simulation)+`_reps_`+(number of repetitions)+`.pkl`. Similarlty it produces a file with statistics, named similarly, but with `Statistics` instead of `Mean_Distribution`. This includes the average $N$ (number of introduced breaks), as the average number of DNA bases in the simulated fragments. 

3. The Jupyter Notebook `Results.ipynb` combines the output files for each gel, to produce .csv files with the results.

4. We have formatted and combined the .csv files with Statistics into a single Excel file `DNAs_per_Ribo.xlsx'


# Explaining the algorithm
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

  $$\texttt{minsize} = \min(\ell_i)$$

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

   $$\texttt{InitialMedian} = \ell_k \quad \text{where } k \text{ is the maximal index such that } \frac{\sum_{i=1}^k c_i}{\sum_{i=1}^n c_i} \leq 0.5$$

2. `Goal_Median`: Median of the `Treated` table.

3. `Fragments`: A dataframe with the data from `NonTreated` normalized to contain approximately `Total` elements:

   $$F(\ell_i) = \left\lfloor \frac{\texttt{NonTreated}(\ell_i)}{\sum_j \texttt{NonTreated}(\ell_j)} \times \texttt{Total} \right\rfloor$$

Where $F(\ell_i)$ is the count corresponding to length $\ell_i$ in the normalized table.

### Main Simulation

**Inputs**: $F = \texttt{Fragments}$, $N$ = number of cuts, $\mu = \texttt{minsize}$.

Steps:

1. Compute weights: $w_i = \ell_i \times F(\ell_i)$.
2. Sample $N$ elements from $\{\ell_1, \ldots, \ell_n\}$ with probabilities:

   $$\mathbb{P}(\ell_i) = \frac{w_i}{\sum_{j=1}^n w_j}$$

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

         $$F(\ell) \leftarrow F(\ell) - 1; \quad F(l_j) \leftarrow F(l_j) + 1$$

**Returns**: $F$, `total_cuts`, `new_lengths`

### `cut_given_length`

**Input**: $\ell$ = fragment length, $\mu$ = min size, $m$ = number of cuts

Steps:

1. Compute:

   $$\texttt{extraspace} = \ell - \mu(m+1)$$

2. Sample $m$ integers $x_1, \dots, x_m$ (with replacement) from $\{0, \dots, \texttt{extraspace}\}$.
3. Order $x_i$ so that $x_1 \leq \cdots \leq x_m$, and define $x_0 = 0$, $x_{m+1} = \texttt{extraspace}$.
4. Compute lengths:

   $$l_i = \mu + (x_i - x_{i-1})$$

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

     $$\texttt{AverageMedian} = \text{mean}(\texttt{medianstry})$$

   - Decision:
     - If $|\texttt{AverageMedian} - \texttt{GoalMedian}| \leq 1$: break
     - If $\texttt{AverageMedian} < \texttt{GoalMedian}$:
       - $N_{\max} = N$
       - $N_{\text{new}} = \frac{N_{\min} + N_{\max}}{2}$
     - If $\texttt{AverageMedian} > \texttt{GoalMedian}$:
       - $N_{\min} = N$
       - $N_{\text{new}} = \frac{N_{\min} + N_{\max}}{2}$
     - If $\frac{|N_{\text{new}} - N|}{\texttt{Total}} \leq \texttt{thresholdN}$: break

**Return**: `N`, `medians_try`, `Fragments_try`, `total_cuts_try`, `new_lengths_try`
