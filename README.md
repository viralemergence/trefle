# A data-inflated host-virus association database

C'est quoi, `trefle`? It is a data product derived from the [`clover`][clover]
database of mammals-virus association. Specifically, `trefle` was produced using
LF-SVD imputation, a two-step algorithm where novel host-virus associations are
recommended based on truncated singular value decomposition applied to initial
values based on a linear filter.

[clover]: https://github.com/viralemergence/clover

## LF-SVD

Associations in `trefle` are recommended based on the output of a two-step
process. First [linear filtering][LF] is used to generate an initial value based
on network properties. The linear filter has four hyper-parameters (the four
weights assigned to the initial association, the connectance, and the in and out
degree of the nodes), constrained as their values must sum to one.

[LF]: https://www.nature.com/articles/srep45908

Second, we apply truncated SVD to the modified `clover` wherein the missing
association we impute get its initial value from to the linear filter. The rank
of truncation for the low-rank approximation is a fifth hyper-parameter in this
model.

In short, `trefle` is a giant LOOCV dataset. This has consequences for how much
computational resources are required to *produce* it, which we will approximate
as: hella. We will discuss the computational requirements more below.

## Hyper-parameters tuning

In practice, we can get away with removing the first hyper-parameter of the
linear filter, as we have reasons to suspect that negative associations can
often be false negatives. This leaves us with four hyper-parameters to tune.

Because exploring the grid of linear filter parameters would be prohibitive in
terms of computing time (but also would lead to less interpretable model
inputs), we picked three initial models: the initial value is the same for all
associations and determined by the connectance of `clover` (`connectance`); the
initial value is given by the averaged relative degree of the host and the virus
(`degree`); the initial value is given by the average of the previous two models
(`hybrid`).

We applied each model at various depth of low-rank approximation, *i.e.* by
truncating the SVD to its 1st to 20th singular value. Within each model-rank
combination, we imputed the value of 780 positive interactions (which we should
assume are true positive given the nature of the `clover` data), and of 780
negative interactions (about which we will refrain from making assumptions),
using LOOCV.

The performance of each model-rank combination was measured using ROC-AUC,
assuming that negative interactions are true negatives. Note that owing to the
dimensions of `clover`, the training sample represents less than 1/1000 of the
entire dataset. Further, for each model we decided on a threshold of evidence
above which the pseudo-probability should be indicative of an actual association
by picking the value of evidence which maximizes Youden's J statistic. In the
overwheling majority of cases, this value of evidence *also* maximized the
accuracy of the model.

## Output values

The output value in `trefle` is akin to an association probability (but it is
not a probability of association in the sense of [probabilistic ecological
networks][pen]). The final value after imputation is divided by the initial
value before imputation. If the association "score" does not change, this gives
a value of 1. We transform this by substracting one from the result, yielding an
*evidence* value for the association: positive evidence makes the association
more likely. To convert the evidence into a pseudo-probability, we put it
through the logistic function. This returns values in [0;1]. In practice, owing
to the numerical imprecisions involved in measuring the logistic on even
moderately large floating-point numbers on 64 bits, it is common to have final
pseudo-probability values of 1, and we rely on the *evidence* for ranking.

[pen]: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12468

The following figure is an illustration of the resulting probabilities in an ensemble model of all of the model candidates used during tuning - the little bump in values that are `false` around 1 are candidate false negatives:

![proba ensemble](model_performance/probabilities.png)

## Model performance

### Top 10 models

The following table has the 10 best models ranked from first to last, as well as
the usual measures of model performance derived from the confusion table. In
addition to the AUC and cutoff (expressed as a *pseudo-probability*), we report
the true positive and true negative rates (TPR, TNR), the positive and negative
predictive values (PPV, NPV), the false negative and positive rates (FNR, FPR),
the false discovery and false omission rates (FDR, FOR), the critical success
index (CSI), accuracy (ACC), and Youden's J.

| model         | rank | AUC   | cutoff | TPR   | TNR   | PPV   | NPV   | FNR   | FPR   | FDR   | FOR   | CSI   | ACC   | J     |
|---------------|------|-------|--------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| `connectance` | 12   | 0.849 | 0.846  | 0.720 | 0.925 | 0.906 | 0.769 | 0.28  | 0.074 | 0.093 | 0.230 | 0.669 | 0.823 | 0.645 |
| `connectance` | 11   | 0.846 | 0.908  | 0.684 | 0.936 | 0.914 | 0.75  | 0.315 | 0.063 | 0.085 | 0.25  | 0.643 | 0.811 | 0.621 |
| `connectance` | 17   | 0.844 | 0.929  | 0.692 | 0.935 | 0.913 | 0.754 | 0.307 | 0.064 | 0.086 | 0.245 | 0.649 | 0.814 | 0.627 |
| `connectance` | 8    | 0.842 | 0.705  | 0.701 | 0.895 | 0.868 | 0.751 | 0.298 | 0.104 | 0.131 | 0.248 | 0.634 | 0.798 | 0.596 |
| `hybrid`      | 12   | 0.841 | 0.707  | 0.703 | 0.877 | 0.851 | 0.748 | 0.296 | 0.122 | 0.148 | 0.251 | 0.626 | 0.790 | 0.581 |
| `connectance` | 14   | 0.839 | 0.902  | 0.700 | 0.929 | 0.907 | 0.758 | 0.299 | 0.070 | 0.092 | 0.241 | 0.653 | 0.815 | 0.629 |
| `hybrid`      | 11   | 0.837 | 0.820  | 0.647 | 0.918 | 0.888 | 0.723 | 0.352 | 0.081 | 0.111 | 0.276 | 0.598 | 0.783 | 0.566 |
| `connectance` | 5    | 0.836 | 0.931  | 0.660 | 0.940 | 0.916 | 0.735 | 0.339 | 0.059 | 0.083 | 0.264 | 0.623 | 0.800 | 0.600 |
| `connectance` | 7    | 0.836 | 0.948  | 0.655 | 0.957 | 0.939 | 0.735 | 0.344 | 0.042 | 0.060 | 0.264 | 0.628 | 0.806 | 0.613 |
| `connectance` | 16   | 0.835 | 0.961  | 0.667 | 0.945 | 0.923 | 0.741 | 0.332 | 0.054 | 0.076 | 0.258 | 0.632 | 0.807 | 0.613 |

Following these results, we have conducted the imputation with on the model
based on connectance and a rank 12 approximation. Visualisations of all these
metrics are provided in `model_performance/metrics`.

### Overview of the best model

The following figure is the ROC AUC, with a depiction of the point maximizing
Youden's J and the probability cutoff associated:

![ROC-AUC](model_performance/roc/rank-12-model-connectance.png)

Visualisations of the same curve for all model-rank combinations are in
`model_performance/roc`.
## Computational resources

We assembled `trefle` on the [beluga][beluga] supercomputer, operated by *Calcul
Québec*, using a pipeline built entirely in [Julia][jl] (1.5.2).

[beluga]: https://www.computecanada.ca/featured/beluga-the-latest-supercomputer-for-canadian-researchers/
[jl]: https://julialang.org/

Tuning the hyper-parameters required about 2400 core hours, and imputation took
approximately 59500 core hours. Rounding up, using recent ARC hardware, the
assembly of `trefle` takes 62000 core hours, or just above 7 core years.
Assuming a cost of $0.051 per hour (equivalent to what a commercial cloud
computing provider would charge), the entire `trefle` production process costs
about $3200.

## How to use `trefle`

⚠️ `trefle` should not be incorporated into your own databases. The associations
are predictions, and we can estimate how many of them are false positives, and
how many are missing. In addition, the probability score is not a biologically
meaningful probability. Unless your database is able to accomodate these
subtlelties and convey them clearly to the user, we advise you against consuming
`trefle`. ⚠️

Contact: `timothee.poisot@umontreal.🇨🇦`

## Repository content

- `hpc` contains all the code used to run the tuning and simulation using `slurm`
    - `inputs` is the main location for the bash scripts and helper functions
    - `outputs` is where the output files are located -- note that they are not written here by default, this was us doing some post-processing
        - `tuning.csv` is the file for model selection (about 6MB)
        - `predictions.csv` is the output of imputation (about 85MB)
    - as a side-note, each thread is responsible for its own files (and works on its own copy of the data, so think about memory use)
    - as an additional side-note, not all species pairs in `clover` are in `trefle`, because some proportion (<1%) of runs fail for reasons that always mean that the association is [almost surely][almost_surely] not happening 
- `data` is storing all the data that are *not* directly generated by `trefle`
- `model_performance` has the file for model selection *and* the figures generated as part of this process
    - `roc` has all the plots of ROC-AUCs
    - `metrics` has the plots of all metrics presented in the table above
- `imputation` has the files to read the data from `hpc/outputs` and do the analyses
- `artifacts` has derived data tables
    - `modelselection.csv` is the list of all models considered during hyper-parameters tuning

[almost_surely]: https://en.wikipedia.org/wiki/Almost_surely
## Main results


## Get involved
