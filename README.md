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

## Model performance

### Top 10 models

| model        | rank  |  AUC       |  cutoff    |  TPR      | TNR      |  PPV       | NPV       | FNR      | FPR       | FDR       | FOR      | CSI      | ACC      | J         |
| ------------ |------ |----------- | ---------- | --------- |--------- |----------- |---------- |----------|-----------|-----------|----------|----------|----------|-------- |
| connectance  |    12 |  0.849131  |  0.846847  |  0.72     | 0.925926 |  0.90604   | 0.769231  | 0.28     | 0.0740741 | 0.0939597 | 0.230769 | 0.669975 | 0.823373 | 0.645926 |
| connectance  |    11 |  0.846133  |  0.908909  |  0.684421 | 0.936759 |  0.914591  | 0.75      | 0.315579 | 0.0632411 | 0.0854093 | 0.25     | 0.643304 | 0.811258 | 0.62118 |
| connectance  |    17 |  0.844302  |  0.92993   |  0.692102 | 0.935271 |  0.913428  | 0.754797  | 0.307898 | 0.0647292 | 0.0865724 | 0.245203 | 0.649497 | 0.814495 | 0.627373 |
| connectance  |     8 |  0.842319  |  0.705706  |  0.701473 | 0.895225 |  0.868988  | 0.75167   | 0.298527 | 0.104775  | 0.131012  | 0.24833  | 0.634383 | 0.798801 | 0.596698 |
| hybrid       |    12 |  0.841677  |  0.707708  |  0.703209 | 0.877822 |  0.851133  | 0.748584  | 0.296791 | 0.122178  | 0.148867  | 0.251416 | 0.62619  | 0.790806 | 0.581031 |
| connectance  |    14 |  0.839073  |  0.902903  |  0.700269 | 0.929708 |  0.907666  | 0.758658  | 0.299731 | 0.0702918 | 0.0923345 | 0.241342 | 0.653701 | 0.815754 | 0.629977 |
| hybrid       |    11 |  0.837223  |  0.820821  |  0.647137 | 0.918991 |  0.888483  | 0.723093  | 0.352863 | 0.0810093 | 0.111517  | 0.276907 | 0.598522 | 0.783245 | 0.566128 |
| connectance  |     5 |  0.83667   |  0.931932  |  0.660453 | 0.940318 |  0.916821  | 0.735477  | 0.339547 | 0.0596817 | 0.0831793 | 0.264523 | 0.623116 | 0.800664 | 0.600771 |
| connectance  |     7 |  0.83667   |  0.948949  |  0.655218 | 0.957784 |  0.939394  | 0.735562  | 0.344782 | 0.0422164 | 0.0606061 | 0.264438 | 0.628644 | 0.806601 | 0.613002 |
| connectance  |    16 |  0.835764  |  0.961962  |  0.66756  | 0.945623 |  0.923933  | 0.741935  | 0.33244  | 0.0543767 | 0.0760668 | 0.258065 | 0.632783 | 0.807333 | 0.613184 |

### Overview of the best model
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

## Main results

## Get involved
