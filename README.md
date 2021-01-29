# A data-inflated host-virus association database

C'est quoi, `trefle`? It is a data product derived from the [`clover`][clover]
database of mammals-virus interaction. Specifically, `trefle` was produced using
LF-SVD imputation, a two-step algorithm where novel host-virus associations are
recommended based on truncated singular value decomposition applied to initial
values based on a linear filter.

[clover]: https://github.com/viralemergence/clover