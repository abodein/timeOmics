Changes in version 0.99.1 (2020-04-03)
+ Submitted to Bioconductor

Changes in version 1.1.2 (2020-06-03)
+ fix `left_join` issue

Changes in version 1.1.3 (2020-06-26)
+ fix `tuneCluster.spls` columns drop

Changes in version 1.1.4 (2020-08-19)
+ add block renaming option in `plotLong`

Changes in version 1.1.5 - 1.1.8 (2020-10-15)
+ fix `lmms` removed from CRAN
+ temporary fix before including `lmmSpline()` in timeOmics

Changes in version 1.3.1 (2020-12-08)
+ fix factor level in silhouette (change to R 4.1)

Changes in version 1.3.2 (2021-04-12)
+ add filters to `getCluster()`

Changes in version 1.3.4 (2021-04-13)
+ add *UpDown clustering* (based on the sign of variation between 2 timepoints) 

Changes in version 1.5.2 (2021-08-09)
+ update vignette

Changes in version 1.7.1 (2021-10-06)
+ Add title parameter to `plot.ncomp.tune.silhouette()`

Changes in version 1.7.2 (2022-03-01)
+ Update citation info
+ fix bug with `tuneCluster.spca()`: empty choice.keepX returned when no features selected on the component (pos or neg).
    Now, if any feature is selected in pos or neg, the function returns the minimum value.
    
Changes in version 1.9.1 (2022-10-04)
+ fix is() multiple tested classes replaced by class(...); bioc 3.16

Changes in version 1.11.1 (2023-03-16)
+ Remove `propr` dependencies (removed from CRAN)
