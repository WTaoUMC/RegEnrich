# RegEnrich 0.99.17
* Fix bugs of regenrich_diffExpr function.
* Add an example for %>%.
* Use Authors@R [cre] designation.
* Update R version dependency from 4.0 to 4.1.

# RegEnrich 0.99.16
* Import magrittr

# RegEnrich 0.99.15
* Fix a bug in `regenrich_rankScore` function.
* Reexport pipe `%>%`.

# RegEnrich 0.99.14
* Add `\donttest` tags.
* Optimize `COEN` function.
* Replace `bplapply` by `lapply` in `pickSoftThreshold2` function.

# RegEnrich 0.99.13
* Fix bugs of `plotSoftPower` on Windows.
* Hide documentations of `COEN` and `GRN` functions.
* Remove `doParallel` package from imports.

# RegEnrich 0.99.12
* Remove .Rpoj file.

# RegEnrich 0.99.11
* Remove `\dontrun` tags from all examples.
* Remove lazyData: true from `DESCRIPTION`.
* Formatting vignettes using `BiocStyle`.
* Replace `NEWS` file by `NEWS.md` file to track changes to the package.
* Some packages in `Imports` are moved to `Depends` in `DESCRIPTION` file.
* `DeaSet` class inherits `SummarizedExperiment` class.
* `RegenrichSet` class inherits `DeaSet` class.
* `TopNetwork` class inherits `BiocSet` class.
* `topResult` and `allResult` slots of `Enrich` object are `tibble` rather than `data frame`.
* `Score` class inherits `tibble`.
* The `show` methods for `DeaSet`, `TopNetwork`, `Enrich`, `Score`, and `RegenrichSet` object have been optimized.


# RegEnrich 0.99.1
* Submission to Bioconductor. 
