# RegEnrich 0.99.10
* Remve `\dontrun` tags from all examples.
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
