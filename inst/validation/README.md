# Public estimator numerical validations

The user successfully executed these comparisons in RStudio on
moderate and strongest scenarios (scalar: nine metrics; semi-competing:
six metrics, with zero displayed differences).
To rerun after installing the package:

```r
source(system.file("validation", "validate_winoddsEIF_scalar.R",
                   package = "winoddsEIF", mustWork = TRUE))
source(system.file("validation", "validate_winoddsEIF_semicomp.R",
                   package = "winoddsEIF", mustWork = TRUE))
```

These scripts execute immediately on sourcing. They compare against the
archived simulations using the same datasets, fold assignments, and settings.
They do not verify the standalone tamoxifen spline script.
