# Reproducing the tamoxifen trial results

Both files are complete, standalone scripts using public `collett::tamoxifen`
data, where A=1 is tamoxifen + radiotherapy and A=0 is tamoxifen alone.
They prioritize death over the first observed local, axillary, or distant
relapse, with death first for same-day ties. Both use five-fold
cross-fitting (seed 2027), baseline age10, log_size, and HR, and full
estimated censoring augmentation q.

## Main manuscript analysis

`tamoxifen_penalized_spline_fullq.R`: penalized cubic B-spline
baseline log-hazards on a 16-cell grid, with 1 interior knot,
spline penalty lambda=8, M_q=256, at 4- and 5-year horizons.

The spline script contains preliminary definitions for the exponential
working-hazard analysis; the **later spline section overrides those
functions** before the final driver runs. Always source the complete
file to reproduce spline results.

## Supporting analysis

`tamoxifen_exponential_fullq.R`: exponential working-hazard models,
not the manuscript's main spline specification. Retained as supporting
reproducibility code.

## Usage after package installation

```r
install.packages("collett")
library(winoddsEIF)
file <- system.file("analysis", "tamoxifen_penalized_spline_fullq.R",
                    package = "winoddsEIF", mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file, envir = env)
print(env$spline_results$table)
```

The script contains top-level commands and **runs the full analysis on
sourcing**; results (CSV and RDS) are saved in the current working directory.
It is not run on package installation or attachment. Check `getwd()` to
choose where outputs go.

For the exponential working-hazard script, use a *fresh environment*:

```r
file <- system.file("analysis", "tamoxifen_exponential_fullq.R",
                    package = "winoddsEIF", mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file, envir = env)
print(env$results)
```

The scalar and semi-competing-risk numerical validations for the public
package functions do not independently verify the spline analysis.
The real-data spline results require separate numerical checks.
