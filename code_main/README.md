# `code_main`

This directory contains the main .R functions that execute,

The following four scripts execute the simulations and the case study. They should all be executable out-of-the-box:

-   run_simulation1.R: execute simulation exercise 1 (exploration of parameter space)

-   run_simulation2.R: execute simulation exercise 2 (changing detection windows, gaps)

-   run_simulation34.R: execute simulation exercises 3 and 4 (goodness-of-fit and evaluating the clustered model)

-   run_snapshot_casestudy.R: execute the case study based on SnapshotUSA 2020 data (see input_data/README to learn about the input data provided)

The following two functions process, summarize, and visualize the output of the preceding four scripts. These will only work if the preceding four scripts have been executed locally. They must be run in order (i.e. vis depends on summaries):

-   summarize_output.R: processes the raw simulation output into useful data formats

-   vis.R: creates the plots from the manuscript
