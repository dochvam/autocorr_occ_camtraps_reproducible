# input_data

This directory contains several sources of data needed for running the code. It contains the following files:

-   subdirectory `snapshot2020/`
    -   This subdirectory contains a number of .RDS files containing the Snapshot USA data used in the case study. These data have been processed and filtered according to the criteria in the manuscript. **Any use of these data should cite the Snapshot 2020 dataset under the DOI: 10.1002/ecy.3775**
    -   Each .RDS file contains an R object of the class `unmarkedFrameOccu` for use with the unmarked `occu()` function. To learn more about this object type, see `?unmarked::unmarkedFrameOccu`. Its elements are:
        -   y: a detection history matrix (row = camera deployment, col = 1-day sampling occasion)
        -   siteCovs: a matrix of scaled occupancy covariates (forest cover and maximum temperature)
        -   obsCovs: a matrix of scaled detection covariates (canopy height and log distance to road)
-   litreview_results.csv: a .csv file containing the results of the literature review of detection windows described in the manuscript. Each row of the dataset represents one of the 100 randomly selected papers returned by a literature search. Its columns are:
    -    Title: the title of the paper
    -   DOI: the paper's DOI
    -   year: the year in which the paper was published
    -   Qualified: whether or not the paper qualified for review (it was a paper using a hierarchical occupancy model or occupancy model extension to analyze camera trap data)
    -   survey_duration_days: for qualifying papers, roughly how long was the average camera deployment?
    -   detection_window_days: for qualifying papers, what length of time (in days) was used to define sampling occasions?
    -   detection_gap_days: for qualifying papers, what length of time (in days) was used to set gaps between sampling occasions?
-   casestudy_target_specs.csv: a file listing the species used in the case study containing some summary data regarding each species. Each species listed has an associated .RDS data file in the `snapshot2020` subdirectory. Its columns are:
    -   sciname: the scientific name of the species
    -   common_name: the common name of the species
