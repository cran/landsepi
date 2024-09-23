## landsepi 1.5.1 (2024-09-23)

* fixed bugs

## landsepi 1.5.0 (2024-09-02)

* Redesigned shiny app: simplified parameterization, include global/local inoculum, include more outputs (graphics & tables)
* Updated webpage

## landsepi 1.4.0 (2024-03-12)

* [**breaking change**] Changed parameterization of inoculum to allow specific spatial location and genetic composition of inoculum (thanks to Manon)
* Changed parameterization of pathogen survival probability to allow different probabilities for different croptypes and years
* Added 2 parameters (`adaptation_cost` and `relative_advantage`) to allow two types of fitness costs paid by pathogens adapted to resistance (on all hosts, or only for unnecessary virulence) 
* Help parameterization of treatments with `loadTreatment()`
* Updated documentation and vignettes

## landsepi 1.3.0 (2023-07-19)

* [**breaking change**] Added parameterization of grapevine mildew (thanks to Marta)
* [**breaking change**] Added parameterization of black sigatoka of banana
* [**breaking change**] Added contact chemical treatments
* Changed the dynamics of host growth
* Added the output `contribution`
* Updated documentation and vignettes


## landsepi 1.2.2 (2022-10-06)

* [**breaking change**] Added sexual reproduction with gene recombination (thanks to Marta)
* [**breaking change**] `loadDispersalPathogen()` now return a list of two matrix, the first for clonal and the second for sexual reproductions
* Changed the mutation process
* Updated Shiny app
* Added a `NEWS.md` file to track changes to the package
* Updated documentation and vignettes


## landsepi 1.1.0 (2021-07-26)

* Updated Shiny app
* Updated documentation & vignettes
* Added economic parameters and outputs


## landsepi 1.0.1 (2020-07-02)

* [**breaking change**] Added simulation of real landscapes using shapefiles
* [**breaking change**] Created a Shiny app (thanks to Jean-Loup)
* Created vignettes (tutorials)
* Added multiple epidemiological and evolutionary outputs
* Added Adult Plant Resistance (APR) genes


## landsepi 0.0.2 (2018-03-26)

* Creation of the package
* Parameterization for rusts of cereal crops
