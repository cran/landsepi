# landsepi 1.4.0

* Added Inoculum strategies (Thanks to Manon)

# landsepi 1.2.0

* Added Sexual Reproduction model (Thanks to Marta)

* Update Shiny interface
* Added a `NEWS.md` file to track changes to the package.
* [**breaking change**] `loadDispersalPathogen()` now return a list of two matrix, the first for clonal and the second for sexual reproductions.


# landsepi 1.3.1

* [**breaking change**] Changed parameterisation of inoculum to allow specific spatial location and genetic composition of inoculum
* Changed parameterisation of pathogen survival probability to allow different probabilities for different croptypes and years
* Added two types of fitness costs paid by pathogens adapted to resistance (on all hosts, or only for unnecessary virulence)
* Help parameterisation of treatments with `loadTreatment()`