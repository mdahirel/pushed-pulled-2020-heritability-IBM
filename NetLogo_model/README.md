# A one-dimensional range expansion model with dispersal evolution

(this README is a copy of the `Info` pane of the NetLogo model)

written by Maxime Dahirel and Chloé Guicharnaud, building on and expanding the basic simulation model in Dahirel et al. (2021)

Please note that this model is made to be used and analyzed through R (versions 4.1.0 and higher, R Core Team 2021), using the `nlrx` package (Salecker et al. 2019). As a result, while code from this NetLogo model can be run as a standalone, it is not fully independent from the rest of its repository.

## WHAT IS IT?

This is a general spatial range expansion model, designed to study how trait distribution at the start of an expansion can shape evolutionary dynamics, and therefore influencing expansion velocity on "small" population sizes (small here is meant by contrast to the functionally infinite sizes in many theoretical models).

We specifically focus on traits driving the position of expansions on the pushed/pulled continuum (Birzu et al., 2018, 2019), namely traits describing the shape of the density-dispersal function. 

The model allows for highly flexible density-dependency as well as individual variation in dispersal. Growth rates decline with increased population sizes (Allee effects are *not* supported). Both dispersal and growth are stochastic. Traits heritability and phenotypic variances values can also be changed. Reproduction can be set to be clonal or sexual.

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

### Set-up phase
The model is initiated by introducing *K* adult individuals that are yet to reproduce or disperse in the patch of coordinates `pxcor` 0 and `pycor` 0, and by drawing their allelic values at the neutral locus from a Bernoulli distribution with *p* = 0.5. Their three phenotypic dispersal traits are initialized by summing a genetic and a noise component both drawn from Normal distributions.

### Go phase
Individuals then live the following life cycle:

- once adults, they disperse or not with a probability *d*, depending on their individual traits and current patch population size. This dispersal function is directly based on Kun and Scheuring (2006). It is basically a three parameter logistic function with two parameters (here termed `slope` and `midpoint`) influencing the shape of the function, and one parameter (here `dmax`) being the asymptote and describing the maximum possible dispersal rate for the individual.

- they reproduce clonally or sexually and transmit their neutral allele(s) and genetic values to offspring without any mutation. If they reproduce sexually, they choose a partner at random among individuals in the same patch (all are hermaphrodite). The noise components for each traits are drawn again from Normal distributions. The number of offspring is drawn from a Poisson distribution, with the mean based on a Ricker model (and thus depending on local population density). As a simplification, the result is directly the number of offspring post-competition in one step (no explicit juvenile competition phase), to avoid wasting computing power by creating individuals that would then be killed.

- they die

Meanwhile the model records various patch-level summaries (about genetic diversity, mean and variances of trait values, population sizes), to be exported and used in further analyses.


## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

### Global parameters

The following global parameters can be set up in the Interface tab:

- `K`: (average) carrying capacity/ equilibrium population density
- `duration`: duration of a run, in generations
- `fecundity`: average growth rate at the start of the run
- `heritability`: **initial** traits heritability. For simplicity, we assume all three dispersal traits have the same
- `start_dmax`: **initial** average maximal dispersal probability
- `start_slope`: **initial** average relative slope parameter of the dispersal-density function
- `start_midpoint`: **initial** average midpoint (density at which dispersal probability *d* is 50% of `dmax`) of the dispersal-density function
- `VP_logit_dmax`: **initial** phenotypic variance for `dmax` (on a logit scale)
- `VP_slope`: **initial** phenotypic variance for `slope` parameter
- `VP_midpoint`: **initial** phenotypic variance for `midpoint`
- `dispersal mortality`: mortality during dispersal (conditional on the individual actually dispersing)
- `reproduction`: either clonal or sexual (assume hermaphrodite species) reproduction

### Landscape window

The landscape window shows the progression of the current wave through the 1D landscape. The whiter the patch, the larger the population, empty patches are black.

### Report graphs

Because the model is designed to be analysed through R (package `nlrx`), there is only one graph in the Interface tab, showing the total metapopulation size through time. But more can easily be added (see suggestions in "extending the model")


## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

There is not many things to notice in the interface, as the model is not meant to be analyzed through it, but it is always possible to manually vary the initial parameter combinations and observe a qualitative expansion success or failure.


## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

Even if the model is made to be analyzed with R, and we invite you to run it through the `nlrx` package, you can try to see what happen when reproduction is set to clonal or sexual, or compare expansion with different levels of heritability


## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

- More complex evolutionary dynamics, dispersal/growth trade-offs

- Landscape heterogeneity in space and or time

- Invasions in 2D space

- Add graphs: 
	- front position (or maximal *x* coordinate of patch with at least 1 individual) through time
	- Genetic diversity (formula for expected neutral heterozygosity when 2 alleles = 2*pq* where *p* and *q* are allelic proportions) in front and/or core patches through time

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

- nothing especially unusual, but see the workarounds used around the fact there is no native random-binomial and logit fonctions in Netlogo (look at how dispersal rates are set)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)

GitHub repository for the model: https://github.com/mdahirel/pushed-pulled-2020-heritability-IBM

Birzu, G., Hallatschek, O., & Korolev, K. S. (2018). Fluctuations uncover a distinct class of traveling waves. *Proceedings of the National Academy of Sciences*, 115: E3645–E3654. https://doi.org/10.1073/pnas.1715737115

Birzu, G., Matin, S., Hallatschek, O., & Korolev, K. S. (2019). Genetic drift in range expansions is very sensitive to density dependence in dispersal and growth. *Ecology Letters*, 22: 1817–1827. https://doi.org/10.1111/ele.13364

Dahirel, M., Bertin, A., Haond, M., Blin, A., Lombaert, E., Calcagno, V., Fellous, S., Mailleret, L., Malausa, T. and Vercken, E. (2021), Shifts from pulled to pushed range expansions caused by reduction of landscape connectivity. *Oikos*, 130: 708-724. https://doi.org/10.1111/oik.08278

Kun, Á., & Scheuring, I. (2006). The evolution of density-dependent dispersal in a noisy spatial population model. *Oikos*, 115: 308–320. https://doi.org/10.1111/j.2006.0030-1299.15061.x

R Core Team. (2021). *R: a language and environment for statistical computing* [Computer software]. R Foundation for Statistical Computing. https://www.R-project.org/

Salecker, J., Sciaini, M., Meyer, K. M., & Wiegand, K. (2019). The nlrx R package: a next-generation framework for reproducible NetLogo model analyses. *Methods in Ecology and Evolution*, 10: 1854–1863. https://doi.org/10.1111/2041-210X.13286
