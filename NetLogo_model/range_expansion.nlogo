;; Individual-based model aiming to study population dynamics, dispersal evolution and neutral evolution during pushed vs pulled range expansions
;;
;; authors: Maxime Dahirel and Chloé Guicharnaud
;; building on a previous model by Maxime Dahirel and Marjorie Haond (see https://doi.org/10.5281/zenodo.3969988 or https://doi.org/10.5281/zenodo.3702252)
;; written using Netlogo versions: 6.2.0 and 6.1.1

;; -999 is used as a placeholder for NAs throughout the code

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; PART 1: define patch and individual variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals[
  ;; fixed once set up

  ;; parameters shaping the initial distributions of traits
  logit_start_dmax     ;; a placeholder global variable used to translate the initial maximal dispersal rate (entered on the probability scale as it's easier) to the logit scale
  VA_logit_dmax         ;; initial genetic additive variance of logit(dmax)
  VR_logit_dmax         ;; residual (i.e. environmental) variance of logit(dmax)
  VA_slope              ;; initial genetic additive variance of the dispersal-density reaction norm slope
  VR_slope              ;; residual (i.e. environmental) variance of the dispersal-density reaction norm slope
  VA_midpoint           ;; same for midpoint
  VR_midpoint           ;; same for midpoint

  ;; mean initial values (t=0) for some patch-level summaries (see turtles-own for meaning of each)
  start_d0
  start_dK
  start_davg_0_K
  start_slopeA_0_K
  start_slopeA_0_avg
]

turtles-own[
  ;; dispersal
  disp                    ;; individual dispersal probability
  x_birth                 ;; x coordinate of birth site; useful to know if individual has dispersed (and how far in cases dispersal>1 patch is allowed)
  dmax                    ;; maximal dispersal probability
  logit_dmax              ;; logit of dmax
  midpoint                ;; density at which dispersal probability is 50% of dmax
  slope                   ;; "slope" of the dispersal-density reaction norm
  available_moves         ;; a list of all possible movements if one disperses (in 1D space, typically -1 or +1 unless landscape boundary is reached)
  alleles_logit_dmax      ;; allele(s) of logit_dmax
  genotype_logit_dmax     ;; the genetic value of logit_dmax (average of the alleles)
  alleles_slope           ;; same for slope, but no logit this time
  genotype_slope          ;; ...
  alleles_midpoint        ;; and same for midpoint
  genotype_midpoint       ;; ...
  ;; dispersal, summary/derived traits
  d0                      ;; dispersal probability at N = 0
  davg_0_K                ;; average dispersal probability over the range 0-K (approximated by taking it every 0.1K from 0 to 1K)
  dK                      ;; dispersal probability at N = K
  slopeA_0_K              ;; difference between d0 and dK
  slopeA_0_avg            ;; difference between d0 and davg_0_K

  ;; growth and reproduction
  adult          ;; a 0/1 flag indicating if the individual is adult (used during reproductive phase)
  has_reproduced ;; a 0/1 flag indicating if the individual has reproduced (individuals can only reproduce once)
  ind_fecundity  ;; actual individual fecundity

  ;; neutral genetic diversity
  neutral_locus ;; two possible allele values (0 ; 1). Inherited with no selection; used for analyses of changes in neutral genetic diversity

  ;; misc.
  momID      ;; [sexual] unique Netlogo-created ID of the mother, useful for pedigree. [clonal] momID set and kept to the unique ID of the starting individual of the clonal line
  PgrandmaID ;; [sexual] unique ID of the paternal grandmother of the turtle, set to -999 for the first generation (because unknown pedigree). [clonal] Set and kept equal to momID
  MgrandmaID ;; [sexual] unique ID of the maternal grandmother of the turtle, set to -999 for the first generation (because unknown pedigree). [clonal] Set and kept equal to momID
]


patches-own [

  founding ;; generation/tick in which the patch was *first* populated

  ;; _pre or _predispersal suffix refers to metrics estimated immediately before the dispersal phase
  ;; _post or _postdispersal suffix refers to metrics estimated immediately after the dispersal phase

  ;; counts of individuals
  carrying_capacity ;; carrying capacity/ equilibrium pop size of the patch,
                    ;; regulation happens by making growth rate fall below 0 if population size > carrying capacity
  population_size   ;; current population size (refreshed several times through life cycle, as up to date values needed quite often)
  N_predispersal    ;; adult population size right before the dispersal phase
  N_postdispersal   ;; adult population size after the dispersal phase
  N_sedentary       ;; number of turtles who did not move
  N_disp_dead       ;; number of turtles born in this patch dead during dispersal

  ;; counts of neutral alleles
  N_allele0_pre     ;; (measured during pre-dispersal phase) number of adults with neutral_locus = 0
  N_allele1_pre     ;; (measured during pre-dispersal phase) number of adults with neutral_locus = 1
  N_allele0_post    ;; (measured during post-dispersal phase) number of adults with neutral_locus = 0
  N_allele1_post    ;; (measured during post-dispersal phase) number of adults with neutral_locus = 1

  ;; patch summaries of dispersal traits (measured pre-dispersal only, to be able to compare with %dispersal + to help ensure a variance exists)

  ;; genetic components
  mean_genotype_logit_dmax   ;; average genotypic value for logit(dmax)
  var_genotype_logit_dmax    ;; variance of genotypic values
  mean_genotype_slope        ;; you get the idea...
  var_genotype_slope
  mean_genotype_midpoint
  var_genotype_midpoint

  ;; overall trait on observed scale
  mean_dmax
  ;var_dmax
  mean_midpoint
  ;var_midpoint
  mean_slope
  ;var_slope

  ;; other summary traits
  mean_d0                    ;; mean dispersal probability at N = 0
  ;var_d0
  mean_dK                    ;; mean dispersal probability at N = K
  ;var_dK
  mean_davg_0_K
  ;var_davg_0_K
  mean_slopeA_0_K            ;; average absolute slope over the range 0-K
  ;var_slopeA_0_K
  mean_slopeA_0_avg          ;;
  ;var_slopeA_0_avg
  ]


to-report invlogit [number]
  ifelse (number <= -710)
    [ report 0 ]  ;; because exp(- -710 or more) would trigger overflow error, and 1/(1+exp(- -710 or more)) is virtually 0 anyway
    [ report (1 / (1 + exp( - number))) ]
end

to-report logit [number]
  ifelse ((number >= 0) and (number <= 1))
    [report ln (number / (1 - number) )]
    [error "input should be a proportion (>= 0 and <= 1)"]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;PART 2: SETUP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to setup
  clear-all
  define-landscape
  setup-initial-variances
  setup-patches
  setup-turtles
  setup-initial-summaries
  reset-ticks
end


to define-landscape
  set-patch-size 3       ;; graphical argument if NetLogo GUI used
  resize-world 0 300 0 0 ;; generate the correct landscape size ( xmin xmax ymin ymax) ymax = ymin = 0 for 1D landscapes
  ;; (0 = initial patch; see move_turtles for how we block individuals from trying to move left if x=0)
  ask patches [set pcolor black] ;; graphical argument if NetLogo GUI used (patches will become whiter as population size increases)
end

to setup-initial-variances
  set logit_start_dmax logit start_dmax ;; get starting logit(dmax) from input dmax (more natural to input dispersal rate on observed scale rather than logit)
  set VA_logit_dmax heritability * VP_logit_dmax       ;; set initial additive genetic variance VA from input heritability and total phenotypic variance
  set VR_logit_dmax (1 - heritability) * VP_logit_dmax ;; same with residual variance VR
  set VA_slope heritability * VP_slope
  set VR_slope (1 - heritability) * VP_slope
  set VA_midpoint heritability * VP_midpoint
  set VR_midpoint (1 - heritability) * VP_midpoint
end

to setup-patches
  ask patches [set carrying_capacity K]                   ;; all patches have the same constant carrying_capacity; extension possible to stochastic spatial variation by using (random-poisson K) instead of K
  ask patches [set founding -999]
  ask patches with [pxcor = 0] [sprout carrying_capacity] ;; we generate individuals only in one initial patch with coordinate (pxcor = 0)
  ask patches with [pxcor = 0] [set founding 0]
  ask patches [check_population_size]
end

to setup-turtles
  ask turtles[
    hide-turtle  ;; we don't visualise individuals on the GUI; we only show the patch-level summary, saves memory
    set adult 1
    set has_reproduced 0
    set ind_fecundity 0
    set x_birth xcor

    ifelse reproduction = "clonal"   ;; for clonal reproduction, only one allele for each trait and neutral locus is drawn
      [set neutral_locus (list random 2) ;; NB: important: random 2 reports 0 or 1, not 1 or 2

       set alleles_logit_dmax random-normal logit_start_dmax sqrt(VA_logit_dmax)
       set alleles_slope random-normal start_slope sqrt(VA_slope)
       set alleles_midpoint random-normal start_midpoint sqrt(VA_midpoint)
       ;; assign genotypic trait values from the global means and genetic variance

       set genotype_logit_dmax alleles_logit_dmax  ;; no need to average anything if clonal
       set genotype_slope alleles_slope
       set genotype_midpoint alleles_midpoint

       ;; at the next generation, individuals get genotypic values from parent(s), and redraw the residual noise from the random normal distribution

       ;; assigning a line ID ; for clonal reproduction, momID correspond to the unique ID of a clonal line. MgrandmaID and PgrandmaID are equal to momID
       set momID who
       set MgrandmaID who
       set PgrandmaID who
      ]
      [set neutral_locus list (random 2) (random 2)   ;; for sexual reproduction, 2 alleles are drawn, otherwise, same as clonal reproduction concerning traits values

       set alleles_logit_dmax list (random-normal logit_start_dmax sqrt(VA_logit_dmax)) (random-normal logit_start_dmax sqrt(VA_logit_dmax))
       set alleles_slope list (random-normal start_slope sqrt(VA_slope)) (random-normal start_slope sqrt(VA_slope))
       set alleles_midpoint list (random-normal start_midpoint sqrt(VA_midpoint)) (random-normal start_midpoint sqrt(VA_midpoint))

       set genotype_logit_dmax mean (alleles_logit_dmax)
       set genotype_slope mean (alleles_slope)
       set genotype_midpoint mean (alleles_midpoint)

       ;; assigning an unknown pedigree for sexual reproduction ; -999 is used as a placeholder for missing values at the start of the experiment
       set momID -999
       set MgrandmaID -999
       set PgrandmaID -999
       ]

   set_individual_traits
  ]
end

to set_individual_traits
        set logit_dmax  genotype_logit_dmax + (random-normal 0 sqrt(VR_logit_dmax))
        set dmax  invlogit logit_dmax  ;; trait = inverse logit of trait on latent scale
        set slope genotype_slope + (random-normal 0 sqrt(VR_slope))
        set midpoint genotype_midpoint + (random-normal 0 sqrt(VR_midpoint))
;; for each trait we add residual noise ; if VR = 0, there is no residual noise and the final trait value correspond to the genotypic trait value

        ;; set "secondary" traits
        set d0 (dmax * invlogit ( slope * (0 - midpoint)))
        set dK (dmax * invlogit ( slope * (1 - midpoint)))
        set slopeA_0_K (dK - d0)

        let dNs map [x -> (dmax * invlogit ( slope * (x - midpoint))) ] (range 0 1.01 0.1)
        ;; we use a sequence going from 0 to 1.01 by steps of 0.1 because weirdly the steps Netlogo adds aren't exactly 0.1 but a bit more
        ;; (ask "show map [x -> x ] (range 0 1.01 0.1)" in the command center to check)
        ;; so if we asked "from 0 to 1 by steps of 0.1", the "1" step would actually be 1.000000...001 or something, so out of boundaries, so ignored
        ;; unless the boundary is a bit above 1

        set davg_0_K mean (dNs)

        set slopeA_0_avg davg_0_K - d0

end

to check_population_size
  set population_size count turtles-here
  set pcolor scale-color green population_size 0 (carrying_capacity * 1.1) ;; graphical argument if NetLogo GUI used, updates the patch colour based on current population size
end

to setup-initial-summaries
  set start_d0 mean ([d0] of turtles)
  set start_dK mean ([dK] of turtles)
  set start_davg_0_K mean ([davg_0_K] of turtles)
  set start_slopeA_0_K mean ([slopeA_0_K] of turtles)
  set start_slopeA_0_avg mean ([slopeA_0_avg] of turtles)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;PART 3 : GO ! (a.k.a. the life cycle)
;;;;;;;;;;;;;;;;;;;;;;;;;;


to go

  if (ticks > duration) [stop] ;;ending condition

  ;; predispersal count here, used to shape density-dependent dispersal
  ask patches[check_population_size]
  ask patches[set N_predispersal population_size]

  ;; reset patch-level summaries for traits and allelic frequencies
  ask patches[reset_summaries]

  ;; record patch-level summaries for traits and predispersal allele distribution
  ask patches with [N_predispersal > 0] [update_alleles_pre]

  ask patches with [N_predispersal > 0] [update_summaries_means]

  ask patches with [N_predispersal > 1] [update_summaries_variances]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; dispersal step
  ask turtles[move_turtles]
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; postdispersal count here, now that individuals have moved, used to shape reproduction + used for census
  ask patches[check_population_size]
  ask patches[
    set N_postdispersal population_size
    set N_sedentary count turtles-here with [x_birth = pxcor]
  ]

  ;; record patch-level allelic frequencies for the *neutral* locus post-dispersal; easier to store for export
  ask patches with [N_postdispersal > 0] [update_alleles_post]

  ask patches with [N_postdispersal > 0 and founding < 0][set founding (ticks + 1)]
  ;; if patch is populated now, but wasn't yet (founding = -999), then set founding date
  ;; (which is ticks + 1 because the ticks number is updated to n only **at the end** of the nth generation

    ;; reproduction step
    ( ifelse reproduction = "clonal"
      [ask turtles[reproduce_clonal]]
      [ask turtles[reproduce_sexual]]
      )
    ;; NB the reproduction formulas implicitly include competition, i.e. we only "materialize" offspring that will reach the adult stage
    ask turtles[check_death] ;; enforce non-overlapping generations: kill all adults and leave only juveniles produced in previous round
    ask turtles[set adult 1] ;; once previous adults are removed, juveniles can become adults

  tick
end



to move_turtles ;; dispersal
  set available_moves [1 -1] ;; individuals can only move left or right

  ;; but if animal is at landscape border can only disperse in one direction (but still with the same total dispersal probability overall)
  ;; only relevant for the starting patch in normal use, as usually landscape size is set to be larger than the total movements possible for the duration of the run :
  if xcor = max-pxcor [set available_moves [-1]]
  if xcor = min-pxcor [set available_moves [1]]

  set disp (dmax * invlogit ( slope * ((population_size / carrying_capacity) - midpoint)))
  ;; sets the individual dispersal probability based on its trait and pre-dispersal population size
  ;; based on kun and scheuring 2006 doi:10.1111/j.2006.0030-1299.15061.x

  if ( (random-float 1) < (disp) ) [    ;; sets whether the individual moves (this should be/approximate a Bernoulli distribution with probability of success disp)
    if (random-float 1) < (dispersal_mortality) [
      set N_disp_dead (N_disp_dead + 1)
      die]                                      ;; mortality during dispersal
    set xcor xcor + one-of available_moves]     ;; sets where it moves if it decides to move and survives dispersal

end


to reproduce_sexual  ;; sexual reproduction, no mutation
  if (has_reproduced = 0 and adult = 1 ) [ ;; safety to avoid multiple reproductions.
    ;; avoid multiple matings as partners of several "focals"
    let mom self
    let mate nobody
    set mate one-of ( (other turtles-here) with [has_reproduced = 0 and adult = 1]) ;; choose one mate among the unmated individuals of the patch


    if mate != nobody [ ;; if they find a mate
      set ind_fecundity random-poisson exp(ln(fecundity) * (1 - population_size / carrying_capacity) )
      ;; ricker function ; for each individual, fecundity is Poisson-distributed around its mean fecundity
      hatch ind_fecundity [ ;; mom reproduce
        hide-turtle         ;; needs to hide again newborn individuals ;; we don't visualise individuals on the GUI; we only show the patch-level summary, saves memory

        set adult 0
        set has_reproduced 0
        set ind_fecundity 0
        set x_birth xcor

        ;; neutral alleles
        set neutral_locus list (one-of [neutral_locus] of mom) (one-of [neutral_locus] of mate)
        ;; in case of sexual reproduction, one allele comes from mom, another from mate

        ;; misc.
        set momID [who] of mom
        set PgrandmaID [momID] of mate
        set MgrandmaID [momID] of mom

        ;; trait determination
        set alleles_logit_dmax list (one-of [alleles_logit_dmax] of mom) (one-of [alleles_logit_dmax] of mate)
        set alleles_slope list (one-of [alleles_slope] of mom) (one-of [alleles_slope] of mate)
        set alleles_midpoint list (one-of [alleles_midpoint] of mom) (one-of [alleles_midpoint] of mate)

        set genotype_logit_dmax mean alleles_logit_dmax
        set genotype_slope mean alleles_slope
        set genotype_midpoint mean alleles_midpoint
        ;; genotypic value inherited from parent(s)

        set_individual_traits
      ]

      ask mate [set ind_fecundity random-poisson exp(ln(fecundity) * (1 - population_size / carrying_capacity) )]

      hatch [ind_fecundity] of mate [ ;; the mate also reproduces
        hide-turtle                   ;; needs to hide again newborn individuals ;; we don't visualise individuals on the GUI; we only show the patch-level summary, saves memory

        set adult 0
        set has_reproduced 0
        set ind_fecundity 0
        set x_birth xcor

        ;; neutral alleles
        set neutral_locus list (one-of [neutral_locus] of mom) (one-of [neutral_locus] of mate)
        ;; in case of sexual reproduction, one allele comes from mom, another from mate

        ;; misc.
        set momID [who] of mate
        set PgrandmaID [momID] of mom
        set MgrandmaID [momID] of mate

        ;; trait determination
        set alleles_logit_dmax list (one-of [alleles_logit_dmax] of mom) (one-of [alleles_logit_dmax] of mate)
        set alleles_slope list (one-of [alleles_slope] of mom) (one-of [alleles_slope] of mate)
        set alleles_midpoint list (one-of [alleles_midpoint] of mom) (one-of [alleles_midpoint] of mate)

        set genotype_logit_dmax mean alleles_logit_dmax
        set genotype_slope mean alleles_slope
        set genotype_midpoint mean alleles_midpoint
        ;; genotypic value inherited from parent(s)

        set_individual_traits
      ]

      set has_reproduced 1
      ask mate [set has_reproduced 1] ;; the mate cannot be chosen/reproduce again
    ]
  ]
end

to reproduce_clonal  ;; clonal reproduction, no mutation
  if (has_reproduced = 0 and adult = 1 ) [
    let mom self
    set ind_fecundity random-poisson exp(ln(fecundity) * (1 - population_size / carrying_capacity) )
      ;; ricker function ; for each individual, fecundity is Poisson-distributed around its mean fecundity

      hatch ind_fecundity [
        hide-turtle  ;; needs to hide again newborn individuals
                     ;; we don't visualise individuals on the GUI; we only show the patch-level summary, saves memory

        set adult 0
        set has_reproduced 0
        set ind_fecundity 0
        set x_birth xcor

        ;;neutral alleles
        set neutral_locus [neutral_locus] of mom

        ;;misc.
        set momID [momID] of mom
        set MgrandmaID [momID] of mom
        set PgrandmaID [momID] of mom


        ;;inheritance and trait determination
        set alleles_logit_dmax [alleles_logit_dmax] of mom
        set alleles_slope  [alleles_slope] of mom
        set alleles_midpoint [alleles_midpoint] of mom

        set genotype_logit_dmax alleles_logit_dmax
        set genotype_slope alleles_slope
        set genotype_midpoint alleles_midpoint

        set_individual_traits
      ]
      set has_reproduced 1
    ]
end

to reset_summaries
;; reset all patches with our NA placeholder to avoid carryovers from previous gen
  set mean_genotype_logit_dmax -999
  set var_genotype_logit_dmax  -999

  set mean_genotype_slope -999
  set var_genotype_slope  -999

  set mean_genotype_midpoint -999
  set var_genotype_midpoint  -999

  set mean_dmax -999
  set mean_midpoint -999
  set mean_slope -999

  set mean_d0 -999
  set mean_dK -999
  set mean_davg_0_K -999
  set mean_slopeA_0_K -999
  set mean_slopeA_0_avg -999

  ;set var_dmax -999
  ;set var_midpoint -999
  ;set var_slope -999

  ;set var_d0 -999
  ;set var_dK -999
  ;set var_davg_0_K -999
  ;set var_slopeA_0_K -999
  ;set var_slopeA_0_avg -999

  set N_disp_dead 0
  set N_allele0_pre 0
  set N_allele1_pre 0
  set N_allele0_post 0
  set N_allele1_post 0
end

to update_alleles_pre
    set N_allele1_pre sum (reduce sentence ([neutral_locus] of turtles-here)) ;; reduce sentence used to "collapse" all the alleles stored in individual lists in one big list we can sum
    ( ifelse reproduction = "clonal"
      [set N_allele0_pre population_size - N_allele1_pre]
      [set N_allele0_pre (2 * population_size) - N_allele1_pre]               ;; for sexual reproduction, twice as many alleles than individuals
    )
end

to update_summaries_means
 ;; mean trait values
  set mean_genotype_logit_dmax mean ([genotype_logit_dmax] of turtles-here)
  set mean_genotype_slope mean ([genotype_slope] of turtles-here)
  set mean_genotype_midpoint mean ([genotype_midpoint] of turtles-here)

  set mean_dmax mean ([dmax] of turtles-here)
  set mean_midpoint mean ([midpoint] of turtles-here)
  set mean_slope mean ([slope] of turtles-here)

  set mean_d0 mean ([d0] of turtles-here)
  set mean_dK mean ([dK] of turtles-here)
  set mean_davg_0_K mean ([davg_0_K] of turtles-here)
  set mean_slopeA_0_K mean ([slopeA_0_K] of turtles-here)
  set mean_slopeA_0_avg mean ([slopeA_0_avg] of turtles-here)
end

to update_summaries_variances
 ;; trait variances; sample variances only meaningful for N > 1
  set var_genotype_logit_dmax  variance ([genotype_logit_dmax] of turtles-here)
  set var_genotype_slope  variance ([genotype_slope] of turtles-here)
  set var_genotype_midpoint  variance ([genotype_midpoint] of turtles-here)

  ;set var_dmax variance ([dmax] of turtles-here)
  ;set var_midpoint variance ([midpoint] of turtles-here)
  ;set var_slope variance ([slope] of turtles-here)

  ;set var_d0 variance ([d0] of turtles-here)
  ;set var_dK variance ([dK] of turtles-here)
  ;set var_davg_0_K variance ([davg_0_K] of turtles-here)
  ;set var_slopeA_0_K variance ([slopeA_0_K] of turtles-here)
  ;set var_slopeA_0_avg variance ([slopeA_0_avg] of turtles-here)

end

to update_alleles_post
    set N_allele1_post sum (reduce sentence ([neutral_locus] of turtles-here)) ;; reduce sentence used to "collapse" all the alleles stored in individual lists in one big list we can sum
    ( ifelse reproduction = "clonal"
      [set N_allele0_post population_size - N_allele1_post]
      [set N_allele0_post (2 * population_size) - N_allele1_post]              ;; for sexual reproduction, twice as many alleles than individuals
    )
end

to check_death
    if adult = 1 [ die ]
  ;; Note: we check death based on the "adult" flag rather than the "has_reproduced" flag to accomodate the scenario with sexual reproduction
  ;; (as not all adults may be able to find a mate)
end
@#$#@#$#@
GRAPHICS-WINDOW
13
29
924
41
-1
-1
3.0
1
10
1
1
1
0
1
1
1
0
300
0
0
1
1
1
ticks
30.0

BUTTON
30
85
93
118
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
29
129
92
162
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
38
281
210
314
K
K
10
1000
250.0
10
1
NIL
HORIZONTAL

SLIDER
239
324
410
357
VP_logit_dmax
VP_logit_dmax
0
3
1.52
0.01
1
NIL
HORIZONTAL

SLIDER
237
174
409
207
start_dmax
start_dmax
0.01
0.99
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
441
174
613
207
fecundity
fecundity
0
20
1.5
0.1
1
NIL
HORIZONTAL

SLIDER
237
221
409
254
start_slope
start_slope
-10
10
0.0
0.5
1
NIL
HORIZONTAL

SLIDER
239
362
412
395
VP_slope
VP_slope
0
50
4.9
0.1
1
NIL
HORIZONTAL

SLIDER
38
328
210
361
duration
duration
0
1000
100.0
10
1
NIL
HORIZONTAL

PLOT
659
325
859
475
count individuals
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles"

SLIDER
441
220
613
253
heritability
heritability
0
1
0.32
0.01
1
NIL
HORIZONTAL

CHOOSER
29
175
167
220
reproduction
reproduction
"clonal" "sexual"
0

SLIDER
238
471
410
504
dispersal_mortality
dispersal_mortality
0
1
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
237
265
409
298
start_midpoint
start_midpoint
-0.5
1.5
0.51
0.01
1
NIL
HORIZONTAL

SLIDER
237
404
411
437
VP_midpoint
VP_midpoint
0
5
0.27
0.01
1
NIL
HORIZONTAL

@#$#@#$#@
# Range expansion model

## WHAT IS IT?

This is a general spatial range expansion model, designed to study what leads to pushed versus pushed range expansions, and the possible evolutionary consequences.

This model is designed to operate at "low" densities (equilibrium population size < 1000), at least several orders of magnitude lower than the one classically used in the theoretical literature, which are adapted to microbial species, but may not be applicable to macroscopic species with lower population sizes.

From an initial population, haploid individuals reproduce, compete and disperse and as a result the species colonizes the landscape.

The model allows for density-dependency in dispersal and growth, as well as individual variation in dispersal and growth.

When there is individual variation in traits, evolution can be activated (individuals inherit their parental allele) or turned down (trait values are drawn at random every generation).

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

### Set-up phase
The model is initiated by placing K adult individuals in the patch of coordinates pxcor 0 and pycor 0, and setting their phenotypic traits from the distribution described by parameters, and by drawing the allele at the neutral locus from a Bernoulli distribution with p = 0.5.

### Go phase
Individuals then live the following life cycle:

-once adults, they disperse or not based on their dispersal-density reaction norm and current patch population size

-they reproduce clonally and transmit their neutral allele (and in the evolutionary setting trait alleles to offspring
the fecundity formula directly gives the number of offspring post-competition in one step, to avoid wasting computing power by creating individuals that would then be killed

-they die


## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

### Global parameters

The following global parameters can be set up in the Interface tab:

- *trait_variation*: two settings; whether or not phenotypic traits are re-drawn at random at each generation, blocking evolutionary change
- *K*: (average) carrying capacity/ equilibrium population density
- *duration*: duration of a run, in generations
- *disp0_mean*: average dispersal rate at the start of the run
- *slope_disp_mean*: slope of the relationship between relative density (population size/K) and dispersal rate **on the logit scale**
- *logit_disp0_sd*: standard deviation of the distribution of dispersal rates at the start of the run **on the logit scale**
- *slope_disp_sd*: standard deviation of the distribution of dispersal reaction norm slopes at the start of the run **on the logit scale**
- *fecundity*: average fecundity at the start of the run

### Landscape window

The landscape window shows the progression of the current wave through the landscape. The whiter the patch, the larger the population, empty patches are black.

### Report graphs

Because the model is designed to be analysed through R (package nlrx), there is only one graph in the Interface tab, showing the total metapopulation size through time. But more can easily be added
(suggestions: front position (or maximal x coordinate of patch with at least 1 individual) through time
Genetic diversity (formula for expected heterozygosity when 2 alleles = 2pq where p and q are allelic proportions) in front patch through time)



## THINGS TO NOTICE

(suggested things for the user to notice while running the model)


- to do

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

- to do

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

More complex evolutionary dynamics:
- dispersal/growth trade-offs
- heritability values different from 0 (implied in current "reshuffled" mode) or 1 (implied in current "evolutionary" mode)

Landscape heterogeneity in space and or time

Invasions in 2D space

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

- nothing especially unusual, but see the workarounds used around the fact there is no native random-binomial and logit fonctions in Netlogo (look at how dispersal rates are set)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)

written by Maxime Dahirel, from an initial Matlab model by Marjorie Haond
DOI: https://doi.org/10.5281/zenodo.3702252
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
