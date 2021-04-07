;; This model is on the interaction of salps and krill. Salps can occurr in three life stages Oozoids, blastozoids and chains and Krill can occurr as adult krill or larval krill
;; larval inclduing juvenile krill < 30 mm will stay active all year round while adult krill will hibernate. The larval stage is not actively represented and only relevant in the growing process
;; To model larvae cohorts are modelled called clutches to avoid computational problems


breed [oozoids oozoid ]        ; one life from of salps that reproduces asexually by releasing chains
breed [blastozoids blastozoid] ; one life form of salps - here the male blastozoids are called blastozoids in the model
breed [chains chain]           ; chains are chains of femlae blastozoids
breed [debkrill sdebkrill]     ; krill modelled individually
breed [clutches clutch]        ; eggs and early stages modelled as cohorts to reduce computational effort

patches-own [chla]

globals [
  gr-o
  gr-b
  maxn
  monthlycount
  rc
  rc_max
  ratio-list
  fitness
  maxlo
  maxlb
  temp_factor
  chlaK
  peakabundance
  peakabundancelist
  max_gen
  n_oozoids
  n_blasto
  resolution
  immigrationtime
  immigrationtimelist
  survivalevents
  regenerationtimelist
  peakimmigrationtime
  migrationevent?
  T
  debsurvival
  maxspawn
]

oozoids-own [
  number                   ; always 1
  l                        ; body length in cm
  age                      ; age in days
  oozoid-repro-events      ; counts number of chain releases
  starvation-days          ; counts days where respiration cannot be fullfilled
  generation               ; counts generation in the season
  grown                    ; not used anymore?
  regenerationtime         ; measures the time from first chain release to last chain release
  cw                       ; carbon weight
  reprobuffer              ; carbon storage for reproduction
]

blastozoids-own [
  number                   ; always 1
  l                        ; body length in cm
  age                      ; age in days
  starvation-days          ; counts days where respiration cannot be fullfilled
  generation               ; counts generation in the season
  sex                      ; 0 = female, 1 = male
  cw                       ; carbon weight
  reprobuffer              ; carbon storage for reproduction
]

chains-own [
  number                   ; number of buds (aggregates) in the chain
  l                        ; body length in cm
  age                      ; age in days
  starvation-days          ; counts days where respiration cannot be fullfilled
  generation               ; counts generation in the season
  sex                      ; 0 = female, 1 = male
  cw                       ; carbon weight
  reprobuffer              ; carbon storage for reproduction
]

debkrill-own [             ; follows the notation of Jager et al. 2015
  WV
  WB
  WR
  L                        ; structural body lengh in mm
  JA
  JM
  JV
  JR
  spawningevent
  age
  ageoffirstrepro
  starvation-days
]

; reproduction will be computational problematic. This is why I model that as clutches
; this needs to be changed if direct interaction between krill and salps will be modelled
clutches-own [
  WV
  WB
  WR
  L                       ; structural body length in mm
  JA
  JM
  JV
  JR
  spawningevent
  age
  ageoffirstrepro
  number
  starvation-days
]


; model setup ---------------------------------------------------------------------------------------------------
to setup

  ca; resets everything

  file-close; make sure there is no file still open

  set immigrationtime 0; this variable stores the time when salps occur in the season

  set max_gen 0; this variable stores the number of reproductive cycles of salps in one season

  set resolution 16; this defines how much m^3 are represented by one patch - this is important for the amount of food one animal can access per day

  if (not file-exists? "chla.txt")[stop]; chla.txt stores a time series of chla values from a lognormal distribution based on the amlr data from crhistian reiss and colleagues

  let nstar 0; nstar is the max chla value in the season

  ifelse (PPMode = "Lognorm") ; if this path is taken, than the chla (food) amount will vary in time and is read in from a text file
  [
    file-open "chla.txt"; text file with time series of chla values

    ifelse (file-at-end?); if file is at the end the simulation stops
    [
      stop
    ]
    [ ; if not a new value is read in - the calculation below makes sure that this value will be reached - deltachla is the decay rate of chla per day and rchla is the growth rate of chla
      set nstar (read-from-string file-read-line) / 100
      set chlaK nstar / ( 1 - deltachla / rchla)
    ]
  ]
  [
    set chlaK const_food / ( 1 - deltachla / rchla); this is the path with const max chla each year - the maximum would be than const_food each each
    set nstar const_food
  ]

  ask patches [ set chla nstar * resolution]; each patch will contain this amount of chla, this depends on the resolution, i.e. how much m^3 are represented by one patch

  set peakabundance 0; this variable should store the maximum abundance during a season
  set peakabundancelist []; this list stores the abundance peak in each season
  set immigrationtimelist []; this list stores the time when salps enter the arena
  set regenerationtimelist []; this list stores the time between the first chain release and the last chain release of oozoids

  set temp_factor 1; based on the temperature physiological processes are reduced

  set rc 0 ; generation time counter
  set rc_max 0; maximum number of reproductive cycles in a season

  set monthlycount [0 0 0 0 0 0 0 0 0 0 0 0]; here the intraannual distribution of salp abundances will be stored

  set maxn 0; maximum number of salps
  set maxlo 0 ; maximum length of oozoids
  set maxlb 0 ; maximum length of blastozoids

  set gr-o []; list where growth rates of oozoids are stored
  set gr-b []; list where growth rates of blastozoids are stored
  set ratio-list [] ; list where ratio of oozoids and blastozoids are stored


  ask patches [set pcolor scale-color  green chla 0 1 ]; set color of patches depending on chla content
  set-default-shape oozoids "default"; shape of oozoids
  set-default-shape blastozoids "default"; shape of blastozoids

  ; creat two oozoids wiht body length l = 2 cm, the carbon weight cw is calculated according to an equation provided by Henschke et al 2018, and an initial value is attributed to the reproductive storage reprobuffer
  create-oozoids 2 [set number 1 set l 2 set size 2 set color 47 + random-float 2 - 1 set oozoid-repro-events 0 set generation 0 set grown true set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4]

  set n_oozoids count oozoids; number of oozoids
  set n_blasto 0; number of blastozoids

  set migrationevent? true; this boolean is true if in that season a migration event has happened


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; N_krill krill will be initialized, with initial Length of 1.7 mm value is given by Ikeda 1984 J. Exp. Mar. Biol. Ecol., which is C1 stage with 30 days of age

  create-debkrill N_krill
  [

   set L 0.34  ; structural length in mm, to get to the real length one has to divide by 0.2
   set WV 0.22 * ( L) ^ 3 ; structural body mass
   set WB 0.0; Assimilate buffer in egg - this is from Jager et al. 2015
   set WR 0 ; Build-up of reproduction buffer
   set spawningevent 0; number of spawning events
   set age 30; age in days
   set ageoffirstrepro 0; age of first reproduction
   setxy random-xcor random-ycor; random location
   set starvation-days 0
  ]

  set debsurvival 1; krill survival is dependent on the biomass krill assuming that predators will be attracted by krill - the survival is updated in the update patch routine
  set maxspawn 0; maximum spawn events of single krill
  reset-ticks; starting the netlogo clock

end



; in the go procuedure the main procedures are called

to go

  if (ticks >= 100000 OR count debkrill <= 0) [stop]; simulation stops either after it has reached its time span or all krill have died

  if (any? (turtle-set oozoids chains blastozoids)); if any salps are around determine the highest generation time reflecting the number of life cycles
  [
    set max_gen max list max_gen (max [generation] of (turtle-set oozoids blastozoids chains))
  ]

  if ( ticks mod 365 = 180) ; in the middle of the winter determine the largest generation time of the last season - end of june
  [
    detreprocyclesalp
  ]


  grow; determine growth in body length
  asexual_repro; asexual reproduction of salps (oozoids)
  sexual-repro; sexual reproduction of krill and salps (blastozoids)
  death; mortality

  if salps? = true ; if salps should be simulated than the immigration of salps into the simulation arena is modelled
  [
    immigration
  ]
  update-patches; patches will be updated, chla content also the density dependent survival of krill is determined here
  move; random walk of krill and salps

  do_graphs; do graphs

  tick ; advance time step by one

end

; in the grow procedure the change in body length will be calculated for all animals, this is done for each patch checking whether there is enough food available or not

to grow

   set T cos ((ticks) / 365 * 360) * 2 + 273; determining a temperature based on the time of the year

   set temp_factor exp (8000 / (2 + 273) - 8000 / T); Arrhenius temperature factor - matches well with the reported q10 of 2.8 -

   let tempdep exp ( 8000 / 273 - 8000 / ( T)) ; Arrhenius relation for krill, that can be potentially different than for salps

   let dayt ticks mod 365; day of the year - used for krill below

   ask patches ; the following procedures are conducted per patch - assuming this is the local area where individuals compete for ressources

   [

    let need 0

    let chla_control chla

    let chla_rescaled chla

    ; this is assuming that all animals have the same Holling Type II functional response, meaning having the same half saturation curve, which is certainly not true
    ; in essence it means more food is better for krill and salps which is debated
    ; below the need in a patch is calculated and if need is larger than the available food the uptake is restricted, i.e. it cannot be more food consumed than it is present

    let fr (chla / resolution) / ( (chla / resolution) + halfsat); this is the Holling type II functional response that animals would have given enough ressources - this depends on the density of chla per m^3 and not the total amount

    let fr-old fr

    ask oozoids-here

    [
      set need need + temp_factor * number * (grazing_factor * fr * l ^ 2 ); for oozoids number is always one, the uptake of food depends on the squared body length
    ]

    ask blastozoids-here

    [
      set need need + temp_factor * number * (grazing_factor * fr * l ^ 2 ); for blastozoids number is always one, the uptake of food depends on the squared body length
    ]

    ask chains-here

    [
      set need need + temp_factor * number * (grazing_factor * fr * l ^ 2 ); for chains number can vary, since number represents the number of blastozoids in the chain
    ]

    ; potential food uptake for krill

    let JAtemp 0

    ask debkrill-here
    [

      ifelse (dayt > 90 and dayt < 270 and l > 30 * 0.2); this is unique to krill, if individuals are adults here l > 30 mm, than they can reduce metabolic activity during winter
      [

        let JaAm HibernationFactor * 0.044; the hibernationfactor changes the food uptake

        set JAtemp fr * JaAM * L ^ 2 * tempdep; this would be the potential food uptake

        set need need +  JAtemp / 48 ; converting JA into food requires to divide by Yax = 0.8 and to convert mg C into mg chla c:chla = 60 and therefore 48 = 0.8 * 60
      ]
      [
        let JaAm 0.044

        set JAtemp fr * JaAM * L ^ 2 * tempdep; seems o.k.


        set need need + JAtemp / 48 ;converting JA into food requires to divide by Yax = 0.8 and to convert mg C into mg chla c:chla = 60 and therefore 48 = 0.8 * 60
      ]
    ]

    ask clutches-here
    [

      if (age > 30)
      [
        let JaAm 0.044

        ; here uptake is multiplied by the number of small krill in the cohort

        set JAtemp fr * JaAM * L ^ 2 * tempdep * number; seems o.k.


        set need need + JAtemp / 48 ;
      ]

    ]

    if chla < need ; if more food is requested than available, than the functional response has to be updated
    [
      set fr fr * chla / need;
    ]

    ; growth of oozoids

    ask oozoids-here
    [

      let carbon_assi temp_factor * fr * grazing_factor * l ^ 2 * 0.64 * 50; this is assumed to be the assimilated carbon assuming a C:Chla ratio of 50
      let carbon_growth 0.85 * carbon_assi; 85% of the carbon is allocated to growth
      let carbon_repro 0.15 * carbon_assi; 15% of the carbon is allocated to reproduction
      let oldl l
      let starvationlocal false

      set chla chla - temp_factor * number * (grazing_factor * fr * l ^ 2); reduction in the overall available chla in the patch

      let growth_resp temp_factor * oozoid_resp * cw; respiration costs

      ifelse (growth_resp > carbon_growth); if respiration is higher than carbon allocated to growth

      [
        ifelse growth_resp < (carbon_growth + carbon_repro); if respiration can be covered by assimilated carbon
        [
          set carbon_repro (carbon_repro + carbon_growth - growth_resp)
          set carbon_growth 0
        ]
        [; if respiration can be covered by assimilated carbon and the reproduction storage
          ifelse growth_resp < (carbon_growth + carbon_repro + reprobuffer)
          [
            set carbon_growth 0
            set carbon_repro 0
            set reprobuffer (reprobuffer + carbon_growth + carbon_repro - growth_resp)

          ]
          [; respiration loss cannot be covered at all - animal starves
            set carbon_growth 0
            set carbon_repro 0
            set reprobuffer 0
            set starvationlocal true
          ]
      ]
      ]
      ; if respiration is less than carbon allocation to growth
      [set carbon_growth carbon_growth - growth_resp if (carbon_growth < 0) [user-message "carbon-growth"]]

      let oldcw cw; this is for testing

      set cw cw + carbon_growth; new growth

      if (cw - oldcw < 0 ) [user-message "cw shrinkage"]; salps are not allowed to shrink

      set reprobuffer (1 - oozoid_resp) * reprobuffer + carbon_repro ; reprobuffer is updated - here are respiration costs for the reproductive tissue calculated - this is not done in the DEB

      set l (17 * cw ^ 0.4) / 10; conversion from carbon weight to body length - from an empirical relationship cited in Henschke et al. 2018

      if (l - oldl < -0.00001) [user-message word "l smaller than old l" (l - oldl) user-message word "l " l user-message word "oldl " oldl user-message word "old cw" oldcw user-message word "cw " cw]

      ifelse (starvationlocal = true) [ set starvation-days starvation-days + 1 set grown false]; increasing starvationday counter if animal could not fullfill respiration requirements

      [set grown true]; I do not use that anymore

      if MeasureInc? [if (ticks < 5000 ) [ set gr-o lput (l - oldl) gr-o]]; measuring the increment in growth - this cannot be done through the complete simulation, because that would need to much memory space

      set size 2; this is only how the agent is displayed.


    ]

; same idea as above
  ask blastozoids-here
[


      let carbon_assi temp_factor * fr * grazing_factor * l ^ 2 * 0.64 * 50
      let carbon_growth 0.8 * carbon_assi
      let carbon_repro 0.2 * carbon_assi
      let oldl l
      let starvationlocal false

      set chla chla - temp_factor * number * (grazing_factor * fr * l ^ 2)

      let growth_resp temp_factor * blasto_resp * cw

      ifelse (growth_resp > carbon_growth); if respiration is higher than carbon allocated to growth

      [
        ifelse growth_resp < (carbon_growth + carbon_repro); if respiration is higher than the carbon allocated to growth and the reprobuffer the animal starves
        [
          set carbon_repro (carbon_repro + carbon_growth - growth_resp)
          set carbon_growth 0
        ]
        [; otherwise the animal does not grow and pays from the reprobuffer
          ifelse growth_resp < (carbon_growth + carbon_repro + reprobuffer)
          [
            set carbon_growth 0
            set carbon_repro 0
            set reprobuffer (reprobuffer + carbon_growth + carbon_repro - growth_resp)

          ]
          [
            set carbon_growth 0
            set carbon_repro 0
            set reprobuffer 0
            set starvationlocal true
          ]
        ]
      ]
      ; if respiration is less than carbon allocation to growth
      [set carbon_growth carbon_growth - growth_resp if (carbon_growth < 0) [user-message "carbon-grwoth"]]

      set cw cw + carbon_growth


      set reprobuffer (1 - blasto_resp) * reprobuffer + carbon_repro
      set l (17 * cw ^ 0.4) / 10

      ifelse (starvationlocal = true) [ set starvation-days starvation-days + 1 ]

      [];[set starvation-days starvation-days - 1 set grown true]

      if MeasureInc? [if (ticks < 5000 ) [ set gr-b lput (l - oldl) gr-b]]

      set size 2


]

; same idea as above
  ask chains-here
[


      let carbon_assi temp_factor * fr * grazing_factor * l ^ 2 * 0.64 * 50
      let carbon_growth 0.85 * carbon_assi
      let carbon_repro 0.15 * carbon_assi
      let oldl l
      let starvationlocal false

      set chla chla - temp_factor * number * (grazing_factor * fr * l ^ 2)

      let growth_resp temp_factor * blasto_resp * cw

      ifelse (growth_resp > carbon_growth); if respiration is higher than carbon allocated to growth

      [
        ifelse growth_resp < (carbon_growth + carbon_repro); if respiration is higher than the carbon allocated to growth and the reprobuffer the animal starves
        [


          set carbon_repro (carbon_repro + carbon_growth - growth_resp)
          set carbon_growth 0


        ]
        [; otherwise the animal does not grow and pays from the reprobuffer
          ifelse growth_resp < (carbon_growth + carbon_repro + reprobuffer)
          [
            set carbon_growth 0
            set carbon_repro 0
            set reprobuffer (reprobuffer + carbon_growth + carbon_repro - growth_resp)

          ]
          [
            set carbon_growth 0
            set carbon_repro 0
            set reprobuffer 0
            set starvationlocal true
          ]
        ]
      ]
      ; if respiration is less than carbon allocation to growth
      [set carbon_growth carbon_growth - growth_resp if (carbon_growth < 0) [user-message "carbon-grwoth"]]

      set cw cw + carbon_growth


      set reprobuffer (1 - blasto_resp) * reprobuffer + carbon_repro
      set l (17 * cw ^ 0.4) / 10

      ifelse (starvationlocal = true) [ set starvation-days starvation-days + 1 ]

      [];[set starvation-days starvation-days - 1 set grown true]

      if MeasureInc? [if (ticks < 5000 ) [ set gr-b lput (l - oldl) gr-b]]
      set size 2

    ]

    ask debkrill-here
    [
      ifelse (l > 11 * 0.2)
      [
        determine_fluxes fr; fluxes deterimend after Jager et al. 2015
      ]
      [
        determine_fluxes_larvae fr
      ]

      set WV WV + JV; new WV is WV + JV - shrinkage happens already in the deterime_fluxes procedure



      let dV 0.22

      set L (WV / dV) ^ (1 / 3) ;

      set chla chla - (JA / 48)

    ]


     ask clutches-here
    [
      if age > 30
      [
       ifelse (l > 11 * 0.2)
        [
          determine_fluxes fr; fluxes deterimend after Jager et al. 2015
        ]
        [
          determine_fluxes_larvae fr
        ]

        set WV WV + JV; new WV is WV + JV - shrinkage happens already in the deterime_fluxes procedure
        set WR WR + JR


        let dV 0.22

        set L (WV / dV) ^ (1 / 3) ;

        set chla chla - number * (JA / 48)

      ]
    ]

    if chla < -0.00001 [user-message "Negative chla values"] ; this is due to numerical issues that zero may be not exactly zero
    if (chla <= 0.005) [  set chla 0.005] ; It is important to set a minimum amount of chla at the end. This is due to diffusion from ohter patches or deeper layser. chla > 0 is important to model the ressource.

  ]



end

; chain release is the asexual part of reproduction - here it is length dependent; if salps reach a given length the release a chain of with a given number of bud (aggregates), but not more
; than it has been observed

to asexual_repro

  ask oozoids [
  if (l > 6.5 and grown = true and oozoid-repro-events = 0)
  [
      let budnumber floor (reprobuffer / 0.0329)
      set budnumber min list budnumber 150



      hatch-chains 1 [

      set generation generation + 1 set number budnumber set l 0.5 set size 2 set age 0 set sex 0 set starvation-days 0 set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4] ;; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
      set oozoid-repro-events oozoid-repro-events + 1
      set regenerationtime age
      set reprobuffer reprobuffer - budnumber * 0.0329


   ]

  if (l > 8 and grown = true and oozoid-repro-events = 1)
  [
      let budnumber floor (reprobuffer / 0.0329)
      set budnumber min list budnumber 180

      hatch-chains 1 [set generation generation + 1 set number 180  set l 0.5 set size 2 set age 0 set sex 0 set starvation-days 0 set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4] ;; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
      set oozoid-repro-events oozoid-repro-events + 1
      set reprobuffer reprobuffer - budnumber * 0.0329

    ]

    if (l > 9.5 and grown = true and oozoid-repro-events = 2)
    [
      let budnumber floor (reprobuffer / 0.0329)
      set budnumber min list budnumber 210



      hatch-chains 1 [ set generation generation + 1  set number budnumber set l 0.5 set size 2 set age 0 set sex 0 set starvation-days 0 set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4] ;; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
      set oozoid-repro-events oozoid-repro-events + 1
      set reprobuffer reprobuffer - budnumber * 0.0329
    ]

    if (l > 10.5 and grown = true and oozoid-repro-events = 3)
    [


      let budnumber floor (reprobuffer / 0.0329)
      set budnumber min list budnumber 240

      hatch-chains 1 [ set generation generation + 1 set number budnumber set l 0.5 set size 2 set age 0 set sex 0 set starvation-days 0 set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4] ;; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
      set oozoid-repro-events oozoid-repro-events + 1
      set regenerationtimelist lput (age - regenerationtime) regenerationtimelist
      die
    ]

  ]

end


; sexual part of the salp reproduction. At a given size and sufficient carbon allocated in the reprobuffer a embryo is released with 0.7 survival. After that the blastozoids change sex and the chain dies,
; since the male blastozoids do not form chains
to sexual-repro

  ask chains [
    if (l > 2.5 and reprobuffer > 0.0269) [
      repeat number [
        ifelse (random-float 1 < 0.7)
        [  hatch-oozoids 1 [ set number 1 set generation generation + 1  set l 0.4 set size 2  set age 0 set color (color + (random-float 1 - 0.5) ) set starvation-days 0 set oozoid-repro-events 0 set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4]  hatch-blastozoids 1 [set sex 1 set number 1 ] ]
        [ hatch-blastozoids 1 [set sex 1 set number 1] ]
      ]
      die; here the chain dies and only individuals survive
    ] ;; embryio release is following Pakhomov and Hunt 2017, Deep- Sea Research II
  ]

  ; for krill we follow the ideas of Jager et al. 2015 - clutch size is abitrary and also the buffer that remains with the krill 35 is open for further manipulation
  ask debkrill
  [

   set WR WR + JR

    if (WR > 56 + 35)       ; the number 56 is following Jager et al. 2015 value for WB0 (dry weight of an egg) and a clutch size of 2000, i.e. 2000 * 0.038 mgdwt = 56 - the number 35 is an ad hoc buffer to buffer starvation
    [
      if spawningevent = 0 [set ageoffirstrepro age]

      set spawningevent spawningevent + 1

      set WR WR - 56

      hatch-clutches 1
      [

        set l (1.72 * 0.2) ; this is already the length of C1 stage, but we assume that the first 30 days krill will not feed, the growth will be fed by the energy from the eggs
        set WV 0.22 * ( L) ^ 3 ; structural body mass
        set WB 0.0; Assimilate buffer in egg - this is from Jager et al. 2015
        set WR 0 ; Build-up of reproduction buffer
        set JA 0
        set JV 0
        set JR 0
        set JM 0
        set spawningevent 0; number of spawning events
        set age 0; age in days
        set ageoffirstrepro 0; age of first reproduction
        setxy random-xcor random-ycor; random location
        set number 2000 set age 0
        set starvation-days 0


      ]

    ]


  ]

end

; salps die after 500 days or due to daily mortality or due to starvation (30 days)

to death

ask oozoids [ set age age + 1 if (age > 500 or oozoid-repro-events = 5) [die]];
ask blastozoids [ set age age + 1 if (age > 500) [die]]
ask chains [ set age age + 1 if (age > 500) [die]]
ask (turtle-set oozoids blastozoids chains) [if (starvation-days > starvation) [ die]]

ask oozoids [ if random-float 1 < DailyMort [die]]

ask blastozoids [ if random-float 1 < DailyMort [die]]

ask chains [
    let counter 0
    repeat number [
      if random-float 1 < DailyMort [
        set counter counter + 1 ]
      ]
    set number number - counter
    if number <= 0 [die]
  ]

  ask debkrill [ set age age + 1
    if (age > (8 * 365)) or (random-float 1 < (1 - debsurvival) )[die]
  ]

  ; To reduce computational effort in the beginning clutches are computed as cohorts. If the number drops below 10 they are treated as individuals

  ask clutches
  [
    set age age + 1

    set number number * 0.95

    if number < 10
    [

      hatch-debkrill (round number)
      [
        set WV 0.22 * L ^ 3 ; structural body mass
        ;set WB 0.028; Assimilate buffer in egg;
        ;set WR 0 ; Build-up of reproduction buffer
        set spawningevent 0
        ;set age 0
        set ageoffirstrepro 0
        ;setxy random-xcor random-ycor

      ]
      die
     ]

  ]
end

; only random walk for salps and krill, but certainly krill may move further and even avoid salps?
to move

  ask turtles [ move-to one-of neighbors]

end

to update-patches

  if (PPMode = "Lognorm")
  [
    if (ticks mod 365 = 180 - vegetation_delay)
    [
      let nstar (read-from-string file-read-line) / 100
      set chlaK nstar / ( 1 - deltachla / rchla)

    ]
  ]

  let algae_growth rchla * (0.5 * cos ((ticks + vegetation_delay) / 365 * 360) + 0.5)

  ask patches
  [

    set chla max list (chla + resolution * (algae_growth * (chla / resolution) * ( 1 - (chla / resolution) / chlaK) - deltachla * (chla / resolution))) (0.005 * resolution)
    set pcolor scale-color  green chla 10 0
    if (chla < 0) [user-message "Warning! Chla < 0 in update patches" ]
  ]

  ;; the survival for krill should be density dependent, assuming one large adult krill has a body dry weight of 220 the following decline in survival would result in a survival of 0 at 100,000 krill

  let WVtemp sum [WV] of debkrill

  let tempsurvival ( 1 - WVtemp / 22000000)
  if tempsurvival < 0 [set tempsurvival 0]


   set debsurvival tempsurvival

  ;print debsurvival

end

; there is no wintersurvival in the model for salps. Thus they migrate into the simulatin arena at some point in the season. This is modelled here

to immigration


if (count (turtle-set oozoids chains blastozoids) < 1) [

    if( random-float 1 < immiprob)
    [

      create-oozoids ni [set number 1 set l sizeofmigra set size 2  (set color 47 + random-float 2 - 1) setxy random-xcor random-ycor set oozoid-repro-events 0 set generation 0 set starvation-days 0 set cw (l * 10 / 17) ^(1 / 0.4) set reprobuffer cw / 4]

      set immigrationtime ticks
      set migrationevent? true

    ]
  ]

end


to plot-histo

set-current-plot "Histogram Oozoids lengths (cm)"

set-plot-pen-mode 1

set-plot-x-range 0 15

set-histogram-num-bars 50

histogram [l] of oozoids


set-current-plot "Histogram Blastozoids lengths (cm)"

set-plot-pen-mode 1

set-plot-x-range 0 5

set-histogram-num-bars 50

histogram [l] of (turtle-set blastozoids chains)



set-current-plot "Age of Oozoids"

set-plot-pen-mode 1

set-plot-x-range 0 500

set-histogram-num-bars 50

histogram [age] of oozoids



set-current-plot "Age of Blastozoids"

set-plot-pen-mode 1

set-plot-x-range 0 500

set-histogram-num-bars 50

histogram [age] of (turtle-set blastozoids chains)

if (ticks mod 365 = 50) ; to reduce computational effort this is only computed once a year
  [

    if MeasureInc?
    [

      set-current-plot "Growth Rates Oozoids"

      clear-plot

      set-plot-pen-mode 1

      set-plot-x-range 0 (max (gr-o) + 0.5)

      set-histogram-num-bars 50

      histogram gr-o



      set-current-plot "Growth Rates Blastozoids"

      set-plot-pen-mode 1

      let upper-value 1

      if (length gr-b > 0)
      [
        set upper-value max(gr-b)
      ]

      set-plot-x-range 0 (upper-value + 0.5)

      set-histogram-num-bars 50

      histogram gr-b

    ]
    set-current-plot "Annual distribution"

    clear-plot

    set-plot-pen-mode 1

    set-histogram-num-bars 12

    set-plot-x-range 0 13

    let summ sum monthlycount

    let monthlycount_temp monthlycount

    (foreach [0 1 2 3 4 5 6 7 8 9 10 11 ] monthlycount_temp[[a b] -> set monthlycount_temp (replace-item a monthlycount_temp ( ((item a monthlycount_temp / summ)  ))  )])
    (foreach [0 1 2 3 4 5 6 7 8 9 10 11 ] monthlycount_temp  [[a b] -> set monthlycount_temp (replace-item a monthlycount_temp ( precision (item a monthlycount_temp) 3   ))] )

    set-plot-y-range 0 (max (monthlycount_temp) + 0.1  )

    (foreach [1 2 3 4 5 6 7 8 9 10 11 12] monthlycount_temp[[a b] -> plotxy a b])

  ]

  set-current-plot "Chain releases"

  set-plot-pen-mode 1

  set-plot-x-range 0 6

  set-histogram-num-bars 6

  histogram [oozoid-repro-events] of oozoids


  set-current-plot "Bud/Oozoid"

  set-plot-pen-mode 0


  let temp n_blasto
  let temp2 n_oozoids

  if (temp2 > 0)
  [ plotxy ticks temp / temp2
  set ratio-list lput (temp / temp2) ratio-list
  ]


  set-current-plot "Krill Growth Curves"

  foreach [who] of debkrill [ x -> plotxy (ticks ) ([l / 0.2] of sdebkrill x)]


end


to determine_fluxes [frac]; fluxes determined following Jager et al. 2015
  ; important in winter adult krill is descreasing its metabolic activity! This is expressed by the case differentiation below.



let dayt ticks mod 365
let tempdep exp ( 8000 / 273 - 8000 / ( T)) ; Arrhenius relation from appendix of jager paper

ifelse (dayt > 90 and dayt < 270 and l > 30 * 0.2); adult krill may have a reduced metabolic activity during winter
  [
    ; assimilation
    let JaAm HibernationFactor * 0.044
    set JA frac * JaAM * L ^ 2 * tempdep; seems o.k.

    ; respiration
    let JVM HibernationFactor * 0.0032
    set JM JVM * L ^ 3 * tempdep ; seems o.k.

    ; growth
    let YVA 0.8
    let kappa 0.8
    set JV (YVA * (kappa * JA - JM))
    ; reproduction
    set JR (1 - kappa) * JA


    ;; if JM > JA

    if (JV < 0)
    [

      let delta JV

      ifelse (JR > abs (delta))
      [

        set JV 0
        set JR JR + delta



      ]
      [


        set JV JV + JR
        set JR 0


        ifelse ( WR > abs (JV))
        [


          set WR WR + JV
          set JV 0

        ]
        [; respiration cannot be covered by assimilated carbon and reproduction buffer, than krill shrink (WV is getting smaller since JV is negative



          set JV JV + WR
          set WR 0


          set WV WV + JV
          set JV 0
          set starvation-days starvation-days + 1

        ]

      ]

    ]
  ]

  [

    let JaAm 0.044
    set JA frac * JaAM * L ^ 2 * tempdep; seems o.k.


    let JVM 0.0032
    set JM JVM * L ^ 3 * tempdep ; seems o.k.

    let YVA 0.8
    let kappa 0.8
    set JV (YVA * (kappa * JA - JM))

    set JR (1 - kappa) * JA



    ;; if JM > JA

    if (JV < 0)
    [

      let delta JV
      ifelse (JR > abs (delta))
      [
        set JV 0
        set JR JR + delta
      ]
      [
        set JV JV + JR
        set JR 0

        ifelse ( WR > abs (JV))
        [
          set WR WR + JV
          set JV 0
        ]
        [
          set JV JV + WR
          set WR 0
          set WV WV + JV
          set JV 0

          set starvation-days starvation-days + 1

        ]

      ]

    ]
  ]
end

to determine_fluxes_larvae [frac]; frac is the functional response (Holling II)  fluxes determined following Jager et al. 2015
  ; but how are the fluxes different for small krill (larvae and juvenile) from the one of adult krill
  ; here following Jager et al. 2015 - there is only the growth part stored. The flux into reproduction JR is burned



  let dayt ticks mod 365
  let tempdep exp ( 8000 / 273 - 8000 / ( T)) ; Arrhenius relation from appendix of jager paper



  let JaAm 0.044
  set JA frac * JaAM * L ^ 2 * tempdep; seems o.k.


  let JVM 0.0032
  set JM JVM * L ^ 3 * tempdep ; seems o.k.

  let YVA 0.8
  let kappa 0.8
  set JV (YVA * (kappa * JA - JM))
  set JR 0

  ;; if JM > JA

  if (JV < 0)
  [
    set WV WV + JV
    set JV 0
    set starvation-days starvation-days + 1
  ]

end

to detreprocyclesalp; how many reproductive cycles have happened in the seaon

    set rc max_gen
    if (rc > rc_max) [set rc_max rc]
    set max_gen 0
    ask (turtle-set blastozoids oozoids chains) [set generation 0]

    set-current-plot "Seasonal repoduction cycles"
    plot rc / 2

end

to do_graphs

  let index3 (floor ((ticks mod 365) / 30.5))
  set monthlycount (replace-item index3 monthlycount (item index3 monthlycount + sum [number] of (turtle-set oozoids blastozoids chains))); storing the number of salps for different months for the intraannual distribution later on
  set n_oozoids count oozoids
  set n_blasto (count blastozoids + sum [number] of chains); counting all blastozoids, males and females


  plot-histo

  set maxn (max (list maxn (n_oozoids + n_blasto))); maximum salp abundance during the simulation run

  if n_oozoids > 0; storing the largest body length of oozoids
  [
    set maxlo (max list maxlo max [l] of oozoids)
  ]

  if n_blasto > 0; storing the largest body length of blastozoids, males and females (chains)
  [
    set maxlb (max list maxlb max [l] of (turtle-set blastozoids chains))
  ]

  if (n_oozoids + n_blasto > peakabundance) ; storing the highes abundance during the season
  [
    set peakabundance (n_oozoids + n_blasto)
    set peakimmigrationtime immigrationtime
  ]

  if (ticks mod 365 = 210); end of july, at the end of the Southern winter when abundances should be really low
  [

    set peakabundancelist lput peakabundance peakabundancelist; in peakabundancelist the highest salp abundance is stored that occured in the season (here the season ends after 210 days, i.e. end of July)

    ifelse migrationevent? = true; only if in that season migration has happened it make sense to store the time
    [
      set immigrationtimelist lput peakimmigrationtime immigrationtimelist
      set migrationevent? false
    ]
    [; -1 indicates that no migration event has happened in that particular season
      set immigrationtimelist lput -1 immigrationtimelist
    ]

    set-current-plot "Salp peaks"
    plot peakabundance
    set peakabundance 0; resetting peakabundance
  ]

  if (count debkrill > 0) ; maxspawn is the highest number of spawning events that one krill had in its life time
  [
    set maxspawn max (list maxspawn (max [spawningevent] of debkrill))
  ]
end

to ref; idea of this function is to set all parameters back to their reference values

  set starvation 30
  set dailymort 0.02
  set grazing_factor 0.003
  set vegetation_delay 45
  set immiprob 0.0085
  set rchla 0.25
  set deltachla 0.05
  set halfsat 0.2
  set HibernationFactor 0.2
  set ni 10
  set const_food 1
  set sizeofmigra 3
  set N_krill 10
  set oozoid_resp 0.05
  set blasto_resp 0.15

end
@#$#@#$#@
GRAPHICS-WINDOW
108
14
520
427
-1
-1
4.0
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
100
0
100
1
1
1
ticks
30.0

BUTTON
12
30
75
63
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
12
87
75
120
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
1

PLOT
535
15
765
300
Abundances
NIL
NIL
0.0
2.0
0.0
1.5
true
true
"" ""
PENS
"Aggregates" 1.0 0 -16777216 true "" "plotxy (ticks) n_blasto"
"Solitaries" 1.0 0 -2674135 true "" "plotxy (ticks)  n_oozoids"

PLOT
770
15
930
300
Mean Chla (scaled from 1 to 0)
NIL
NIL
0.0
10.0
0.0
1.5
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [chla] of patches / resolution"
"pen-1" 1.0 0 -13840069 true "" "plot max [chla] of patches / resolution"

PLOT
247
460
507
748
Histogram Oozoids lengths (cm)
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
"default" 1.0 0 -16777216 true "" ""

PLOT
517
462
768
748
Histogram Blastozoids lengths (cm)
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
"default" 1.0 0 -16777216 true "" ""

PLOT
763
906
960
1044
Age of Oozoids
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
"default" 1.0 0 -16777216 true "" ""

PLOT
764
754
953
900
Age of Blastozoids
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
"default" 1.0 0 -16777216 true "" ""

PLOT
1257
535
1546
740
Growth Rates Oozoids
NIL
NIL
-0.01
0.2
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

PLOT
1258
848
1552
1052
Growth Rates Blastozoids
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
"default" 1.0 0 -16777216 true "" ""

MONITOR
1259
482
1432
527
NIL
mean gr-o
17
1
11

MONITOR
1268
767
1480
812
NIL
mean gr-b
17
1
11

SLIDER
9
485
181
518
starvation
starvation
0
1000
30.0
1
1
NIL
HORIZONTAL

PLOT
251
755
754
970
Annual distribution
NIL
NIL
0.0
10.0
0.0
10.0
false
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

SLIDER
10
530
182
563
DailyMort
DailyMort
0
1
0.02
0.001
1
NIL
HORIZONTAL

SLIDER
10
575
182
608
grazing_factor
grazing_factor
0
2
0.003
0.001
1
NIL
HORIZONTAL

SLIDER
10
620
182
653
vegetation_delay
vegetation_delay
0
180
45.0
1
1
NIL
HORIZONTAL

BUTTON
13
136
76
169
NIL
Ref
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
976
756
1242
1056
Chain releases
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
"default" 1.0 0 -16777216 true "" ""

PLOT
794
464
1179
744
Seasonal repoduction cycles
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
"default" 1.0 0 -16777216 true "" ""

PLOT
912
306
1083
458
Bud/Oozoid
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
"default" 1.0 0 -16777216 true "" ""

MONITOR
6
327
87
372
Max Density
precision (maxn /  (101 * 101 * 16)) 2
17
1
11

PLOT
545
307
705
457
Temp_Factor
NIL
NIL
0.0
365.0
0.0
3.0
false
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot temp_factor"

CHOOSER
2
184
94
229
PPMode
PPMode
"Const" "Lognorm"
0

PLOT
711
307
911
457
Salp peaks
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
"default" 1.0 0 -16777216 true "" ""

SLIDER
10
665
182
698
immiprob
immiprob
0
1
0.0085
0.001
1
NIL
HORIZONTAL

SWITCH
1
237
104
270
Chains?
Chains?
0
1
-1000

SLIDER
10
705
182
738
rchla
rchla
0
1
0.25
0.01
1
NIL
HORIZONTAL

SLIDER
10
750
182
783
deltachla
deltachla
0
1
0.05
0.01
1
NIL
HORIZONTAL

SWITCH
10
439
136
472
MeasureInc?
MeasureInc?
0
1
-1000

SLIDER
10
795
182
828
halfsat
halfsat
0
1
0.2
0.01
1
NIL
HORIZONTAL

SLIDER
10
885
182
918
ni
ni
0
100
10.0
1
1
NIL
HORIZONTAL

SLIDER
10
970
182
1003
sizeofmigra
sizeofmigra
0
5
3.0
0.1
1
NIL
HORIZONTAL

PLOT
935
15
1134
296
Krill Growth Curves
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
"default" 1.0 2 -16777216 true "" ""

SLIDER
10
930
182
963
const_food
const_food
0
5
1.0
0.001
1
NIL
HORIZONTAL

SLIDER
10
1010
182
1043
N_krill
N_krill
0
2000
1.0
1
1
NIL
HORIZONTAL

SWITCH
0
282
103
315
salps?
salps?
0
1
-1000

SLIDER
10
840
182
873
HibernationFactor
HibernationFactor
0
1
0.2
0.001
1
NIL
HORIZONTAL

PLOT
1140
15
1310
295
Krill abundance
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
"default" 1.0 0 -16777216 true "" "plot count debkrill + sum [number] of clutches"

PLOT
1087
308
1287
458
Debsurvival
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot debsurvival"

PLOT
1290
308
1490
458
Maximum krill life time spawn events 
NIL
NIL
0.0
10.0
0.0
12.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot maxspawn"

SLIDER
245
985
417
1018
oozoid_resp
oozoid_resp
0
1
0.05
0.001
1
NIL
HORIZONTAL

SLIDER
245
1030
417
1063
blasto_resp
blasto_resp
0
1
0.15
0.001
1
NIL
HORIZONTAL

PLOT
1315
15
1490
295
Abundance adult krill
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
"default" 1.0 0 -16777216 true "" "plot count debkrill"

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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

caterpillar
true
0
Polygon -7500403 true true 165 210 165 225 135 255 105 270 90 270 75 255 75 240 90 210 120 195 135 165 165 135 165 105 150 75 150 60 135 60 120 45 120 30 135 15 150 15 180 30 180 45 195 45 210 60 225 105 225 135 210 150 210 165 195 195 180 210
Line -16777216 false 135 255 90 210
Line -16777216 false 165 225 120 195
Line -16777216 false 135 165 180 210
Line -16777216 false 150 150 201 186
Line -16777216 false 165 135 210 150
Line -16777216 false 165 120 225 120
Line -16777216 false 165 106 221 90
Line -16777216 false 157 91 210 60
Line -16777216 false 150 60 180 45
Line -16777216 false 120 30 96 26
Line -16777216 false 124 0 135 15

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

monster
false
0
Polygon -7500403 true true 75 150 90 195 210 195 225 150 255 120 255 45 180 0 120 0 45 45 45 120
Circle -16777216 true false 165 60 60
Circle -16777216 true false 75 60 60
Polygon -7500403 true true 225 150 285 195 285 285 255 300 255 210 180 165
Polygon -7500403 true true 75 150 15 195 15 285 45 300 45 210 120 165
Polygon -7500403 true true 210 210 225 285 195 285 165 165
Polygon -7500403 true true 90 210 75 285 105 285 135 165
Rectangle -7500403 true true 135 165 165 270

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
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2910"/>
    <metric>mean [l] of debkrill</metric>
    <metric>min [l] of debkrill</metric>
    <metric>max [l] of debkrill</metric>
    <metric>mean [spawningevent] of debkrill</metric>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ni">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rchla">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sizeofmigra">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_krill">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="const_food">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chains?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grazing_factor">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PPMode">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-o">
      <value value="0.027"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DailyMort">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-b">
      <value value="0.0155"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immiprob">
      <value value="0.0085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltachla">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment2" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2910"/>
    <metric>mean [l] of debkrill</metric>
    <metric>mean [spawningevent] of debkrill</metric>
    <metric>mean peakabundancelist</metric>
    <metric>max  peakabundancelist</metric>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ni">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rchla">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sizeofmigra">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_krill">
      <value value="1"/>
      <value value="50"/>
      <value value="100"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="const_food">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chains?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grazing_factor">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salps?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PPMode">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-o">
      <value value="0.027"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DailyMort">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-b">
      <value value="0.0155"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immiprob">
      <value value="0.0085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltachla">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentfood05" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2910"/>
    <metric>mean [l] of debkrill</metric>
    <metric>mean [spawningevent] of debkrill</metric>
    <metric>mean peakabundancelist</metric>
    <metric>max  peakabundancelist</metric>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ni">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rchla">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sizeofmigra">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_krill">
      <value value="1"/>
      <value value="50"/>
      <value value="100"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="const_food">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chains?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grazing_factor">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salps?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PPMode">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-o">
      <value value="0.027"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DailyMort">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-b">
      <value value="0.0155"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immiprob">
      <value value="0.0085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltachla">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentfood15" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2910"/>
    <metric>mean [l] of debkrill</metric>
    <metric>mean [spawningevent] of debkrill</metric>
    <metric>mean peakabundancelist</metric>
    <metric>max  peakabundancelist</metric>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ni">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rchla">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sizeofmigra">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_krill">
      <value value="1"/>
      <value value="50"/>
      <value value="100"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="const_food">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chains?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grazing_factor">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salps?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PPMode">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-o">
      <value value="0.027"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DailyMort">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-b">
      <value value="0.0155"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immiprob">
      <value value="0.0085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltachla">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentvegdelay90" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2910"/>
    <metric>mean [l] of debkrill</metric>
    <metric>mean [spawningevent] of debkrill</metric>
    <metric>mean peakabundancelist</metric>
    <metric>max  peakabundancelist</metric>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ni">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rchla">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sizeofmigra">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_krill">
      <value value="1"/>
      <value value="50"/>
      <value value="100"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="const_food">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chains?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grazing_factor">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salps?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PPMode">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-o">
      <value value="0.027"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DailyMort">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-b">
      <value value="0.0155"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immiprob">
      <value value="0.0085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltachla">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentfullfactorial" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2910"/>
    <metric>mean [l] of debkrill</metric>
    <metric>mean [spawningevent] of debkrill</metric>
    <metric>mean peakabundancelist</metric>
    <metric>max  peakabundancelist</metric>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="30"/>
      <value value="45"/>
      <value value="60"/>
      <value value="75"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ni">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rchla">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sizeofmigra">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N_krill">
      <value value="1"/>
      <value value="50"/>
      <value value="100"/>
      <value value="200"/>
      <value value="300"/>
      <value value="400"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="const_food">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
      <value value="1.25"/>
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chains?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="grazing_factor">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salps?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PPMode">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-o">
      <value value="0.027"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DailyMort">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rb-b">
      <value value="0.0155"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immiprob">
      <value value="0.0085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltachla">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hibernationfactor">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
1
@#$#@#$#@
