;; This model is on the interaction of salps and krill. Salps can occurr in three life stages Oozoids, blastozoids and chains and Krill can occurr as adult krill or larval krill
;; larval inclduing juvenile krill < 35 mm will stay active all year round while adult krill will hibernate. The larval stage is not actively represented and only relevant in the growing process
;; To model larvae cohorts are modelled called clutches to avoid computational problems


breed [oozoids oozoid ]        ; one life from of salps that reproduces asexually by releasing chains
breed [blastozoids blastozoid] ; one life form of salps - here the male blastozoids are called blastozoids in the model
breed [chains chain]           ; chains are chains of femlae blastozoids
breed [krills krill]           ; krill modelled individually
breed [clutches clutch]        ; eggs and early stages modelled as cohorts to reduce computational effort

patches-own [chla]

globals [
  gr_o                      ; list where growth rates of oozoids are stored
  gr_b                      ; list where growth rates of blastozoids are stored
  n_s_max_total             ; maximum number of salps
  monthlycount              ; list that stores the intra-annual distribution of salp abundances
  reg_cycles                ; number of regeneration cycles each season
  reg_cycles_max            ; maximum number of regeneration cycles in a season
  ratio_list                ; list where ratio of oozoids and blastozoids are stored
  l_o_max                   ; maximum length of oozoids
  l_b_max                   ; maximum length of blastozoids
  T_factor_s                ; based on the temperature physiological processes are reduced
  chla_max                  ; maximum chla density in the season
  n_s_max                   ; maximum abundance of salps in the season
  n_s_max_list              ; list of peak abundances
  n_s_med                   ; median of n_s_max (Salp peakabundances)
  rep_cycles_max            ; max number of reproductive cycles of salps in one season
  n_oozoids                 ; number of oozoids in the world
  n_blasto                  ; number of blastozoids in the world
  resolution                ; how many cubic meters are represented by a netlogo patch
  t_immigration             ; time during the year when a small cohort of oozoids entered the simulation arena
  t_immigration_list        ; list of the immigration times during the simulation
  t_regeneration_list       ; this measures the time between first chain release and last
  t_immigration_max         ; this variable store the maximum immigration time during a season
  migrationevent?           ; true if a migration event happened during this season
  T                         ; daily temperature (Kelvin)
  spawn_max                 ; maximum spawn events of single krill
  l_k_max                   ; maximum size (l) of single krill
  n_k_max                   ; maximum abundance of krill
  n_k_max_eggs              ; maximum amaount of eggs produced by one individual krill
  d_k_firstrepro            ; day of first krill reproduction
]

oozoids-own [
  number                   ; always 1
  l                        ; body length in cm
  d_age                    ; age in days
  n_repro                  ; counts number of chain releases
  d_starvation             ; counts days where respiration cannot be fullfilled
  generation               ; counts generation in the season
  t_regeneration           ; measures the time from first chain release to last chain release
  c_weight                 ; carbon weight
  c_reproduction           ; carbon storage for reproduction
]

blastozoids-own [
  number                   ; always 1
  l                        ; body length in cm
  d_age                    ; age in days
  d_starvation             ; counts days where respiration cannot be fullfilled
  generation               ; counts generation in the season
  sex                      ; 0 = female, 1 = male
  c_weight                 ; carbon weight
  c_reproduction           ; carbon storage for reproduction
]

chains-own [
  number                   ; number of buds (aggregates) in the chain
  l                        ; body length in cm
  d_age                    ; age in days
  d_starvation             ; counts days where respiration cannot be fullfilled
  generation               ; counts generation in the season
  sex                      ; 0 = female, 1 = male
  c_weight                 ; carbon weight
  c_reproduction           ; carbon storage for reproduction
]

krills-own [                ; follows the notation of Jager et al. 2015
  W_V                       ; structural body mass (mg dry weight)
  W_B                       ; assimilate buffer in egg (mg dry weight)
  W_R                       ; build-up of reproduction buffer (mg dry weight)
  L                         ; structural body length (mm), divide by 0.2 to get real length
  J_A                       ; assimilation (mg dry weight per day)
  J_M                       ; somativ maintenance (mg dry weight per day)
  J_V                       ; structural growth (mg dry weight per day)
  J_R                       ; investment reproduction buffer (mg dry weight per day)
  n_spawn                   ; number of spawning events
  n_eggs                    ; number of eggs produced
  d_age                     ; age (days)
  d_first_reproduction      ; age of first reproduction (days)
  d_starvation              ; counts days where respiration cannot be fullfilled
]

clutches-own [
  ; reproduction will be computational problematic. This is why I model that as clutches
  ; this needs to be changed if direct interaction between krill and salps will be modelled
  W_V                       ; structural body mass (mg dry weight)
  W_B                       ; assimilate buffer in egg (mg dry weight)
  W_R                       ; build-up of reproduction buffer (mg dry weight)
  L                         ; structural body length (mm)
  J_A                       ; assimilation (mg dry weight per day)
  J_M                       ; somativ maintenance (mg dry weight per day)
  J_V                       ; structural growth (mg dry weight per day)
  J_R                       ; investment reproduction buffer (mg dry weight per day)
  n_spawn                   ; number of spawning events
  d_age                     ; age in days
  d_first_reproduction      ; age of first reproduction (days)
  number                    ; number of juveniles in the clutch
  d_starvation              ; counts days where respiration cannot be fullfilled
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; setup procedure                     ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup

  ca                                        ; resets everything

  set n_s_max_list []                       ; this list stores the abundance peak in each season
  set t_immigration_list []                 ; this list stores the time when salps enter the arena
  set t_regeneration_list []                ; this list stores the time between the first chain release and the last chain release of oozoids
  set migrationevent? false

  ; section for random parameter values during sensitivity analysis
  if SA? [
    let SArange 0.5          ; define parameter deviation for SA
    let SAmin (1 - SArange)  ; define minimum level for SA
    let SAspan (2 * SArange) ; define parameter range for SA

    ; for each parameter, select a random value between min (SAmin) and max (min + SAspan)
    set chla_growth precision (SAmin * 0.25 + random-float (SAspan * 0.25)) 3           ; round to 3 digits
    set chla_decay precision (SAmin * 0.05 + random-float (SAspan * 0.05)) 4            ; round to 4 digits
    set vegetation_delay (precision (SAmin * 45) 0 + random precision (SAspan * 45) 0)  ; round to 0 digits
    set salp_halfsat precision (SAmin * 0.20 + random-float (SAspan * 0.20)) 3          ; round to 3 digits
    set salp_immiprob (precision (SAmin * 0.85 + random-float (SAspan * 0.85)) 3)       ; round to 3 digits
    set salp_amount (SAmin * 10 + random (SAspan * 10 + 1))                             ; no rounding
    set salp_length (precision (SAmin * 3 + random-float (SAspan * 3)) 1)               ; round to 1 digit
    set salp_starvation (SAmin * 30 + random (SAspan * 30 + 1))                         ; no rounding
    set salp_mortality (precision (SAmin * 2.5 + random-float (SAspan * 2.5)) 2)        ; round to 2 digits
    set oozoid_resp (precision (SAmin * 5.0 + random-float (SAspan * 5.0)) 2)           ; round to 2 digits
    set blasto_resp (precision (SAmin * 15.0 + random-float (SAspan * 15.0)) 2)         ; round to 2 digits
    set krill_halfsat precision (SAmin * 0.106 + random-float (SAspan * 0.106)) 4       ; round to 4 digits
    set krill_amount (SAmin * 30 + random (SAspan * 30 + 1))                            ; no rounding
    set krill_mortality (precision (SAmin * 0.07 + random-float (SAspan * 0.07)) 3)     ; round to 3 digits
    set krill_hibernation (precision (SAmin * 20 + random-float (SAspan * 20)) 1)       ; round to 1 digit
  ]

  set monthlycount [0 0 0 0 0 0 0 0 0 0 0 0] ; here the intraannual distribution of salp abundances will be stored

  set gr_o []                                ; list for oozoid growth rates
  set gr_b []                                ; list for blastozoid growth rates
  set ratio_list []                          ; list for ratio of oozoids and blastozoids

  set resolution 16                          ; one patch resembles 16 m^3 of water

  ; calculation of chlorophyl a values ---------------------------------------------------------------------------------------------------------------;

  ; using a lognormal distribution based on the amlr data from christion reiss and colleagues
  ifelse (chla_supply = "Lognorm") [
    random-seed 0    ; set fixed seed from interface input
    let mu 3.83      ; meanlog value from amlr data
    let sd 0.58      ; sdlog from amlr data

    ; calculate chla_max using lognormal distribution and mu and sd from amlr data
    set chla_max (e ^ (mu + sd * (random-normal 0 1)) / 100)

    ; randomize seed again
    random-seed new-seed

    ; renormalization using the logistic equation term with loss term
    set chla_max (chla_max / (1 - chla_decay / chla_growth))

  ][
  ; const max chla each year using the logmean from the amlr data
    set chla_max (log 3.83 10 / ( 1 - chla_decay / chla_growth))
  ]

  ask patches [
    ; calculates amount of chla for each patch (correction with decay and growth
    ; is done to achieve the same maximum throughout the simulation if constant food supply is chosen
    set chla (chla_max * resolution * (1 - chla_decay / chla_growth))

    ; scale green coloring in relation to maximum possible chla value (chla_max * resolution * (1 - chla_decay / chla_growth))
    ; and stretch it over a range of 4 color values
    set pcolor (59 - chla / (chla_max * resolution * (1 - chla_decay / chla_growth)) * 4)
  ]

  ; creation of oozoids ----------------------------------------------------------------------------------------------------------------------------- ;

  if species != "krill" [

    create-oozoids 2 [                             ; creates 2 oozoids
      set number 1                                 ; 1 individual
      set l 2                                      ; body length l = 2 cm
      set size 2                                   ; size as depicted in NetLogo
      set color 45                                 ; salps are displayed as yellow
      set n_repro 0                                ; no chains released so far
      set generation 0                             ; no generation this season so far
      set c_weight ((l * 10 / 17) ^(1 / 0.4))      ; carbon weight c_weight is calculated according to an equation provided by Henschke et al 2018
      set c_reproduction (c_weight / 4)            ; an initial value is attributed to the reproductive storage c_reproduction
      setxy random-xcor random-ycor                ; random location
    ]
    set migrationevent? true                       ; store that immigration of salps took place
    set n_oozoids (count oozoids)                  ; number of oozoids
  ]


  ; creation of krill ------------------------------------------------------------------------------------------------------------------------------- ;

  if species != "salps" [            ; krill is produced only if chosen so
    create-krills krill_amount [     ; create N_krill
      set L 0.34                     ; structural length of 1.7 mm - C1 stage with 30 days of age (Ikeda 1984 J. Exp. Mar. Biol. Ecol.)
      set W_V 0.22 * ( L) ^ 3        ; structural body mass
      set W_B 0.0                    ; no assimilate buffer in egg so far
      set W_R 0                      ; no reproduction buffer so far
      set n_spawn 0                  ; no spawning events so far
      set d_age 30                   ; age in days
      set d_first_reproduction 0     ; no reproduction so far
      setxy random-xcor random-ycor  ; random location
      set d_starvation 0             ; no starvation so far
      set color black                ; krill are displayed as black
    ]
  ]

  reset-ticks                      ; update all plots and start the NetLogo clock

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; go procedure                        ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go


  ;if (ticks >= 100000) [stop]                                 ; simulation stops after it has reached its time span
  if (ticks >= (20 * 365 - 30)) [stop]                        ; simulation stops after it has reached its time span

  if SA? and (ticks >= (6 * 365 - 30)) [stop]                 ; stop simulation during SA if one krill lifecycle finished (6 years minus 30 days)

  let salp-set (turtle-set oozoids chains blastozoids)        ; create turtle-set of all present salp life stations
  if (any? salp-set) [                                        ; check if any salps are around
    set rep_cycles_max (max [generation] of (salp-set ))      ; store the highest generation time reflecting number of life cycles
  ]

  if (ticks mod 365 = 180) [                           ; in the middle of the winter - end of june
    calculate_salp_cycle                               ; determine the largest generation time of the last season
  ]

  grow                                                 ; determine growth in body length
  asexual_repro                                        ; asexual reproduction of salps (oozoids)
  sexual_repro                                         ; sexual reproduction of krill and salps (blastozoids)
  death                                                ; mortality

  if species != "krill" [                              ; should salps be simulated?
    immigration                                        ; immigration of salps into the simulation arena
  ]

  update-patches                                       ; patches will be updated (e.g. chla content, density dependent krill survival)
  move                                                 ; random walk of krill and salps
  calc_globals                                         ; calculate globals and results
  if plots? [
    do_graphs                                          ; update plots
  ]

  tick                                                 ; advance time step by one

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; store amount of reproductive cycles ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calculate_salp_cycle

  set reg_cycles rep_cycles_max                                            ; store max amount of repro cycles
  if (reg_cycles > reg_cycles_max) [set reg_cycles_max reg_cycles]         ; correct amount of repro cycles to maximum possible
  set rep_cycles_max 0                                                     ; set max-gen to 0

  ; set generation of all life cycles of salp to 0
  ask (turtle-set blastozoids oozoids chains) [
    set generation 0
  ]

  if plots?[
    ; update Plot on reproduction cycles
    set-current-plot "Seasonal repoduction cycles"
    plot reg_cycles / 2
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; grow procedure                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to grow

  ; in the grow procedure the change in body length will be calculated for all animals, this is done for each patch checking whether there is enough food available or not

  set T (cos ((ticks) / 365 * 360) * 2 + 273)         ; set temperature based on time of year
  set T_factor_s (exp (8000 / 275 - 8000 / T))        ; Arrhenius relation for salps taken from Groenefeld et al. (2020)
  ;let T_factor_k (exp (7421 / 275 - 7421 / T))     ; Arrhenius relation for krill taken from Bahlburg et al. (2021)
  let T_factor_k (exp (8000 / 275 - 8000 / T))     ; Arrhenius relation for krill taken from Bahlburg et al. (2021)
  let dayt ticks mod 365                              ; day of the year - used for krill below


  ; the following procedures are conducted per patch - assuming this is the local area where individuals compete for ressources ------------------
  ask patches [

    let need 0               ; need is the amount of food that all individuals on the patch want to consume

    ; this is assuming that all animals have the same Holling Type II functional response, meaning having the same half saturation curve, which is certainly not true
    ; in essence it means more food is better for krill and salps which is debated
    ; below the need in a patch is calculated and if need is larger than the available food the uptake is restricted, i.e. it cannot be more food consumed than it is present

    ; this is the Holling type II functional response (fr) that animals would have given enough ressources - this depends on the density of chla per m^3 and not the total amount
    ; halfsat (!slider!, standard value: 0.2)
    let fr_s (chla / resolution) / ( (chla / resolution) + salp_halfsat)
    let fr_k (chla / resolution) / ( (chla / resolution) + krill_halfsat)

    ask oozoids-here [
      set need (need + T_factor_s * number * (0.0025 * fr_s * l ^ 2)) ; for oozoids number is always one, the uptake of food depends on the squared body length
    ]

    ask blastozoids-here [
      set need (need + T_factor_s * number * (0.0025 * fr_s * l ^ 2)) ; for blastozoids number is always one, the uptake of food depends on the squared body length
    ]

    ask chains-here [
      set need (need + T_factor_s * number * (0.0025 * fr_s * l ^ 2)) ; for chains number can vary, since number represents the number of blastozoids in the chain
    ]

    ; potential food uptake for krill
    let J_A_local 0

    ask krills-here [
      ifelse (dayt > 90 and dayt < 270 and l > 35 * 0.2) [        ; this is unique to krill, if individuals are adults here l > 35 mm, than they can reduce metabolic activity during winter
        let J_a_Am (krill_hibernation / 100 * 0.044)              ; J_a_Am: maximum area-specific assimilation rate (Jager & Ravagnan 2015)
        set J_A_local (fr_k * J_a_Am * L ^ 2 * T_factor_k)          ; potential food uptake
        set need (need + J_A_local / 48)                          ; converting J_A into food requires to divide by Yax = 0.8 and to convert mg C into mg chla c:chla = 60 and therefore 48 = 0.8 * 60
      ][
        let J_a_Am 0.044                                          ; J_a_Am: maximum area-specific assimilation rate (Jager & Ravagnan 2015)
        set J_A_local (fr_k * J_a_Am * L ^ 2 * T_factor_k)          ; potential food uptake
        set need (need + J_A_local / 48)                          ; converting J_A into food requires to divide by Yax = 0.8 and to convert mg C into mg chla c:chla = 60 and therefore 48 = 0.8 * 60
      ]
    ]

    ask clutches-here [
      if (d_age > 30) [
        let J_a_Am 0.044                                            ; J_a_Am: maximum area-specific assimilation rate (Jager & Ravagnan 2015)
        set J_A_local (fr_k * J_a_Am * L ^ 2 * T_factor_k * number)   ; here uptake is multiplied by the number of small krill in the cohort
        set need (need + J_A_local / 48)                            ; converting J_A into food requires to divide by Yax = 0.8 and to convert mg C into mg chla c:chla = 60 and therefore 48 = 0.8 * 60
      ]
    ]

    if chla < need [ ; if more food is requested than available, than the functional response has to be updated
      set fr_s (fr_s * chla / need)
      set fr_k (fr_k * chla / need)
    ]

    ; growth of oozoids
    ask oozoids-here [
      let c_assimilated (T_factor_s * fr_s * 0.0025 * l ^ 2 * 0.64 * 50)      ; this is assumed to be the assimilated carbon assuming a C:Chla ratio of 50
      let c_growth (0.85 * c_assimilated)                                   ; 85% of the carbon is allocated to growth
      let c_repro_local (0.15 * c_assimilated)                              ; 15% of the carbon is allocated to reproduction
      let oldl l
      let starvationlocal false
      set chla (chla - T_factor_s * number * (0.0025 * fr_s * l ^ 2))         ; reduction of the overall available chla in the patch
      let c_respiration (T_factor_s * oozoid_resp / 100 * c_weight)         ; respiration costs based on carbon weight (Iguchi 2004)
      ifelse (c_respiration > c_growth) [                                   ; if respiration is higher than carbon allocated to growth
        ifelse c_respiration < (c_growth + c_repro_local) [                 ; if respiration can be covered by assimilated carbon
          set c_repro_local (c_repro_local + c_growth - c_respiration)
          set c_growth 0
        ][; if respiration can be covered by assimilated carbon and the reproduction storage
          ifelse c_respiration < (c_growth + c_repro_local + c_reproduction) [
            set c_growth 0
            set c_repro_local 0
            set c_reproduction (c_reproduction + c_growth + c_repro_local - c_respiration)
          ][; respiration loss cannot be covered at all - animal starves
            set c_growth 0
            set c_repro_local 0
            set c_reproduction 0
            set starvationlocal true
          ]
        ]
      ][; if respiration is less than carbon allocation to growth
        set c_growth (c_growth - c_respiration)
        if (c_growth < 0) [user-message "carbon-growth"]
      ]
      let oldc_weight c_weight                                                        ; this is for testing
      set c_weight (c_weight + c_growth)                                              ; new growth
      if (c_weight - oldc_weight < 0 ) [user-message "c_weight shrinkage"]            ; salps are not allowed to shrink
      set c_reproduction ((1 - oozoid_resp / 100) * c_reproduction + c_repro_local)   ; c_reproduction is updated - here are respiration costs for the reproductive tissue calculated - this is not done in the DEB
      set l ((17 * c_weight ^ 0.4) / 10)                                              ; conversion from carbon weight to body length - from an empirical relationship cited in Henschke et al. 2018
      if (l - oldl < -0.00001) [
        user-message word "l smaller than old l" (l - oldl)
        user-message word "l " l
        user-message word "oldl " oldl
        user-message word "old c_weight" oldc_weight
        user-message word "c_weight " c_weight
      ]
      ifelse (starvationlocal = true) [
        set d_starvation d_starvation + 1 ; increasing starvationday counter if animal could not fullfill respiration requirements
      ][]

      if MeasureInc? [
        if (ticks < 5000 ) [             ; limited duration of measurement because of memory issues
         set gr_o lput (l - oldl) gr_o   ; measuring the increment in growth
        ]
      ]
      set size 2 ; displayed agent size
    ]

    ; growth of blastozoids
    ask blastozoids-here [
      let c_assimilated T_factor_s * fr_s * 0.0025 * l ^ 2 * 0.64 * 50       ; this is assumed to be the assimilated carbon assuming a C:Chla ratio of 50
      let c_growth (0.8 * c_assimilated)                                   ; 85% of the carbon is allocated to growth
      let c_repro_local (0.2 * c_assimilated)                              ; 15% of the carbon is allocated to reproduction
      let oldl l
      let starvationlocal false
      set chla (chla - T_factor_s * number * (0.0025 * fr_s * l ^ 2))        ; reduction of the overall available chla in the patch
      let c_respiration (T_factor_s * blasto_resp / 100 * c_weight)        ; respiration costs based on carbon weight (Iguchi 2004)
      ifelse (c_respiration > c_growth) [                                  ; if respiration is higher than carbon allocated to growth
        ifelse c_respiration < (c_growth + c_repro_local) [                ; if respiration is higher than the carbon allocated to growth and the c_reproduction the animal starves
          set c_repro_local (c_repro_local + c_growth - c_respiration)
          set c_growth 0
        ][                                                                   ; otherwise the animal does not grow and pays from the c_reproduction
          ifelse c_respiration < (c_growth + c_repro_local + c_reproduction) [
            set c_growth 0
            set c_repro_local 0
            set c_reproduction (c_reproduction + c_growth + c_repro_local - c_respiration)
          ][
            set c_growth 0
            set c_repro_local 0
            set c_reproduction 0
            set starvationlocal true
          ]
        ]
      ][ ; if respiration is less than carbon allocation to growth
        set c_growth (c_growth - c_respiration )
        if (c_growth < 0) [user-message "carbon-growth"]
      ]
      set c_weight (c_weight + c_growth)
      set c_reproduction ((1 - blasto_resp / 100) * c_reproduction + c_repro_local)
      set l ((17 * c_weight ^ 0.4) / 10)
      ifelse (starvationlocal = true) [
        set d_starvation (d_starvation + 1)
      ][]
      if MeasureInc? [
        if (ticks < 5000 ) [
          set gr_b (lput (l - oldl) gr_b)
        ]
      ]
      set size 2
    ]

    ; growth of chains
    ask chains-here [
      let c_assimilated (T_factor_s * fr_s * 0.0025 * l ^ 2 * 0.64 * 50)
      let c_growth (0.85 * c_assimilated)
      let c_repro_local (0.15 * c_assimilated)
      let oldl l
      let starvationlocal false
      set chla (chla - T_factor_s * number * (0.0025 * fr_s * l ^ 2))
      let c_respiration (T_factor_s * blasto_resp / 100 * c_weight)
      ifelse (c_respiration > c_growth) [                                       ; if respiration is higher than carbon allocated to growth
        ifelse c_respiration < (c_growth + c_repro_local) [                     ; if respiration is higher than the carbon allocated to growth and the c_reproduction the animal starves
          set c_repro_local (c_repro_local + c_growth - c_respiration)
          set c_growth 0
        ][                                                                     ; otherwise the animal does not grow and pays from the c_reproduction
          ifelse c_respiration < (c_growth + c_repro_local + c_reproduction) [
            set c_growth 0
            set c_repro_local 0
            set c_reproduction (c_reproduction + c_growth + c_repro_local - c_respiration)
          ][
            set c_growth 0
            set c_repro_local 0
            set c_reproduction 0
            set starvationlocal true
          ]
        ]
      ][                                                                       ; if respiration is less than carbon allocation to growth
        set c_growth (c_growth - c_respiration)
        if (c_growth < 0) [user-message "carbon-growth"]
      ]
      set c_weight (c_weight + c_growth)
      set c_reproduction ((1 - blasto_resp / 100) * c_reproduction + c_repro_local)
      set l ((17 * c_weight ^ 0.4) / 10)
      ifelse (starvationlocal = true) [
        set d_starvation (d_starvation + 1)
      ][]
      if MeasureInc? [
        if (ticks < 5000 ) [
          set gr_b (lput (l - oldl) gr_b)
        ]
      ]
      set size 2
    ]

    ; growth of adult krill
    ask krills-here [
      ifelse (l > 11 * 0.2) [        ; check if krill is juvenile or adult (Jager & Ravagnan 2015)
        determine_fluxes fr_k          ; fluxes deterimend after Jager & Ravagnan (2015)
      ][
        determine_fluxes_larvae fr_k
      ]
      set W_V (W_V + J_V)               ; new W_V is W_V + J_V, shrinkage happens in deterime_fluxes
      let dV 0.22                       ; dry weight density (Jager & Ravagnan 2015)
      set L ((W_V / dV) ^ (1 / 3))
      set chla (chla - (J_A / 48))
    ]

    ; growth of larvae krill
    ask clutches-here [
      if d_age > 30 [
        ifelse (l > 11 * 0.2) [      ; check if kirll is juvenile or adult (Jager & Ravagnan 2015)
          determine_fluxes fr_k        ; fluxes deterimend after Jager & Ravagnan (2015)
        ][
          determine_fluxes_larvae fr_k
        ]
        set W_V (W_V + J_V)             ; new W_V is W_V + J_V, shrinkage happens in deterime_fluxes
        set W_R (W_R + J_R)
        let dV 0.22                     ; dry weight density (Jager & Ravagnan 2015)
        set L ((W_V / dV) ^ (1 / 3))
        set chla (chla - number * (J_A / 48))
      ]
    ]

    ; It is important to set a minimum amount of chla at the end. This is due to diffusion from other patches or deeper layers.
    ; chla > 0 is important to model the ressource.
    if (chla <= 0.005) [set chla 0.005]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; determine fluxes                    ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to determine_fluxes [frac]
  ; fluxes determined following Jager & Ravagnan (2015)

  let dayt (ticks mod 365) ; get day of year based on ticks

  let T_factor_k (exp (7421 / 274.15 - 7421 / T)) ; Arrhenius relation from Bahlburg et al. (2021)

  ; adult krill may have a reduced metabolic activity during winter
  ; length for adults from Jager & Ravagnan (2015)
  ifelse (dayt > 90 and dayt < 270 and l > 35 * 0.2) [

    ; assimilation
    let J_a_Am (krill_hibernation / 100 * 0.044)   ; maximum area-specific assimilation rate (Jager & Ravagnan 2015)
    set J_A frac * J_a_Am * L ^ 2 * T_factor_k     ; mass flux for assimilation

    ; respiration
    let J_VM (krill_hibernation / 100 * 0.0032)   ; volume-specific maintenance cost (Jager & Ravagnan 2015)
    set J_M J_VM * L ^ 3 * T_factor_k             ; mass flux for maintenance

    ; growth
    let Y_VA 0.8                             ; yield of assimilates on structure (Jager & Ravagnan 2015)
    let kappa 0.8                            ; fraction of assimilation flux for soma (Jager & Ravagnan 2015)
    set J_V (Y_VA * (kappa * J_A - J_M))     ; mass flux for structure

    ; reproduction
    set J_R (1 - kappa) * J_A                ; mass flux to reproduction buffer

    if (J_V < 0) [
      let delta J_V
      ifelse (J_R > abs (delta)) [
        set J_V 0                            ; mass flux for structure
        set J_R (J_R + delta)                ; mass flux to reproduction buffer
      ][
        set J_V (J_V + J_R)                  ; mass flux for structure
        set J_R 0                            ; mass flux to reproduction buffer
        ifelse ( W_R > abs (J_V)) [
          set W_R (W_R + J_V)                          ; mass of reproduction buffer in adult
          set J_V 0                                    ; mass flux for structure
        ][; respiration cannot be covered by assimilated carbon and reproduction buffer,
          ; than krill shrink (W_V is getting smaller since J_V is negative)
          set J_V (J_V + W_R)                          ; mass flux for structure
          set W_R 0                                    ; mass of reproduction buffer in adult
          set W_V (W_V + J_V)                          ; mass of sturctural body
          set J_V 0                                    ; mass flux for structure
          set d_starvation (d_starvation + 1)          ; increase counter for starvation
        ]
      ]
    ]
  ][; fluxes during summer
    let J_a_Am 0.044                               ; maximum area-specific assimilation rate (Jager & Ravagnan (2015)
    set J_A (frac * J_a_Am * L ^ 2 * T_factor_k)   ; assimilation flux
    let J_VM 0.0032                                ; volume specific maintenance cost (Jager & Ravagnan 2015)
    set J_M (J_VM * L ^ 3 * T_factor_k)            ; maintenance flux
    let Y_VA 0.8                                   ; yield of assimilates on structure (Jager & Ravagnan 2015)
    let kappa 0.8                                  ; fraction of assimilation flux for soma (Jager & Ravagnan 2015)
    set J_V (Y_VA * (kappa * J_A - J_M))           ; mass flux for structure
    set J_R ((1 - kappa) * J_A)                    ; mass flux for reproduction buffer
    if (J_V < 0) [
      let delta J_V
      ifelse (J_R > abs (delta)) [
        set J_V 0                                     ; mass flux for structure
        set J_R (J_R + delta)                         ; mass flux to reproduction buffer
      ][
        set J_V (J_V + J_R)                           ; mass flux for structure
        set J_R 0                                     ; mass flux for reproduction buffer
        ifelse (W_R > abs (J_V)) [
          set W_R (W_R + J_V)                         ; mass of reproduction buffer in adult
          set J_V 0                                   ; mass flux for structure
        ][
          set J_V (J_V + W_R)                         ; mass flux for structure
          set W_R 0                                   ; mass of reproduction buffer in adult
          set W_V (W_V + J_V)                         ; mass of structural body
          set J_V 0                                   ; mass flux for structure
          set d_starvation (d_starvation + 1)         ; increase starvation counter
        ]
      ]
    ]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; determine fluxes larvae             ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to determine_fluxes_larvae [frac]
  ; frac is the functional response (Holling II)  fluxes determined following Jager et al. 2015
  ; but how are the fluxes different for small krill (larvae and juvenile) from the one of adult krill
  ; here following Jager et al. 2015 - there is only the growth part stored.
  ; the flux into reproduction J_R is burned

  let dayt (ticks mod 365)                             ; get day of year from tick number
  let T_factor_k (exp (7421 / 274.15 - 7421 / (T)))    ; Arrhenius relation for krill from Bahlburg et al. (2021)

  let J_a_Am 0.044                                     ; maximum area-specific assimilation rate (Jager & Ravagnan 2015)
  set J_A (frac * J_a_Am * L ^ 2 * T_factor_k)         ; mass flux for assimilation

  let J_VM 0.0032                                      ; volume-specific maintenance cost (Jager & Ravagnan 2015)
  set J_M (J_VM * L ^ 3 * T_factor_k)                  ; mass flux for maintenance

  let Y_VA 0.8                                        ; yield of assimilates on structure (Jager et al. 2015)
  let kappa 0.8                                       ; fraction of assimilation flux for soma (Jager et al. 2015)
  set J_V (Y_VA * (kappa * J_A - J_M))                ; mass flux for structure
  set J_R 0                                           ; mass flux to reproduction buffer

  if (J_V < 0) [
    set W_V (W_V + J_V)                               ; mass of structural body
    set J_V 0                                         ; mass flux for storage
    set d_starvation (d_starvation + 1)               ; increase counter for days starving
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; asexual reproduction                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to asexual_repro
  ; chain release is the asexual part of reproduction - here it is length dependent
  ; if salps reach a given length the release a chain with a given number of bud (aggregates), but not more
  ; than it has been observed

  ask oozoids [
    if (l > 6.5 and n_repro = 0) [
      let budnumber (floor (c_reproduction / 0.0329))
      set budnumber (min list budnumber 150)

      hatch-chains 1 [ ;; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
        set generation (generation + 1)
        set number budnumber
        set l 0.5
        set size 2
        set d_age 0
        set sex 0
        set d_starvation 0
        set c_weight ((l * 10 / 17) ^(1 / 0.4))
        set c_reproduction (c_weight / 4)
      ]
      set n_repro (n_repro + 1)
      set t_regeneration d_age
      set c_reproduction (c_reproduction - budnumber * 0.0329)
    ]
    if (l > 8 and n_repro = 1) [
      let budnumber (floor (c_reproduction / 0.0329))
      set budnumber (min list budnumber 180)
      hatch-chains 1 [ ; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
        set generation (generation + 1)
        set number budnumber
        set l 0.5
        set size 2
        set d_age 0
        set sex 0
        set d_starvation 0
        set c_weight ((l * 10 / 17) ^(1 / 0.4))
        set c_reproduction (c_weight / 4)
      ]
      set n_repro (n_repro + 1)
      set c_reproduction (c_reproduction - budnumber * 0.0329)
    ]

    if (l > 9.5 and n_repro = 2) [
      let budnumber (floor (c_reproduction / 0.0329))
      set budnumber (min list budnumber 210)
      hatch-chains 1 [ ; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
        set generation (generation + 1)
        set number budnumber
        set l 0.5
        set size 2
        set d_age 0
        set sex 0
        set d_starvation 0
        set c_weight ((l * 10 / 17) ^(1 / 0.4))
        set c_reproduction (c_weight / 4)
      ]
      set n_repro (n_repro + 1)
      set c_reproduction (c_reproduction - budnumber * 0.0329)
    ]

    if (l > 10.5 and n_repro = 3) [
      let budnumber (floor (c_reproduction / 0.0329))
      set budnumber (min list budnumber 240)
      hatch-chains 1 [ ; initial size following Pakhomov and Hunt 2017 release size 3 - 5 mm
        set generation (generation + 1)
        set number budnumber
        set l 0.5
        set size 2
        set d_age 0
        set sex 0
        set d_starvation 0
        set c_weight ((l * 10 / 17) ^(1 / 0.4))
        set c_reproduction (c_weight / 4)
      ]
      set n_repro (n_repro + 1)
      set t_regeneration_list (lput (d_age - t_regeneration) t_regeneration_list)
      die
    ]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sexual reproduction                 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to sexual_repro
  ; sexual part of the salp reproduction. At a given size and sufficient carbon allocated in the c_reproduction
  ; a embryo is released with 0.7 survival. After that the blastozoids change sex and the chain dies,
  ; since the male blastozoids do not form chains

  ask chains [ ; embryo release is following Pakhomov and Hunt 2017, Deep- Sea Research II
    if (l > 2.5 and c_reproduction > 0.0269) [
      repeat number [
        ifelse (random-float 1 < 0.7) [
          hatch-oozoids 1 [
            set number 1
            set generation (generation + 1)
            set l 0.4
            set size 2
            set d_age 0
            set color 45
            set d_starvation 0
            set n_repro 0
            set c_weight ((l * 10 / 17) ^(1 / 0.4))
            set c_reproduction (c_weight / 4)
          ]
          hatch-blastozoids 1 [
            set sex 1
            set number 1
          ]
        ][
          hatch-blastozoids 1 [
            set sex 1
            set number 1
          ]
        ]
      ]
      die ; here the chain dies and only individuals survive
    ]
  ]

  ask krills [
    set W_R (W_R + J_R)
    ; the threshold for production of eggs is calculated based on the individual length
    let egg_threshold ((150.83 * L / 0.2 - 3027) * 0.028) ; equation from Bahlburg et al. (2021)
    ; if W_R exceeds egg threshold, length exceeds 35 mm and its a day of October to March
    ; (which was taken from Bahlburg et al. (2021), eggs are produced
    if (ticks mod 365 < 91) or (ticks mod 365 > 272 ) [
      if (W_R >= egg_threshold and L >= 7) [
        if n_spawn = 0 [
          set d_first_reproduction d_age       ; store age of first reproduction event
          if d_first_reproduction < d_k_firstrepro or d_k_firstrepro = 0 [
            set d_k_firstrepro d_first_reproduction ; store day of first repro event if its the earliest
          ]
        ]
        set n_spawn (n_spawn + 1)
        set n_eggs round (n_eggs + (egg_threshold / 0.028))
        if n_eggs > n_k_max_eggs [
          set n_k_max_eggs n_eggs
        ]
        set W_R (W_R - egg_threshold)          ; update reproduction buffer
        hatch-clutches 1 [
          set l (1.72 * 0.2)                   ; this is already the length of C1 stage, but we assume that the first 30 days krill will not feed, the growth will be fed by the energy from the eggs
          set W_V (0.22 * (L) ^ 3)             ; structural body mass
          set W_B 0                            ; assimilate buffer in egg
          set W_R 0                            ; reproduction buffer
          set J_A 0
          set J_V 0
          set J_R 0
          set J_M 0
          set n_spawn 0                        ; number of spawning events
          set d_age 0                          ; age in days
          set d_first_reproduction 0           ; age of first reproduction
          setxy random-xcor random-ycor        ; random location
          set number (egg_threshold / 0.028)   ; number of eggs in clutch
          set d_starvation 0                   ; days without food
        ]
      ]
    ]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; aging and dying                     ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to death
  ; salps die after half a year (Siegel 2017 mentiones a few month as maximum age)
  ; or due to daily mortality or due to starvation (30 days)

  ask oozoids [
    set d_age (d_age + 1)
    if (d_age > 182 or n_repro = 5) [die]
    if (d_starvation > salp_starvation) [die]
    if (random-float 1 < (salp_mortality / 100)) [die]
  ]

  ask blastozoids [
    set d_age (d_age + 1)
    if (d_age > 182) [die]
    if (d_starvation > salp_starvation) [die]
    if (random-float 1 < (salp_mortality / 100)) [die]
  ]

  ask chains [
    set d_age (d_age + 1)
    if (d_age > 182) [die]
    if (d_starvation > salp_starvation) [die]
    let counter 0
    if (number > 0) [
      repeat number [
        if (random-float 1 < (salp_mortality / 100)) [
          set counter (counter + 1)
        ]
      ]
    ]
    set number (number - counter)
    if (number <= 0) [die]
  ]

  ; krill dies after 6 years or due to daily mortality
  ask krills [
    set d_age (d_age + 1)
    if (d_age > (6 * 365)) or (random-float 1 < (krill_mortality / 100))[die]
    ; daily mortality from Auerswald et al. (2015)
  ]

  ; To reduce computational effort in the beginning clutches are computed as cohorts.
  ; If the number drops below 10 they are treated as individuals and inherit variables such as age and so on.
  ask clutches [
    set d_age (d_age + 1)
    set number (number * 0.95)
    if number < 10 [
      hatch-krills (round number) [
        set W_V (0.22 * L ^ 3) ; structural body mass
        set n_spawn 0
        set d_first_reproduction 0
      ]
      die
    ]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; immigration of salps                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to immigration
  ; there is no wintersurvival in the model for salps. Thus they migrate into the simulation area at some point in the season.
  if (count (turtle-set oozoids chains blastozoids) < 1) [
    if (random-float 1 < (salp_immiprob / 100)) [
      create-oozoids salp_amount [
        set number 1
        set l salp_length
        set size 2
        set color 45
        setxy random-xcor random-ycor
        set n_repro 0
        set generation 0
        set d_starvation 0
        set c_weight ((l * 10 / 17) ^(1 / 0.4))
        set c_reproduction (c_weight / 4)
      ]
      set t_immigration ticks
      set migrationevent? true
    ]
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; update patches                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-patches

  if (chla_supply = "Lognorm") [
    if (ticks mod 365 = 180 - vegetation_delay) [
      ; set fixed random seed
      random-seed 0
      let mu 3.83 ; meanlog from amlr data
      let sd 0.58 ; sdlog from amlr data

      ; repeat calclulation to match given year number
      repeat (ceiling (ticks / 365) + 1) [
        ; calculate chla_max using meanlog and sdlog from the amlr data
        set chla_max (e ^ (mu + sd * (random-normal 0 1)) / 100)
      ]
      ; randomized seed
      random-seed new-seed

      ; renormalization using the logistic equation term with loss term
      set chla_max (chla_max / (1 - chla_decay / chla_growth))
    ]
  ]


  ; ?
  let algae_growth (chla_growth * (0.5 * cos ((ticks + vegetation_delay) / 365 * 360) + 0.5))

  ask patches [
    ; update chla value: if chla is below 0.005 per m³, it is set to 0.005 per m³ to avoid chla depletion due to diffusion for example
    set chla (max list (chla + resolution * (algae_growth * (chla / resolution) * (1 - (chla / resolution) / chla_max) - chla_decay * (chla / resolution))) (0.005 * resolution))

    ; scale green coloring in relation to maximum possible chla value (chla_max * resolution * 0.8)
    ; and stretch it over a range of 4 color values
    set pcolor (59 - chla / (chla_max * resolution * (1 - chla_decay / chla_growth)) * 4)
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; movement                            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to move
  ; only random walk for salps and krill, but certainly krill may move further and even avoid salps?
  ask turtles [move-to one-of neighbors]

  ; Individual swimming speeds have been measured at between 15 and 30 cm/s (Kils 1981; Murphy et al. 2011).

  ; [...] In particular, the direction of swarm movement deviated from the background flow most when fluorescence levels were high,
  ; which is likely to be a response to retaining a favourable feeding environment (Tarling & Fielding 2016).

  ; Tarling and Thorpe (2014):
  ; [...] In the natural environment, both individuals and swarms have to contend with one of the world's strongest ocean currents, the
  ; Antarctic Circumpolar Current, with average flow rates of between 10 and 20 cm/s.
  ; [...] The overall cost of transport was estimated to be around 73% of total daily metabolic expenditure, at least during early summer. This
  ; compares to an estimate of 60% of total daily metabolic expenditure by Kils (1981), which he calculated through measuring the cost of hovering and assuming
  ; that forward propulsion at speeds of up to 25 cm s?1 could be attained at no extra cost through changing body angle.
  ; [...] At the individual level, Antarctic krill has been reported to alter activity levels according to time of year.

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate globals                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to calc_globals

  let n_month (floor ((ticks mod 365) / 30.5)) ; calculate month number based on ticks

  ; storing the number of salps for different months for the intraannual distribution later on
  set monthlycount (replace-item n_month monthlycount (item n_month monthlycount + sum [number] of (turtle-set oozoids blastozoids chains)))

  set n_oozoids (count oozoids) ; counting oozoids
  set n_blasto (count blastozoids + sum [number] of chains); counting all blastozoids, males and females

  set n_s_max_total (max (list n_s_max_total (n_oozoids + n_blasto))); maximum salp abundance during the simulation run

  if n_oozoids > 0 [ ; storing the largest body length of oozoids
    set l_o_max (max list l_o_max max [l] of oozoids)
  ]

  if n_blasto > 0 [ ; storing the largest body length of blastozoids, males and females (chains)
    set l_b_max (max list l_b_max max [l] of (turtle-set blastozoids chains))
  ]

  if (n_oozoids + n_blasto > n_s_max) [ ; storing the highest abundance during the season
    set n_s_max (n_oozoids + n_blasto)
    set t_immigration_max t_immigration
  ]

  if (count krills > 0) [
    ; calculate max spawning events
    set spawn_max max (list spawn_max (max [n_spawn] of krills))

    ; calculate max abundance
    set n_k_max max (list n_k_max (precision (count krills + sum [number] of clutches) 0))

    ; calc max size
    set l_k_max max (list l_k_max (precision (max [l] of krills) 2))
  ]

  ; end of july, at the end of the Southern winter when abundances should be really low
  if (ticks mod 365 = 210) [
    ; in n_s_max_list the highest salp abundance is stored that occured in the season (here the season ends after 210 days, i.e. end of July)
    set n_s_max_list (lput n_s_max n_s_max_list)
    ; only if in that season migration has happened it make sense to store the time
    ifelse migrationevent? = true [
      set t_immigration_list lput t_immigration_max t_immigration_list
      set migrationevent? false
    ][
      ; -1 indicates that no migration event has happened in that particular season
      set t_immigration_list lput -1 t_immigration_list
    ]

    if plots? [
      set-current-plot "Salp peaks"
      plot n_s_max
    ]

    ; resetting n_s_max
    set n_s_max 0

    set n_s_med (median n_s_max_list)
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; update all plots                    ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do_graphs

  set-current-plot "Abundances"
  set-current-plot-pen "Aggregates"
  plotxy (ticks) n_blasto
  set-current-plot-pen "Solitaries"
  plotxy (ticks)  n_oozoids

  set-current-plot "Max Chla"
  set-current-plot-pen "pen-1"
  plotxy (ticks) (max [chla] of patches / resolution)

  set-current-plot "Maximum krill life time spawn events"
  set-current-plot-pen "default"
  plotxy (ticks) spawn_max

  set-current-plot "Krill abundance"
  set-current-plot-pen "default"
  plotxy (ticks) (count krills + sum [number] of clutches)

  set-current-plot "Abundance adult Krill"
  set-current-plot-pen "default"
  plotxy (ticks) (count krills)

  set-current-plot "T_factor_s"
  set-current-plot-pen "default"
  plot T_factor_s

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
  histogram [d_age] of oozoids

  set-current-plot "Age of Blastozoids"
  set-plot-pen-mode 1
  set-plot-x-range 0 500
  set-histogram-num-bars 50
  histogram [d_age] of (turtle-set blastozoids chains)

  ; to reduce computational effort this is only computed once a year
  if (ticks mod 365 = 50) and species != "krill" [
    if MeasureInc? [

      set-current-plot "Growth Rates Oozoids"
      clear-plot
      set-plot-pen-mode 1
      set-plot-x-range 0 (max (gr_o) + 0.5)
      set-histogram-num-bars 50
      histogram gr_o

      set-current-plot "Growth Rates Blastozoids"
      set-plot-pen-mode 1
      let upper-value 1
      if (length gr_b > 0) [
        set upper-value max(gr_b)
      ]
      set-plot-x-range 0 (upper-value + 0.5)
      set-histogram-num-bars 50
      histogram gr_b
    ]

    set-current-plot "Annual distribution"
    clear-plot
    set-plot-pen-mode 1
    set-histogram-num-bars 12
    set-plot-x-range 0 13
    let summ sum monthlycount
    let monthlycount_temp monthlycount
    (foreach [0 1 2 3 4 5 6 7 8 9 10 11 ] monthlycount_temp [[a b] -> set monthlycount_temp (replace-item a monthlycount_temp (((item a monthlycount_temp / summ))))])
    (foreach [0 1 2 3 4 5 6 7 8 9 10 11 ] monthlycount_temp [[a b] -> set monthlycount_temp (replace-item a monthlycount_temp (precision (item a monthlycount_temp) 3))])
    set-plot-y-range 0 (precision (max (monthlycount_temp) + 0.1) 2)
    (foreach [1 2 3 4 5 6 7 8 9 10 11 12] monthlycount_temp [[a b] -> plotxy a b])
  ]

  set-current-plot "Chain releases"
  set-plot-pen-mode 1
  set-plot-x-range 0 6
  set-histogram-num-bars 6
  histogram [n_repro] of oozoids

  set-current-plot "Bud/Oozoid"
  set-plot-pen-mode 0
  let temp n_blasto
  let temp2 n_oozoids
  if (temp2 > 0) [
    plotxy ticks (temp / temp2)
    set ratio_list (lput (temp / temp2) ratio_list)
  ]

  set-current-plot "Krill Growth Curves"
  foreach [who] of krills [x -> plotxy (ticks) ([l / 0.2] of krill x)]

end
@#$#@#$#@
GRAPHICS-WINDOW
495
10
907
423
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
15
20
95
53
NIL
setup
NIL
1
T
OBSERVER
NIL
S
NIL
NIL
1

BUTTON
105
20
185
53
NIL
go
T
1
T
OBSERVER
NIL
W
NIL
NIL
1

PLOT
910
10
1140
295
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
"Aggregates" 1.0 0 -16777216 true "" ";plotxy (ticks) n_blasto"
"Solitaries" 1.0 0 -2674135 true "" ";plotxy (ticks)  n_oozoids"

PLOT
1145
10
1305
295
Max Chla
NIL
NIL
0.0
10.0
0.0
1.1
true
false
"" ""
PENS
"pen-1" 1.0 0 -13840069 true "" ";plot max [chla] of patches / resolution"

PLOT
655
465
920
745
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
925
465
1176
745
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
1080
905
1277
1043
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
1079
754
1279
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
1572
465
1861
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
1565
760
1859
964
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
1675
300
1750
345
NIL
mean gr_o
5
1
11

MONITOR
1755
300
1835
345
NIL
mean gr_b
5
1
11

SLIDER
15
540
185
573
salp_starvation
salp_starvation
0
1000
21.0
1
1
d
HORIZONTAL

PLOT
654
755
1069
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
16
585
188
618
salp_mortality
salp_mortality
0
100
2.1
0.1
1
% / d
HORIZONTAL

SLIDER
15
315
185
348
vegetation_delay
vegetation_delay
0
180
47.0
1
1
d
HORIZONTAL

BUTTON
195
20
275
53
ref-values
set chla_growth 0.25\nset chla_decay 0.05\nset vegetation_delay 45\nset salp_halfsat 0.20\nset salp_immiprob 0.85\nset salp_amount 10\nset salp_length 3.0\nset oozoid_resp 5.0\nset blasto_resp 15.0\nset salp_starvation 30\nset salp_mortality 2.5\nset krill_halfsat 0.106\nset krill_hibernation 20\nset krill_amount 30\nset krill_mortality 0.07
NIL
1
T
OBSERVER
NIL
A
NIL
NIL
1

PLOT
1291
756
1557
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
1180
465
1565
745
Seasonal repoduction cycles
NIL
NIL
0.0
10.0
0.0
3.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

PLOT
1287
301
1458
453
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
1755
350
1835
395
Max Density
precision (n_s_max /  (101 * 101 * 16)) 2
17
1
11

PLOT
920
302
1080
452
T_factor_s
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
"default" 1.0 0 -16777216 true "" ";plot temp_factor"

CHOOSER
15
70
115
115
chla_supply
chla_supply
"Const" "Lognorm"
0

PLOT
1086
302
1286
452
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
15
420
185
453
salp_immiprob
salp_immiprob
0
100
0.969
0.01
1
%
HORIZONTAL

SLIDER
15
230
185
263
chla_growth
chla_growth
0
1
0.318
0.01
1
NIL
HORIZONTAL

SLIDER
15
275
187
308
chla_decay
chla_decay
0
1
0.0457
0.01
1
NIL
HORIZONTAL

SWITCH
16
169
142
202
MeasureInc?
MeasureInc?
1
1
-1000

SLIDER
15
380
257
413
salp_halfsat
salp_halfsat
0
1
0.208
0.01
1
mg Chla / m³
HORIZONTAL

SLIDER
15
460
187
493
salp_amount
salp_amount
0
100
14.0
1
1
n
HORIZONTAL

SLIDER
15
500
187
533
salp_length
salp_length
0
5
2.6
0.1
1
cm
HORIZONTAL

PLOT
1310
10
1509
291
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
15
765
185
798
krill_amount
krill_amount
0
2000
20.0
1
1
n
HORIZONTAL

SLIDER
15
845
185
878
krill_hibernation
krill_hibernation
0
100
13.0
0.1
1
%
HORIZONTAL

PLOT
1515
10
1685
290
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
"default" 1.0 0 -16777216 true "" ";plot count debkrill + sum [number] of clutches"

PLOT
1465
300
1665
450
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
"default" 1.0 0 -16777216 true "" ";plot maxspawn"

SLIDER
15
625
185
658
oozoid_resp
oozoid_resp
0
100
5.27
0.1
1
% / d
HORIZONTAL

SLIDER
15
665
185
698
blasto_resp
blasto_resp
0
100
8.35
0.1
1
% / d
HORIZONTAL

PLOT
1690
10
1865
290
Abundance adult Krill
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
"default" 1.0 0 -16777216 true "" ";plot count debkrill"

TEXTBOX
190
540
350
581
starvation capability (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
195
585
350
626
daily mortality (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
195
325
355
351
(Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
190
420
370
450
immigration probability (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
190
230
360
271
rate of primary production (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
190
275
350
306
(Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
265
380
420
415
half saturarion salps (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
195
845
460
875
reduced metabolism of adult Krill during winter; Atkinson et al. (2002): < 30 %
12
0.0
1

TEXTBOX
190
460
355
495
immigration group size (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
190
500
360
530
length of individuals (Groeneveld et al. 2020)
12
0.0
1

TEXTBOX
190
765
365
795
initial population size
12
0.0
1

MONITOR
1675
350
1750
395
# krill
(count krills) + (round (sum [number] of clutches))
0
1
11

TEXTBOX
195
625
345
660
C loss for respiration (Iguchi 2004): 2.5-4.9
12
0.0
1

TEXTBOX
195
665
345
695
C loss for respiration (Iguchi 2004): 0.8-6.1
12
0.0
1

SLIDER
15
805
187
838
krill_mortality
krill_mortality
0
100
0.102
0.01
1
% / d
HORIZONTAL

TEXTBOX
20
360
235
386
----------- Salp parameter -----------
12
0.0
1

TEXTBOX
15
705
255
731
----------- Krill parameter -----------
12
0.0
1

TEXTBOX
195
805
350
846
daily mortality (Auerswald et al., 2015)
12
0.0
1

TEXTBOX
120
70
290
111
constant or varying Chla (based on AMLR data)
12
0.0
1

TEXTBOX
155
170
305
200
on: measure Salp growth
12
0.0
1

TEXTBOX
15
210
260
228
----------- world parameter -----------
12
0.0
1

SWITCH
510
450
613
483
SA?
SA?
0
1
-1000

MONITOR
1675
400
1750
445
year
floor (ticks / 365) + 1
0
1
11

MONITOR
1755
400
1835
445
NIL
chla_max
3
1
11

SWITCH
510
490
613
523
plots?
plots?
1
1
-1000

SLIDER
15
725
247
758
krill_halfsat
krill_halfsat
0
1
0.0641
0.01
1
mg Chla / m³
HORIZONTAL

TEXTBOX
255
725
415
760
half saturation krill (Atkinson et al. 2006)
12
0.0
1

CHOOSER
15
120
115
165
species
species
"krill" "salps" "both"
2

TEXTBOX
120
115
370
161
\"krill\" - only krill individuals, \"salps\" - only salp individuals, \"both\" - both species simultaneously
12
0.0
1

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
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="SA" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>vegetation_delay</metric>
    <metric>oozoid_resp</metric>
    <metric>krill_amount</metric>
    <metric>chla_growth</metric>
    <metric>salp_immiprob</metric>
    <metric>salp_mortality</metric>
    <metric>krill_mortality</metric>
    <metric>blasto_resp</metric>
    <metric>salp_amount</metric>
    <metric>chla_decay</metric>
    <metric>salp_length</metric>
    <metric>salp_starvation</metric>
    <metric>krill_hibernation</metric>
    <metric>salp_halfsat</metric>
    <metric>krill_halfsat</metric>
    <metric>n_s_max_total</metric>
    <metric>n_s_med</metric>
    <metric>d_k_firstrepro</metric>
    <metric>l_k_max</metric>
    <metric>n_k_max_eggs</metric>
  </experiment>
  <experiment name="k-halfsat" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2160"/>
    <metric>krill_halfsat</metric>
    <metric>l_k_max</metric>
    <metric>n_k_max_eggs</metric>
    <metric>d_k_firstrepro</metric>
    <enumeratedValueSet variable="SA?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vegetation_delay">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="oozoid_resp">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salp_immiprob">
      <value value="0.85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chla_growth">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salp_starvation">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salp_amount">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salp_mortality">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="krill_amount">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="plots?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="krill_hibernation">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="blasto_resp">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salps?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chla_decay">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salp_length">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="salp_halfsat">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeasureInc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chla_supply">
      <value value="&quot;Const&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="krill_mortality">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="population" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="4380"/>
    <metric>n_s_max_total</metric>
    <metric>n_s_med</metric>
    <metric>l_o_max</metric>
    <metric>l_b_max</metric>
    <metric>rep_cycles_max</metric>
    <metric>n_k_max</metric>
    <metric>l_k_max</metric>
    <metric>n_k_max_eggs</metric>
    <metric>d_k_firstrepro</metric>
    <metric>spawn_max</metric>
    <enumeratedValueSet variable="species">
      <value value="&quot;krill&quot;"/>
      <value value="&quot;salps&quot;"/>
      <value value="&quot;both&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chla_supply">
      <value value="&quot;Const&quot;"/>
      <value value="&quot;Lognorm&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="population2" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2160"/>
    <metric>n_s_max_total</metric>
    <metric>n_s_med</metric>
    <metric>l_o_max</metric>
    <metric>l_b_max</metric>
    <metric>rep_cycles_max</metric>
    <metric>n_k_max</metric>
    <metric>l_k_max</metric>
    <metric>n_k_max_eggs</metric>
    <metric>d_k_firstrepro</metric>
    <metric>spawn_max</metric>
    <enumeratedValueSet variable="species">
      <value value="&quot;krill&quot;"/>
      <value value="&quot;salps&quot;"/>
      <value value="&quot;both&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chla_supply">
      <value value="&quot;Const&quot;"/>
      <value value="&quot;Lognorm&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ranges" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>chla_growth</metric>
    <metric>chla_decay</metric>
    <metric>vegetation_delay</metric>
    <metric>salp_halfsat</metric>
    <metric>salp_immiprob</metric>
    <metric>salp_amount</metric>
    <metric>salp_length</metric>
    <metric>salp_starvation</metric>
    <metric>salp_mortality</metric>
    <metric>oozoid_resp</metric>
    <metric>blasto_resp</metric>
    <metric>krill_halfsat</metric>
    <metric>krill_amount</metric>
    <metric>krill_mortality</metric>
    <metric>krill_hibernation</metric>
    <metric>min [L] of krills</metric>
    <metric>max [L] of krills</metric>
    <metric>min [W_R] of krills</metric>
    <metric>max [W_R] of krills</metric>
    <metric>min [J_A] of krills</metric>
    <metric>max [J_A] of krills</metric>
    <metric>min [J_M] of krills</metric>
    <metric>max [J_M] of krills</metric>
    <metric>min [J_V] of krills</metric>
    <metric>max [J_V] of krills</metric>
    <metric>min [J_R] of krills</metric>
    <metric>max [J_R] of krills</metric>
    <metric>max [d_age] of krills</metric>
    <metric>min [W_V] of krills</metric>
    <metric>max [W_V] of krills</metric>
    <metric>max [n_spawn] of krills</metric>
    <metric>max[d_starvation] of krills</metric>
    <metric>max[d_first_reproduction] of krills</metric>
    <metric>min [W_B] of krills</metric>
    <metric>max[W_B] of krills</metric>
    <metric>min [L] of clutches</metric>
    <metric>max [L] of clutches</metric>
    <metric>min [W_R] of clutches</metric>
    <metric>max [W_R] of clutches</metric>
    <metric>min [J_A] of clutches</metric>
    <metric>max [J_A] of clutches</metric>
    <metric>min [J_M] of clutches</metric>
    <metric>max [J_M] of clutches</metric>
    <metric>min [J_V] of clutches</metric>
    <metric>max [J_V] of clutches</metric>
    <metric>min [J_R] of clutches</metric>
    <metric>max [J_R] of clutches</metric>
    <metric>max [d_age] of clutches</metric>
    <metric>min [W_V] of clutches</metric>
    <metric>max [W_V] of clutches</metric>
    <metric>min [number] of clutches</metric>
    <metric>max [number] of clutches</metric>
    <metric>max[d_starvation] of clutches</metric>
    <metric>max[d_first_reproduction] of clutches</metric>
    <metric>min [W_B] of clutches</metric>
    <metric>max[W_B] of clutches</metric>
    <metric>min [l] of oozoids</metric>
    <metric>max [l] of oozoids</metric>
    <metric>max [d_age] of oozoids</metric>
    <metric>max [n_repro] of oozoids</metric>
    <metric>max [generation] of oozoids</metric>
    <metric>min [t_regeneration] of oozoids</metric>
    <metric>max [t_regeneration] of oozoids</metric>
    <metric>min [c_weight] of oozoids</metric>
    <metric>max [c_weight] of oozoids</metric>
    <metric>min [c_reproduction] of oozoids</metric>
    <metric>max [c_reproduction] of oozoids</metric>
    <metric>max[d_starvation] of oozoids</metric>
    <metric>min [l] of blastozoids</metric>
    <metric>max [l] of blastozoids</metric>
    <metric>max [d_age] of blastozoids</metric>
    <metric>max [generation] of blastozoids</metric>
    <metric>min [c_weight] of blastozoids</metric>
    <metric>max [c_weight] of blastozoids</metric>
    <metric>min [c_reproduction] of blastozoids</metric>
    <metric>max [c_reproduction] of blastozoids</metric>
    <metric>max[d_starvation] of blastozoids</metric>
    <metric>min [l] of chains</metric>
    <metric>max [l] of chains</metric>
    <metric>max [d_age] of chains</metric>
    <metric>max [generation] of chains</metric>
    <metric>min [c_weight] of chains</metric>
    <metric>max [c_weight] of chains</metric>
    <metric>min [c_reproduction] of chains</metric>
    <metric>max [c_reproduction] of chains</metric>
    <metric>max[d_starvation] of chains</metric>
    <metric>min [number] of chains</metric>
    <metric>max [number] of chains</metric>
    <metric>min [chla] of patches</metric>
    <metric>max [chla] of patches</metric>
    <enumeratedValueSet variable="SA?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species">
      <value value="&quot;both&quot;"/>
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
