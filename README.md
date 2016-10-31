# origins_of_plasticity_f16
Code/scripts for evolutionary origins of phenotypic plasticity experiments. Relevant for ALife Journal article.

## The Question
What are the stepping stones for the evolution of phenotypic plasticity? Do these conditions seem to be
robust under a few different conditions? Can we identify alternative strategies to phenotypic plasticity
in changing environments?

## Hypotheses:
  * 'Stepping stones': Unconditional trait expression prior to conditional trait expression. Imprecise,
    but still beneficial plasticity before perfect plasticity (this is more useful for EC).
  * 'Alternative strategies': At higher mutation rates or longer environmental change rates, bet-hedging
    will emerge as a viable strategy as an alternative to phenotypic plasticity.

## Experimental Design


### Environments
  * Constant Environments (controls):
    * T1(+)
    * T2(+)
    * T1(+)T2(+)
  * Changing Environments
    * [T1(+), T1(-)]: NAND+, NAND-
    * [T2(+), T2(-)]: NOT+, NOT-
    * [T1(+)T2(-), T1(-)T2(+)]
    * [T1(+)T2(+), T1(+)T2(-), T1(-)T2(+), T1(-)T2(-)]
    * [T2(+)T1(+), T2(+)T1(-), T2(-)T1(+), T2(-)T1(-)]
  * Point Mutation Rates (Not yet)
    * Low: 0.0025
    * Baseline: 0.0075
    * High: 0.0125
  * Change Rates (in updates)
    * 1     (50  reps)
    * 25    (50  reps)
    * 50    (200 reps)
    * 100   (200 reps)
    * 200   (200 reps)
    * 400   (200 reps)
    * 800   (50  reps)
    * 1600  (50  reps)
    * 25k   (50  reps)
  * Birth Method
    * Well Mixed (4)
  * Sensors
    * T1 Sensor (1 if rewarded, -1 if punished, 0 otherwise)
    * T2 Sensor (1 if rewarded, -1 if punished, 0 otherwise)
  * Traits:
    * T1: NAND
    * T2: NOT

### Data Collection
  * Population slicing: 50K Updates

### Naming scheme:
  * <config file type>___<treatment arg>_<arg val>__<treatment arg>_<arg val>.... .cfg
    * arg vals may be a list: <arg val>_<arg val>..._<arg val>
  * For 2 task, 4 environment files:
    * Order that tasks are listed indicates which task has a long cycle and which has short cycle
      * <longest cycle>_..._<shortest cycle>
