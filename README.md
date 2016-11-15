# origins_of_plasticity_f16
Code/scripts for evolutionary origins of phenotypic plasticity experiments. Relevant for ALife Journal article.

## The Question
  * What are the stepping stones for the evolution of phenotypic plasticity?
    * Unconditional expression before conditional expression?
      - Unconditional expression -- express the same phenotype across all environments.
      - Conditional expression -- phenotypic expression is a function of the environment (in some way).
      - To test: Check lineages for cases where unconditional expression evolves before conditional expression.
    * Imprecise/sub-optimal plasticity before optimal plasticity?
      - Optimal plasticity -- phenotype is able to be regulated to perfectly match the environment.
      - Sub-optimal plasticity -- conditional expression, but not optimal.
        * Adaptive sub-optimal plasticity -- in total, plasticity results in better fitness than could have been achieved
          via unconditional trait expression.
        * Neutral sub-optimal plasticity -- in total, plasticity results in fitness equal to what could have been
          achieved via unconditional trait expression.
        * Mal-adaptive sub-optimal plasticity -- in total, plasticity results in fitness lower to what could have been
          achieved via unconditional trait expression.
        * In the above definitions, I have several options for fitness: actual fitness scores, simplified fitness
          calculation based only on traits expressed (reward: +1, punish: -1, nothing: 0)
        * Another note: I should be careful to also take into consideration the environmental background.
          Some things may look mal-adaptive as a whole, but at long env. cycle rates, only the current environment
          actually matters. Is there some way to quantify selective pressure given rate of environmental change? How much
          pressure on average is there for traits?
  * Do these conditions seem to be robust under a few different conditions?
  * Can we identify alternative strategies to phenotypic plasticity in changing environments?
  * When do these alternative strategies emerge?

## Hypotheses:
  * Stepping stones:
    * Unconditional trait expression prior to conditional trait expression.
    * Imprecise, but still beneficial plasticity before perfect plasticity (this is more useful for EC).
  * Alternative strategies:
    * Bet-hedging will emerge as a viable strategy as an alternative to phenotypic plasticity at certain cycle rates.
  * At higher environmental change rates, plastic machinery is harder to maintain because of relaxed selection
    pressure.

## Experimental Design


### Environments
  * Constant Environments (controls):
    * T1(+):      [NAND+]
    * T2(+):      [NOT+]
    * T1(+)T2(+): [NAND+NOT+]
  * Changing Environments
    * Traits: 1, Sensors: 1, Actuators: 1, Environments: 2
      * [T1(+), T1(-)]:                                   [NAND+, NAND-]
      * [T2(+), T2(-)]:                                   [NOT+, NOT-]

    * Traits: 2, Sensors: 2 (only need 1), Actuators: 2, Environments: 2
      * [T1(+)T2(-), T1(-)T2(+)]:                         [NAND+NOT-, NOT-NAND+]

    * Traits: 2, Sensors: 2, Actuators: 2, Environments: 4
      * STYLE -- 1
        * [T1(+)T2(+), T1(+)T2(-), T1(-)T2(+), T1(-)T2(-)]: [NAND+NOT+, NAND+NOT-, NAND-NOT+, NAND-NOT-]
        * [T2(+)T1(+), T2(+)T1(-), T2(-)T1(+), T2(-)T1(-)]: [NOT+NAND+, NOT+NAND-, NOT-NAND+, NOT-NAND-]
      * STYLE -- 2
        * [T1(+)T2(-), T1(+)T2(+), T1(-)T2(+), T1(-)T2(-)]
        * [T2(+)T1(-), T2(+)T2(+), T2(-)T1(+), T2(-)T1(-)]

  * Point Mutation Rates (Not yet)
    * Low: 0.0025
    * Baseline: 0.0075
    * High: 0.0125
  * Change Rates (in updates)
    * 1     (50  reps)
    * 25    (50  reps)
    * 50    (50  reps)
    * 100   (50  reps)
    * 200   (50  reps)
    * 400   (50  reps)
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
  * [config file type]___[treatment arg]_[arg val]__[treatment arg]_[arg val].... .cfg
    * arg vals may be a list: [arg val]_[arg val]..._[arg val]
  * For 2 task, 4 environment files:
    * Order that tasks are listed indicates which task has a long cycle and which has short cycle
      * [longest cycle]_..._[shortest cycle]

## Avida Analyses
  * Final Dominant:
    * Details on lineage
    * Trace on final dominant

## Post Avida Analyses
  * Phenotype encoding (Bit ordering top to bottom left to right):
    * Ordering:
      * T1-T2-:T1, T1-T2-:T2
      * T1-T2+:T1, T1-T2+:T2
      * T1+T2-:T1, T1+T2-:T2
      * T1+T2+:T1, T1+T2+:T2
      * T1-:T1
      * T1+:T1
      * T2-:T2
      * T2+:T2
    * 1: organism expressed trait in that environment.
    * 0: organism did not express trait in environment.

## Scripts
  * get_fdom_lineage_results.py
    * REQ: Analysis data from avida's analyze mode (analyze_final_dominants.cfg).
    * Given analysis data, extract fdom lineage results.
  * get_fdom_results.py
    * REQ: Analysis data from avida's analyze mode (analyze_final_dominants.cfg).
    * Given analysis data, extract fdom results.
  * get_lineage_stepping_stones.py
    * REQ: Output from get_fdom_lineage_results.py
    * Given output of fdom_lineage_results script, get information on stepping stones:
      * Unconditional plasticity before conditional plasticity?
      * Sub-optimal plasticity before optimal plasticity?
  * run_avida_analyses.py
    * REQ: Data from run.
    * Run Avida's analyze mode on data specified in settings.json

## To do list:
  * Visualize lineage phenotypes
  * Visualize: unconditional expression before conditional expression? (trait by trait)
  * Visualize: sub-optimal plasticity before optimal plasticity?
  * Characterize plasticity in genome (how tight is regulation with expression)?
    * Physical modularity; functional modularity
  * Alternative strategies? -- Dependency: visualization.
  * Seeds for division of labor experiments: how to set this up?
  * (2) More change rates between (1:25)? Yes: 3, 6, 12
    * Submit jobs
  * (3) Block stepping runs for env 2, traits 2
  * Visualize 'update_x' data as violin plots (should be interesting)
