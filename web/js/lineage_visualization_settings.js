var DEFAULT_TREATMENT = "ntasks_2__envs_2__cr_800__tasks_NAND_NOT__mr_p0075";
var SEQ_DELIMITER = "*";
var ENV_DELIMITER = "|";
var data_fpath = "data/processed_fdom_lineage__2016-11-16.csv";
var exp_settings = {
  "experiment_base_location": "/mnt/home/lalejini/Data/oop_f16",
  "avida_analysis_scripts_location": "../avida_analysis_configs",
  "analysis_scripts_to_run": ["analyze_final_dominants.cfg"],
  "final_update": 200000,
  "run_avida_analyses_settings": {
    "treatments_to_analyze": [
      "ntasks_2__envs_4__cr_1__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_25__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_50__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_100__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_200__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_400__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_800__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_1600__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_25000__tasks_NAND_NOT__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_1__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_25__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_50__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_100__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_200__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_400__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_800__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_1600__tasks_NOT_NAND__style_2__mr_p0075",
      "ntasks_2__envs_4__cr_25000__tasks_NOT_NAND__style_2__mr_p0075"
    ]
  },
  "get_fdom_results_settings": {
    "treatments_to_process": "all"
  },
  "treatments": {
      "ntasks_1__envs_1__cr_0__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_1__tasks_NAND"},
      "ntasks_1__envs_1__cr_0__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_1__tasks_NOT"},
      "ntasks_2__envs_1__cr_0__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_1__tasks_NAND_NOT"},

      "ntasks_1__envs_2__cr_1__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_25__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_50__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_100__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_200__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_400__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_800__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_1600__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_25000__tasks_NAND__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NAND"},
      "ntasks_1__envs_2__cr_1__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_25__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_50__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_100__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_200__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_400__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_800__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_1600__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},
      "ntasks_1__envs_2__cr_25000__tasks_NOT__mr_p0075": {"settings": "ntasks_1__envs_2__tasks_NOT"},

      "ntasks_2__envs_2__cr_1__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_25__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_50__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_100__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_200__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_400__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_800__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_1600__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},
      "ntasks_2__envs_2__cr_25000__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_2__tasks_NAND_NOT"},

      "ntasks_2__envs_4__cr_1__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_25__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_50__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_100__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_200__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_400__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_800__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_1600__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_25000__tasks_NAND_NOT__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT"},
      "ntasks_2__envs_4__cr_1__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_25__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_50__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_100__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_200__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_400__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_800__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_1600__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},
      "ntasks_2__envs_4__cr_25000__tasks_NOT_NAND__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND"},

      "ntasks_2__envs_4__cr_1__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_25__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_50__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_100__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_200__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_400__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_800__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_1600__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_25000__tasks_NAND_NOT__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NAND_NOT__style_2"},
      "ntasks_2__envs_4__cr_1__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_25__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_50__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_100__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_200__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_400__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_800__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_1600__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"},
      "ntasks_2__envs_4__cr_25000__tasks_NOT_NAND__style_2__mr_p0075": {"settings": "ntasks_2__envs_4__tasks_NOT_NAND__style_2"}

    },
  "treatment_settings": {
    "ntasks_1__envs_1__tasks_NAND": {
      "test_environments": ["REWARD_NAND"],
      "experienced_environments": ["REWARD_NAND"],
      "replicates": [1, 50]
    },
    "ntasks_1__envs_1__tasks_NOT": {
      "test_environments": ["REWARD_NOT"],
      "experienced_environments": ["REWARD_NOT"],
      "replicates": [1, 20]
    },
    "ntasks_2__envs_1__tasks_NAND_NOT": {
      "test_environments": ["REWARD_NAND_NOT"],
      "experienced_environments": ["REWARD_NAND_NOT"],
      "replicates": [1, 50]
    },
    "ntasks_1__envs_2__tasks_NAND": {
      "test_environments": ["REWARD_NAND", "PUNISH_NAND"],
      "experienced_environments": ["REWARD_NAND", "PUNISH_NAND"],
      "replicates": [1, 50]
    },
    "ntasks_1__envs_2__tasks_NOT": {
      "test_environments": ["REWARD_NOT", "PUNISH_NOT"],
      "experienced_environments": ["REWARD_NOT", "PUNISH_NOT"],
      "replicates": [1, 20]
    },
    "ntasks_2__envs_2__tasks_NAND_NOT": {
      "test_environments": ["REWARD_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND"],
      "experienced_environments": ["REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND"],
      "replicates": [1, 50]
    },
    "ntasks_2__envs_4__tasks_NAND_NOT": {
      "test_environments": ["REWARD_NAND_NOT", "PUNISH_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND"],
      "experienced_environments": ["REWARD_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND", "PUNISH_NAND_NOT"],
      "replicates": [1, 50]
    },
    "ntasks_2__envs_4__tasks_NOT_NAND": {
      "test_environments": ["REWARD_NAND_NOT", "PUNISH_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND"],
      "experienced_environments": ["REWARD_NAND_NOT", "REWARD_NOT__PUNISH_NAND", "REWARD_NAND__PUNISH_NOT", "PUNISH_NAND_NOT"],
      "replicates": [1, 20]
    },
    "ntasks_2__envs_4__tasks_NAND_NOT__style_2": {
      "test_environments": ["REWARD_NAND_NOT", "PUNISH_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND"],
      "experienced_environments": ["REWARD_NAND__PUNISH_NOT", "REWARD_NAND_NOT", "REWARD_NOT__PUNISH_NAND", "PUNISH_NAND_NOT"],
      "replicates": [1, 50]
    },
    "ntasks_2__envs_4__tasks_NOT_NAND__style_2": {
      "test_environments": ["REWARD_NAND_NOT", "PUNISH_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "REWARD_NOT__PUNISH_NAND"],
      "experienced_environments": ["REWARD_NOT__PUNISH_NAND", "REWARD_NAND_NOT", "REWARD_NAND__PUNISH_NOT", "PUNISH_NAND_NOT"],
      "replicates": [1, 20]
    }
  },
  "phenotype_encoding": {
    "environment_order": ["PUNISH_NAND_NOT", "REWARD_NOT__PUNISH_NAND", "REWARD_NAND__PUNISH_NOT", "REWARD_NAND_NOT", "PUNISH_NAND", "REWARD_NAND", "PUNISH_NOT", "REWARD_NOT"],
    "trait_order": ["NAND", "NOT"]
  }
};
