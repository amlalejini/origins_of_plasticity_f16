"""
This script pulls out relevant results on the final dominant organisms (just the fdoms, not their lineages).

Things this will output:
 * By rep:
    * Plastic?
    * Phenotype
    * Generation length x environment
    * Fitness x environment
    * Tree depth

Data files:
    * treatment, rep, plastic, phenotype, generation length, fitness, tree depth
"""

import json, sys, os, datetime
from utilities.utilities import *
from utilities.parse_avida_output import *

def main():
    settings = None
    try:
        settings_fn = sys.argv[1]
        with open(settings_fn) as fp:
            settings = json.load(fp)
    except:
        print "Failed to load specified settings file."
        exit(-1)

    def GetEnvironments(treatment):
        """
        Given a treatment, look in settings for relevant environments.
        """
        return settings["treatment_settings"][settings["treatments"][treatment]["settings"]]["test_environments"]
    def GetEnvironmentTasks(environment):
        """
        Given an environment, return the tasks that were in the environment.
        Return tasks in rank order (provided by phenotype encoding settings).
        """
        tasks = set([])
        attrs = environment.split("__")
        for attr in attrs:
            tasks = tasks.union(set(attr.split("_")[1:]))
        return [task for task in settings["phenotype_encoding"]["trait_order"] if task in tasks]
    CalcMaxPhenotypeScoreMemo = {}
    def CalcMaxPhenotypeScore(environment):
        """
        Given an environment, return the max score possible in this environment.
        This function is memoized.
        """
        # Have we seen this input?
        if environment in CalcMaxPhenotypeScoreMemo:
            return CalcMaxPhenotypeScoreMemo[environment]
        # If not, calculate the max score:
        something = {thing.split("_")[0]:set(thing.split("_")[1:]) for thing in env.split("__")}
        if "REWARD" in something: return len(something["REWARD"])
        return 0
    def GetPhenotypeScore(phenotype_by_env):
        """
        Given phenotypes by environment, calculate a score using the following:
        If (trait is expressed):
            If (trait is punished):
                score += (-)1
            If (trait is rewarded):
                score += (+)1
        Else: score += 0
        """
        score_by_env = {env:0 for env in phenotype_by_env}
        for env in phenotype_by_env:
            something = {thing.split("_")[0]:set(thing.split("_")[1:]) for thing in env.split("__")}
            rewarded_tasks = set([]) if not "REWARD" in something else something["REWARD"]
            punished_tasks = set([]) if not "PUNISH" in something else something["PUNISH"]
            env_traits = GetEnvironmentTasks(env)
            for trait in env_traits:
                if phenotype_by_env[env][trait] == "1" and trait in rewarded_tasks:
                    score_by_env[env] += 1
                elif phenotype_by_env[env][trait] == "1" and trait in punished_tasks:
                    score_by_env[env] += -1
        return score_by_env


    # Pull out locations of interest.
    experiment_base_loc = settings["experiment_base_location"]
    experiment_analysis_loc = os.path.join(experiment_base_loc, "analysis")
    experiment_processed_loc = os.path.join(experiment_base_loc, "processed")
    experiment_configs_loc = os.path.join(experiment_base_loc, "configs")

    # If it doesn't already exist, make the processed data location.
    mkdir_p(experiment_processed_loc)

    # What treatments should I process?
    treatments_to_process = None
    if settings["get_fdom_results_settings"]["treatments_to_process"] == "all":
        treatments_to_process = settings["treatments"].keys()
    else:
        treatments_to_process = settings["get_fdom_results_settings"]["treatments_to_process"]
    # Collect available reps for each treatment to process.
    reps_by_treatment = {t:[r for r in os.listdir(experiment_analysis_loc) if t in r] for t in treatments_to_process}
    # Prepare data file.
    header_attrs = ["treatment", "replicate", "tree_depth", "is_plastic", "is_optimal", "phenotype_score", "max_phenotype_score", "phenotype_signature", "phenotype_by_env", "generation_length_by_env", "fitness_by_env", "env_order"]
    data_fpath = os.path.join(experiment_processed_loc, "processed_fdom__%s.csv" % str(datetime.datetime.now()).split(" ")[0])
    with open(data_fpath, "w") as fp:
        fp.write(",".join(header_attrs) + "\n")
    for treatment in treatments_to_process:
        content = "" # variable to keep track of data content for all reps of this treatment
        print "\n\nProcessing %s" % treatment
        # Extract info x environment
        # What environments are relevant?
        environments = GetEnvironments(treatment)
        print "  " + str(environments)
        # Create environment ordering for this treatment.
        environment_order = [env for env in settings["phenotype_encoding"]["environment_order"] if env in environments]
        print "  Environment order: " + str(environment_order)
        for rep in reps_by_treatment[treatment]:
            # Extract replicate ID and replicate attributes from rep name.
            rep_fn = rep
            rep = rep.split("_")
            rep_id = rep[-1]
            rep = "_".join(rep[:len(rep) - 1])
            rep = rep.split("__")
            rep_attrs = {ra.split("_")[0]:ra.split("_")[1:] for ra in rep}
            print "  processing rep: " + str(rep_id)
            # Some dictionaries of relevant information.
            details_by_env = {e:{} for e in environments}
            phenotype_by_env = {e:{} for e in environments}
            verbose_phenotype_by_env = {e: {} for e in environments}
            max_phenotype_scores_by_env = {e:None for e in environments}
            for env in environments:
                # Extract relevant information from data files.
                detail_fn = "fdom_details.dat"
                detail_fpath = os.path.join(experiment_analysis_loc, rep_fn, "final_dominant", "ENV__" + env, detail_fn)
                details = ParseDetailFile(detail_fpath)[0] # Only one organism per fdom detail file.
                details_by_env[env] = details
                env_tasks = GetEnvironmentTasks(env)
                phenotype_by_env[env] = {t:details[t.lower()] for t in env_tasks}
                phenotype_by_env[env]["encoding"] = "".join([details[t.lower()] for t in env_tasks])
                max_phenotype_scores_by_env[env] = CalcMaxPhenotypeScore(env)
            ##### THINGS TO RECORD: ####
            #   * Does the organism express different traits in different environments?
            # * Is plastic?
            is_plastic = len(set([phenotype_by_env[env]["encoding"] for env in environments])) > 1
            # * Score the phenotype
            #   * Is it optimal? How good is it?
            score_by_env = GetPhenotypeScore(phenotype_by_env)
            total_score = sum(score_by_env.values())
            max_score = sum(max_phenotype_scores_by_env.values())
            is_optimal = total_score == max_score
            # * Environment order key
            env_order = "|".join(environment_order)
            # * Phenotype
            # print "\tPhenotype x env: " + str(phenotype_by_env)
            full_encoding = "".join([phenotype_by_env[env]["encoding"] for env in environment_order])
            # * Verbose phenotype by environment
            verbose_phenotype_by_env = "|".join([phenotype_by_env[env]["encoding"] for env in environment_order])
            # * Gen len X Env
            genlen_by_env = "|".join([details_by_env[env]["gestation time"] for env in environment_order])
            # * Fitness X Env
            fitness_by_env = "|".join([details_by_env[env]["fitness"] for env in environment_order])
            # * Tree depth
            tree_depth = details_by_env[environments[0]]["tree depth"]
            #### Add rep info to content ####
            # header_attrs = ["treatment", "replicate", "tree_depth", "is_plastic", "is_optimal", "phenotype_score", "phenotype_signature", "phenotype_by_env", "generation_length_by_env", "fitness_by_env", "env_order"]
            fdom_output = map(str, [treatment, rep_id, tree_depth, int(is_plastic), int(is_optimal), total_score, max_score, full_encoding, verbose_phenotype_by_env, genlen_by_env, fitness_by_env, env_order])
            content += ",".join(fdom_output) + "\n"
        #### Add treatment content to data file ####
        with open(data_fpath, "a") as fp:
            fp.write(content)
    print "Done"


if __name__ == "__main__":
    main()
