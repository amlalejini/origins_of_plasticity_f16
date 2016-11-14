"""
This scirpt pulls out relevant results on the lineages of final dominant organisms.

Things to output:
    * Fdom Plastic?
    * Phenotype Signature key
    * Env Signature key
    * Score sequence
    * Phenotype Signature sequence
Data file:
    * treatment, rep, final_is_plastic, final_is_optimal, final_phenotype_score, final_phenotype_signature, phenotype_signature_sequence, ...keys...
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

    # Some useful functions for extracting information from settings.
    def GetEnvironments(treatment):
        """
        Given a treatment, look in settings for relevant environments.
        """
        return settings["treatment_settings"][settings["treatments"][treatment]["settings"]]["experienced_environments"]
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
        something = {thing.split("_")[0]:set(thing.split("_")[1:]) for thing in environment.split("__")}
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

    # Pull out some locations of interest.
    experiment_base_loc = settings["experiment_base_location"]
    experiment_analysis_loc = os.path.join(experiment_base_loc, "analysis")
    experiment_processed_loc = os.path.join(experiment_base_loc, "processed")

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
    # * treatment, rep, final_is_plastic, final_is_optimal, final_phenotype_score, final_phenotype_signature, phenotype_signature_sequence, ...keys...
    header_attrs = ["treatment", "replicate", "final_is_plastic", "final_is_optimal", "final_phenotype_score", "max_phenotype_score", "enviroment_key", "lineage_phenotype_signature_sequence_abbrev", "lineage_phenotype_verbose_signature_sequence_abbrev", "lineage_phenotype_score_sequence_abbrev", "lineage_is_plastic_sequence_abbrev", "lineage_is_optimal_sequence_abbrev", "lineage_phenotype_start_updates_abbrev", "lineage_phenotype_duration_updates_abbrev"]
    data_fpath = os.path.join(experiment_processed_loc, "processed_fdom_lineage__%s.csv" % str(datetime.datetime.now()).split(" ")[0])
    with open(data_fpath, "w") as fp:
        fp.write(",".join(header_attrs) + "\n")

    for treatment in treatments_to_process:
        content = "" # Variable to keep track of data content for all reps of this treatment.
        print "\n\nProcessing %s" % treatment
        # Extract info x environment.
        # What environments are relevant?
        environments = GetEnvironments(treatment)
        max_phenotype_scores_by_env = {e:CalcMaxPhenotypeScore(e) for e in environments}
        # Create an environment ordering for this treatment.
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
            lineage_details_by_env = {} # Will store lineage details by environment
            lineage_len = set()
            for env in environments:
                # Extract information from data files.
                detail_fn = "fdom_lineage_details.dat"
                detail_fpath = os.path.join(experiment_analysis_loc, rep_fn, "final_dominant", "ENV__" + env, detail_fn)
                lineage_details = ParseDetailFile(detail_fpath)
                lineage_details_by_env[env] = lineage_details
                lineage_len.add(len(lineage_details))
            assert len(lineage_len) == 1
            lineage_len = next(iter(lineage_len)) # Set length should be == 1
            # Re-structure how lineage is stored.
            lineage_analysis_details = [{env:lineage_details_by_env[env][i] for env in environments} for i in range(0, lineage_len)]
            # Do some more processing on each ancestor x environment
            # Things to collect per lineage:
            # * Full sequences
            full_signature_sequence = [None for _ in range(0, lineage_len)]
            full_signature_verb_sequence = [None for _ in range(0, lineage_len)]
            full_total_score_sequence = [None for _ in range(0, lineage_len)]
            full_start_updates_sequence = [None for _ in range(0, lineage_len)]
            full_duration_updates_sequence = [None for _ in range(0, lineage_len)]
            full_is_plastic_sequence = [None for _ in range(0, lineage_len)]
            full_is_optimal_sequence = [None for _ in range(0, lineage_len)]
            # * Abbreviated sequences (only add new element on change)
            abbrev_signature_sequence = []
            abbrev_signature_verb_sequence = []
            abbrev_total_score_sequence = []
            abbrev_start_updates_sequence = []
            abbrev_duration_updates_sequence = []
            abbrev_is_plastic_sequence = []
            abbrev_is_optimal_sequence = []
            # Book keeping for abbreviated sequence tracking:
            cur_phen_state = None
            cur_signature_verb = None
            cur_is_plastic = None
            cur_is_optimal = None
            cur_total_score = None
            cur_start = 0
            cur_duration = 0
            for li in range(0, lineage_len):
                ancestor_by_env = lineage_analysis_details[li]
                phenotype_by_env = {}
                # Get phenotype x environment.
                for env in environments:
                    # Grab the environment tasks (in order for encoding)
                    env_tasks = GetEnvironmentTasks(env)
                    # Characterize phenotype of ancestor in environment.
                    phenotype = {t:ancestor_by_env[env][t.lower()] for t in env_tasks}
                    phenotype["encoding"] = "".join([ancestor_by_env[env][t.lower()] for t in env_tasks])
                    phenotype_by_env[env] = phenotype

                # Get info for this ancestor organism.
                # * Start update
                start_update = set([ancestor_by_env[env]["update born"] for env in environments])
                assert len(start_update) == 1
                start_update = int(next(iter(start_update)))
                if start_update == -1: start_update = 0 # Correct that ancestor starts at update -1
                # * Duration in updates (next start - this start) or (total updates - this start)
                next_start_update = None
                if li + 1 < lineage_len:
                    next_start_update = set([lineage_analysis_details[li + 1][env]["update born"] for env in environments])
                    assert len(next_start_update) == 1
                    next_start_update = int(next(iter(next_start_update)))
                else:
                    next_start_update = settings["final_update"]
                duration_update = next_start_update - start_update
                # * Phenotype signature
                signature = "".join([phenotype_by_env[env]["encoding"] for env in environment_order])
                # * Verbose phenotype by environment.
                signature_verb = "|".join([phenotype_by_env[env]["encoding"] for env in environment_order])
                # * Score by env
                score_by_env = GetPhenotypeScore(phenotype_by_env)
                # * is_plastic
                is_plastic = len(set([phenotype_by_env[env]["encoding"] for env in environments])) > 1
                # * max score?
                max_score = sum(max_phenotype_scores_by_env.values())
                # * total score?
                total_score = sum(score_by_env.values())
                # * is_optimal?
                is_optimal = total_score == max_score

                # Record relevant values.
                full_signature_sequence[li] = signature
                full_signature_verb_sequence[li] = signature_verb
                full_total_score_sequence[li] = total_score
                full_start_updates_sequence[li] = start_update
                full_duration_updates_sequence[li] = duration_update
                full_is_plastic_sequence[li] = int(is_plastic)
                full_is_optimal_sequence[li] = int(is_optimal)

                # Update abbrev. state, record abbreviated values if necessary.
                if cur_phen_state == None:
                    # First phenotype.
                    cur_phen_state = signature
                    cur_start = start_update
                    cur_duration = duration_update
                    cur_signature_verb = signature_verb
                    cur_is_plastic = int(is_plastic)
                    cur_is_optimal = int(is_optimal)
                    cur_total_score = total_score
                elif (li == lineage_len - 1):
                    # At end of lineage, unconditionally clip state.
                    # Clip previous state.
                    abbrev_signature_sequence.append(cur_phen_state)
                    abbrev_signature_verb_sequence.append(cur_signature_verb)
                    abbrev_total_score_sequence.append(cur_total_score)
                    abbrev_start_updates_sequence.append(cur_start)
                    abbrev_duration_updates_sequence.append(cur_duration)
                    abbrev_is_plastic_sequence.append(cur_is_plastic)
                    abbrev_is_optimal_sequence.append(cur_is_optimal)
                elif cur_phen_state == signature:
                    # This phenotype collapses.
                    cur_duration += duration_update
                else:
                    # This is a new phenotype and not at the end of the lineage.
                    # Clip previous state.
                    abbrev_signature_sequence.append(cur_phen_state)
                    abbrev_signature_verb_sequence.append(cur_signature_verb)
                    abbrev_total_score_sequence.append(cur_total_score)
                    abbrev_start_updates_sequence.append(cur_start)
                    abbrev_duration_updates_sequence.append(cur_duration)
                    abbrev_is_plastic_sequence.append(cur_is_plastic)
                    abbrev_is_optimal_sequence.append(cur_is_optimal)
                    # Reset book-keeping variables.
                    cur_phen_state = signature
                    cur_start = start_update
                    cur_duration = duration_update
                    cur_signature_verb = signature_verb
                    cur_is_plastic = int(is_plastic)
                    cur_is_optimal = int(is_optimal)
                    cur_total_score = total_score

            env_order = "|".join(environment_order)
            max_score = sum(max_phenotype_scores_by_env.values())
            final_is_plastic = abbrev_is_plastic_sequence[-1]
            final_is_optimal = abbrev_is_optimal_sequence[-1]
            final_score = abbrev_total_score_sequence[-1]
            lineage_phenotype_signature_sequence_abbrev = "*".join(abbrev_signature_sequence)
            lineage_phenotype_verbose_signature_sequence_abbrev = "*".join(abbrev_signature_verb_sequence)
            lineage_phenotype_score_sequence_abbrev = "*".join(map(str, abbrev_total_score_sequence))
            lineage_is_plastic_sequence_abbrev = "*".join(map(str, abbrev_is_plastic_sequence))
            lineage_is_optimal_sequence_abbrev = "*".join(map(str, abbrev_is_optimal_sequence))
            lineage_phenotype_start_updates_abbrev = "*".join(map(str, abbrev_start_updates_sequence))
            lineage_phenotype_duration_updates_abbrev = "*".join(map(str, abbrev_duration_updates_sequence))

            #       ["treatment", "replicate", "final_is_plastic", "final_is_optimal", "final_phenotype_score", "max_phenotype_score", "enviroment_key", "lineage_phenotype_signature_sequence_abbrev", "lineage_phenotype_verbose_signature_sequence_abbrev", "lineage_phenotype_score_sequence_abbrev", "lineage_is_plastic_sequence_abbrev", "lineage_is_optimal_sequence_abbrev", "lineage_phenotype_start_updates_abbrev", "lineage_phenotype_duration_updates_abbrev"]
            lineage_attrs = map(str,[treatment,    rep_id,      final_is_plastic,   final_is_optimal,   final_score,            max_score,              env_order,        lineage_phenotype_signature_sequence_abbrev,   lineage_phenotype_verbose_signature_sequence_abbrev,   lineage_phenotype_score_sequence_abbrev,   lineage_is_plastic_sequence_abbrev,   lineage_is_optimal_sequence_abbrev,   lineage_phenotype_start_updates_abbrev,   lineage_phenotype_duration_updates_abbrev])
            content += ",".join(lineage_attrs) + "\n"
        #### Add treatment content to file ####
        with open(data_fpath, "a") as fp:
            fp.write(content)
    print "Done"

if __name__ == "__main__":
    main()
