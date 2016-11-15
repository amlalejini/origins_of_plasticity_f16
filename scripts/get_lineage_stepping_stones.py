"""
This script pulls out information on the evolutionary stepping stones for plasticity:
    * Does unconditional expression precede conditional expression?
    * Does sub-optimal plasticity precede optimal plasticity?

Input:
    * settings file (.json)
    * fdom_lineage file (.csv)
Details Output:
    * treatment, rep,
        - (x) any_plasticity -- does this lineage have any plastic organisms?
        - (x) any_adaptive_plasticity -- does this lineage have any adaptively plastic organisms?
            -- ( ) Is adaptive plasticity lost? How many times? (low change rate = hard to maintain plasticity?)
        - (x) any_optimality -- does this lineage have any optimal organisms?
        - (x) update_plastic -- what update is the first plastic organism on this lineage?
        - (x) update_adaptive_plastic -- what update is the first adaptively plastic organism on this lineage?
        - (x) update_optimal -- what update is the first optimal organism on this lineage?
        - (x) nand_uncon_before_nand_con -- is NAND expressed unconditionally before being expressed conditionally?
            -- ( ) time between nand uncon and nand con
            -- should I also look at *adaptive* conditional expression? Or is just conditional in general fine?
        - (x) not_uncon_before_not_con -- is NOT expressed unconditionally before being expressed conditionally?
            -- ( ) time between not uncon and not con
        - (x) any_uncon_before_con -- is any form of conditional expression preceded by any form of unconditional expression?
            -- ( ) time between uncon and con
        - (x) nand_subopt_before_nand_opt -- is nand expressed suboptimally before being expressed optimally?
            -- ( ) time between nand subopt before nand opt
        - (x) not_subopt_before_not_opt -- is not expressed suboptimally before being expressed optimally?
            -- ( ) time between not subopt before not opt
        - (x) any_subopt_before_opt -- is completely optimal expression preceded by any form of suboptimal expression?
            -- ( ) time between subopt before opt

Overview Output:
    * treatment
        - (x) total_reps -- How many replicates of this treatment are there?
        - (x) any_plasticity_cnt -- How many reps have any plasticity?
        - (x) any_adaptive_plasticity_cnt -- How many reps have any adaptive plasticity?
        - (x) any_optimality_cnt -- How many reps have an optimal organism in them?
        - (x) nand_uncon_before_nand_con_cnt -- In how many reps is nand expressed unconditionally before it is expressed conditionally?
        - (x) nand_con_cnt -- In how many reps is nand expressed conditionally at some point?
        - (x) not_uncon_before_not_con_cnt -- In how many reps is not expressed unconditionally before it is expressed conditionally?
        - (x) not_con_cnt -- In how many reps is not expressed conditionally at some point?
        - (x) any_uncon_before_con_cnt -- In how many reps is something expressed unconditionally before anything is expressed conditionally?
        - (x) con_cnt -- In how many reps is something expressed conditionally?
        - (x) nand_subopt_plastic_before_nand_opt_plastic_cnt -- In how many reps is nand expressed suboptimally plastic prior to being expressed optimally plastic?
        - (x) nand_opt_plastic_cnt -- In how many reps is nand expressed optimally plastic?
        - (x) not_subopt_plastic_before_not_opt_plastic_cnt -- In how many reps is not expressed suboptimally plastic prior to being expressed optimally plastic?
        - (x) not_opt_plastic_cnt -- In how many reps is not expressed optimally plastic?
        - (x) any_subopt_plastic_before_opt_plastic_cnt -- In how many reps is *anything* expressed sub-optimally plastic prior to being expressed optimally plastic?
        - (x) any_opt_plastic_cnt -- In how many reps is there optimal plasticity?
"""

import json, sys, os, datetime
from utilities.utilities import *

SEQ_DELIMITER = "*"
ENV_DELIMITER = "|"

def main():
    settings = None
    try:
        settings_fn = sys.argv[1]
        with open(settings_fn) as fp:
            settings = json.load(fp)
    except:
        print "Failed to load specified settings file."
        exit(-1)

    # Some useful functions for accessing information in settings.
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

    # Pull out some locations of interest.
    experiment_base_loc = settings["experiment_base_location"]
    experiment_analysis_loc = os.path.join(experiment_base_loc, "analysis")
    experiment_processed_loc = os.path.join(experiment_base_loc, "processed")

    # If it doesn't already exist, make the processed data location.
    mkdir_p(experiment_processed_loc)

    content = None
    try:
        lineage_results_fpath = sys.argv[2]
        with open(lineage_results_fpath) as fp:
            content = fp.read().strip("\n")
    except:
        print "Failed to load specified data file."
        exit(-1)

    content = content.split("\n")
    header = [attr for attr in content[0].split(",")]
    content = content[1:]
    lineage_data = [{header[i]:line.split(",")[i] for i in range(0, len(header))} for line in content]
    results_by_treatment = {}
    for rep_data in lineage_data:
        treatment = rep_data["treatment"]
        # Add new dictionary to results_by_treatment for this treatment if we haven't seen it before.
        if treatment not in results_by_treatment:
            results_by_treatment[treatment] = [] # Will be sorted on rep id

        environments = GetEnvironments(treatment)
        tasks = list(set([task for env in environments for task in GetEnvironmentTasks(env)]))
        rep_id = rep_data["replicate"]
        print "Processing: " + str(treatment) + " - " + str(rep_id)
        print "Environments: " + str(environments)
        print "Tasks: " + str(tasks)
        max_phen_score = rep_data["max_phenotype_score"]

        env_key = rep_data["environment_key"].split(ENV_DELIMITER)
        lineage_sig_seq = rep_data["lineage_phenotype_verbose_signature_sequence_abbrev"].split(SEQ_DELIMITER)
        lineage_phen_by_env_seq = [{env_key[ei]:phen.split(ENV_DELIMITER)[ei] for ei in range(0, len(env_key))} for phen in lineage_sig_seq]

        lineage_score_seq = rep_data["lineage_phenotype_score_sequence_abbrev"].split(SEQ_DELIMITER)
        lineage_is_plastic_seq = rep_data["lineage_is_plastic_sequence_abbrev"].split(SEQ_DELIMITER)
        lineage_is_optimal_seq = rep_data["lineage_is_optimal_sequence_abbrev"].split(SEQ_DELIMITER)
        lineage_start_update_seq = rep_data["lineage_phenotype_start_updates_abbrev"].split(SEQ_DELIMITER)
        lineage_duration_seq = rep_data["lineage_phenotype_duration_updates_abbrev"].split(SEQ_DELIMITER)
        lineage_len = len(lineage_sig_seq)

        # Does this lineage have *any* plasticity in it? And, what update is the first instance of adaptive plasticity?
        any_plasticity = False
        update_plastic = -1
        # Does this lineage have *any* adaptive plasticity in it? And, what update is the first instance of adaptive plasticity?
        any_adaptive_plasticity = False
        update_adaptive_plastic = -1
        # Does this lineage have *any* optimal phenotypes in it? And what update is the first instance of optimal plasticity?
        any_optimal = False
        update_optimal = -1
        # Is NAND expressed unconditionally before being expressed conditionally?
        first_nand_uncon_update = -1
        first_nand_con_update = -1
        # Is NOT expressed unconditionally before being expressed conditionally?
        first_not_uncon_update = -1
        first_not_con_update = -1
        # Is anything expressed unconditionally before anything is expressed conditionally?
        first_uncon_update = -1
        first_con_update = -1
        # Is nand expressed sub-optimally plastically before it is expressed optimally plastic?
        first_suboptimal_nand_plasticity = -1
        first_optimal_nand_plasticity = -1
        # Is not expressed sub-optimally plastically before it is expressed optimally plastic?
        first_suboptimal_not_plasticity = -1
        first_optimal_not_plasticity = -1
        # Does some form of sub-optimal plasticity precede fully optimal expression (across all tasks?)
        first_suboptimal_plasticity = -1
        first_optimal_plasticity = -1
        for i in range(0, lineage_len):
            # Any adaptive plasticity and when?
            any_adaptive_plasticity = any_adaptive_plasticity or (lineage_is_plastic_seq[i] == "1" and int(lineage_score_seq[i]) > 0)
            if any_adaptive_plasticity and update_adaptive_plastic == -1:
                update_adaptive_plastic = int(lineage_start_update_seq[i]) # What update is the first instance of adaptive plasticity?
            # Any plasticity and when?
            any_plasticity = any_plasticity or (lineage_is_plastic_seq[i] == "1")
            if any_plasticity and update_plastic == -1:
                update_plastic = int(lineage_start_update_seq[i]) # What update is the first instance of plasticity?
            # Any optimal phenotype?
            any_optimal = any_optimal or (lineage_is_optimal_seq[i] == "1")
            if any_optimal and update_optimal == -1:
                update_optimal = int(lineage_start_update_seq[i])

            phenotype_by_env = lineage_phen_by_env_seq[i]
            decoded_phenotype_by_env = {env:{GetEnvironmentTasks(env)[ti]:phenotype_by_env[env][ti] for ti in range(0, len(GetEnvironmentTasks(env)))} for env in phenotype_by_env}
            env_pressures = {env:{thing.split("_")[0]:set(thing.split("_")[1:]) for thing in env.split("__")} for env in phenotype_by_env}
            ## Treatment-specific questions: ##
            if "NAND" in tasks:
                # Is NAND expressed unconditionally before being expressed conditionally?
                # * Find first instance of unconditional nand expression.
                # * Find first instance of conditional nand expression
                nand_expression = {env:decoded_phenotype_by_env[env]["NAND"] for env in decoded_phenotype_by_env}
                nand_unconditional = len(set(nand_expression.values())) == 1 and "1" in set(nand_expression.values())
                nand_conditional = len(set(nand_expression.values())) > 1
                if nand_unconditional and first_nand_uncon_update == -1:
                    first_nand_uncon_update = int(lineage_start_update_seq[i])
                if nand_conditional and first_nand_con_update == -1:
                    first_nand_con_update = int(lineage_start_update_seq[i])
                # Is NAND expressed sub-optimally plastic before being expressed optimally plastic?
                # * Find max NAND score: for each environment, how many times does nand show up as rewarded?
                max_nand_score = 0
                for env in env_pressures:
                    pressures = env_pressures[env]
                    if "REWARD" in pressures:
                        if "NAND" in pressures["REWARD"]: max_nand_score += 1
                # * Calc achieved NAND score
                nand_score = 0
                for env in env_pressures:
                    pressures = env_pressures[env]
                    rewarded_tasks = set([]) if not "REWARD" in pressures else pressures["REWARD"]
                    punished_tasks = set([]) if not "PUNISH" in pressures else pressures["PUNISH"]
                    if "NAND" in rewarded_tasks and nand_expression[env] == "1":
                        nand_score += 1
                    elif "NAND" in punished_tasks and nand_expression[env] == "1":
                        nand_score -= 1
                # * sub-optimal plasticity: nand_is_plastic and achieved nand score < max nand score
                if nand_conditional and nand_score < max_nand_score and first_suboptimal_nand_plasticity == -1:
                    first_suboptimal_nand_plasticity = int(lineage_start_update_seq[i])
                # * optimal plasticity: is_plastic and achieved nand score == max nand score
                if nand_conditional and nand_score == max_nand_score and first_optimal_nand_plasticity == -1:
                    first_optimal_nand_plasticity = int(lineage_start_update_seq[i])

            if "NOT" in tasks:
                not_expression = {env:decoded_phenotype_by_env[env]["NOT"] for env in decoded_phenotype_by_env}
                not_unconditional = len(set(not_expression.values())) == 1 and "1" in set(not_expression.values())
                not_conditional = len(set(not_expression.values())) > 1
                if not_unconditional and first_not_uncon_update == -1:
                    first_not_uncon_update = int(lineage_start_update_seq[i])
                if not_conditional and first_not_con_update == -1:
                    first_not_con_update = int(lineage_start_update_seq[i])
                # Is NOT expressed sub-optimally plastic before being expressed optimally plastic?
                # * Find max NOT score
                max_not_score = 0
                for env in env_pressures:
                    pressures = env_pressures[env]
                    if "REWARD" in pressures:
                        if "NOT" in pressures["REWARD"]: max_not_score += 1
                # * Calc achieved NOT score
                not_score = 0
                for env in env_pressures:
                    pressures = env_pressures[env]
                    rewarded_tasks = set([]) if not "REWARD" in pressures else pressures["REWARD"]
                    punished_tasks = set([]) if not "PUNISH" in pressures else pressures["PUNISH"]
                    if "NOT" in rewarded_tasks and not_expression[env] == "1":
                        not_score += 1
                    elif "NOT" in punished_tasks and not_expression[env] == "1":
                        not_score -= 1
                # * sub-optimal plasticity: is_plastic and achieved not score < max not score
                if not_conditional and not_score < max_not_score and first_suboptimal_not_plasticity == -1:
                    first_suboptimal_not_plasticity = int(lineage_start_update_seq[i])
                # * optimal plasticity: is_plastic and achieved not score == max not score
                if not_conditional and not_score == max_not_score and first_optimal_not_plasticity == -1:
                    first_optimal_not_plasticity = int(lineage_start_update_seq[i])

            # When is the first unconditional expression of anything?
            if first_uncon_update == -1 and first_not_uncon_update != -1:
                first_uncon_update = first_not_uncon_update
            elif first_uncon_update == -1 and first_nand_uncon_update != -1:
                first_uncon_update = first_nand_uncon_update
            # When is the first conditional expression of anything?
            if first_con_update == -1 and first_not_con_update != -1:
                first_con_update = first_not_con_update
            elif first_con_update == -1 and first_nand_con_update != -1:
                first_con_update = first_nand_con_update
            # Does any sub-optimal plasticity precede optimal plasticity?
            # * sub-optimal plasticity: is_plastic and not is_optimal
            if first_suboptimal_plasticity == -1 and lineage_is_plastic_seq[i] == "1" and lineage_is_optimal_seq[i] != "1":
                first_suboptimal_plasticity = int(lineage_start_update_seq[i])
            # * optimal plasticity: is_plastic and is_optimal
            if first_optimal_plasticity == -1 and lineage_is_plastic_seq[i] == "1" and lineage_is_optimal_seq[i] == "1":
                first_optimal_plasticity = int(lineage_start_update_seq[i])


        # Is nand expressed unconditionally before being expressed conditionally?
        nand_uncon_before_nand_con = None
        if first_nand_con_update != -1:
            # Nand is conditionally expressed at some point
            if first_nand_uncon_update == -1 or first_nand_con_update < first_nand_uncon_update:
                # Conditional expression came first.
                nand_uncon_before_nand_con = False
            elif first_nand_con_update > first_nand_uncon_update:
                # Unconditional expression came first.
                nand_uncon_before_nand_con = True
            else:
                print "This is bad, and we shouldn't be here."
                exit(-1)

        # Is not expressed unconditionally before being expressed conditionally?
        not_uncon_before_not_con = None
        if first_not_con_update != -1:
            # Not is conditionally expressed at some point.
            if first_not_uncon_update == -1 or first_not_con_update < first_not_uncon_update:
                # Conditional expression came first.
                not_uncon_before_not_con = False
            elif first_not_con_update > first_not_uncon_update:
                # Unconditional expression came first.
                not_uncon_before_not_con = True
            else:
                print "This is bad, and we shouldn't be here."
                exit(-2)

        # Is anything expressed unconditionally before anything is expressed conditionally?
        any_uncon_before_any_con = None
        if first_con_update != -1:
            # Something is conditionally expressed at some point.
            if first_uncon_update == -1 or first_con_update < first_uncon_update:
                # Conditional expression came first.
                any_uncon_before_any_con = False
            elif first_con_update > first_uncon_update:
                # Unconditional expression came first.
                any_uncon_before_any_con = True
            else:
                print "This is bad, and we shouldn't be here."
                exit(-3)

        # Is nand expressed sub-optimally plastic before being expressed optimally plastic?
        nand_subopt_before_nand_opt = None
        if first_optimal_nand_plasticity != -1:
            # Something is optimally plastic about nand at some point.
            if first_suboptimal_nand_plasticity == -1 or first_optimal_nand_plasticity < first_suboptimal_nand_plasticity:
                nand_subopt_before_nand_opt = False
            elif first_optimal_nand_plasticity > first_suboptimal_nand_plasticity:
                nand_subopt_before_nand_opt = True
            else:
                print "This is bad, and we shouldn't be here."
                exit(-4)
        # Is not expressed sub-optimally plastic before being expressed optimally plastic?
        not_subopt_before_not_opt = None
        if first_optimal_not_plasticity != -1:
            # Something is optimally plastic about not at some point.
            if first_suboptimal_not_plasticity == -1 or first_optimal_not_plasticity < first_suboptimal_not_plasticity:
                not_subopt_before_not_opt = False
            elif first_optimal_not_plasticity > first_suboptimal_not_plasticity:
                not_subopt_before_not_opt = True
            else:
                print "This is bad, and we shouldn't be here."
                exit(-5)
        # Does sub-optimal plasticity (of any form) precede optimal plasticity?
        any_subopt_plasticity_before_opt_plasticity = None
        if first_optimal_plasticity != -1:
            # Something is optimally plastic at some point.
            if first_suboptimal_plasticity == -1 or first_optimal_plasticity < first_suboptimal_plasticity:
                any_subopt_plasticity_before_opt_plasticity = False
            elif first_suboptimal_plasticity < first_optimal_plasticity:
                any_subopt_plasticity_before_opt_plasticity = True
            else:
                print "This is bad, and we shouldn't be here."
                exit(-6)

        # Is the extant organism plastic?
        final_plastic = bool(int(rep_data["final_is_plastic"]))
        # Is the extant organism optimal?
        final_optimal = bool(int(rep_data["final_is_optimal"]))

        results_by_treatment[treatment].append({"treatment": treatment,
                                                "replicate": rep_id,
                                                "any_plasticity": any_plasticity,
                                                "any_adaptive_plasticity": any_adaptive_plasticity,
                                                "any_optimal": any_optimal,
                                                "first_plastic_update": update_plastic,
                                                "first_adaptively_plastic_update": update_adaptive_plastic,
                                                "first_optimal_update": update_optimal,
                                                "nand_uncon_before_nand_con": nand_uncon_before_nand_con,
                                                "not_uncon_before_not_con": not_uncon_before_not_con,
                                                "any_uncon_before_any_con": any_uncon_before_any_con,
                                                "nand_subopt_plastic_before_nand_opt_plastic": nand_subopt_before_nand_opt,
                                                "not_subopt_plastic_before_not_opt_plastic": not_subopt_before_not_opt,
                                                "any_subopt_plastic_before_opt_plastic": any_subopt_plasticity_before_opt_plasticity,
                                                "final_plastic": final_plastic,
                                                "final_optimal": final_optimal
                                                })


    ### Generate overview stats ###
    overview_header = ["treatment",
                       "total_reps",
                       "any_plasticity_cnt",
                       "any_adaptive_plasticity_cnt",
                       "any_optimality_cnt",
                       "nand_uncon_before_nand_con_cnt",
                       "nand_con_cnt",
                       "not_uncon_before_not_con_cnt",
                       "not_con_cnt",
                       "any_uncon_before_con_cnt",
                       "con_cnt",
                       "nand_subopt_plastic_before_nand_opt_plastic_cnt",
                       "nand_opt_plastic_cnt",
                       "not_subopt_plastic_before_not_opt_plastic_cnt",
                       "not_opt_plastic_cnt",
                       "any_subopt_plastic_before_opt_plastic_cnt",
                       "any_opt_plastic_cnt"]
    overview_fpath = os.path.join(experiment_processed_loc, "lineage_stepping_stones_overview__%s.csv" % str(datetime.datetime.now()).split(" ")[0])
    with open(overview_fpath, "w") as fp:
        fp.write(",".join(overview_header) + "\n")
    treatment_overview_content = ""
    ### Write details out to csv file ###
    details_header = ["treatment",
                      "replicate",
                      "any_plasticity",
                      "any_adaptive_plasticity",
                      "any_optimal",
                      "first_plastic_update",
                      "first_adaptively_plastic_update",
                      "first_optimal_update",
                      "nand_uncon_before_nand_con",
                      "not_uncon_before_not_con",
                      "any_uncon_before_any_con",
                      "nand_subopt_plastic_before_nand_opt_plastic",
                      "not_subopt_plastic_before_not_opt_plastic",
                      "any_subopt_plastic_before_opt_plastic",
                      "final_plastic",
                      "final_optimal"]
    details_fpath = os.path.join(experiment_processed_loc, "lineage_stepping_stones_details__%s.csv" % str(datetime.datetime.now()).split(" ")[0])
    with open(details_fpath, "w") as fp:
        fp.write(",".join(details_header) + "\n")

    for treatment in results_by_treatment:
        treatment_details_content = ""

        total_reps = 0
        any_plasticity_cnt = 0
        any_adaptive_plasticity_cnt = 0
        any_optimality_cnt = 0

        nand_uncon_before_nand_con_cnt = 0
        nand_con_cnt = 0
        not_uncon_before_not_con_cnt = 0
        not_con_cnt = 0
        any_uncon_before_con_cnt = 0
        con_cnt = 0

        nand_subopt_plastic_before_nand_opt_plastic_cnt = 0
        nand_opt_plastic_cnt = 0
        not_subopt_plastic_before_not_opt_plastic_cnt = 0
        not_opt_plastic_cnt = 0
        any_subopt_plastic_before_opt_plastic_cnt = 0
        any_opt_plastic_cnt = 0

        for rep_results in results_by_treatment[treatment]:
            # make details line for this replicate.
            treatment_details_content += ",".join([str(rep_results[attr]) for attr in details_header]) + "\n"
            total_reps += 1
            if rep_results["any_plasticity"]:
                any_plasticity_cnt += 1
            if rep_results["any_adaptive_plasticity"]:
                any_adaptive_plasticity_cnt += 1
            if rep_results["any_optimal"]:
                any_optimality_cnt += 1
            # Uncon before con
            if rep_results["nand_uncon_before_nand_con"] != None:
                nand_con_cnt += 1
                if rep_results["nand_uncon_before_nand_con"]: nand_uncon_before_nand_con_cnt += 1
            if rep_results["not_uncon_before_not_con"] != None:
                not_con_cnt += 1
                if rep_results["not_uncon_before_not_con"]: not_uncon_before_not_con_cnt += 1
            if rep_results["any_uncon_before_any_con"] != None:
                con_cnt += 1
                if rep_results["any_uncon_before_any_con"]: any_uncon_before_con_cnt += 1
            # Sub-opt before opt
            if rep_results["nand_subopt_plastic_before_nand_opt_plastic"] != None:
                nand_opt_plastic_cnt += 1
                if rep_results["nand_subopt_plastic_before_nand_opt_plastic"]: nand_subopt_plastic_before_nand_opt_plastic_cnt += 1
            if rep_results["not_subopt_plastic_before_not_opt_plastic"] != None:
                not_opt_plastic_cnt += 1
                if rep_results["not_subopt_plastic_before_not_opt_plastic"]: not_subopt_plastic_before_not_opt_plastic_cnt += 1
            if rep_results["any_subopt_plastic_before_opt_plastic"] != None:
                any_opt_plastic_cnt += 1
                if rep_results["any_subopt_plastic_before_opt_plastic"]: any_subopt_plastic_before_opt_plastic_cnt += 1
        # Update details file.
        with open(details_fpath, "a") as fp:
            fp.write(treatment_details_content)

        treatment_overview_content += ",".join(map(str, [
                                      treatment,
                                      total_reps,
                                      any_plasticity_cnt,
                                      any_adaptive_plasticity_cnt,
                                      any_optimality_cnt,
                                      nand_uncon_before_nand_con_cnt,
                                      nand_con_cnt,
                                      not_uncon_before_not_con_cnt,
                                      not_con_cnt,
                                      any_uncon_before_con_cnt,
                                      con_cnt,
                                      nand_subopt_plastic_before_nand_opt_plastic_cnt,
                                      nand_opt_plastic_cnt,
                                      not_subopt_plastic_before_not_opt_plastic_cnt,
                                      not_opt_plastic_cnt,
                                      any_subopt_plastic_before_opt_plastic_cnt,
                                      any_opt_plastic_cnt])) + "\n"
    with open(overview_fpath, "a") as fp:
        fp.write(treatment_overview_content)
    print "Done"


if __name__ == "__main__":
    main()
