#!/usr/bin/python2

"""
Manage the running of Avida's analyze mode given the settings specified in the given settings file.
"""

import json, os, subprocess, sys

def main():
    """
    Main script
    """
    settings_fn = None
    try:
        settings_fn = sys.argv[1]
    except:
        print "Failed to load specified settings file."
        exit(-1)
    # Load settings file.
    with open(settings_fn) as fp:
        settings = json.load(fp)["avida_analyses"]
    # Pull out locations of interest.
    experiment_base_loc = settings["experiment_base_location"]
    data_loc = os.path.join(experiment_base_loc, "data")
    analysis_scripts_loc = settings["avida_analysis_scripts_location"]
    configs_loc = os.path.join(experiment_base_loc, "configs")
    #
    avida_cmd_by_treatment = {}
    with open(os.path.join(configs_loc, "run_list"), "r") as fp:
        for line in fp:
            if "./avida" in line:
                mline = line.split(" ./avida ")
                treatment = mline[0].split(" ")[-1]
                args = mline[1].strip()
                avida_cmd_by_treatment[treatment] = args.replace("-s $seed", "")
    #
    analyses = settings["analysis_scripts_to_run"]
    # Find all things to analyze.
    treatments = settings["treatments"]
    # Hack for testing
    #testing_treatments = ["_".join(name.split("_")[:len(name.split("_")) - 1]) for name in os.listdir(data_loc) if "ntasks_" in name and os.path.isdir(os.path.join(data_loc, name))]
    #print set(testing_treatments)
    #for treatment in testing_treatments:
    for treatment in treatments:
        print "Analyzing %s" % treatment
        treatment_settings = settings["treatment_settings"][treatments[treatment]["settings"]]

        start_rep = treatment_settings["replicates"][0]
        end_rep = treatment_settings["replicates"][1]
        final_update = settings["final_update"]
        test_environments = treatment_settings["test_environments"]
        env_tests = "\n".join(["TEST_ENV___%s $t $i" % test for test in test_environments])

        for ascript in analyses:
            print "\tRunning %s" % ascript
            ascript_fpath = os.path.join(analysis_scripts_loc, ascript)
            ## Build temporary analysis file. ##
            # <start_replicate>
            # <end_replicate>
            # <final_update>
            # <base_experiment_directory>
            # <treatments>
            # <environment_tests>
            temp_ascript_content = ""
            with open(ascript_fpath, "r") as fp:
                temp_ascript_content = fp.read()
            temp_ascript_content = temp_ascript_content.replace("<start_replicate>", str(start_rep))
            temp_ascript_content = temp_ascript_content.replace("<end_replicate>", str(end_rep))
            temp_ascript_content = temp_ascript_content.replace("<final_update>", str(final_update))
            temp_ascript_content = temp_ascript_content.replace("<base_experiment_directory>", experiment_base_loc)
            temp_ascript_content = temp_ascript_content.replace("<treatments>", treatment)
            temp_ascript_content = temp_ascript_content.replace("<environment_tests>", env_tests)
            ## Write out the temporary analysis file to the run location. ##
            temp_ascript = os.path.join(configs_loc, "temp_" + ascript)
            with open(temp_ascript, "w") as fp:
                fp.write(temp_ascript_content)
            ## Build analysis comand. ##
            #  cmd = "./avida %s -a -set ANALYZE_FILE %s" % (avida_cmd_args[treatment], ascript)
            cmd = "./avida %s -a -set ANALYZE_FILE %s" % (avida_cmd_by_treatment[treatment], temp_ascript)
            ## Run Avida analysis ##
            return_code = subprocess.call(cmd, shell = True, cwd = configs_loc)
            ## Clean up analysis script ##
            return_code = subprocess.call("rm %s" % "temp_" + ascript, shell = True, cwd = configs_loc)


if __name__ == "__main__":
    main()
