// Globals
var currentTreatment = DEFAULT_TREATMENT;

var margin = {top: 20, right: 40, bottom: 20, left: 100};
var frameWidth = 940;
var frameHeight = 1500;
var canvasWidth = frameWidth - margin.left - margin.right;
var canvasHeight = frameHeight - margin.top - margin.bottom;

var zoomMult = 0.20;
var zoomRate = 0.05;
var minZoom = 0.05;
var maxZoom = 1.00;

var slicedRanges = [[0, 500], [95000, 105000], [195000, 200000]];
var updateTickInterval = 500;

var lineageSpacer = 3;
var lineageWidth = 20;
var sliceSpacer = 50;
var envWidth = 10;
var envSpacer = 5;
var maxReps = 75;

var repsByTreatment = {};
var envSequenceByTreatment = {};
var treatments = [];
var tasksByEnvironmentMemo = {};
var maxScoreByTreatment = {};

var getEnvTooltipHTML = function(env_data) {
  var tt_content = "<div class='well'>" +
                    "Environment: " + env_data.state +
                    "<br/>Start Update: " + env_data.true_start +
                    "<br/>Duration: " + env_data.true_duration +
                    "</div>";
  return tt_content;
}

var getLineageStateTooltipHTML = function(state_data) {
  var phenotypeTable = "<table class='table table-striped'>";
  var traitSet = new Set();
  // Build header.
  phenotypeTable += "<tr><th>Environment</th>";
  for (env in state_data.phenotype) {
    for (trait in state_data.phenotype[env]) {
      traitSet.add(trait);
    }
  }
  traitSet = [...traitSet];
  for (var ti = 0; ti < traitSet.length; ti++) {
    phenotypeTable += "<th>Trait: " + traitSet[ti] + "</th>";
  }
  phenotypeTable += "</tr>";
  for (env in state_data.phenotype) {
    phenotypeTable += "<tr><td>" + env + "</td>";
    for (var ti = 0; ti < traitSet.length; ti++) {
      phenotypeTable += "<td>" + state_data.phenotype[env][traitSet[ti]] + "</td>";
    }
    phenotypeTable += "</tr>";
  }

  phenotypeTable += "</table>"
  var tt_content = "<div class='well'>" +
                   "<br/>Start update: " + state_data.true_start  +
                   "<br/>Duration: " + state_data.true_duration +
                   "<br/>Phenotype Signature: " + state_data.state +
                   "<br/>Phenotype: " + phenotypeTable +
                   "</div>";
  return tt_content;
}

var getTraitsExpressed = function(phenotype) {
  var traits_expressed = new Set();
  for (var env in phenotype) {
    for (var trait in phenotype[env]) {
      if (phenotype[env][trait] == "1") {
        traits_expressed.add(trait);
      }
    }
  }
  return traits_expressed;
}

var tooltip = d3.select("body")
                        .append("div")
                        .attr({"class": "t-tip"})
                        .style("position", "absolute")
                        .style("z-index", "10")
                        .style("visibility", "hidden");

var getEnvironmentTasks = function(environment) {
  /*
  Given an environment, return the tasks that were in the environment.
  */
  // Have we seen this environment before?
  if (!(environment in tasksByEnvironmentMemo)) {
    var env_tasks = [];
    var attrs = environment.split("__");
    for (var i = 0; i < attrs.length; i++) {
      var tasks = attrs[i].split("_").slice(1);
      for (var ti = 0; ti < tasks.length; ti++) {
        if (!(tasks[ti] in env_tasks)) {
          env_tasks.push(tasks[ti])
        }
      }
    }
    // Put tasks in rank order.
    tasksByEnvironmentMemo[environment] = env_tasks;
  }
  return tasksByEnvironmentMemo[environment];
}

var getPhenotype = function(treatment, encoded_phenotype, environment_key) {
  /*
    verb_phenotype: xx|xx|xx|xx
    environment_key: env1|env2|env3|env4
    decoded_phenotype_by_env = {env:{GetEnvironmentTasks(env)[ti]:phenotype_by_env[env][ti] for ti in range(0, len(GetEnvironmentTasks(env)))} for env in phenotype_by_env}
  */
  var decoded_phen_by_env = {};
  var encoded_phens = encoded_phenotype.split(ENV_DELIMITER);
  var env_key = environment_key.split(ENV_DELIMITER);
  var trait_order = exp_settings["phenotype_encoding"]["trait_order"];
  for (var ei = 0; ei < env_key.length; ei++) {
    var env = env_key[ei];
    var enc_env_phenotype = encoded_phens[ei];
    var env_traits = getEnvironmentTasks(env);
    var env_traits_in_order = [];
    for (var ti = 0; ti < trait_order.length; ti++) {
      if (env_traits.indexOf(trait_order[ti]) != -1) {
        env_traits_in_order.push(trait_order[ti]);
      }
    }
    decoded_phen_by_env[env] = {};
    for (var o = 0; o < env_traits_in_order.length; o++) {
      decoded_phen_by_env[env][env_traits_in_order[o]] = enc_env_phenotype[o];
    }
  }
  return decoded_phen_by_env;
}

var generateEnvironmentSequence = function(treatment) {
  /*
    Relevant things:
      * settings[settings["treatments"][treatment]["settings"]]["experienced_environments"]
      * change rate
      rep_attrs = {ra.split("_")[0]:ra.split("_")[1:] for ra in rep}
  */
  // Get change rate:
  var attributes = getTreatmentAttributes(treatment);
  var change_rate = Number(attributes["cr"][0]);
  var total_updates = exp_settings["final_update"];
  var environments = exp_settings["treatment_settings"][exp_settings["treatments"][treatment]["settings"]]["experienced_environments"];
  var cur_env = 0;
  var collapsed_seq = [];
  if (change_rate == 0) {
    collapsed_seq.push({start:0, duration:total_updates + 1, state:environments[cur_env]});
    return collapsed_seq;
  }
  for (var i = 0; i < total_updates; i += change_rate) {
    var start_update = i;
    var duration = change_rate;
    // Clip duration if it would go above max update.
    if ((start_update + duration) > total_updates) {
      duration = total_updates - start_update;
    }
    // Add env to sequence.
    collapsed_seq.push({"state":environments[cur_env], "start":start_update, "duration":duration});
    // Get the current environment.
    cur_env = (cur_env + 1) % environments.length;
  }
  return collapsed_seq;
}

var getTreatmentAttributes = function(treatment) {
  /*
   * Given a treatment name, return dictionary of attribute-value pairs.
  */
  var treatment_attrs = treatment.split("__");
  var attributes = {};
  for (var i = 0; i < treatment_attrs.length; i++) {
    var attribute = treatment_attrs[i].split("_")[0];
    var values = treatment_attrs[i].split("_").slice(1);
    attributes[attribute] = values;
  }
  return attributes;
}

var lineageDataAccessor = function(row) {
  /*
    This function is the data accessor for the fdom lineage .csv data file.
    Attributes:
      treatment,
      replicate,
      final_is_plastic,
      final_is_optimal,
      final_phenotype_score,
      max_phenotype_score,
      environment_key,
      lineage_phenotype_signature_sequence_abbrev,
      lineage_phenotype_verbose_signature_sequence_abbrev,
      lineage_phenotype_score_sequence_abbrev,
      lineage_is_plastic_sequence_abbrev,
      lineage_is_optimal_sequence_abbrev,
      lineage_phenotype_start_updates_abbrev,
      lineage_phenotype_duration_updates_abbrev
  */
  // Get everything of interest from the data file.
  var treatment = row.treatment;
  var replicate = row.replicate;
  var env_key = row.environment_key;
  var final_is_plastic = row.final_is_plastic;
  var final_is_optimal = row.final_is_optimal;
  var final_score = row.final_phenotype_score;
  var max_score = row.max_phenotype_score;
  var state_seq = row.lineage_phenotype_signature_sequence_abbrev.split(SEQ_DELIMITER);
  var state_verb_seq = row.lineage_phenotype_verbose_signature_sequence_abbrev.split(SEQ_DELIMITER);
  var score_seq = row.lineage_phenotype_score_sequence_abbrev.split(SEQ_DELIMITER).map(Number);
  var is_plastic_seq = row.lineage_is_plastic_sequence_abbrev.split(SEQ_DELIMITER).map(Number);
  var is_optimal_seq = row.lineage_is_optimal_sequence_abbrev.split(SEQ_DELIMITER).map(Number);
  var start_updates_seq = row.lineage_phenotype_start_updates_abbrev.split(SEQ_DELIMITER).map(Number);
  var duration_updates_seq = row.lineage_phenotype_duration_updates_abbrev.split(SEQ_DELIMITER).map(Number);
  // Build state sequence: {state:phen_sig, duration:duration_updates, start:start_update}
  var seq = [];
  for (i = 0; i < state_seq.length; i++) {
    seq.push({state: state_seq[i],
              duration: duration_updates_seq[i],
              start: start_updates_seq[i],
              is_plastic: is_plastic_seq[i],
              is_optimal: is_optimal_seq[i],
              score: score_seq[i],
              phenotype: getPhenotype(treatment, state_verb_seq[i], env_key),
            });
  }
  // Update treatments, reps x treatment, and env seq by treatment.
  if (treatments.indexOf(treatment) == -1) {
    treatments.push(treatment);
    // Generate environment sequence for treatment.
    envSequenceByTreatment[treatment] = generateEnvironmentSequence(treatment);
    maxScoreByTreatment[treatment] = Number(max_score);
    repsByTreatment[treatment] = [];
  }
  repsByTreatment[treatment].push(replicate);
  return {
    treatment: treatment,
    replicate: replicate,
    final_is_plastic: final_is_plastic,
    final_is_optimal: final_is_optimal,
    final_score: final_score,
    max_score: max_score,
    state_sequence: seq
  };
}

var sliceSequence = function(sequence, ranges) {
  /* Given a sequence and a range, return a sliced version of the sequence. */
  var slicedSeq = [];
  for (var si = 0; si < sequence.length; si++) {
    for (var ri = 0; ri < ranges.length; ri++) {
      var start = sequence[si]["start"];
      var end = sequence[si]["duration"] + start;
      if ((start <= ranges[ri][1] && start >= ranges[ri][0]) ||
          (end <= ranges[ri][1] && end >= ranges[ri][0]) ||
          start <= ranges[ri][0] && end >= ranges[ri][1]) {
          var clipped_state = {};
          for (var key in sequence[si]) {
            clipped_state[key] = sequence[si][key];
          }
          clipped_state["true_duration"] = sequence[si]["duration"];
          clipped_state["true_start"] = sequence[si]["start"];
          if (start < ranges[ri][0]) {
            clipped_state.start = ranges[ri][0];
            clipped_state.duration = end - clipped_state.start;
          }
          if (end > ranges[ri][1]) {
            clipped_state.duration = ranges[ri][1] - clipped_state.start;
          }
          slicedSeq.push(clipped_state);
      }
    }
  }
  return slicedSeq;
}

var clipToRange = function(val, minVal, maxVal) {
  if (val < minVal) { return minVal; }
  else if (val > maxVal) { return maxVal; }
  else { return val; }
}

// var transformDataSequences = function(data, fun, ...fun_args) {
//   console.log("transformDataSequences");
//   var transformedData = [];
//   for (var i = 0; i < data.length; i++) {
//
//     transformedData.push(fun(data[i].state_sequence, fun_args));
//   }
//   return transformedData;
// }

var runVisualization = function(data) {
  /*
    This function is called when data is loaded.
    It initializes and runs the visualization.
  */
  //////// Setup the canvas ////////
  var chartArea = d3.select("#chart_area");
  var frame = chartArea.append("svg");
  var canvas = frame.append("g");
  var envCanvas = canvas.append("g").attr({"class": "env_canvas"});
  var dataCanvas = canvas.append("g").attr({"class": "data_canvas"});

  //////// Collect some user input ////////
  var displayFull = $("#slice_toggle").prop("checked");
  var display = $('input[name="display"]:checked').val();

  //////// Visualization functions ////////
  var refreshDashboard = function() {
    /* The function refreshes the visualization dashboard. */
    // Update treatment dropdown text.
    var treatmentDropDownButton = $("#treatment_selector").text(currentTreatment);
  }

  var treatmentDropdownCallback = function() {
    /* Called on treatment drop down click. */
    var selection = $(this).attr("value");
    // Update the current treatment
    currentTreatment = selection;
    // Refresh the dashboard.
    refreshDashboard();
    // Redraw the visualization.
    update();
  }

  var update = function() {
    /* Call this function to redraw visualization on screen. */
    // What data are we chugging?
    console.log("=== Update ====");
    var treatmentData = data.filter(function(d) { return d.treatment == currentTreatment; })
                            .filter(function(d) { if (display == "plastic") {
                                                    return d.final_is_plastic == 1;
                                                  } else if (display == "nonplastic") {
                                                    return d.final_is_plastic == 0;
                                                  } else if (display == "all") {
                                                    return d;
                                                  }
                                                });
    // Now that we have only relevant data, slice and dice.
    var dataRanges = null;
    if (displayFull) {
      dataRanges = [[0, exp_settings["final_update"]]];
    } else {
      dataRanges = slicedRanges;
    }
    // Calculate magnitude of display range.
    var totalRange = 0;
    for (var i = 0; i < dataRanges.length; i++) {
      totalRange += dataRanges[i][1] - dataRanges[i][0];
    }

    var getRangeID = function(state_seq_obj) {
      /* Helper function:
        Given state sequence object, determine which range it belongs to.
      */
      var start = state_seq_obj.start;
      for (var i = 0; i < dataRanges.length; i++) {
        if (start <= dataRanges[i][1] && start >= dataRanges[i][0]) {
          return i;
        }
      }
      // Failure
      return -1;
    }

    ///////////////////////////////////////////////////////
    // Setup frame/canvas parameters.
    ///////////////////////////////////////////////////////
    frameWidth = $("#vis_panel").width() - 20;
    frameHeight = totalRange;
    canvasWidth = frameWidth - margin.left - margin.right;
    canvasHeight = frameHeight - margin.top - margin.bottom;
    frame.attr({"width": frameWidth, "height": frameHeight});
    canvas.attr({"transform": "translate(" + margin.left + "," + margin.top + ")"});
    xDomain = [0, maxReps * (lineageWidth + lineageSpacer)];
    yDomain = [0, exp_settings["final_update"]];
    // Setup color scale.
    minScore = -1 * maxScoreByTreatment[currentTreatment];
    maxScore = maxScoreByTreatment[currentTreatment];
    colorDomain = [];
    for (var i = minScore; i <= maxScore; i++) {
      colorDomain.push(i);
    }
    var colorScale = d3.scale.ordinal().domain(colorDomain).range(colorbrewer.RdYlGn[colorDomain.length]);
    ///////////////////////////////////////////////////////
    // Setup x axis
    ///////////////////////////////////////////////////////
    // clean up old x axis
    canvas.selectAll("g.x_axis").remove();
    canvas.selectAll("text#x_axis_label").remove();
    // make new x axis
    var xScale = d3.scale.linear();
    xScale.domain(xDomain).range([0, canvasWidth]);
    var xAxis = d3.svg.axis().scale(xScale).tickValues([]).orient("top");
    canvas.append("g").attr({"class": "x_axis"}).call(xAxis);
    // axis labels
    canvas.append("text").attr({"id": "x_axis_label", "class": "axis_label", "x":xScale(xDomain[1] / 2), "y": 0 - margin.top / 2})
                        .style("text-anchor", "middle")
                        .text("Hello there2!");
    ///////////////////////////////////////////////////////
    // Setup y axis
    ///////////////////////////////////////////////////////
    // clean up old y axis
    canvas.selectAll("g.y_axis").remove();
    canvas.selectAll("text#y_axis_label").remove();
    // make new scales/axes
    var yScales = [];
    var prev_range_end = 0;
    for (var i = 0; i < dataRanges.length; i++) {
      var current_range_end = prev_range_end + (((dataRanges[i][1] - dataRanges[i][0]) / totalRange) * (canvasHeight - (sliceSpacer * (dataRanges.length - 1))))
      var yScale = d3.scale.linear();
      yScale.domain(dataRanges[i]).range([prev_range_end, current_range_end]);
      yScales.push(yScale);
      prev_range_end = current_range_end + sliceSpacer;
      var tick_num = (dataRanges[i][1] - dataRanges[i][0]) / updateTickInterval;
      var yAxis = d3.svg.axis().scale(yScale).ticks(tick_num).orient("left");
      canvas.append("g").attr({"class": "y_axis", "id": "y_axis-" + i}).call(yAxis);
    }
    canvas.append("text").attr({"id": "y_axis_label", "class": "axis_label", "x":0 - (canvasHeight / 2), "y":0 - (margin.left / 2), "transform": "rotate(-90)"})
                    .style("text-anchor", "middle")
                    .text("Updates");
    ///////////////////////////////////////////////////////
    // Draw environment indicator
    ///////////////////////////////////////////////////////
    // Slice the environment.
    var slicedEnv = sliceSequence(envSequenceByTreatment[currentTreatment], dataRanges);
    // Draw!
    var envBlocks = envCanvas.selectAll("rect").data(slicedEnv);
    envBlocks.enter().append("rect");
    envBlocks.exit().remove();
    envBlocks.attr({"y": function(d) { return yScales[getRangeID(d)](d.start); },
                    "x": function(d) { return xScale(0); },
                    "width": xScale(envWidth),
                    "height": function(d) {
                                var si = getRangeID(d);
                                return yScales[si](dataRanges[si][0] + d.duration) - yScales[si](dataRanges[si][0]);
                              },
                    "class": function(d) { return d.state; }
                   });
    envBlocks.on("mouseover", function(d) {
                               return tooltip.style("visibility", "visible")
                                             .html(getEnvTooltipHTML(d));
                             })
              .on("mousemove", function() { return tooltip.style("top", (event.pageY-10)+"px").style("left",(event.pageX+10)+"px"); })
              .on("mouseout", function() { return tooltip.style("visibility", "hidden"); });
    ///////////////////////////////////////////////////////
    // Draw data
    ///////////////////////////////////////////////////////
    // Add a group for each lineage.
    var lineages = dataCanvas.selectAll("g").data(treatmentData, function(d) { return d.replicate; });
    lineages.enter().append("g");
    lineages.exit().remove();
    lineages.attr({"id": function(d) { return d.replicate; }});
    // Add a rectangle for each state in the sequence.
    console.log(lineages);
    lineages.each(function(lin_d, i) {
                    // var slicedData = lin_d.state_sequence.filter(function(d) {
                    //                                                 var start = d.start;
                    //                                                 var end = d.start + d.duration;
                    //                                                 var inc = false;
                    //                                                 for (var ri = 0; ri < dataRanges.length; ri++) {
                    //                                                         inc = inc || (start <= dataRanges[ri][1] && start >= dataRanges[ri][0]) ||
                    //                                                               (end <= dataRanges[ri][1] && end >= dataRanges[ri][0]) ||
                    //                                                               (start <= dataRanges[ri][0] && end >= dataRanges[ri][1]);
                    //                                                 }
                    //                                                 return inc;
                    //                                               });
                    slicedData = sliceSequence(lin_d.state_sequence, dataRanges);
                    var rep_id = Number(lin_d.replicate);
                    var stateBlocks = d3.select(this)
                                        .selectAll("rect").data(slicedData);

                    stateBlocks.enter().append("rect");
                    stateBlocks.exit().remove();
                    stateBlocks.attr({"y": function(d) {
                                                        var ri = getRangeID(d);
                                                        return yScales[ri](d.start); },
                                      "x": function(d) { return xScale((i * (lineageSpacer + lineageWidth)) + envWidth + envSpacer); },
                                      "height": function(d) {
                                                  var ri = getRangeID(d);
                                                  return yScales[ri](dataRanges[ri][0] + d.duration) - yScales[ri](dataRanges[ri][0]); },
                                      "width": xScale(lineageWidth),
                                      "class": function(d) { return "C" + d.state; },
                                      "fill": function(d) {
                                        if (d.is_plastic) {
                                            return colorScale(d.score);
                                        } else {
                                          var traitsExpressed = getTraitsExpressed(d.phenotype);
                                          if (traitsExpressed.has("NAND") && traitsExpressed.has("NOT")) {
                                            return "white";
                                          } else if (traitsExpressed.has("NAND")) {
                                            return "cyan";
                                          } else if (traitsExpressed.has("NOT")) {
                                            return "blue";
                                          } else {
                                            return "grey";
                                          }
                                        }
                                      }
                    });
                    stateBlocks.on("mouseover", function(d) {
                                  return tooltip.style("visibility", "visible")
                                                .html(getLineageStateTooltipHTML(d)); })
                                .on("mousemove", function() { return tooltip.style("top", (event.pageY-10)+"px").style("left",(event.pageX+10)+"px"); })
                                .on("mouseout", function() { return tooltip.style("visibility", "hidden"); });
    });

  }

  //////// Build dashboard ////////
  // Populate the treatment dropdown.
  var treatmentDropdown = $("#treatment-selection-dropdown");
  treatmentDropdown.empty();
  $.each(treatments, function(i, p) {
    var li = $("<li/>")
              .appendTo(treatmentDropdown);
    var a = $("<a/>")
              .attr({"value": this, "href": "#"})
              .text(this)
              .appendTo(li);
  });
  var treatmentDropdownButton = $("#treatment_selector").text(currentTreatment);
  $("<span/>").attr({"class": "caret"}).appendTo(treatmentDropdownButton);
  // Setup component listeners.
  $(document).ready(function() {
    // Treatment dropdown listener.
    $("#treatment-selection-dropdown li a").click(treatmentDropdownCallback);
    // Lineage type filtering.
    $("input[type='radio']").on("change", function(){
        display = $('input[name="display"]:checked').val();
        console.log(display);
        refreshDashboard();
        update();
    });
    // Zoom
    // $("#vis_zoom_in").click(function() {
    //   // zoom in!
    //   zoomMult = zoomMult + zoomRate;
    //   if (zoomMult > maxZoom) {
    //     zoomMult = maxZoom;
    //   }
    //   update();
    // });
    // $("#vis_zoom_out").click(function() {
    //   // zoom out!
    //   zoomMult = zoomMult - zoomRate;
    //   if (zoomMult < minZoom) {
    //     zoomMult = minZoom;
    //   }
    //   update();
    // });
    // slice vs. full
    $("#slice_toggle").change(function() {
      displayFull = $(this).prop("checked");
      refreshDashboard();
      update();
    });
    $(window).resize(function() {
      update();
    });
  });
  update();
}

var main = function() {
  // Load lineage data.
  d3.csv(data_fpath, lineageDataAccessor, runVisualization);
}

main();
