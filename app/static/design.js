$(document).ready(function() {
  // bind event handler to change observer frame
  $("#display-observer").change(function() {
    if($(this).val() == "inertial") {
      viewer.scene.postUpdate.addEventListener(icrf);
    } else {
      viewer.scene.postUpdate.removeEventListener(icrf);
    }
  });

  $('#display-start').datetimepicker({
    timeZone: "UTC",
    defaultDate: moment.now(),
    icons: {
      time: 'far fa-clock'
    }
  });
  $("#display-duration").val(60);
  $("#display-step").val(30);

  // function simulate an inertial frame
  function icrf(scene, time) {
    if (scene.mode !== Cesium.SceneMode.SCENE3D) {
      return;
    }
    var icrfToFixed = Cesium.Transforms.computeIcrfToFixedMatrix(time);
    if (Cesium.defined(icrfToFixed)) {
      var camera = viewer.camera;
      var offset = Cesium.Cartesian3.clone(camera.position);
      var transform = Cesium.Matrix4.fromRotationTranslation(icrfToFixed);
      camera.lookAtTransform(transform, offset);
    }
  }

  function getOrbit() {
    if($("#satellite-orbit").val()=="tle") {
      return {
        type: "tle",
        tle: $("#orbit-tle").val().split("\n")
      }
    } else if($("#satellite-orbit").val()=="circular") {
      return {
        type: "circular",
        epoch: $("#orbit-circular-epoch").datetimepicker('viewDate').toISOString(),
        altitude: $("#orbit-circular-altitude").val()*1000,
        inclination: $("#orbit-circular-inclination").val(),
        right_ascension_ascending_node: $("#orbit-circular-raan").val(),
        true_anomaly: $("#orbit-circular-ta").val()
      }
    } else if($("#satellite-orbit").val()=="sso") {
      return {
        type: "sso",
        epoch: $("#orbit-sso-epoch").datetimepicker('viewDate').toISOString(),
        equator_crossing_time: $("#orbit-sso-ect").datetimepicker('viewDate').format("HH:mm"),
        equator_crossing_ascending: $("#orbit-sso-direction").val()=="asc",
        altitude: $("#orbit-sso-altitude").val()*1000,
        true_anomaly: $("#orbit-sso-ta").val()
      }
    } else if($("#satellite-orbit").val()=="keplerian") {
      return {
        type: "keplerian",
        epoch: $("#orbit-keplerian-epoch").datetimepicker('viewDate').toISOString(),
        altitude: $("#orbit-keplerian-altitude").val()*1000,
        inclination: $("#orbit-keplerian-inclination").val(),
        eccentricity: $("#orbit-keplerian-eccentricity").val(),
        perigee_argument: $("#orbit-keplerian-pa").val(),
        right_ascension_ascending_node: $("#orbit-keplerian-raan").val(),
        true_anomaly: $("#orbit-keplerian-ta").val()
      }
    }
  }

  function queryCelestrak() {
    var query = $("#celestrak-name").val();
    if(query.length > 2) {
      $("#celestrak-objects").val([]);
      $("#celestrak-objects").change();
      $("#celestrak-tle").val("");
      $.ajax({
        url: "celestrak/tle?name=" + query,
        type: "GET",
        crossDomain: true,
        success: function(response) {
          var lines = response.split(/\r?\n/);
          if(lines.length >= 3) {
            $("#celestrak-objects").empty();
            for(var i = 0; i < lines.length - 2; i += 3) {
              $("#celestrak-objects").append(
                $('<option label="' + lines[i] + '" value="' + lines[i] + '">')
                .data("tle", lines[i+1] + "\n" + lines[i+2])
              );
            }
          }
        }
      });
    }
  }

  function changeSatellite() {
    const satellite = $("#satellites option:selected").data("satellite");
    const features = $("#satellites option:selected").data("features");
    const color = $("#satellites option:selected").data("color");
    $("#satellite-edit").prop("disabled", satellite?false:true);
    $("#satellite-copy").prop("disabled", satellite?false:true);
    $("#satellite-delete").prop("disabled", satellite?false:true);
    $("#orbit-data-json").prop("disabled", features?false:true);
    $("#orbit-data-csv").prop("disabled", features?false:true);
    if(satellite) {
      $(".satellite-group").show();
      $("#satellite-name").val(satellite.name);
      if(!satellite.type) {
        $("#satellite-type").val("satellite");
      }
      if(satellite.type == "train") {
        $("#satellite-type").val("train");
        $(".constellation-train-group").show();
        $("#train-count").val(satellite.number_satellites);
        $("#train-interval").val(satellite.interval/60);
      } else {
        $(".constellation-train-group").hide();
      }
      if(satellite.type == "walker") {
        $("#satellite-type").val("walker");
        $(".constellation-walker-group").show();
        $("#walker-count").val(satellite.number_satellites);
        $("#walker-planes").val(satellite.number_planes);
        $("#walker-config").val(satellite.configuration);
        $("#walker-spacing").val(satellite.relative_spacing);
      } else {
        $(".constellation-walker-group").hide();
      }
      $("#satellite-orbit").val(satellite.orbit.type);
      $("#satellite-color").val(color);
      if(satellite.orbit.type == "tle") {
        $("#orbit-tle").val(satellite.orbit.tle.join("\n"));
      } else if(satellite.orbit.type == "circular") {
        var listener = $("#orbit-circular-epoch").off("change.datetimepicker");
        $("#orbit-circular-epoch").datetimepicker('date', satellite.orbit.epoch);
        $("#orbit-circular-epoch").on("change.datetimepicker", listener);
        $("#orbit-circular-altitude").val(satellite.orbit.altitude/1000);
        $("#orbit-circular-inclination").val(satellite.orbit.inclination);
        $("#orbit-circular-raan").val(satellite.orbit.right_ascension_ascending_node);
        $("#orbit-circular-ta").val(satellite.orbit.true_anomaly);
      } else if(satellite.orbit.type == "sso") {
        var listener = $("#orbit-sso-epoch").off("change.datetimepicker");
        $("#orbit-sso-epoch").datetimepicker('date', satellite.orbit.epoch);
        $("#orbit-sso-epoch").on("change.datetimepicker", listener);
        $("#orbit-sso-ect").datetimepicker('date', satellite.orbit.equator_crossing_time);
        $("#orbit-sso-direction").val(satellite.orbit.equator_crossing_ascending?"asc":"desc");
        $("#orbit-sso-altitude").val(satellite.orbit.altitude/1000);
        $("#orbit-sso-ta").val(satellite.orbit.true_anomaly);
      } else if(satellite.orbit.type == "keplerian") {
        var listener = $("#orbit-keplerian-epoch").off("change.datetimepicker");
        $("#orbit-keplerian-epoch").datetimepicker('date', satellite.orbit.epoch);
        $("#orbit-keplerian-epoch").on("change.datetimepicker", listener);
        $("#orbit-keplerian-altitude").val(satellite.orbit.altitude/1000);
        $("#orbit-keplerian-inclination").val(satellite.orbit.inclination);
        $("#orbit-keplerian-eccentricity").val(satellite.orbit.eccentricity);
        $("#orbit-keplerian-pa").val(satellite.orbit.perigee_argument);
        $("#orbit-keplerian-raan").val(satellite.orbit.right_ascension_ascending_node);
        $("#orbit-keplerian-ta").val(satellite.orbit.true_anomaly);
      }
      $("#satellite-orbit").val()=="tle" ? $(".orbit-tle-group").show() : $(".orbit-tle-group").hide();
      $("#satellite-orbit").val()=="circular" ? $(".orbit-circular-group").show() : $(".orbit-circular-group").hide();
      $("#satellite-orbit").val()=="sso" ? $(".orbit-sso-group").show() : $(".orbit-sso-group").hide();
      $("#satellite-orbit").val()=="keplerian" ? $(".orbit-keplerian-group").show() : $(".orbit-keplerian-group").hide();
      $("#instruments").empty();
      satellite.instruments.forEach(function(instrument) {
        var option = $("<option></option>")
          .data("instrument", instrument)
          .text(instrument.name);
        if(instrument===$("#satellites option:selected").data("instrument")) {
          option.attr("selected", "selected");
        }
        $("#instruments").append(option);
        changeInstruments();
      });
    } else {
      $(".satellite-group").hide();
    }
  }

  function changeInstruments() {
    const instrument = $("#instruments option:selected").data("instrument");
    $("#instrument-copy").prop("disabled", instrument?false:true);
    $("#instrument-delete").prop("disabled", (instrument?false:true) || $("#instruments option").length <= 1);
    if(instrument) {
      $(".instrument-group").show();
      $("#instrument-name").val(instrument.name);
      $("#instrument-field").val(instrument.field_of_regard);
      $("#instrument-access").val(instrument.min_access_time);
      $("#instrument-target-sunlit").val(
        instrument.req_target_sunlit?"sunlit":
          (instrument.req_target_sunlit===false?"eclipse":"none")
      );
      $("#instrument-self-sunlit").val(
        instrument.req_self_sunlit?"sunlit":
          (instrument.req_self_sunlit===false?"eclipse":"none")
      );
      $("#satellites option:selected").data("instrument", instrument);
    } else {
      $(".instrument-group").hide();
    }
  }

  $('#orbit-circular-epoch').datetimepicker({
    timeZone: "UTC",
    defaultDate: moment.now(),
    icons: {
      time: 'far fa-clock'
    }
  });
  $('#orbit-sso-epoch').datetimepicker({
    timeZone: "UTC",
    defaultDate: moment.now(),
    icons: {
      time: 'far fa-clock'
    }
  });
  $('#orbit-sso-ect').datetimepicker({
    timeZone: "UTC",
    defaultDate: moment.now(),
    format: 'LT'
  });
  $('#orbit-keplerian-epoch').datetimepicker({
    timeZone: "UTC",
    defaultDate: moment.now(),
    icons: {
      time: 'far fa-clock'
    }
  });

  function getSunlitState(value) {
    if(value=="sunlit") {
      return true;
    } else if(value=="eclipse") {
      return false;
    } else {
      return null;
    }
  }

  function collectOrbitTrack() {
    clearOrbitTrack();
    $("#orbit-analyze").prop("disabled", true);
    $("#orbit-analyze .spinner-border").show();
    $("#satellites > option").each(function() {
      const option = $(this);
      $.ajax({
        url: "/analyze/orbit-track",
        type: "POST",
        contentType: "application/json",
        dataType: "json",
        data: JSON.stringify({
          satellite: option.data("satellite"),
          instrument: option.data("instrument"),
          times: {
            start: $('#display-start').datetimepicker('viewDate').toISOString(),
            end: $('#display-start').datetimepicker('viewDate').clone().add(
              $("#display-duration").val(), "minutes").toISOString(),
            delta: $("#display-step").val()
          }
        }),
        success: function(response) {
          checkOrbitTrack(response.task_id, option);
          if(response.group_id) {
            checkProgress(response.group_id, option);
          }
        }
      });
    });
  };

  function checkProgress(groupId, option) {
    $.get(
      "/tasks/"+groupId+"/progress",
      function(progress) {
        if(progress.task_count == progress.completed_count) {
          $("#orbit-progress").hide();
          $("#orbit-progress .progress-bar").css('width', "0%");
          $("#orbit-progress .progress-bar").prop('aria-valuenow', 0);
        } else {
          $("#orbit-progress").show();
          var percent = Math.round(100*(progress.completed_count/progress.task_count));
          $("#orbit-progress .progress-bar").css('width', percent + "%");
          $("#orbit-progress .progress-bar").prop('aria-valuenow', percent);
          setTimeout(
            function() {
              checkProgress(groupId);
            }, 1000
          );
        }
      }
    );
  };
  function checkOrbitTrack(taskId, option) {
    $.get(
      "/tasks/"+taskId+"/status",
      function(status) {
        if(status.ready) {
          $.get(
            "/analyze/orbit-track/"+taskId,
            function(results) {
              option.data("features", results.features);
              $("#orbit-analyze .spinner-border").hide();
              $("#orbit-analyze").prop("disabled", false);
              $("#orbit-data").prop("disabled", false);
              $("#orbit-data-json").prop("disabled", false);
              $("#orbit-data-csv").prop("disabled", false);
              processOrbitTrack(option);
            }
          );
        } else {
          setTimeout(
            function() {
              checkOrbitTrack(taskId, option);
            }, 1000
          );
        }
      }
    );
  };

  var orbitTrack = viewer.entities.add(new Cesium.Entity());
  function processOrbitTrack(option) {
    var start = $('#display-start').datetimepicker('viewDate').toISOString();
    var end = $('#display-start').datetimepicker('viewDate').clone().add($("#display-duration").val(), "minutes").toISOString();
    var satellites = {};
    _.forEach(option.data("features"), function(feature) {
      var satellite = feature.properties.satellite;
      if(!satellites.hasOwnProperty(satellite)) {
        satellites[satellite] = {
          position: new Cesium.SampledPositionProperty(),
          radius: new Cesium.SampledProperty(Number),
          visibility: new Cesium.SampledProperty(Cesium.Color)
        };
        satellites[satellite].position.setInterpolationOptions({
          interpolationDegree: 2,
          interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
        });
      }
      satellites[satellite].position.addSample(
        Cesium.JulianDate.fromIso8601(feature.properties.time),
        Cesium.Cartesian3.fromDegrees(
          feature.geometry.coordinates[0],
          feature.geometry.coordinates[1],
          feature.geometry.coordinates[2]
        )
      );
      satellites[satellite].radius.addSample(
        Cesium.JulianDate.fromIso8601(feature.properties.time),
        feature.properties.swath_width/2
      );
      satellites[satellite].visibility.addSample(
        Cesium.JulianDate.fromIso8601(feature.properties.time),
        Cesium.Color.fromCssColorString(option.data("color")).withAlpha(feature.properties.valid_obs*0.5)
      );
    });
    for(var satellite in satellites) {
      viewer.entities.add({
        parent: orbitTrack,
        availability: new Cesium.TimeIntervalCollection([
          new Cesium.TimeInterval({
            start: Cesium.JulianDate.fromIso8601(start),
            stop: Cesium.JulianDate.fromIso8601(end)
          })
        ]),
        position: satellites[satellite].position,
        orientation: new Cesium.VelocityOrientationProperty(satellites[satellite].position),
        point: {
          pixelSize: 8,
          color: Cesium.Color.fromCssColorString(option.data("color"))
        },
        path: {
          leadTime: $("#display-path").is(":checked") ? $("#display-duration").val()*60 : 0,
          trailTime: $("#display-path").is(":checked") ? $("#display-duration").val()*60 : 0,
          material: new Cesium.ColorMaterialProperty(Cesium.Color.fromCssColorString(option.data("color")))
        },
        ellipse: {
          semiMajorAxis: satellites[satellite].radius,
          semiMinorAxis: satellites[satellite].radius,
          height: 0,
          heightReference: Cesium.HeightReference.CLAMP_TO_GROUND,
          material: new Cesium.ColorMaterialProperty(satellites[satellite].visibility),
          outline: true,
          outlineColor: Cesium.Color.BLACK
        }
      });
    }
    viewer.clock.startTime = Cesium.JulianDate.fromIso8601(start);
    viewer.clock.currentTime = Cesium.JulianDate.fromIso8601(start);
    viewer.clock.stopTime = Cesium.JulianDate.fromIso8601(end);
    viewer.clock.clockRange = Cesium.ClockRange.LOOP_STOP;
    viewer.timeline.zoomTo(
      Cesium.JulianDate.fromIso8601(start),
      Cesium.JulianDate.fromIso8601(end)
    );
  }
  function clearOrbitTrack() {
    clearEntities(orbitTrack);
  }

  $("#display-path").change(function() {
    _.forEach(viewer.entities.values, function(entity) {
      if(entity.parent === orbitTrack) {
        entity.path.leadTime = $("#display-path").is(":checked") ? $("#display-duration").val()*60 : 0;
        entity.path.trailTime = $("#display-path").is(":checked") ? $("#display-duration").val()*60 : 0;
      }
    });
  });
  $("#display-lighting").change(function() {
    viewer.scene.globe.enableLighting = $("#display-lighting").is(":checked");
  });
  $("#satellites").dblclick(function() {
    $("#satellite-dialog").modal('show');
  });
  $("#satellite-add").click(function() {
    const satellite = {
      name: "New Satellite",
      orbit: {
        type: "circular",
        altitude: 400000,
        inclination: 0,
        true_anomaly: 0,
        right_ascension_ascending_node: 0,
        epoch: $("#display-start").datetimepicker("viewDate").toISOString()
      },
      instruments: [
        {
          name: "Default",
          field_of_regard: 180,
          min_access_time: 0,
          req_target_sunlit: null,
          req_self_sunlit: null
        }
      ]
    };
    $("#satellites").append(
      $("<option></option>")
        .data("satellite", satellite)
        .data("instrument", satellite.instruments[0])
        .data("color", "#ff0000")
        .text(satellite.name)
        .attr("selected", "selected")
    );
    changeSatellite();
  });
  $("#satellite-copy").click(function() {
    const selected = $("#satellites option:selected");
    const satellite = $.extend({}, selected.data("satellite"));
    satellite.name += " (Copy)"
    $("#satellites").append(
      $("<option></option>")
        .data("satellite", satellite)
        .data("instrument", satellite.instruments[0])
        .data("color", selected.data("color"))
        .text(satellite.name)
        .attr("selected", "selected")
    );
    changeSatellite();
  });
  $("#satellite-delete").click(function() {
    if(confirm("Delete satellite " + $("#satellites :selected").data("satellite").name + "?")) {
      const index = $("#satellites").prop("selectedIndex");
      $("#satellites :selected").remove();
      if($("#satellites option").length > 0) {
        $("#satellites").prop("selectedIndex", Math.min(index, $("#satellites option").length - 1));
      }
      changeSatellite();
    }
  });

  var orbitTable;
  $("#orbit-dialog").on("show.bs.modal", function() {
    $("#orbit-dialog").find("tbody").empty();
    var features = [];
    $("#satellites option").each(function(index) {
      _.forEach($(this).data("features"), function(feature) {
          features.push({
            date: feature.properties.time,
            satellite: feature.properties.satellite,
            instrument: feature.properties.instrument,
            latitude: feature.geometry.coordinates[1],
            longitude: feature.geometry.coordinates[0],
            altitude: feature.geometry.coordinates[2]/1000,
            obsRadius: feature.properties.swath_width/2000,
            obsValid: feature.properties.valid_obs
          });
      });
    });

    if(orbitTable) {
      orbitTable.destroy();
    }
    orbitTable = $('#orbit-table').DataTable({
      'data': _.map(features, function(feature) {
        return [
          feature.date,
          feature.satellite,
          feature.instrument,
          feature.latitude,
          feature.longitude,
          feature.altitude,
          feature.obsRadius,
          feature.obsValid
        ];
      }),
      'columns': [
        { title: 'Date' },
        { title: 'Satellite' },
        { title: 'Instrument' },
        {
          title: 'Latitude (deg)',
          render: $.fn.dataTable.render.number(',', '.', 6)
        },
        {
          title: 'Longitude (deg)',
          render: $.fn.dataTable.render.number(',', '.', 6)
        },
        {
          title: 'Altitude (km)',
          render: $.fn.dataTable.render.number(',', '.', 2)
        },
        { title: 'Obs. Radius (km)',
          render: $.fn.dataTable.render.number(',', '.', 2)
        },
        { title: 'Obs. Valid' }
      ],
      'order': [[ 0, 'asc'], [ 1, 'asc' ], [2, 'asc']],
      'searching': false
    });
    $("#orbit-data-json").attr(
      "href",
      "data:text/json;charset=utf-8,"
      + encodeURIComponent(
        JSON.stringify(
          _.map(features, function(feature) {
            return {
              date: feature.date,
              satellite: feature.satellite,
              instrument: feature.instrument,
              latitude_deg: feature.latitude,
              longitude_deg: feature.longitude,
              altitude_km: feature.altitude,
              obsRadius_km: feature.obsRadius,
              obsValid: feature.obsValid
            }
          })
        )
      )
    );
    $("#orbit-data-json").attr("download", "orbit.json");
    $("#orbit-data-csv").attr(
      "href",
      "data:text/csv;charset=utf-8,"
      + encodeURIComponent(
        [
          [
            "date",
            "satellite",
            "instrument",
            "latitude_deg",
            "longitude_deg",
            "altitude_km",
            "obsRadius_km",
            "obsValid"
          ].join()
        ].concat(
          _.map(features, function(feature) {
            return [
              feature.date,
              feature.satellite,
              feature.instrument,
              feature.latitude,
              feature.longitude,
              feature.altitude,
              feature.obsRadius,
              feature.obsValid
            ].join();
          })
        ).join("\n")
      )
    );
    $("#orbit-data-csv").attr("download", "orbit.csv");
  });

  $("#satellite-name").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.name = $("#satellite-name").val();
    $("#satellites option:selected").text($("#satellite-name").val());
    $("#orbit-data").prop("disabled", true);
  });
  $("#satellite-color").change(function() {
    $("#satellites option:selected").data("color", $("#satellite-color").val());
  });
  $("#satellite-type").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    if($("#satellite-type").val()=="train") {
      satellite.type = "train";
      satellite.number_satellites = $("#train-count").val();
      satellite.interval = $("#train-interval").val()*60;
      $(".constellation-train-group").show();
    } else {
      delete satellite.interval;
      $(".constellation-train-group").hide();
    }
    if($("#satellite-type").val()=="walker") {
      satellite.type = "walker";
      satellite.number_satellites = $("#walker-count").val();
      satellite.number_planes = $("#walker-planes").val();
      satellite.configuration = $("#walker-config").val();
      satellite.relative_spacing = $("#walker-spacing").val();
      $(".constellation-walker-group").show();
    } else {
      delete satellite.number_planes;
      delete satellite.configuration;
      delete satellite.relative_spacing;
      $(".constellation-walker-group").hide();
    }
    if($("#satellite-type").val()=="satellite") {
      delete satellite.type;
      delete satellite.number_planes;
    }
    $("#orbit-data").prop("disabled", true);
  });
  $("#train-count").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.number_satellites = $("#train-count").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#train-interval").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.interval = $("#train-interval").val()*60;
    $("#orbit-data").prop("disabled", true);
  });
  $("#walker-count").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.number_satellites = $("#walker-count").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#walker-planes").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.number_planes = $("#walker-planes").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#walker-config").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.configuration = $("#walker-config").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#walker-spacing").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.relative_spacing = $("#walker-spacing").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#satellite-orbit").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    $("#satellite-orbit").val()=="tle" ? $(".orbit-tle-group").show() : $(".orbit-tle-group").hide();
    $("#satellite-orbit").val()=="circular" ? $(".orbit-circular-group").show() : $(".orbit-circular-group").hide();
    $("#satellite-orbit").val()=="sso" ? $(".orbit-sso-group").show() : $(".orbit-sso-group").hide();
    $("#satellite-orbit").val()=="keplerian" ? $(".orbit-keplerian-group").show() : $(".orbit-keplerian-group").hide();
    satellite.orbit = getOrbit();
    $("#orbit-data").prop("disabled", true);
  });
  $("#celestrak-dialog").on("shown.bs.modal", queryCelestrak);
  $("#celestrak-name").change(queryCelestrak);
  $("#celestrak-objects").change(function() {
    if($("#celestrak-objects option:selected").length > 0) {
      $("#celestrak-tle").val($("#celestrak-objects option:selected").data("tle"));
      $("#celestrak-tle-ok").prop('disabled', false);
    } else {
      $("#celestrak-tle-ok").prop('disabled', true);
    }
  });
  $("#celestrak-tle-ok").click(function() {
    $("#orbit-tle").val($("#celestrak-tle").val());
    $("#orbit-tle").change();
  })
  $("#orbit-tle").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.tle = $("#orbit-tle").val().split("\n");
    $("#orbit-data").prop("disabled", true);
  });
  $("#orbit-circular-epoch").on("change.datetimepicker", function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.epoch = $("#orbit-circular-epoch").datetimepicker("viewDate").toISOString();
    $("#orbit-data").prop("disabled", true);
  });
  $(".orbit-circular-group").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.altitude = $("#orbit-circular-altitude").val()*1000;
    satellite.orbit.inclination = $("#orbit-circular-inclination").val();
    satellite.orbit.right_ascension_ascending_node = $("#orbit-circular-raan").val();
    satellite.orbit.true_anomaly = $("#orbit-circular-ta").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#orbit-sso-epoch").on("change.datetimepicker", function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.epoch = $("#orbit-sso-epoch").datetimepicker("viewDate").toISOString();
    $("#orbit-data").prop("disabled", true);
  });
  $("#orbit-sso-ect").on("change.datetimepicker", function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.equator_crossing_time = $("#orbit-sso-ect").datetimepicker("viewDate").format("HH:mm");
    $("#orbit-data").prop("disabled", true);
  });
  $(".orbit-sso-group").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.equator_crossing_ascending = $("#orbit-sso-direction").val()=="asc";
    satellite.orbit.altitude = $("#orbit-sso-altitude").val()*1000;
    satellite.orbit.true_anomaly = $("#orbit-sso-ta").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#orbit-keplerian-epoch").on("change.datetimepicker", function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.epoch = $("#orbit-keplerian-epoch").datetimepicker("viewDate").toISOString();
    $("#orbit-data").prop("disabled", true);
  });
  $(".orbit-keplerian-group").change(function() {
    var satellite = $("#satellites option:selected").data("satellite");
    satellite.orbit.altitude = $("#orbit-keplerian-altitude").val()*1000;
    satellite.orbit.inclination = $("#orbit-keplerian-inclination").val();
    satellite.orbit.eccentricity = $("#orbit-keplerian-eccentricity").val();
    satellite.orbit.perigee_argument = $("#orbit-keplerian-pa").val();
    satellite.orbit.right_ascension_ascending_node = $("#orbit-keplerian-raan").val();
    satellite.orbit.true_anomaly = $("#orbit-keplerian-ta").val();
    $("#orbit-data").prop("disabled", true);
  });
  $("#instrument-add").click(function() {
    const satellite = $("#satellites option:selected").data("satellite");
    const instrument = {
      name: "New Instrument",
      field_of_regard: 180,
      min_access_time: 0,
      req_target_sunlit: null,
      req_self_sunlit: null
    };
    satellite.instruments.push(instrument);
    $("#satellites option:selected").data("instrument", instrument);
    $("#instruments").append(
      $("<option></option>")
        .data("instrument", instrument)
        .text(instrument.name)
        .attr("selected", "selected")
    );
    changeInstruments();
    $("#orbit-data").prop("disabled", true);
  });
  $("#instrument-copy").click(function() {
    const satellite = $("#satellites option:selected").data("satellite");
    const selected = $("#instruments option:selected");
    const instrument = $.extend({}, selected.data("instrument"));
    instrument.name += " (Copy)";
    satellite.instruments.push(instrument);
    $("#satellites option:selected").data("instrument", instrument);
    $("#instruments").append(
      $("<option></option>")
        .data("instrument", instrument)
        .text(instrument.name)
        .attr("selected", "selected")
    );
    changeInstruments();
    $("#orbit-data").prop("disabled", true);
  });
  $("#instrument-delete").click(function() {
    if(confirm("Delete instrument " + $("#instruments :selected").data("instrument").name + "?")) {
      const index = $("#instruments").prop("selectedIndex");
      $("#instruments :selected").remove();
      $("#instruments").prop("selectedIndex", Math.min(index, $("#instruments option").length - 1));
      changeInstruments();
      $("#orbit-data").prop("disabled", true);
    }
  });
  $("#instrument-name").change(function() {
    var instrument = $("#instruments option:selected").data("instrument");
    instrument.name = $("#instrument-name").val();
    $("#instruments option:selected").text($("#instrument-name").val());
    $("#orbit-data").prop("disabled", true);
  });
  $(".instrument-group").change(function() {
    var instrument = $("#instruments option:selected").data("instrument");
    instrument.field_of_regard = $("#instrument-field").val();
    instrument.min_access_time = $("#instrument-access").val();
    instrument.req_self_sunlit = getSunlitState($("#instrument-self-sunlit").val());
    instrument.req_target_sunlit = getSunlitState($("#instrument-target-sunlit").val());
    $("#orbit-data").prop("disabled", true);
  });

  $("#orbit-analyze").click(collectOrbitTrack);
  $("#satellite-dialog").on('hidden.bs.modal', collectOrbitTrack);
  $("#nav-design-tab").on("show.bs.tab", changeSatellite);
  $("#nav-design-tab").on("show.bs.tab", collectOrbitTrack);
  $("#nav-design-tab").on("hide.bs.tab", clearOrbitTrack);

  // load initial data
  var satellites = [
    {
      name: "International Space Station",
      orbit: {
        type: "tle",
        tle: [
          "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
          "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754"
        ]
      },
      instruments: [
        {
          name: "Default",
          field_of_regard: 180,
          min_access_time: 0,
          req_target_sunlit: null,
          req_self_sunlit: null
        }
      ]
    }
  ];
  satellites.forEach(function(satellite) {
    $("#satellites").append(
      $("<option></option>")
        .data("satellite", satellite)
        .data("instrument", satellite.instruments[0])
        .data("color", "#ff0000")
        .text(satellite.name)
    );
  });
  $("#satellites").change(changeSatellite);
  $("#instruments").change(changeInstruments);
});
