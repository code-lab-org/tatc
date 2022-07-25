$(document).ready(function() {
  $('#analysis-start').datetimepicker({
    timeZone: "UTC",
    defaultDate: moment.now(),
    icons: {
      time: 'far fa-clock'
    }
  });
  $("#analysis-duration").val(30);
  function collectCoverage() {
    clearCoverage();
    coveragePoints.show = true;
    coverageCells.show = true;
    $("#coverage-analyze").prop("disabled", true);
    $("#coverage-analyze .spinner-border").show();
    $.ajax({
      url: "/analyze/coverage",
      type: "POST",
      contentType: "application/json",
      dataType: "json",
      data: JSON.stringify({
        points: {
            method: $("#point-method").val(),
            distance: $("#point-distance").val()*1000,
            mask: getRegionMask()
        }, // TODO replace with list of points for better responsiveness
        cells: {
            method: $("#cell-method").val(),
            distance: $("#cell-distance").val()*1000,
            mask: getRegionMask()
        }, // TODO replace with list of cells for better responsiveness
        satellites: _.map($("#satellites > option"), function(option) {
          return $(option).data("satellite");
        }),
        start: $('#analysis-start').datetimepicker('viewDate').toISOString(),
        end: $('#analysis-start').datetimepicker('viewDate').clone().add(
          $("#analysis-duration").val(), "days").toISOString(),
      }),
      success: function(response) {
        checkCoverage(response.task_id);
        checkProgress(response.group_id);
      }
    });
  }

  function checkProgress(groupId) {
    $.get(
      "/tasks/"+groupId+"/progress",
      function(progress) {
        if(progress.task_count == progress.completed_count) {
          $("#coverage-progress").hide();
          $("#coverage-progress .progress-bar").css('width', "0%");
          $("#coverage-progress .progress-bar").prop('aria-valuenow', 0);
        } else {
          $("#coverage-progress").show();
          var percent = Math.round(100*(progress.completed_count/progress.task_count));
          $("#coverage-progress .progress-bar").css('width', percent + "%");
          $("#coverage-progress .progress-bar").prop('aria-valuenow', percent);
          setTimeout(
            function() {
              checkProgress(groupId);
            }, 1000
          );
        }
      }
    );
  }

  function checkCoverage(taskId) {
    $.get(
      "/tasks/"+taskId+"/status",
      function(status) {
        if(status.ready) {
          $.get(
            "/analyze/coverage/"+taskId,
            function(results) {
              $("#coverage-analyze").data("points", results.points.features);
              $("#coverage-analyze").data("cells", results.cells.features);
              $("#coverage-analyze .spinner-border").hide();
              $("#coverage-analyze").prop("disabled", false);
              $("#coverage-data").prop("disabled", false);
              processCoverage();
            }
          );
        } else {
          setTimeout(
            function() {
              checkCoverage(taskId);
            }, 1000
          );
        }
      }
    );
  }

  function scaleColorbar(points) {
    var cmap = createColormap({
      colormap: "viridis",
      nlevels: 72
    });
    var minRevisit;
    var maxRevisit;
    if($("#revisit-custom").is(":checked")) {
      minRevisit = $("#revisit-min").val()*60*60;
      maxRevisit = $("#revisit-max").val()*60*60;
    } else if(points.length == 0) {
      minRevisit = 0;
      maxRevisit = 1;
    } else {
      minRevisit = Math.min(
        ..._.map(
          _.filter(points, function(feature) {
            return !isNaN(feature.properties.revisit);
          }),
          function(feature) {
            return feature.properties.revisit;
          }
        )
      );
      maxRevisit = Math.max(
        ..._.map(
          _.filter(points, function(feature) {
            return !isNaN(feature.properties.revisit);
          }),
          function(feature) {
            return feature.properties.revisit;
          }
        )
      );
    }
    $("#revisit-colorbar").empty();
    $("#revisit-colorbar").append(
      $("<span></span>").text(
        Number.parseFloat(minRevisit/60/60).toFixed(1) + " hr"
      ).addClass("m-1")
    );
    _.forEach(_.clone(cmap).reverse(), function(color) {
      $("#revisit-colorbar").append(
        $("<span>&nbsp;</span>").css("background-color", color)
      ).addClass("m-1");
    });
    $("#revisit-colorbar").append(
      $("<span></span>").text(
        Number.parseFloat(maxRevisit/60/60).toFixed(1) + " hr"
      ).addClass("m-1")
    );
    return function(revisit) {
      if(revisit < minRevisit) {
        return cmap[0];
      } else if(revisit > maxRevisit) {
        return cmap[cmap.length-1];
      } else {
        return cmap[
          Math.round(
            (
              (maxRevisit - revisit)
              / (maxRevisit - minRevisit)
            )*(cmap.length-1)
          )
        ]
      }
    };
  }
  scaleColorbar([]);

  function processCoverage() {
    clearCoverage();
    $("#coverage-data").prop("disabled", false);
    var points = $("#coverage-analyze").data("points");
    var cells = $("#coverage-analyze").data("cells");
    var getColor = scaleColorbar(points);
    _.forEach(points, function(feature) {
      var color = getColor(feature.properties.revisit);
      coveragePoints.add({
        position: Cesium.Cartesian3.fromDegrees(
          feature.geometry.coordinates[0],
          feature.geometry.coordinates[1]
        ),
        color: Cesium.Color.fromCssColorString(
          color===undefined ? "#000000" : color
        )
      });
    });
    _.forEach(cells, function(feature) {
      var color = getColor(feature.properties.revisit);
      const vertices = _.map(
        feature.geometry.coordinates[0],
        function(coordinate){
          return Cesium.Cartesian3.fromDegrees(coordinate[0], coordinate[1])
        }
      );
      coverageCells.add(
        new Cesium.GroundPrimitive({
          geometryInstances: new Cesium.GeometryInstance({
            geometry: new Cesium.PolygonGeometry({
              polygonHierarchy: new Cesium.PolygonHierarchy(vertices),
              height: 0,
              arcType: getArcType()
            }),
            attributes: {
              color: Cesium.ColorGeometryInstanceAttribute.fromColor(
                Cesium.Color.fromCssColorString(
                  color===undefined ? "#000000" : color
                ).withAlpha(0.5)
              )
            }
          })
        })
      );
    });
  }
  function clearCoveragePoints() {
    coveragePoints.removeAll();
  }
  function clearCoverageCells() {
    coverageCells.removeAll();
  }

  function clearCoverage() {
    clearCoveragePoints();
    clearCoverageCells();
    $("#coverage-data").prop("disabled", true);
  }

  var coverageTable;
  $("#coverage-dialog").on("show.bs.modal", function() {
    var features = [];
    _.forEach($("#coverage-analyze").data("points"), function(feature) {
        features.push({
          pointId: feature.properties.point_id,
          latitude: feature.geometry.coordinates[1],
          longitude: feature.geometry.coordinates[0],
          revisitMean: feature.properties.revisit/60/60,
          accessMean: feature.properties.access,
          samples: feature.properties.samples
        });
    });

    if(coverageTable) {
      coverageTable.destroy();
    }
    coverageTable = $('#coverage-table').DataTable({
      'data': _.map(features, function(feature) {
        return [
          feature.pointId,
          feature.latitude,
          feature.longitude,
          feature.revisitMean,
          feature.accessMean,
          feature.samples
        ];
      }),
      'columns': [
        { title: 'Point ID' },
        {
          title: 'Latitude (deg)',
          render: $.fn.dataTable.render.number(',', '.', 6)
        },
        {
          title: 'Longitude (deg)',
          render: $.fn.dataTable.render.number(',', '.', 6)
        },
        {
          title: 'Mean Revisit (hr)',
          render: $.fn.dataTable.render.number(',', '.', 1)
        },
        { title: 'Mean Access (s)',
          render: $.fn.dataTable.render.number(',', '.', 1)
        },
        { title: '# Samples' }
      ],
      'order': [[ 0, 'asc']],
      'searching': false
    });

    $("#coverage-data-json").attr(
      "href",
      "data:text/json;charset=utf-8,"
      + encodeURIComponent(
        JSON.stringify(
          _.map(features, function(feature) {
            return {
              pointId: feature.point_id,
              latitude_deg: feature.latitude,
              longitude_deg: feature.longitude,
              revisitMean_hr: feature.revisitMean,
              accessMean_s: feature.accessMean,
              obsCount: feature.obsCount
            }
          })
        )
      )
    );
    $("#coverage-data-json").attr("download", "coverage.json");
    $("#coverage-data-csv").attr(
      "href",
      "data:text/csv;charset=utf-8,"
      + encodeURIComponent(
        [
          [
            "pointId",
            "latitude_deg",
            "longitude_deg",
            "revisitMean_hr",
            "accessMean_s",
            "obsCount"
          ].join()
        ].concat(
          _.map(features, function(feature) {
            return [
              feature.pointId,
              feature.latitude,
              feature.longitude,
              feature.revisitMean,
              feature.accessMean,
              feature.obsCount
            ].join();
          })
        ).join("\n")
      )
    );
    $("#coverage-data-csv").attr("download", "coverage.csv");
  });

  $("#revisit-custom").change(processCoverage);
  $("#revisit-min").change(processCoverage);
  $("#revisit-max").change(processCoverage);
  $("#coverage-analyze").click(collectCoverage);

  $("#region-type").change(clearCoverage);

  $("#nav-analysis-tab").on("hide.bs.tab", clearCoverage);
});
