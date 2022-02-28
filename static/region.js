var regionPoints = [];

function getRegionCoordinates() {
  var coordinates = _.map(
    // close polygon
    regionPoints.concat(regionPoints[0]),
    function(point) {
      const ellipsoid = viewer.scene.globe.ellipsoid;
      const cartographic = Cesium.Cartographic.fromCartesian(point, ellipsoid);
      return [cartographic.longitude*180/Math.PI, cartographic.latitude*180/Math.PI];
    }
  );
  return coordinates;
}

function getRegionMask() {
  if($("#region-type").val()=="custom") {
    return {
      coordinates: [getRegionCoordinates()],
      type: "Polygon"
    }
  } else if($("#region-type").val()=="latitude") {
    return {
      coordinates: [[
        [-180, $("#region-latitude-min").val()],
        [180, $("#region-latitude-min").val()],
        [180, $("#region-latitude-max").val()],
        [-180, $("#region-latitude-max").val()],
        [-180, $("#region-latitude-min").val()]
      ]],
      type: "Polygon"
    }
  } else if($("#region-type").val()=="latlon") {
    return {
      coordinates: [[
        [$("#region-longitude-min").val(), $("#region-latitude-min").val()],
        [$("#region-longitude-max").val(), $("#region-latitude-min").val()],
        [$("#region-longitude-max").val(), $("#region-latitude-max").val()],
        [$("#region-longitude-min").val(), $("#region-latitude-max").val()],
        [$("#region-longitude-min").val(), $("#region-latitude-min").val()]
      ]],
      type: "Polygon"
    }
  } else if($("#region-type").val()=="global") {
    return null;
  } else {
    return $("#region-type").val();
  }
}

$(document).ready(function() {
  const regionEntity = viewer.entities.add({
    polygon: {
      hierarchy: regionPoints,
      arcType: Cesium.ArcType.RHUMB,
      height: 0,
      material: Cesium.Color.WHITE.withAlpha(0.1),
      outline: true,
      outlineColor: Cesium.Color.BLACK
    },
    show: !$("#region-type").val()=="global"
  });

  function updateCoordinates() {
    $("#region-coordinates").val(
      JSON.stringify(
        _.flatten(
          _.map(regionPoints, function(point) {
            const ellipsoid = viewer.scene.globe.ellipsoid;
            const cartographic = Cesium.Cartographic.fromCartesian(point, ellipsoid);
            return [
              parseFloat((cartographic.longitude*180/Math.PI).toFixed(6)),
              parseFloat((cartographic.latitude*180/Math.PI).toFixed(6))
            ]
          })
        )
      )
    );
  }

  var pointsReady;
  var cellsReady;
  function generateRegion() {
    $("#region-generate").prop("disabled", true);
    $("#region-generate .spinner-border").show();
    pointsReady = false;
    cellsReady = false;
    generatePoints();
    generateCells();
  }

  function checkRegionReady(ready) {
    if(ready.points) {
      pointsReady = true;
    }
    if(ready.cells) {
      cellsReady = true;
    }
    if(pointsReady && cellsReady) {
      $("#region-generate .spinner-border").hide();
      $("#region-generate").prop("disabled", false);
    }
  }

  function clearRegion() {
    clearPoints();
    clearCells();
  }

  var editRegionMode = false;

  function editRegionCoordinates() {
    if(editRegionMode) {
      editRegionMode = false;
      $("#region-coordinates").prop("disabled", true);
      $("#edit-coordinates").children("i").removeClass("fa-check").addClass("fa-edit");
      $("#edit-coordinates").removeClass("btn-success").addClass("btn-secondary");
      try {
        regionPoints = Cesium.Cartesian3.fromDegreesArray(
          JSON.parse($("#region-coordinates").val())
        );
        regionEntity.polygon.hierarchy = regionPoints;
        generateRegion();
      } catch(e) {
        updateCoordinates();
      }
    } else {
      editRegionMode = true;
      $("#region-coordinates").prop("disabled", false);
      $("#edit-coordinates").children("i").removeClass("fa-edit").addClass("fa-check");
      $("#edit-coordinates").removeClass("btn-secondary").addClass("btn-success");
    }
  };

  function editRegionListener(e) {
    e.preventDefault();
    const mousePosition = new Cesium.Cartesian2(e.layerX, e.layerY);
    const ellipsoid = viewer.scene.globe.ellipsoid;
    const cartesian = viewer.camera.pickEllipsoid(mousePosition, ellipsoid);
    if (cartesian) {
      regionPoints.push(cartesian);
      regionEntity.polygon.hierarchy = regionPoints;
      updateCoordinates();
    }
  };
  function selectRegionCoordinates() {
    if(editRegionMode) {
      editRegionMode = false;
      $("#select-coordinates").children("i").removeClass("fa-check").addClass("fa-map-marker");
      $("#select-coordinates").removeClass("btn-success").addClass("btn-secondary");
      document.body.style.cursor = 'default';
      viewer.canvas.removeEventListener('click', editRegionListener);
      generateRegion();
    } else {
      editRegionMode = true;
      clearRegion();
      $("#select-coordinates").children("i").removeClass("fa-map-marker").addClass("fa-check");
      $("#select-coordinates").removeClass("btn-secondary").addClass("btn-success");
      document.body.style.cursor = 'crosshair';
      regionPoints = [];
      regionEntity.polygon.hierarchy = regionPoints;
      viewer.canvas.addEventListener('click', editRegionListener);
    }
  }

  function changeRegionType() {
    regionEntity.show = $("#region-type").val()=="custom";
    if($("#region-type").val()=="global") {
      $("#cell-strips").val("none");
    }
    if($("#region-type").val()=="latitude" || $("#region-type").val()=="latlon") {
      $("#region-latitude-group").show();
      if($("#region-type").val()=="latlon") {
        $("#region-longitude-group").show();
      } else {
        $("#region-longitude-group").hide();
      }
    } else {
      $("#region-latitude-group").hide();
      $("#region-longitude-group").hide();
    }
    if($("#region-type").val()=="custom") {
      $("#region-coordinates-group").show();
    } else {
      $("#region-coordinates-group").hide();
    }
    $("#cell-strips").prop("disabled", $("#region-type").val()!="custom");
    $("#select-coordinates").prop("disabled", $("#region-type").val()!="custom");
    $("#edit-coordinates").prop("disabled", $("#region-type").val()!="custom");
    $("#region-latitude-min").prop("disabled", $("#region-type").val()!="latitude" && $("#region-type").val()!="latlon");
    $("#region-latitude-max").prop("disabled", $("#region-type").val()!="latitude" && $("#region-type").val()!="latlon");
    $("#region-longitude-min").prop("disabled", $("#region-type").val()!="latlon");
    $("#region-longitude-max").prop("disabled", $("#region-type").val()!="latlon");
  }

  function changeRegionLatitude() {
    $("#region-latitude-min").attr("max", $("#region-latitude-max").val());
    $("#region-latitude-max").attr("min", $("#region-latitude-min").val());
  };
  function changeRegionLongitude() {
    $("#region-longitude-min").attr("max", $("#region-longitude-max").val());
    $("#region-longitude-max").attr("min", $("#region-longitude-min").val());
  };

  const points = viewer.scene.primitives.add(new Cesium.PointPrimitiveCollection());
  function generatePoints() {
    clearPoints();
    points.show = true;
    $.ajax({
      url: "/generate/points",
      type: "POST",
      contentType: "application/json",
      dataType: "json",
      data: JSON.stringify({
          method: $("#point-method").val(),
          distance: $("#point-distance").val()*1000,
          mask: getRegionMask()
      }),
      success: function(response) {
        _.forEach(response.features, function(feature) {
          points.add({
              position: Cesium.Cartesian3.fromDegrees(
                feature.geometry.coordinates[0],
                feature.geometry.coordinates[1]
              ),
              color: Cesium.Color.BLUE
          });
        });
        checkRegionReady({points: true});
      }
    });
  };
  function clearPoints() {
    points.removeAll();
  }

  const cells = viewer.scene.primitives.add(new Cesium.PrimitiveCollection());
  const cellOutlines = viewer.scene.primitives.add(new Cesium.PolylineCollection());
  function generateCells() {
    clearCells();
    cells.show = true;
    $.ajax({
      url: "/generate/cells",
      type: "POST",
      contentType: "application/json",
      dataType: "json",
      data: JSON.stringify({
          method: $("#cell-method").val(),
          distance: $("#cell-distance").val()*1000,
          mask: getRegionMask(),
          strips: $("#cell-strips").val() == "none" ? null : $("#cell-strips").val()
      }),
      success: function(response) {
        _.forEach(response.features, function(feature) {
          const vertices = _.map(
            feature.geometry.coordinates[0],
            function(coordinate){
              return Cesium.Cartesian3.fromDegrees(coordinate[0], coordinate[1])
            }
          );
          cells.add(
            new Cesium.GroundPrimitive({
              geometryInstances: new Cesium.GeometryInstance({
                geometry: new Cesium.PolygonGeometry({
                  polygonHierarchy: new Cesium.PolygonHierarchy(vertices),
                  height: 0,
                  arcType: getArcType()
                }),
                attributes: {
                  color: Cesium.ColorGeometryInstanceAttribute.fromColor(Cesium.Color.BLUE.withAlpha(0.2))
                }
              })
            })
          );
          cellOutlines.add({
            positions: vertices,
            width: 1.0,
            loop: true,
            material: new Cesium.Material.fromType('Color', { color: Cesium.Color.BLACK })
          });
        });
        checkRegionReady({cells: true});
      }
    });
  }
  function clearCells() {
    cells.removeAll();
    cellOutlines.removeAll();
  }

  $("#region-type").change(changeRegionType);
  $("#region-latitude-min").change(changeRegionLatitude);
  $("#region-latitude-max").change(changeRegionLatitude);
  $("#region-longitude-min").change(changeRegionLongitude);
  $("#region-longitude-max").change(changeRegionLongitude);
  $("#edit-coordinates").click(editRegionCoordinates);
  $("#select-coordinates").click(selectRegionCoordinates);

  $("#coverage-analyze").click(function() {
    points.show = false;
    cells.show = false;
  });

  $("#nav-analysis-tab").on("show.bs.tab", changeRegionType);
  $("#nav-analysis-tab").on("show.bs.tab", generateRegion);
  $("#nav-analysis-tab").on("hide.bs.tab", clearRegion);
  $("#region-dialog").on("hide.bs.modal", generateRegion);
  $("#region-type").change(generateRegion);

  // set initial custom region coordinates
  regionPoints = Cesium.Cartesian3.fromDegreesArray(
    [
      -130, 20,
      -130, 50,
      -60, 50,
      -60, 20
    ]
  );
  regionEntity.polygon.hierarchy = regionPoints;
  updateCoordinates();
});
