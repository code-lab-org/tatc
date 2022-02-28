// Initialize the Cesium Viewer in the HTML element with the `cesiumContainer` ID.
var viewer = new Cesium.Viewer('cesiumContainer', {
  baseLayerPicker: false,
  homeButton: false,
  infoBox: false,
  geocoder: false,
  selectionIndicator: false,
  navigationHelpButton: false,
  navigationInstructionsInitiallyVisible: false,
  timeline: true,
});

function clearEntities(parent) {
  var entitiesToRemove = [];
  _.forEach(viewer.entities.values, function(entity) {
    if(entity.parent === parent) {
      entitiesToRemove.push(entity);
    }
  });
  _.forEach(entitiesToRemove, function(entity) {
    viewer.entities.remove(entity);
  });
}

function getArcType() {
  if($("#region-type").val()=="global"
      || ($("#region-type").val()=="latitude"
          && ($("#region-latitude-min").val()<=-90
            || $("#region-latitude-min").val()>=90))) {
    return Cesium.ArcType.GEODESIC;
  } else {
    return Cesium.ArcType.RHUMB
  }
}

function initializePage(user) {
  $("#login-modal").modal("hide");
  $("#login-link").hide();
  $("#logout-link").show();
  $("#user-email").text(user.email);
  $("#nav-info-tab").off("hide.bs.tab");
}
function getUserInformation() {
  // perform get request
  $.get("/users/me")
  .done(initializePage)
  .fail(function() {
    // if failed, show the login link
    $("#login-link").show();
    $("#logout-link").hide();
    $("#user-email").text("");
    // configure login modal to show when navigating from home tab
    $('#nav-info-tab').on('hide.bs.tab', function (e) {
      e.preventDefault();
      $("#login-modal").modal("show");
    });
  });
}
function initializeCesium() {
  // perform get request
  $.get("/cesium/token")
  .done(function(data) {
    Cesium.Ion.defaultAccessToken = data;
    viewer.entities.removeAll();
    viewer.destroy();
    // Initialize the Cesium Viewer in the HTML element with the `cesiumContainer` ID.
    viewer = new Cesium.Viewer('cesiumContainer', {
      terrainProvider: Cesium.createWorldTerrain(),
      baseLayerPicker: false,
      homeButton: false,
      infoBox: false,
      geocoder: false,
      selectionIndicator: false,
      navigationHelpButton: false,
      navigationInstructionsInitiallyVisible: false,
      timeline: true,
      imageryProvider: new Cesium.IonImageryProvider({ assetId: 3845 })
    });
  });
}

$(document).ready(function() {
  window.Split(['#editor', '#viewer'], { sizes: [25, 75] });
  // behavior when the login form is submitted
  $("#login-form").on("submit", function(e) {
    // reset error messages
    $("#login-email-message").text("");
    $("#login-password-message").text("");
    // try to login
    $.post("/login", $(this).serialize())
    .done(function() {
      getUserInformation();
      initializeCesium();
    })
    .fail(function() {
      // try to register a new account
      $.ajax({
        type: "POST",
        url: "/register",
        data: JSON.stringify({
          email: $("#login-email").val(),
          password: $("#login-password").val(),
          passcode: $("#login-password").val()
        }),
        contentType: "application/json"
      }).done(function() {
        // if successful, login again
        $("#login-form").submit();
      }).fail(function(xhr, text, error) {
        // if failed, display error message
        if(xhr.status == 422) {
          $("#login-email-message").text("Invalid email format.");
        } else if(xhr.status == 400) {
          $("#login-password-message").text("Incorrect password.");
        }
      });
    });
    return false;
  });
  // behavior when the logout link is clicked
  $("#logout-link").click(function() {
    $.post("/logout").done(function(){
      // show home tab and reset user information
      $('#nav-info-tab').tab("show");
      getUserInformation();
    });
  });
  // configure login modal dialog
  $("#login-modal").modal({
    backdrop: "static",
    keyboard: false,
    show: false
  });
  // check if the user is logged in
  getUserInformation();
});
