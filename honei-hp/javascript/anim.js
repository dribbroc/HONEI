$(document).ready(function() {
  function animateUp() {
    $("#logo").animate({"bottom": "+=3px"}, "slow", animateDown);
  }
  function animateDown() {
    $("#logo").animate({"bottom": "-=3px"}, "slow", animateUp);
  }

  animateUp();
});
