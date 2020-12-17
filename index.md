---
title: Home
layout: template
filename: index
---

<head>
<style>
h1 {text-align: center;}
h2 {text-align: center;}

* {
  box-sizing: border-box;
}

img {
  vertical-align: middle;
}

/* Position the image container (needed to position the left and right arrows) */
.container {
  position: relative;
}

/* Hide the images by default */
.mySlides {
  display: none;
}

/* Add a pointer when hovering over the thumbnail images */
.cursor {
  cursor: pointer;
}

/* Next & previous buttons */
.prev,
.next {
  cursor: pointer;
  position: absolute;
  top: 40%;
  width: auto;
  padding: 16px;
  margin-top: -50px;
  color: white;
  background-color: rgba(0,0,0,0.4)
  font-weight: bold;
  font-size: 20px;
  border-radius: 0 3px 3px 0;
  user-select: none;
  -webkit-user-select: none;
}

/* Position the "next button" to the right */
.next {
  right: 0;
  border-radius: 3px 0 0 3px;
}

/* On hover, add a black background color with a little bit see-through */
.prev:hover,
.next:hover {
  background-color: rgba(0, 0, 0, 0.8);
}

/* Number text (1/3 etc) */
.numbertext {
  color: #f2f2f2;
  font-size: 12px;
  padding: 8px 12px;
  position: absolute;
  top: 0;
}

/* Container for image text */
.caption-container {
  text-align: center;
  background-color: #222;
  padding: 2px 16px;
  color: white;
}

.row:after {
  content: "";
  display: table;
  clear: both;
}

/* Six columns side by side */
.column {
  float: left;
  width: 16.66%;
}

/* Add a transparency effect for thumnbail images */
.demo {
  opacity: 0.6;
}

.active,
.demo:hover {
  opacity: 1;
}

</style>
</head>

<h1> Ab-initio transport </h1>
<h2> Simple and fast </h2>

Phoebe is an open-source code for the ab-initio computation of electron and phonon transport properties of crystals using HPC computers.

### Current functionalities

<div class="container">
  <div class="mySlides">
    <div class="numbertext">1 / 6</div>
    <img src="pictures/home/1.png" style="width:100%">
  </div>

  <div class="mySlides">
    <div class="numbertext">2 / 6</div>
    <img src="pictures/home/2.png" style="width:100%">
  </div>

  <div class="mySlides">
    <div class="numbertext">3 / 6</div>
    <img src="pictures/home/3.png" style="width:100%">
  </div>
    
  <div class="mySlides">
    <div class="numbertext">4 / 6</div>
    <img src="pictures/home/4.png" style="width:100%">
  </div>

  <div class="mySlides">
    <div class="numbertext">5 / 6</div>
    <img src="pictures/home/5.png" style="width:100%">
  </div>
    
  <div class="mySlides">
    <div class="numbertext">6 / 6</div>
    <img src="pictures/home/6.jpg" style="width:100%">
  </div>
    
  <a class="prev" style="text-decoration: none" onclick="plusSlides(-1)">❮</a>
  <a class="next" style="text-decoration: none" onclick="plusSlides(1)">❯</a>

  <div class="caption-container">
    <p id="caption"></p>
  </div>

  <div class="row">
    <div class="column">
      <img class="demo cursor" src="pictures/home/1.png" style="width:100%" onclick="currentSlide(1)" alt="Electronic conductivity limited by electron-phonon scattering">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/2.png" style="width:100%" onclick="currentSlide(2)" alt="Thermal conductivity with phonon-phonon scattering">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/3.png" style="width:100%" onclick="currentSlide(3)" alt="Electron and phonon viscosity">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/4.png" style="width:100%" onclick="currentSlide(4)" alt="Band Structure and density of states calculations">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/5.png" style="width:100%" onclick="currentSlide(5)" alt="GPU acceleration with Kokkos">
    </div>    
    <div class="column">
      <img class="demo cursor" src="pictures/home/6.jpg" style="width:100%" onclick="currentSlide(6)" alt="Support for ab-initio data from Quantum ESPRESSO">
    </div>
  </div>
</div>

In details

* Electronic transport coefficients (mobility, conductivity, thermal conductivity, and Seebeck coefficient);

* Electron-phonon scattering with Wannier interpolation;

* Electron-phonon scattering within the EPA approximation;

* Phonon thermal conductivity with 3-phonon scattering;

* Calculation of electron and phonon band structure

* Calculation of electron and phonon linewidths or relaxation times.

* Calculation of the electron or phonon density of states.

and more to come...



<script>
var slideIndex = 1;
showSlides(slideIndex);

function plusSlides(n) {
  showSlides(slideIndex += n);
}

function currentSlide(n) {
  showSlides(slideIndex = n);
}

function showSlides(n) {
  var i;
  var slides = document.getElementsByClassName("mySlides");
  var dots = document.getElementsByClassName("demo");
  var captionText = document.getElementById("caption");
  if (n > slides.length) {slideIndex = 1}
  if (n < 1) {slideIndex = slides.length}
  for (i = 0; i < slides.length; i++) {
      slides[i].style.display = "none";
  }
  for (i = 0; i < dots.length; i++) {
      dots[i].className = dots[i].className.replace(" active", "");
  }
  slides[slideIndex-1].style.display = "block";
  dots[slideIndex-1].className += " active";
  captionText.innerHTML = dots[slideIndex-1].alt;
}
</script>

