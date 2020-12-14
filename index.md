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

<h1>Website under construction!</h1>

<h1> Ab-initio transport </h1>
<h2> Simple and fast </h2>

Phoebe is an open-source code for the ab-initio computation of electron and phonon transport properties of crystals using HPC computers.


### Current functionalities



GPU acceleration
    Electron mobility
    Phonon thermal conductivity
    Electron and phonon Viscosity
    Wigner distribution
    Density of States
    Band Structure
    Quantum-ESPRESSO + Phono3py

With more in development!


<div class="container">
  <div class="mySlides">
    <div class="numbertext">1 / 8</div>
    <img src="pictures/home/1.jpg" style="width:100%">
  </div>

  <div class="mySlides">
    <div class="numbertext">2 / 8</div>
    <img src="pictures/home/2.png" style="width:100%">
  </div>

  <div class="mySlides">
    <div class="numbertext">3 / 8</div>
    <img src="pictures/home/3.png" style="width:100%">
  </div>
    
  <div class="mySlides">
    <div class="numbertext">4 / 8</div>
    <img src="pictures/home/4.png" style="width:100%">
  </div>

  <div class="mySlides">
    <div class="numbertext">5 / 8</div>
    <img src="pictures/home/5.png" style="width:100%">
  </div>
    
  <div class="mySlides">
    <div class="numbertext">6 / 8</div>
    <img src="pictures/home/6.png" style="width:100%">
  </div>
    
  <div class="mySlides">
    <div class="numbertext">7 / 8</div>
    <img src="pictures/home/7.jpg" style="width:100%">
  </div>
    
  <div class="mySlides">
    <div class="numbertext">8 / 8</div>
    <img src="pictures/home/8.jpg" style="width:100%">
  </div>
    
  <a class="prev" onclick="plusSlides(-1)">❮</a>
  <a class="next" onclick="plusSlides(1)">❯</a>

  <div class="caption-container">
    <p id="caption"></p>
  </div>

  <div class="row">
    <div class="column">
      <img class="demo cursor" src="pictures/home/1.png" style="width:100%" onclick="currentSlide(1)" alt="Electron mobility">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/2.png" style="width:100%" onclick="currentSlide(2)" alt="Thermal conductivity">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/3.png" style="width:100%" onclick="currentSlide(3)" alt="Quasiparticle viscosity">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/4.png" style="width:100%" onclick="currentSlide(4)" alt="Wigner distribution">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/5.png" style="width:100%" onclick="currentSlide(5)" alt="Density of States">
    </div>    
    <div class="column">
      <img class="demo cursor" src="pictures/home/6.png" style="width:100%" onclick="currentSlide(6)" alt="Band structure">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/7.jpg" style="width:100%" onclick="currentSlide(7)" alt="GPU acceleration">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/8.jpg" style="width:100%" onclick="currentSlide(8)" alt="Quantum ESPRESSO">
    </div>
  </div>
</div>

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

