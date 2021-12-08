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

/* Next & previous buttons to the left/right of picture */
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
  background-color: #333;
}

/* background color without hover with a little bit see-through */
.prev, .next {
  background-color: #555;
  opacity: 0.7;
}

/* Number text (1/3 etc) */
.numbertext {
  color: #f2f2f2;
  font-size: 12px;
  padding: 8px 12px;
  position: absolute;
  top: 0;
}

/* Container for caption below the images */
.caption-container {
  text-align: center;
  background-color: #555;
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

<h1> Ab-initio transport, simple and fast.</h1>
<p style="text-align:center">
Phoebe is an open-source code for the ab-initio computation of electron and phonon transport properties of crystalline materials.
</p>

<p style="text-align:center;">
It is designed to take advantage of HPC systems via MPI-OpenMP hybrid parallelism, memory-distributed computing via ScaLAPACK, and GPU accelerated calculation of scattering rates.
</p>

<div class="container" style="margin-top:2em;">
  <div class="mySlides">
    <div class="numbertext">1 / 6</div>
    <img src="pictures/home/colorPhdisp.png" style="max-height:450px; margin-left:auto; margin-right:auto; display: block">
  </div>

  <div class="mySlides">
    <div class="numbertext">2 / 6</div>
    <img src="pictures/home/wigner.png" style="max-height:450px; margin-left:auto; margin-right:auto; display: block">
  </div>

  <div class="mySlides">
    <div class="numbertext">3 / 6</div>
    <img src="pictures/home/el.bands.tau.png" style="max-height:450px; margin-left:auto; margin-right:auto; display: block">
  </div>

  <div class="mySlides">
    <div class="numbertext">4 / 6</div>
    <img src="pictures/home/3.png" style="max-height:450px; margin-left:auto; margin-right:auto; display: block">
  </div>

  <div class="mySlides">
    <div class="numbertext">5 / 6</div>
    <img src="pictures/home/gan_rta_ph_relaxation_times.png" style="max-height:450px; margin-left:auto; margin-right:auto; display: block">
  </div>

  <div class="mySlides">
    <div class="numbertext">6 / 6</div>
    <img src="pictures/home/gan-thermal-conductivity.png" style="max-height:450px; margin-left:auto; margin-right:auto; display: block">
  </div>

  <a class="prev" style="text-decoration: none" onclick="plusSlides(-1)">❮</a>
  <a class="next" style="text-decoration: none" onclick="plusSlides(1)">❯</a>

  <div class="caption-container">
    <p id="caption"></p>
  </div>

  <div class="row">
    <div class="column">
      <img class="demo cursor" src="pictures/home/colorPhdisp.png" style="max-height:120px" onclick="currentSlide(1)" alt="Phonon linewidths projected onto dispersion of silicon">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/wigner.png" style="max-height:120px" onclick="currentSlide(2)" alt="Electron-phonon limited conductivity with the Wigner transport equation">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/el.bands.tau.png" style="max-height:120px" onclick="currentSlide(3)" alt="Electron-phonon linewidths for GaN along its band structure">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/3.png" style="max-height:120px" onclick="currentSlide(4)" alt="Resistivity of n-type (CoSb3) skutterudite, for which the RTA (dark blue triangles) and iterative (light blue 'x's) perform similarly.">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/gan_rta_ph_relaxation_times.png" style="max-height:120px" onclick="currentSlide(5)" alt="RTA phonon lifetimes of GaN, where colors from blue to green represent different phonon branches of increasing energy.">
    </div>
    <div class="column">
      <img class="demo cursor" src="pictures/home/gan-thermal-conductivity.png" style="max-height:120px" onclick="currentSlide(6)" alt="Lattice thermal conductivity of GaN along the 'a' crystal axis using the RTA, variational, and relaxons solvers.">
    </div>
  </div>
</div>

<br>

<h2 style="text-align:left; padding-bottom:0.5em; border-bottom:solid"> Current functionalities</h2>

<h3 style="margin:1em 0 1em 0;">Electronic Transport</h3>
  <ul style="padding-left:2em;">
    <li style="margin: 0 0 0.5em 0;"> Electron-phonon scattering by Wannier interpolation</li>
    <li style="margin: 0 0 0.5em 0;"> Electron-phonon scattering within the electron-phonon averaged (EPA) approximation</li>
    <li style="margin: 0 0 0.5em 0;"> Polar correction and boundary scattering contributions to transport</li>
    <li style="margin: 0 0 0em 0;"> Electronic transport coefficients (mobility, conductivity, thermal conductivity, and Seebeck coefficient)</li>
  </ul>
<h3 style="margin:1em 0 1em 0;">Phonon Transport</h3>
  <ul style="padding-left:2em;">
    <li style="margin: 0 0 0.5em 0;"> 3-phonon scattering from thirdOrder.py/ShengBTE or Phono3py force constants </li>
    <li style="margin: 0 0 0.5em 0;"> Boundary and isotope scattering contributions to transport</li>
    <li style="margin: 0 0 0em 0;"> Phonon (lattice) thermal conductivity</li>
  </ul>
<h3 style="margin:1em 0 1em 0;">And more...</h3>
  <ul style="padding-left:2em;">
    <li style="margin: 0 0 0.5em 0;"> BTE solutions by RTA, iterative, variational, and relaxons solvers</li>
    <li style="margin: 0 0 0.5em 0;"> Calculation of electron and phonon linewidths or relaxation times on a path</li>
    <li style="margin: 0 0 0.5em 0;"> Wigner transport equation correction for electrons and phonons (Zener tunneling contribution to electron transport)</li>
    <li style="margin: 0 0 0.5em 0;"> Hydrodynamic transport properties (viscosity) for electrons and phonons</li>
  </ul>

<h3 style="margin:1em 0 1em 0;"> For the full details of the Phoebe 1.0 release, see:</h3>

Phoebe: a collection of Phonon and Electron Boltzmann Equation solvers.<br>
A. Cepellotti, J. Coulter, A. Johansson, N. S. Fedorova, B. Kozinsky.<br>
[arXiv:2111.14999](https://arxiv.org/abs/2111.14999). (2021).


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

