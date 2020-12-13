---
title: Home
layout: template
filename: index
---

<head>
<style>
h1 {text-align: center;}
h2 {text-align: center;}

/*  SECTIONS  */
.section {
	clear: both;
	padding: 0px;
	margin: 0px;
}

/*  COLUMN SETUP  */
.col {
	display: block;
	float:left;
	margin: 1% 0 1% 1.6%;
}
.col:first-child { margin-left: 0; }

/*  GROUPING  */
.group:before,
.group:after { content:""; display:table; }
.group:after { clear:both;}
.group { zoom:1; /* For IE 6/7 */ }

/*  GRID OF TWO  */
.span_2_of_2 {
	width: 100%;
}
.span_1_of_2 {
	width: 49.2%;
}

/*  GO FULL WIDTH AT LESS THAN 480 PIXELS */

@media only screen and (max-width: 480px) {
	.col { 
		margin: 1% 0 1% 0%;
	}
}

@media only screen and (max-width: 480px) {
	.span_2_of_2, .span_1_of_2 { width: 100%; }
}


</style>
</head>

<h1> Ab-initio transport </h1>
<h2> Simple and fast </h2>

Phoebe is a code for the ab-initio compoutation of electron and phonon transport properties of crystals using HPC computers.


### Current functionalities

<div class="section group">
  <div class="col span_1_of_2">
    Electron mobility
  </div>
  <div class="col span_1_of_2">
    Phonon thermal conductivity
  </div>  
</div>

<div class="section group">
  <div class="col span_1_of_2">
    Electron and phonon Viscosity
  </div>
  <div class="col span_1_of_2">
    Wigner distribution
  </div>  
</div>

<div class="section group">
  <div class="col span_1_of_2">
    Density of States
  </div>
  <div class="col span_1_of_2">
    Band Structure
  </div>  
</div>

<div class="section group">
  <div class="col span_1_of_2">
     <img src="pictures/home/gpu.jpg" alt="GPU" style="width:100%">
    GPU acceleration
  </div>
  <div class="col span_1_of_2">
     <img src="pictures/home/qe.png" alt="QE" style="width:100%">
     Quantum-ESPRESSO + Phono3py
  </div>  
</div>

With more in development!

