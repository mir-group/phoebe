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


### Functionalities


Current capabilities include:

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


<ul>
<li> Electron mobility</li>
<li> Phonon thermal conductivity</li>
<li> Electron and phonon viscosity</li>
<li> Wigner distribution corrections to Boltzmann equation</li>
<li> Density of states</li>
<li> Plot of band structure dispersion.</li>
<li> GPU acceleration</li>
<li> Powered by Quantum-ESPRESSO and Phonon3py</li>
</ul>


With more in development!

