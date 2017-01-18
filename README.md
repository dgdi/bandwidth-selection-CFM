# bandwidth-selection-CFM## bandwidth-selection-CFM


***
This repo provides the R functions to evaluate LSCV errors; the first step in a "Conditional Frontier Model" analysis.
The whole "Conditional Frontier Model" methodology has been developed and used in the fourth Chapter of Duccio Gamannossi degl'Innocenti PhD thesis:

Italian municipalities efﬁciency: A conditional frontier model approach.

(A short version of the chapter is about to be submitted to the Journal of Productivity Analysis.)

However, the whole methodology is not publicly available yet. 
For informations or requests write to duccio.gama@gmail.com
 
For Further informations on Conditional Frontier Models see the seminal papers:
* [How to measure the impact of environmental factors in a nonparametric production model](http://www.sciencedirect.com/science/article/pii/S0377221712004833) by Luiza Bădin, Cinzia Daraio, Léopold Simar.
* [Introducing Environmental Variables in Nonparametric Frontier Models: a Probabilistic Approach](http://link.springer.com/article/10.1007/s11123-005-3042-8#page-1) by Cinzia Daraio and Léopold Simar.

***

This Repo provides the functions: 

Ker_LSCV_IN.R 
Ker_LSCV_OUT.R

That implement in the R language the MATLAB function:

Ker_LSCV_OUT.m

Presented in:

* [Optimal bandwidth selection for conditional efficiency measures: A data-driven approach](http://www.sciencedirect.com/science/article/pii/S0377221709002148) by Luiza Bădin, Cinzia Daraio, Léopold Simar.

The R functions proposed extend the MATLAB one allowing to evaluate the bandwidths' LSCV errors also in an input oriented setting.

A thorough validation is provided with respect to simulated datasets reproducing the DGPs in:

* [Optimal bandwidth selection for conditional efficiency measures: A data-driven approach](http://www.sciencedirect.com/science/article/pii/S0377221709002148) by Luiza Bădin, Cinzia Daraio, Léopold Simar.

* [“Two-stage DEA: caveat emptor](http://link.springer.com/article/10.1007/s11123-011-0230-6) by Léopold Simar.

The *reproducible* validation reports ( are:

*validation_assessment - EJOR 2010.pdf* 
*validation_assessment JPA 2011.pdf*

The code to perform the validation is available in Bandwidth Selection\validation\ 

