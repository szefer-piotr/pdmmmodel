This zip file bundles R code and input data to reconstruct abstract
trophic traits and trophic niche space of based on trait and diet
data, as demonstrated by Nagelkerke and Rossberg ("Trophic niche-space
imaging, using resource and consumer traits", Theoretical Ecology,
2014).

It contains the following files:

README.txt                      This file
TraitSelection.R                Batch code for trait selection
ModelFit.R                      Interactive code to analyze model fit
functions.R                     Helper functions                        
nichespace.R                    Functions for fitting quadratic form
Tana_Consumer_Traits.dat        Input data
Tana_Interaction_Strength.dat   Input data
Tana_Resource_Traits.dat        Input data

To run the batch code TraitSelection.R, first adjust the value of
target.dim at the top of the file, and then run the code in R,
e.g. using the command line

     R --no-save < TraitSelection.R

This will run for several hours, and leave a log file.  In the log
file, the last line in the file starting with "XB" gives you the
selected trait combination.

ModelFit.R contains interactive code snippets, to be evaluated,
e.g. with RStudio.  Please adjust the assignments to the variables
dims.to.keep, r.list, c.list according to the desired number of
nichespace dimensions and your selected resource and consumer traits.

The five other files are read by TraitSelection.R and ModelFit.R as
they are needed.

Feedback welcome!

Axel Rossberg and Leo Nagelkerke
