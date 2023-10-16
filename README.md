README

The following readme describes the structure of the R code used for the manuscript titled "Using joint species distribution modelling to predict distributions of seafloor taxa and identify vulnerable marine ecosystems in New Zealand waters"

authored by Fabrice Stephenson, David A Bowden, Ashley A Rowden, Owen F Anderson, Malcolm R Clark, Matthew Bennion, Brittany Finucci, Matt H Pinkerton, Savannah Goode, Caroline Chin, Niki Davey, Alan Hart, Rob Stewart. 

R CODE:

Workflow of JSDM model fitting and spatial predictions

1. S1_read_data_ZBD201901.R  - Read data and prepare model objects (based on code from Ovaskainen and Abrego 2020).
2. S2_define_models_ZBD201901.R - Code for defining the HSMC models (presence-absense and abundance conditional on presence) for 67 invertebrate taxa (based on code from Ovaskainen and Abrego 2020).
3. S3_fit_models_ZBD201901.R - Code for fitting the HSMC models (based on code from Ovaskainen and Abrego 2020).
4. S4_evaluate_convergence_ZBD201901.R - Code for evaluating convergence of the HSMC models (based on code from Ovaskainen and Abrego 2020).
5. S5_compute_model_fit_ZBD201901.R - Code for assessing HSMC model fits using training data and 10-fold cross validated data (based on code from Ovaskainen and Abrego 2020).
6. S6_show_model_fit_ZBD201901.R - Code for visualising and saving HSMC model fits using training data and 10-fold cross validated data (based on code from Ovaskainen and Abrego 2020).
7. S7_show_parameter_estimates_ZBD201901.R - Code for visualising and saving HSMC model parameters (based on code from Ovaskainen and Abrego 2020).
8. S8_HiRes_predictions_ZBD201901.R - Code for generating spatial predictions of occurrence and abundance, and richness for 67 taxa.
9. S9_Mapping_and_VME_distributions_ZBD201901.R - Code for generating spatial estimates of VME indices and estimating most likley location of VME.

