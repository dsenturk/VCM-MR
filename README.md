CONTENTS OF THIS FOLDER 
——————————————

* VCMMR_tutorial.R : A step-by-step implementation of VCM-MR and the associated procedures described in "Modeling Time-Varying Effects of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis".  

* VCMMR_estimation.R : Function for estimation of VCM-MR model parameters described in "Modeling Time-Varying Effects of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis", including estimation of facility-specific fixed effects, time-varying effects of multilevel risk factors, calendar year effect and variance of subject specific random effects.

* VCMMR_simulation.R : Function for simulating one data set under the simulation design described in Web Appendix E of the supplementary materials.

* VCMMR_bootstrap.R : Function for performing bootstrap inference for effects of multilevel factors and calendar year effects as described in "Modeling Time-Varying Effects of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis".

* VCMMR_hypothesistest.R : Function for performing bootstrap hypothesis testing for the significance of the facility-specific fixed effects as described in Web Appendix B of "Modeling Time-Varying Effects of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis". 


INTRODUCTION
——————————————

The contents of this folder allow for implementation of the VCM-MR estimation and inference described in "Modeling Time-Varying Effects of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis". Users can simulate a sample data frame (VCMMR_simulation.R) and apply the proposed estimation algorithm (VCMMR_estimation.R). Also, we include tools to perform inference on multilevel risk factors via a bootstrap procedure (VCMMR_bootstrap.R), allowing users to test whether the effect of a patient-level or a facility-level risk factor is significant. Further, users can test the significance of the facility-specific fixed effects via the proposed hypothesis testing procedure (VCMMR_hypothesistest.R). Detailed instructions on how to perform the aforementioned procedures, make predictions of patient- and facility-level risk trajectories and visualize results are included in VCMMR_tutorial.R.

REQUIREMENTS 
——————————————

The included R programs require R 3.3.3 (R Core Team, 2017) and the packages listed in VCMMR_tutorial.R.

INSTALLATION
——————————————

Load the R program files into the global environment and install required packages using commands in VCMMR_tutorial.R
