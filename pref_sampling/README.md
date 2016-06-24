pref_sampling
===========

Code to implement preferential sampling analyses described in "Confronting preferential sampling in wildlife surveys: diagnosis and model-based triage" by Conn, Thorson, and Johnson using Template Model Builder (TMB).

Example simulation code is provided in ./inst/run_sim_example.R, but only one such simulation is conducted (the full suite of simulations took several days).  The script ./inst/plot_sim_example then constructs plots from this simulation.

To run the bearded seal example, two scripts were used.  The main results reported were from ./inst/run_bearded_models.R; however, one of the models produced high predictions in southern portions of the study area where there was no ice.  This model was rerun in run_bearded_models_ice0.R with zero values in the southern portion of the landscape to remove this feature when providing an apparent abundance estimate.  

TMB code for conducting analyses is in the ./inst/TMB_version directory.





