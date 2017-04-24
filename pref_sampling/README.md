pref_sampling
===========

Code to implement preferential sampling analyses described in "Confronting preferential sampling in wildlife surveys: diagnosis and model-based triage" by Conn, Thorson, and Johnson using Template Model Builder (TMB).

Example simulation code is provided in ./inst/run_sim_example.R, but only one such simulation is conducted (the full suite of simulations took several days).  

The script ./inst/run_bearded_models.R runs the bearded seal analysis, while the script ./inst/run_bearded_models_cross_val.R computes cross validation metrics.  

TMB code for conducting analyses is in the ./src directory.





