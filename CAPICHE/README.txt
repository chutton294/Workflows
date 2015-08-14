Catchment Prediction in Changing Environments (CAPICHE)

Experiment Aims
Calibrate different hydrological models to consecutive time periods that cover a change in catchment forcing (e.g. land cover change) to understand how well the models simulate the change, and whether they simulate the change for the correct hydrological reasons. That is, to understand whether the correct model parameters change to simulate hydrological response to change, and compensatory parameter effects and interactions occur in the model simulations.

Qualitative workflow

A series of hydrological signatures are calculated in each time-window. The models are calibrated to these hydrologic signatures. Instead of applying an optimisation approach, and to account for data uncertainty, an interval is set around each hydrologic signature. When combining these signatures, the intervals combine to produce a hyper-volume in multidimensional signature space. A multi-objective optimisation algorithm is then iteratively applied at multiple locations in order to populate the signature space and identify behavioural models.
Once the behavioural models are identified, sensitivity analysis is conducted to evaluate parameter sensitivity across the different time-periods to see whether the correct model parameters change to simulate change. 

Practical workflow

The qualitative workflow has been executed with a series of separate R-scripts. See relevant scripts for required libraries and 
source files.

window_signatures.R: with provided empirical data, calculates hydrological signatures for each calibration window.

    Inputs: empirical catchment data (see script for format).
    Outputs: window_signatures.txt
    
window_optimisation.R: with provided empiirical data, and window_signatures.txt, the code identifies a set of behavioural models 
in each calibration window

    Inputs: empirical catchment data (see script for format); window_signatures.txt
    Outputs: behavioural_models.txt
    
window_sensitivity_analysis.R: produces sensitivity analysis for different windows to evaluate model performance

   Inputs: behavioural_models.txt
   Outputs: graphs of model sensitivity analysis



