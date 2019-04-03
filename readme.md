Functions and implementations of the wavelet-based verification tools presented in the paper "Assessment of wavelet-based spatial verification by means of a stochastic precipitation model" by Sebastian Buschow, Jakiw Pidstrigach and Petra Friederichs.

Unless noted otherwise, all code was written by Sebastian Buschow.

---
### Primary functions:
**`verification_functions.r`**
Contains the central function `fld2S()` which performs the stationary wavelet transform, takes care of boundaries, etc. also load all necessary packages and contains a few utility functions which are called by `fld2S()`.

**`wavelet_functions.r`**
Contains an altered version of the RDWT-routine cddews from the `LS2W`-package, which does not re-compute the correction matrix at each function call. Also contains the function which generates, saves and loads these correction matrices.

**`stochastic_rain.r`**
Function to simulate the stochastic rain fields from [Hewer 2018: "Stochastisch-physikalische Modelle fuer Windfelder und Niederschlagsextreme"](http://hss.ulb.uni-bonn.de/2018/5122/5122.htm). Original code by Rüdiger Hewer, simplified implementation (for the special case of zero correlations between q, &chi; and &psi;) by Sebastian Buschow.


---
### Example implementation:
**`simulation_example.r`**
Simulates a few realizations of the model-configurations used in the paper and saves them as `example_fields.rdata`

**`verif_example.r`**
Loads the stochastic rain realizations, transforms them and performs a short determinstic verification experiment via the EMD between the mean spectra. 


**`map_of_scales_example.r`**
Calculates and plots the map of central scales for one of the example cases from the `SpatialVx`-package.

---
### Full experiment:
*Take care, these programs are poorly optimized and hungry for RAM, especially when executed on many CPUs in parallel.*

**`simulation_full.r`**
Simulates a large number of realizations of the stochastic rain model and separates them into potential observations and (ensemble) forecasts. *Make sure to select an appropriate path to store the output.*

**`simulation_nonstat.r`**
As `simulation_full.r` but for the non-stationary model.

**`verification_full.r`**
This script can perform all of the syntehtic verification experiments from the paper, including the wavelet transform, aggregation and comparison via the propsed scores, as well as alternatives from the literature. *Make sure to put in the path where you stored the output of `simulation_full.r`.*

---
### Other things:
**`get_centres.r`**
Calculates the centres of mass for the daughter wavelets to which their positions are shifted after the transformation, saves them in `DBcentres.rdata`

**`verification_plots.r`**
Probably the most unnecessarily complicated script here, this is only used to visulaize the output of `verification_full.r`.

**`other_scores`**
Folder containing fortran source code and R-functions for the variogram scores, developed by Franka Nawrath, used with permission. Also contains R-functions for SAL and eSAL, adapted from the `SpatialVx` R-package.

**`Amats`**
Folder in which the bias correction matrices are stored.

**`results`** and **`plots`**
Folders where the output of `verification_full.r` and `verification_plots.r` are stored.


