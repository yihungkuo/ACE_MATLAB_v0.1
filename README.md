# ACE_MATLAB_v0.1
MATLAB Implementation of Anelastic Convective Entity (ACE) Model version 0.1

This ACE_MATLAB_v0.1 code package is prepared for and released with the submission of Kuo & Neelin (2023) to the Journal of Atmospheric Sciences. Detailed model description can be found in the references list below.

The pre-computed nonlocal bases under `ACE_MATLAB_v0.1/data/basis/` are provided as is.

Running the ACE model requires ARM site sounding data (ARMBEATM) available at http://dx.doi.org/10.5439/1333748. 

A copy of a few samples---not part of this code package, provided only for user’s convenience---can be found here: [soundings.zip](https://drive.google.com/file/d/1XQ6rVE7Izc_5xipvHFaswgSCNdquk61T/view?usp=drive_link)


## To run the package
(1) place sounding data under `ACE_MATLAB_v0.1/data/soundings/` and

(2) execute the driver script `ace_v01.m` in MATLAB.


## Experimenting with the ACE model
In the `Set up parameters` section of `ace_v01.m` there are a few parameters users can vary, including: `arm`, `pidx`, `qc_ramp`, `conti`, `tspan`, as well as variables for `Initial mass flux`, `Initial thermal bubble`, and `External buoyancy forcing` \[most variables are documented with in-line comments in the driver script; initiation and forcing options documented in Kuo & Neelin (2023)].

Changing other variables is not recommended unless you are certain what you are doing.


## References:
Kuo, Y.-H., and J. D. Neelin, 2022: Conditions for convective deep inflow. Geophys. Res. Lett., 49, e2022GL100 552, doi:10.1029/2022GL100552.
Kuo, Y.-H., and J. D. Neelin, 2023: Anelastic Convective Entities: formulation and properties. J. Atmos. Sci., submitted. Preprint available: [KN23-ACE-preprint.pdf](https://drive.google.com/file/d/1vdOy_NckxwoU3WVhTsHkunjNt_fQqR7b/view?usp=drive_link)
