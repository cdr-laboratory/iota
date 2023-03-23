#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# iota — A stochastic model of marine particle biogeochemistry
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# --------------------------------------------------------------------------------
# REQUIREMENTS: Written/tested on Matlab R2018a/R2018b/R2019a
#               [but should be backwards-compatible]
# --------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------
# DIRECTORY STRUCTURE
# --------------------------------------------------------------------------------
# iota/
#   ├── input/
#         ├── calibration_final.mat
#         └── psd/
#               └── cepturtis.mat
#               └── renforth.mat
#   ├── src/
#         └── BC_DIC.m
#         └── BC_O2.m
#         └── BC_TA.m
#         └── csys_CO3.m
#         └── csys_H2CO3.m
#         └── csys_pH.m
#         └── ODE_DIC.m
#         └── ODE_O2.m
#         └── ODE_TA.m
#   ├── LICENSE
#   ├── README.txt
#   ├── iota.m
#   ├── psd.m
# --------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------
# SUMMARY
# --------------------------------------------------------------------------------
# iota comprises a stochastic model of marine particle aggregation/disaggregation 
# embedded in a 1-D reaction-transport model of ocean biogeochemistry. The model
# is designed to repressing the impacts of adding alkaline mineral feedstocks to
# the surface oceans as a carbon dioxide removal (CDR) strategy. The model couples
# particle dynamics with vertical sinking, zooplankton grazing, microbial respiration
# and temperature/pH-dependent mineral transformation and alkalinity release. The 
# initial version of the model is published in:
#
#     [ ------------------------------------------------------- ]
#     [ Fakhraee et al. | 2023 | Environmental Research Letters ]
#     [ preprint doi:10.21203/rs.3.rs-1475007/v1                ]
#     [ ------------------------------------------------------- ]
#
# The supplementary information to Fakhraee et al. details model structure, 
# equations, and parameters and validates the base model against observations
# from the modern ocean.
# --------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------
# FILES
# --------------------------------------------------------------------------------
# iota.m — Main model script. Runs particle model and water column biogeochemistry
#          based on assumed feedstock application rates, net primary production,
#          temperature, and surface water chemistry.
# psd.m  — Script for loading and plotting particle size distributions (PSDs) for
#          milled mineral feedstocks. 
# calibration_final.mat — Validation data for 1-D biogeochemistry model. Source
#                         datasets are referenced w/DOIs in iota.m.
# cepurtis.mat — Particle size distribution for milled feedstock from:
#                [ Cepurtis et al. | 2017 | doi:10.1016/j.apt.2016.11.018 ]
# renforth.mat — Particle size distribution for milled feedstock from:
#                [ Renforth et al. | 2012 | doi:10.1016/j.ijggc.2012.06.011 ]
# BC_DIC.m     — Sets boundary conditions for dissolved inorganic carbon [DIC]
# BC_O2.m      — Sets boundary conditions for dissolved oxygen [O2]
# BC_TA.m      — Sets boundary conditions for dissolved alkalinity [TA]
# csys_CO3.m   — Carbonate system solver for carbonate ion [CO3]
# csys_H2CO3.m — Carbonate system solver for carbonic acid [H2CO3]
# csys_pH.m    — Carbonate system solver for pH
# ODE_DIC.m.   — Sets up ordinary differential equation for DIC
# ODE_O2.m.    — Sets up ordinary differential equation for oxygen
# ODE_TA.m.    — Sets up ordinary differential equation for alkalinity
# --------------------------------------------------------------------------------
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
