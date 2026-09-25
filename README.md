# sloshing-matlab-simulink-model

# Sloshing Phenomenon — Multi-Pendulum Model (MATLAB/Simulink)

Lab assignment for the course *Analysis and Simulation of Aerospace Systems*, 
Politecnico di Milano, Aerospace Engineering, A.Y. 2024/2025.

Group project — Emma Mori, Giacomo Trabacchin, Alex Triolo, Chiara Vaccari.

## Overview
This project models lateral liquid sloshing in a partially filled cylindrical 
tank using an equivalent mechanical system of pendulums, derived from the 
linearized potential-flow solution of the fluid domain (antisymmetric slosh 
modes, based on Bessel function roots).

The model is developed in state-space form and used to:
- derive pendulum masses, lengths and natural frequencies from tank geometry
- build the linearized equations of motion, with and without modal damping
- compare analytical, modal-analytical and numerical (`ode45`, `impulse`, `step`) 
  responses to impulse and step inputs
- determine the minimum number of pendulums needed for an accurate, 
  computationally efficient model (convergence study)
- derive and compare the frequency response of the full dynamic model against 
  a simplified "frozen-liquid" model
- validate the MATLAB state-space model against an equivalent Simulink 
  implementation

## Contents
- `report.pdf` — full report: theory, derivations, results and plots
- MATLAB scripts and Simulink model files implementing the multi-pendulum 
  system, response analyses and convergence study

## Key results
- Optimal model order: 4 pendulums, within a 1% tolerance on impulse response amplitude
- Full agreement between analytical, modal-analytical, numerical (MATLAB) and 
  Simulink solutions
- At high frequency, system dynamics converge to the frozen-liquid 
  approximation, governed by the fixed mass component
