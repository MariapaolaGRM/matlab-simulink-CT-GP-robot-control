# Computed Torque Control with Gaussian Processes for a 2R Robot
MATLAB/Simulink implementation of model-based control with Gaussian Process Regression (GPR) for a planar two-link robotic manipulator.

The project investigates whether a Gaussian Process can learn the residual dynamic effects not captured by a nominal robot model and improve trajectory tracking without relying exclusively on high feedback gains.

The implementation is inspired by:

> M. Beckers, S. Umlauft, and S. Hirche,  
> "Stable Model-based Control with Gaussian Process Regression for Robot Manipulators,"  
> IFAC-PapersOnLine, 2017.

This is an educational simulation project with its own robot parameters, dataset-generation procedure, and experimental setup. It is not intended as an exact numerical reproduction of the original paper.

## Requirements
- MATLAB.
- Simulink.
- Statistics and Machine Learning Toolbox, for fitrgp.
- All MATLAB functions and Simulink models included in this repository.

## Repository Structure

The repository contains MATLAB scripts for dataset generation, GP training, validation, and controller comparison, together with the functions and Simulink models required to run the experiments.

```text
.
├── README.md
├── [Report](https://github.com/MariapaolaGRM/matlab-simulink-CT-GP-robot-control/relation.pdf)     # Detailed description of the proposed methods.
├── scripts/
│   ├── generate_gp_dataset.m
│   ├── train_gp.m
│   ├── prepara_simulink.m
│   ├── analisi_gp.m
│   ├── check_dataset.m
│   └── controller_comparison.m
├── simulink/
│   ├── GP_DataCollection_2R.slx
│   └── CT_PD_robot_2R.slx
└── src/
    ├── get_params.m
    ├── H.m
    ├── C.m
    ├── g.m
    ├── dynamics_model.m
    ├── computed_torque_model.m
    └── gp_mean.m
```

The structure above describes the main project components. Additional files, datasets, and result folders may be present depending on the current repository version.

## How to run

| Step | Script | Purpose |
|:---:|---|---|
| 1 | `generate_gp_dataset.m` | Generate training and validation datasets from robot simulations. |
| 2 | `check_dataset.m` | Check dataset integrity and consistency. |
| 3 | `train_gp.m` | Train the Gaussian Process models. |
| 4 | `prepara_simulink.m` | Export GP parameters for Simulink integration. |
| 5 | `analisi_gp.m` | Evaluate prediction accuracy and the effect of training set size. |
| 6 | `controller_comparison.m` | Compare trajectory tracking performance across controllers. |


The scripts expect the generated .mat files to be available in the MATLAB working directory. Check the paths and output locations if the repository structure is changed.


