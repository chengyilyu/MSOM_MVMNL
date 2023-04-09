# Software Package for Assortment Optimization under MVMNL

## Overview

This repository contains numerical experiments for the paper 
[Assortment Optimization with Multi-Item Basket Purchase under the Multivariate MNL Model](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3818886).

## Folders
- `MSOM_Data`: The Instacart Online Grocery Shopping Dataset 2017
- `Fitness_Codes/`: Python files for fitness tests
- `FPTAS_SGS_Uncapacitated/`: Matlab codes for uncapacitated cases
- `FPTAS_SGS_Capacitated/`: Matlab codes for capacitated cases
- `Additional_Numerical/`: Matlab codes for additional numerical tests


## MSOM_Data
- `MSOM_Data/Instacart_data`: Instacart raw data
- `MSOM_Data/Instacart_cleaned`: We have filtered the data in 'orders_unique.csv' and prepared a subset of data in 'orders_unique_test.csv' as inputs for cross-validation.

## Fitness_Codes
Each python file in `Fitness_Codes/` conducts cross validation with 'orders_unique_test.csv' as inputs. File paths in the Python files need to be updated to input data.
- `MNL_CrossValidation.py`: conduct cross validation under MNL model
- `MVMNL_CrossValidation.py`: conduct cross validation under MVMNL model


## FPTAS_SGS_Capacitated
- `Main_Capacitated.m` is used for one run of a certain capacitated numerical experiment, while other files define some functions to be imported by `Main_Capacitated.m`.
- `FPTAS_capacitated`: FPTAS
- `Heuristic_Golden_Section_capacitated.m`: SGS algorithms
- `Optimal_capacitated.m`: find optimal policies
- `GroupwiseMNL_capacitated.m`: a benchmark that treats each group as an MNL model
- `StaticMNL.m`: an algorithm in Rusmevichientong et al. (2010)


## FPTAS_SGS_Uncapacitated
- `Main_Uncapacitated.m` is used for one run of a certain uncapacitated numerical experiment, while other files define some functions to be imported 
by `Main_Uncapacitated.m`. The descriptions of these files are similar to those found in the `FPTAS_Capacitated/`.


## Additional_Numerical
The files are for the additional numerical experiments in Appendix/Online version
- `Main.m` is used for one run of numerical experiment aimed at gaining insights on the interaction term, 
while other files define some functions to be imported by `Main.m`.
