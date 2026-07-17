# BOIN12-PFS
BOIN12-PFS: A Bayesian optimal interval phase I/II design incorporating progression-free survival endpoint

## Overview
This repository provides R code for conducting dose-finding trial simulations based on the BOIN12-PFS design.
The design jointly incorporates toxicity, overall response rate (ORR), and 12-month progression-free survival (PFS) information to identify the Optimal Biological Dose (OBD).
Operating characteristics are evaluated through simulations under a variety of dose-toxicity, efficacy, and PFS scenarios.
This code is intended for biostatisticians involved in the design and evaluation of phase I/II oncology dose-optimization studies.

## Repository Structure

- `boin12_pfs_functions.R`
  - Supporting functions for BOIN12-PFS simulations

- `boin12_pfs_simulation.R`
  - Main simulation program for Operating characteristic evaluation

## Requirements

- R version: 4.6.0 (2026-04-24 ucrt)

- R packages:
  - `mvtnorm`
  - `dplyr`
  - `ReIns`
  - `tidyverse`
  - `Iso`
