# FMskewt-censored

This repository contains the code for "Finite mixture of regression models for censored data based on the skew-t distribution," which was submitted to COST. In this work, we address the problem of modeling censored data using a finite mixture of regression models based on the skew-t distribution. Our approach allows for flexible modeling of the skewness and heavy-tailed often observed in real-world data, while also accounting for the censoring that frequently occurs in many applications.

Our main contributions include a novel approach to modeling data with different censored level using a finite mixture of regression models based  on the skewed-t distribution, as well as a comprehensive implementation of our methodology in R. Our code enable user to get parameter estimation, corresponding standard deviation, and model comparison criteria-AIC, BIC and EDC for four different types of FM-CR model such as FM-NCR, FM-TCR, FM-SNCR and FM-STCR.

We provide numerical experiments to demonstrate the performance of our approach. Specifically, we include two simulation scenarios and a real data analysis, all of which are introduced in the paper.

Our repository includes three source files:
- algEM.fmr.ST.cr.try.R: function for EM algorithm for Finite Mixtures of Cnesored data under N, T, SN and ST. This function allows users to easily determine the type of censoring present in their data using y (for left-censored data) and -y (for right-censored data).
- Integral_nu_float.R: defind the density function of location-scale Student's-T distribution
- Moment_SMSNT.R: Functions for calculating EM algorithm



