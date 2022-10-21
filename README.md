Hi there, all codes used are given in the R script named ‘Code. R’ as well as the input values. Most instructions are given inside and if the results did not match or the code did not run, please track back and see if any line was not run or any values are overwritten accidentally. 
	For the table in the model performance, below are the R variables use:
•	Table 1: MSE_ELC_In (for the Extended LC model), MSE_HU_In (for the Hyndman_Ullah model), and MSE_SDF_In (for the SDF model) where the first row represents a one-factor model, the second row represents two-factor model, and so on.

•	Table 2: MSE_SDF_VAR (for the VAR process), and MSE_SDF (for the ARIMA process) where the first row represents a one-factor model, the second row represents a two-factor model, and so on.


•	Table 3: MSE_SDF (for the SDF model), MSE_HU (for the Hyndman_Ullah model), MSE_SLC (for the Standard LC model), and MSE_ELC (for the Extended LC model)

•	Table 4: The correlation matrix is obtained by the code given in the R script.


•	Table 5: MSE_SDF_In_diff with each value is obtained by code given in R script with a different set of age factors.

•	Table 6: MSE_SDF_diff with each value is obtained by code given in R script with a different set of age factors.

The mortality data used is kept in Mortality.RData with name of FRdata. 
	Please do email if any confusion is made as the code are checked without error and can be run using either a Monash desktop or an HP laptop with Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz 1.99 GHz, 4.00 GB (3.88 GB usable), 64-bit operating system, x64-based processor. And the software used is Rstudio with R version 4.0.3 (2020-10-10) and R version 4.2.0 (2022-04-22 ucrt)




Packages References
•	Villegas AM, Kaishev VK, Millossovich P (2018). “StMoMo: An R Package for Stochastic Mortality Modeling.” _Journal of Statistical Software_, *84*(3), 1-38. doi: 10.18637/jss.v084.i03 (URL: https://doi.org/10.18637/jss.v084.i03).
•	    Rob J Hyndman with contributions from Heather Booth, Leonie Tickle and John Maindonald. (2019). demography: Forecasting Mortality, Fertility, Migration and Population Data. R package version 1.22.  https://CRAN.Rproject.org/package=demography
•	  Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
•	  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
•	  Hans W. Borchers (2022). pracma: Practical Numerical Math Functions. R package version 2.3.8. https://CRAN.R-project.org/package=pracma
•	R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL: https://www.R-project.org/.
•	  Rob Hyndman (2021). fpp3: Data for "Forecasting: Principles and Practice" (3rd Edition). R package version 0.4.0. https://CRAN.R-project.org/package=fpp3
