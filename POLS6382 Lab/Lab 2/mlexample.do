**************************************************************
*POLS6382 Quantitative Method III
* Example: using ml to implement maximum likelihood estimation
* Last update: August 25, 2024
* Ling Zhu
**************************************************************
* Example data
sysuse auto.dta
* Define log-likelihood functionprogram define mylogit          args lnf Xb          replace ‘lnf’ = -ln(1+exp(-‘Xb’)) if $ML_y1==1          replace ‘lnf’ = -‘Xb’ - ln(1+exp(-‘Xb’)) if $ML_y1==0end
* Define model ml model if mylogit(foreign=mpg weight)
* Estimate the logit modelml maximize
* Double check results
logit foreign mpg weight
* End
