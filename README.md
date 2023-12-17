# Hepatocytes-activity-distribution
The proposed set of functions are able to perform a statistical analysis of the data available from The bioIVT(R) database. 
The analysis is performed using the activity ditribution (pmol/min/Mio cells) from CYPs:
  * CYP2A6: coumarin
  * CYP2B6: bupropion
  * CYP2C8: amodiaquine
  * CYP2C9: tolbutamide
  * CYP2C19: mephenytoin
  * CYP2D6: dextromethorphan
  * CYP2E1: chlorzoxazone 
  * CYP3A4(tests.): testosterone
  * CYP3A4(MDZ): midazolam
  
  Additional analysis were reported for UGT1A1, ST (using 7-hydroxycoumarin), and AO  
  
  The analysis needs the input dataset structure reported in the input file "Human_Hepatocytes_data.csv" and the script automatically returns for each isoform:
   * The number of male and female donors investigate for the given isoform
   * The median, mean, SD, the p value for the Shapiro normality assessment, the Barlett's test outcome (homogeneity of the variance between male and female), and the Hartigian's test to evaluate the unimodal distribution
   * The mean and the SD in logarithmic domain
   * The median, mean, SD, best lambda, the p value for the Shapiro normality assessment for male and female, the Barlett's test and the p value from the ANOVA test after the Box-Cox transformation.
  
# Examples cases of the system parameters impact
The code report the simulations reported in the tutorial in order to explore the impact of the:
  * Surface area and volume of the apical side in the transwell of a Gut-Liver OoC.
  * Volume samling and evaporation.

In the codes are available the input data and the output reported as figures in the MS.
