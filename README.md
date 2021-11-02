# Mendelian randomization suggests a bidirectional, causal relationship between physical inactivity and obesity
This is a repository of the Mendelian Randomization analysis ran for the publication "Mendelian randomization suggests a bidirectional, causal relationship between physical inactivity and obesity".

The repository consists of two sections:

1) The folder structure and code to reproduce our analysis in the /R folder.
2) The variants and summary statistics used to produce our results, plus a code to easily navigate through them in /data_&_results

## GWAS summary statistics used:

In this publication we are testing whether there is any causal association between adiposity and physical activity or inactivity. To do so we used the latest and largest GWAS summary statistics for each of the traits. The main analysis are performed for Body Mass Index (BMI) and secondary analysis for three other anthropometric traits.

### Physical activity and inactivity traits

For physical activity we decided to use accelerometer data for vigorous physical activity, moderate physical activity and sedentary time from two different sources.

Vigorous physical activity from the acc425 trait from Klimentidis et al available here: (link)
Moderate physical activity form Doherty et al available here: (link)
And sedentary time from Doherty et al available here: (link)

### Anthropometric traits:

As explained before, the main analysis are for BMI as the main adiposity trait, though we also ran analysis for other traits that inform about central obesity or abdominal fat accumualation: body fat percentage, waist circumference adjusted for BMI and waist-to-hip ratio adjusted for BMI. 

BMI GWAS summary statistics from Pulit et al can be found here: (link)
WHRadjBMI GWAS summary statistics from Pulit et al can be found here: (link)
WCadjBMI GWAS summary statistics from Shungin et al can be found here: (link)
BFP GWAS summary statistics from Elsworth et al can be found here: (link)

## Requirement to run CAUSE:




