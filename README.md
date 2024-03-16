# Yeast-Transcript-Half-Life-Calculation-and-Functional-Enrichment-Analysis
This project aims to analyze the yeast transcript half-lives using time series data and identifies genes with significantly high or low half-lives. 
Then the list of output genes are used to perform a basic functional enrichment analysis to explore potential biological implications.

# 1. Data Loading and Preprocessing: 
After you load the time series data from DecayTimecourse.txt containing multiple columns representing time points and corresponding transcript abundance values.
Three dataframes (tc1, tc2, tc3) are created to represent the three sets of time course data.
To normalize the data and account for exponential decay, a log transformation is applied to each transcript's time series data within each dataframe.

# 2. Half-Life Calculation and Replicate Averaging:

The code performs a linear regression analysis for each transcript on the log-transformed data.
The slope of the regression line represents the transcript's decay rate.
The half-life for each transcript is calculated in each dataframe using the formula: half-life = ln(2) / (absolute value of slope).
Finally,average half-life for each transcript across the three dataframes (tc1, tc2, tc3) is calculated to provide a more robust estimate.

# 3. Identification of High/Low Half-Life Transcripts:

The average half-life for all transcripts is calculated.
Finally we aim to identify the transcripts with the top 10% highest and bottom 10% lowest average half-lives, potentially representing long-lived and short-lived transcripts, respectively.
