# Category_Comparison_Algorithim

This little algorithim finds pairs of individuals from a dataframe (generated from the main file) and comapres then according to a
scoring system.  The algorithim itself came about as my solution to the following question

Assume a cohort of 1000 individuals with ages range from 20 to 40, M/F gender, and 8 chromosomes that may be of African/Asian/European 
ancestry (e.g., 25-M-Af-Eu-Af-As-Eu-Eu-As-As). Your goal is to develop a Matlab/R code that will optimize pairing of the most similar
individuals. The scoring scheme is as follows:
every category/chromosome where individuals match equals one point (age difference of 5<=years is considered a match).
An identical individual pair worth 10 points and the maximum score you can get is 10,000. Leaving individuals unpaired has a 
penalty of 5 points per individual.Include a code to generate this dataset and some plot showing the scoring distribution of your code 
for multiple datasets.

Steps carried out by R script for comaprison algorithim

1. generates dataframe 

2.  creates distance matrix for all combinations of matched pairs according to age difference metric.

3. creates function to calculate score of matched rows.

4. make 1000 x 1000 matrix of where age difference are greater than 5

5.  make  1000 x 1000 matrix of all normalised scores between pairs

6. Add the two 1000 x 1000 matrices together and remove all entries greater than 1 

7. convert matrix of scores back to unnormalised form.

8. Use this matrix to collect all pairs corresponding to the 10, 7, 6, 5, 4, 3,2,1, -5 scores.

9. remove duplicates from 10's data frame 
   remove  duplicates from 7's dataframe
   remove individuals already in 10 from 7 and combine both data frames
   remove duplicates

10. repeat step 9 for 6's 5's, 4's 3' 2's 1's and -5's data frames.


The R script "main.R" takes about 10 minutes to run. It then generates a .xlsx file showing all the matched pairs and
their corresponding score. The "run_results_summary.xlsx" file shows the number of matched pairs, total score, median and mean of 
15 different runs of the "main.R" script.








