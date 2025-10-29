--------------------------------------------------------------------------------
Replication Data for Atsusaka (2021) "A Logical Model for Predicting Minority Representation"
--------------------------------------------------------------------------------

Contents:
   - Citation and Use
   - R Code Summary
   - Data Summary

For further questions, please contact the author at: 
   atsusaka@rice.edu
   yukia887@gmail.com

--------------------------------------------------------------------------------
Citation and Use of this Replication Data
--------------------------------------------------------------------------------

When citing this replication data (either its entirety, code, or datasets), please use the following formats:

For LaTeX Users:

@data{DVN/F2OX6O_2021,
author = {Atsusaka, Yuki},
publisher = {Harvard Dataverse},
title = {{Replication Data for: A Logical Model for Predicting Minority Representation: Application to Redistricting and Voting Rights Cases}},
year = {2021},
version = {DRAFT VERSION},
doi = {10.7910/DVN/F2OX6O},
url = {https://doi.org/10.7910/DVN/F2OX6O}
}

For non-LaTex Users:

Atsusaka, Yuki, 2021, "Replication Data for: A Logical Model for Predicting Minority Representation: Application to Redistricting and Voting Rights Cases", https://doi.org/10.7910/DVN/F2OX6O, Harvard Dataverse, DRAFT VERSION 

--------------------------------------------------------------------------------
R Code Summary
--------------------------------------------------------------------------------

This folder contains nine .R files that are necessary to replicate all the figures and tables in the main text and the Online Appendix.
Each file name indicates what result the code therein replicates:

   Figure2.R
   Figure3.R
   Figure4.R
   Figure5.R
   Table1.R
   Table2.R
   Appendix_A.R
   Appendix_B.R
   Appendix_C.R

Note that Appendix_B.R and Appendix_C.R contain multiple "chuncks" of code for replicating more than one figure or table.

--------------------------------------------------------------------------------
Data Summary
--------------------------------------------------------------------------------

This folder contains two datasets used in the main text and one dataset analyzed in the Online Appendix.

   Data_LAMayoral.csv
   Data_StateLegislative.csv
   Data_LAMayoral_Appendix.csv


"Data_LAMayoral.csv" contains information about outcomes of Black candidate emergence and victory in 303 Louisiana municipalities from 1986 to 2016.

Full title: Louisiana Mayoral Election Dataset
Unit: municipality-election
Number of observations: 2037
Number of variables:    11
Desciprtion of variables:

  entity_name:   string denoting the name of each city
  year:          numeric denoting the year of each election
  run:           binary indicator for Black candidate emergence (1 = At least one Black candidate, 0 = No Black candidate)
  win:           binary indicator for Black candidate victory (1 = the winner is Black, 0 = the winner is not Black)
  M_raw:         numeric representing the RAW racial margin of victory (Top Black candidate's vote share - Top non-Black candidate's vote share [both in %])
  M:             numeric representing the adjusted racial margin of victory (M_raw + 50)
  C:             numeric denoting the percentage of Blacks in the city population
  incumb_ran:    binary variable denoting whether the incumbent ran (1) or not (0)
  unopposed:     binary variable denoting whether the election was unopposed (1) or not (0)
  city_type:     categorical variable for the type of city (Rural, Suburban, or Urban)
  city_council:  categorical variable for the type of city council elections for each city (At Large, Single Member District, or Mixed System)


Data collection procedure:
This dataset is a transformed and reduced version of the Louisisna Mayoral Election Candidate-Level dataset collected under the Local Election in America Project (LEAP) at Rice University (The PI: Melissa Marschall (marschal@rice.edu)). 

For more about the LEAP, please visit at:
   https://leap.rice.edu/
   https://uwmrf.org/technology/ott1366/



"Data_StateLegislative.csv" contains information about outcomes of Black, Hispanic, and Asian candidate emergence and victory in 36 states in 2012 anf 2014.

Full title: State Legislative General Election Dataset
Unit: group-district-election
Number of observations: 1306
Number of variables:    17  
Desciprtion of variables:

   state:         abbreviated state names (CAPITALIZED)
   state.lower:   abbreviated state names (lower-cased)
   sl_chamber:    numeric denoting whether the election is for State Senate (8) or State House (9)
   sl_district:   numeric for the name of the single-member district
   year:          numeric denoting the year of each election 
   group:         categorical variable denoting whether the observation is for Asian ("Asian"), African American ("Black"), or Latin/x ("Hispanic") candidates
   minority_run:  binary variable denoting whether at least one minority candidate from the above group ran (1) or not (0)
   minority_win:  binary variable denoting whether the minority candidate from the above group won (1) or not (0)
   M:             numeric representing the adjusted racial margin of victory (Top Black candidate's vote share - Top non-Black candidate's vote share [both in %] + 50)
   C:             numeric denoting the percentage of minority voters from the above group in Citizen-Voting Age Population (CVAP)
   south:         binary variable denoting whether the election took place in the South:
                  Alabama, Arkansas, Delaware, the District of Columbia, Florida, Georgia, Kentucky, Louisiana, Maryland, Mississippi, North Carolina, Oklahoma, South Carolina, Tennessee, Texas, Virginia, and West Virginia
   deepsouth:     binary variable denoting whether the election took place in the Deep South:
                  Georgia, Alabama, South Carolina, Mississippi, and Louisiana
   rimsouth:      binary variable denoting whether the election took place in the Rim South:
                  Arkansas, Delaware, the District of Columbia, Florida, Kentucky, Maryland, North Carolina, Oklahoma, Tennessee, Texas, Virginia, and West Virginia
   section5:      binary variable denoting whether the district was covered by Sections 4 and 5 of the Voting Rights Act at the time of Shelby County v. Holder (2013) decision
   litigated:     binary variable denoting whether the district was litigated in the 2010 round of redistricting
   white_pct:     numeric for the percentage of White in CVAP
   unusual:       binary indicator for unusual districts (e.g., district numbers changed after redistricting)


Data collection procedure:
This dataset is a transformed and extensively argumented version of the State Legislative Election Dataset complied by Fraga, Juenke, and Shah (2019) (hereafter FJS).
Specifically, the following variables have been directly borrowed from FJS:
   state, state.lower, sl_chamber, sl_district, year, minority_run, C (modified from black_pch, latino_pch, and asian_pch from FJS)

Other variables are coded by the author based on the following methods:

minority_win:
This variable was coded based on the grid search on each election via the Ballotpedia.
For example, to code the Senate State legislative general election in Delaware 19th District, 
(1) I first visited: https://ballotpedia.org/Delaware_State_Senate_elections,_2012#District_19
(2) I then searched if the winner is Asian, African American, Hispanic, or White by looking at the candidate-specific page and/or search outside the Ballotpedia
(3) I finalized my code while coding as "White winner" if no definitive and clear information was available

M:
This variable was coded based on the grid search on each election via the Ballotpedia.
For each election, I coded the vote share of the top Asian, top Black, top Hispanic, and top White candidates in %, respectively.
Here, a top candidate means the most-vote earning candidate for each group.

south, deepsouth, and rimsouth:
To code these variables, I relied on the categorization used in Bullock and Rozell (2018) "The New Politics of the Old South" An Introduction to Southern Politics" Sixth Edition, Lanham, Maryland: Rowman & Little Field

section5:
This variable was coded in the following way:
(1) I first relied on the information about which county and states were covered obtained at: https://www.justice.gov/crt/jurisdictions-previously-covered-section-5
(2) For states with partial coverage, I then coded what districts the covered counties had at the time of Shelby County v. Holder (2013) decision
(3) Finally, I coded 1 if the district is in the set of covered jurisdictions (either states or counties)

litigated:
This variable was coded by relying on the information obtained at: https://www.ncsl.org/research/redistricting/redistricting-case-summaries-2010-present.aspx
(1) I first visited the NCSL website and checked if each state has any redistricting-related litigation on the grounds of minority voting rights
(2) I then scrutinized publicly available documents related to each case and manually recored the name of legislative districts that either plaintiffs, defendants, or other parties cited in the documents

unusual:
This variable was coded by individually scrutinizing each of several seemingly-unreasonable districts (e.g., M and C are both very low, but minority candidates won).
(1) To check if such districts are unusual cases, I first validate the coding by visiting the Ballotpedia 
(2) If the original coding by JFS was incorrect, I simply corrected the error
(3) If the district has experienced any irregular change, for example District 3 became District 5 after redistricting, I coded the distrcit as unusual (1)



For more about the Fraga, Juenke, and Shah (2019) dataset, please visit at:
   "Replication Data for: One Run Leads to Another: Minority Incumbents and the Emergence of Lower Ticket Minority Candidates"
   https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/IUXG2I



"Data_LAMayoral_Appendix.csv" contains additional information on Louisiana mayoral elections for analyses reported in the Online Appendix.


Full title: Louisiana Mayoral Election Dataset for the Online Appendix
Unit: municipality-election
Number of observations: 2037
Number of variables:    21
Desciprtion of variables:

See the above list for the variables contained in "Data_LAMayoral.csv".
Additional variables included to this dataset are as follows:

  woman_run:         binary variable denoting whether at least one woman candidate ran in each election (1) or not (0)
  woman_win:         binary variable denoting whether at least one woman candidate ran in each election (1) or not (0)
  num_black_cand:    numeric denoting the number of Black candidates
  M_t2:              the (adjusted) racial margin of victory at time T-2
  M_t3:              the (adjusted) racial margin of victory at time T-3
  educ_baplus_black: numeric denoting the percentage of African Americans with BAs
  educ_baplus_white: numeric denoting the percentage of Whites with BAs
  new_electiontime:  binary variable denoting whether the election is on-cycle (1) or off-cycle (0)
  white_over65:      numeric denoting the percentage of Whites who are 65 or older
  density:           numeric representing the human density in each municipality



--------------------------------------------------------------------------------
End of This README.txt File
--------------------------------------------------------------------------------