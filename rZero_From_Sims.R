#
# Clear out data
#
rm(list=ls())

#
# Source in the RZero function
#
source("rZero_Function.R")

#
# Do the calculation
#
v <- findRZero(".",
               "ls_1_epidemics",
               2,
               T)
print(v)
