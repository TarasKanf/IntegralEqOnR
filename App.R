library("magrittr") # Pipe operator in R
library("modules")

# Order is important!
source("ProblemSolver.R")

# Landweber
result <- ProblemSolver$start(n=4, regularizationParam=0.1)

print(result)
print("aloha")