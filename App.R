library("magrittr") # Pipe operator in R
library("modules")

# Order is important!
source("ProblemSolver.R")

result <- ProblemSolver$start(n=16, regularizationParam=0.1)

print("solution")
print(result)
print("aloha")