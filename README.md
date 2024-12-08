# data-science
how to run the code : 

## Project Structure
- `scripts/LB.R`: The main R script.
- `data-science/lending_borrowing_data.xlsx`: Input data file.

## How to Run
1. Clone the repository:
   git clone https://github.com/Jeremiedepins/data-science.git
2. Open R or RStudio.
3. Set the working directory to the project root.
   setwd("path/to/data-science")
4. Run the script:
   source("scripts/LB.R")

## what if it does not work ?
1. download the xlsx file and get the path, ex : "C:\Users\moi\Documents\test_python"
2. Open RStudio.
3. copy paste the R code in R-studio for it to work
4. replace the path to the path where you put the xlsx doc and In put : data <- read_excel(" C:\Users\moi\Documents\test_python\lending_borrowing_data ")
