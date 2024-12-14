# data-science
how to run the code : 

## Project Structure
- `scripts/LB.R`: The main R script.
- `data-science/lending_borrowing_data_csv.csv`: Input data file.
- `data-science/lending_borrowing_data.xlsx`: Input data file if csv does not work.

## How to Run #1 way
1. For the first way you just need to download the file in a folder called "project_data_LB"
2. open R-studio and run the code

## How to Run #2 way throught github
1. For the second way way you just need to download the LB.R script
2. open R-studio and run the code


## How to Run #3 way throught local computer
1. download the xlsx or CSV file and get the path, ex : "C:\Users\moi\Documents\test_python" (try with csv first and if it does not work try with xlsx file)
2. Open RStudio.
3. copy paste the R code in R-studio for it to work
4. replace the path to the path where you put the xlsx doc and In put : data <- read_excel(" C:\Users\moi\Documents\test_python\lending_borrowing_data ")


## How to Run with repository (this is another way again to run your code if none of the 3 way arn't working)
1. Clone the repository:
   git clone https://github.com/Jeremiedepins/data-science.git
2. Open R or RStudio.
3. Set the working directory to the project root.
   setwd("path/to/data-science")
4. Run the script:
   source("scripts/LB.R")
