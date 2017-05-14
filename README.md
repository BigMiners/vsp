## How to run this VSP implementation on the IMDb dataset

### Requirements:

- shiny : In a R console, run command 'install.packages("shiny")'
- shinydashboard : In a R console, run command 'install.packages("shinydashboard")'
- leaflet : In a R console, run command 'install.packages("leaflet")'
- igraph : In a R console, run command 'install.packages("igraph")'
- gdata : In a R console, run commend 'install.packages("gdata")'
- slam : In a R console, run command 'install.packages("slam")' 
- networkD3 : In a R console, run command 'install.packages("networkD3")'
- gurobi (temporary : mandatory until SCIP replaces gurobi entirely) :

-----

### In short:

follow instructions on the "Getting started" guide http://www.gurobi.com/documentation/current/quickstart_linux.pdf

-----

### Some details:

- Register with an academic account on http://www.gurobi.com/registration/general-reg#Reg .
- Request a free academic licence on http://user.gurobi.com/download/licenses/free-academic ; you will get on a page with a command of the form "grbgetkey my-licence-key". Run it in your terminal to activate your gurobi installation.
- Set and update the variables GUROBI_HOME, PATH and LD_LIBRARY_PATH according to the installation guide (cf READ_ME)
- Install the Gurobi package in R : run command in a R console : "install.packages('<R-package-file>', repos=NULL)" with <R-package-file> being "$GUROBI_HOME/R/gurobi-version.tar.gz"

### Running the program:

- Make a copy of config.R.DEFAULT and name it config.R
- Edit config.R and set the required variables at the top of the file
- Run command `R < IMDbGraph.R --no-save

Png files for plots will be placed under the temporary directory you have set in config.R.


### Miscelaneous:

- If you want to run the project from rstudio, you must launch it from the terminal as such : 'rstudio &'.
- On a mac, use "install.packages('/Library/gurobi650/mac64/R/gurobi_6.5-0.tgz', repos=NULL)" to install gurobi package in R.
