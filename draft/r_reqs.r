old <- getOption("defaultPackages"); 
r <- getOption("repos")
r["CRAN"] <- "http://cran.stat.auckland.ac.nz"
options(defaultPackages = c(old, "MASS","methods"), repos = r)

packages <- c(
    'gstat'
)

install.packages('tidyr')

for (p in packages) {
    if (!require(p, character.only=TRUE)) {
        install.packages(p)
    }
}

if(!'plotmapbox' %in% installed.packages()[,"Package"]){
    require(devtools)
    install_github('dragonfly-science/plotmapbox')
}

