# Use more CPU cores while building packages
options(Ncpus = 2)

# Use PPM and {pak} with {renv}
# https://github.com/rstudio/renv/
# https://github.com/r-lib/pak
options(renv.config.pak.enabled = TRUE)

# activate renv
source("renv/activate.R")

# Use the most recent snapshot of the BI's PPM
options(repos = c(BI = "https://pm.prod.copernicus.aws.boehringer.com/cran/latest"))
