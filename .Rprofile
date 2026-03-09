# Use more CPU cores while building packages
options(Ncpus = parallel::detectCores())

# Use PPM and {pak} with {renv}
# https://github.com/rstudio/renv/
# https://github.com/r-lib/pak
options(renv.config.pak.enabled = TRUE)

# activate renv
source("renv/activate.R")

# Use the most recent snapshot of the BI's PPM
options(repos = c(RSPM = "https://packagemanager.posit.co/cran/latest"))
