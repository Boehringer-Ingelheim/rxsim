# Use more CPU cores while building packages
options(Ncpus = 2)

# Use PPM and {pak} with {renv}
# https://github.com/rstudio/renv/
# https://github.com/r-lib/pak
options(renv.config.pak.enabled = TRUE)

# Ensure user-local tools are visible to R (e.g. pandoc at ~/local/bin)
local_bin <- path.expand("~/local/bin")
if (dir.exists(local_bin) && !grepl(paste0("(^|", .Platform$path.sep, ")", local_bin, "($|", .Platform$path.sep, ")"), Sys.getenv("PATH"))) {
	Sys.setenv(PATH = paste(local_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
}

pandoc <- Sys.which("pandoc")
if (nzchar(pandoc)) {
	Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc))
}

# activate renv
source("renv/activate.R")

# Use the most recent snapshot of the BI's PPM
options(repos = c(BI = "https://pm.prod.copernicus.aws.boehringer.com/cinema/latest"))
