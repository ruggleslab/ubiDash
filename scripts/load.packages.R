# ++++++++++++++++++++++++++++
# load.packages
# ++++++++++++++++++++++++++++
load.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
  }

  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg)
  }

  sapply(pkg, require, character.only = TRUE, quietly = TRUE)
}
