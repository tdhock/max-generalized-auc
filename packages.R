### Write down what package versions work with your R code, and
### attempt to download and load those packages. The first argument is
### the version of R that you used, e.g. "3.0.2" and then the rest of
### the arguments are package versions. For
### CRAN/Bioconductor/R-Forge/etc packages, write
### e.g. RColorBrewer="1.0.5" and if RColorBrewer is not installed
### then we use install.packages to get the most recent version, and
### warn if the installed version is not the indicated version. For
### GitHub packages, write "user/repo@commit"
### e.g. "tdhock/animint@f877163cd181f390de3ef9a38bb8bdd0396d08a4" and
### we use install_github to get it, if necessary.
works_with_R <- function(Rvers,...){
  local.lib <- file.path(getwd(), "library")
  dir.create(local.lib, showWarnings=FALSE, recursive=TRUE)
  .libPaths(c(local.lib, .libPaths()))
  pkg_ok_have <- function(pkg,ok,have){
    stopifnot(is.character(ok))
    if(!as.character(have) %in% ok){
      warning("works with ",pkg," version ",
              paste(ok,collapse=" or "),
              ", have ",have)
    }
  }
  pkg_ok_have("R",Rvers,getRversion())
  pkg.vers <- list(...)
  for(pkg.i in seq_along(pkg.vers)){
    vers <- pkg.vers[[pkg.i]]
    pkg <- if(is.null(names(pkg.vers))){
      ""
    }else{
      names(pkg.vers)[[pkg.i]]
    }
    if(pkg == ""){# Then it is from GitHub.
      ## suppressWarnings is quieter than quiet.
      if(!suppressWarnings(require(requireGitHub))){
        ## If requireGitHub is not available, then install it using
        ## devtools.
        if(!suppressWarnings(require(devtools))){
          install.packages("devtools")
          require(devtools)
        }
        install_github("tdhock/requireGitHub")
        require(requireGitHub)
      }
      requireGitHub(vers)
    }else{# it is from a CRAN-like repos.
      if(!suppressWarnings(require(pkg, character.only=TRUE))){
        install.packages(pkg)
      }
      pkg_ok_have(pkg, vers, packageVersion(pkg))
      library(pkg, character.only=TRUE)
    }
  }
}
options(repos=c(
  "http://www.bioconductor.org/packages/release/bioc",
  ##"http://r-forge.r-project.org",
  "http://cloud.r-project.org",
  "http://cran.r-project.org"))
works_with_R(
  "4.1.0",
  data.table="1.14.0",
  future="1.21.0",
  future.apply="1.7.0",
  RJSONIO="1.3.1.4",
  R.utils="2.10.1",
  "tdhock/penaltyLearning@558be4c28702396e37477d7fc01bad17d456fd82",
  "tdhock/directlabels@5f8ba1ae60c2a87b96d09a6f9ccaadd6f3859c57",
  jointseg="1.0.2",
  gridExtra="2.3",
  neuroblastoma="1.0",
  tikzDevice="0.12.3.1",
  microbenchmark="1.4.7",
  ##"tdhock/animint2@8a7327e66fa65f31ccc17508f914aaaee41909ef"
  ggplot2="3.3.3")
library(aum)
## requireGitHub::requireGitHub_package(
##   "great-northern-diver",
##   "loon/R",
##   "33442bae2ba2dca66050c36081079f9fb45a9358",
##   "loon")
options(
  tikzDocumentDeclaration=paste(
    "\\documentclass[12pt]{article}\\usepackage{amsmath,amssymb}"),
  tikzMetricsDictionary="tikzMetrics")
