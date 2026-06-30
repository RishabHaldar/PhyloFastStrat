# deploy.R
# Redeploy PhyloFastStrat to shinyapps.io.
#
# Prerequisites (one-time):
#   1. Create a free shinyapps.io account at https://www.shinyapps.io
#   2. Account → Tokens → copy the rsconnect::setAccountInfo(...) snippet
#      and run it once in your R console.
#
# Usage: from R, with this repo as the working directory:
#   source("deploy.R")

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}

rsconnect::deployApp(
  appDir       = ".",
  appName      = "PhyloFastStrat",
  appFiles     = c(
    "app.R",
    list.files("data", full.names = TRUE)
  ),
  forceUpdate  = TRUE,
  launch.browser = TRUE
)
