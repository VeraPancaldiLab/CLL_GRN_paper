# Install system dependencies for the R project
if (Sys.info()["sysname"] == "Linux") {
  message("Detected Linux system. Attempting to install system dependencies.")
  
  # List system dependencies
  dependencies <- c(
    "libcurl4-openssl-dev",
    "libfontconfig1-dev",
    "libfreetype6-dev",
    "libxml2-dev",
    "cmake",
    "libharfbuzz-dev",
    "libfribidi-dev",
    "libpng-dev",
    "libtiff5-dev",
    "libjpeg-dev",
    "libudunits2-dev",
    "libgdal-dev"
  )
  
  # Install system dependencies
  for (dep in dependencies) {
    message(paste("Installing", dep))
    system(paste("sudo apt-get install -y", dep), intern = TRUE)
  }
  
  message("System dependencies installed successfully.")
  
} else {
  message("This script currently only supports Linux systems. Please install system dependencies manually.")
}
