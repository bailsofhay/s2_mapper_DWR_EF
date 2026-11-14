# =============================================================================
# Sentinel-2 Harmonized: Spectral Index Functions
# =============================================================================
#
# Description:
#   Functions for calculating spectral vegetation and water indices from
#   Sentinel-2 Harmonized surface reflectance imagery. These indices are
#   used as predictor variables in aquatic vegetation classification.
#
# Input:
#   A terra SpatRaster containing Sentinel-2 spectral bands, with layer
#   names corresponding to the Sentinel-2 band identifiers (e.g., B2,
#   B3, B4, B5, B8, B8A, B11).
#
# Output:
#   Each function returns a single-layer SpatRaster containing the
#   calculated spectral index, named according to the function.
#
# Notes:
#   - Input bands should contain surface reflectance values.
#   - Reflectance scaling must be consistent with each index equation.
#   - Missing input values propagate to the calculated indices.
#   - Functions can be sourced by other scripts using source().
#
# =============================================================================


# 1. Spectral Index Functions -------------------------------------------------

# Enhanced Vegetation Index (EVI)
# Uses blue, red, and NIR reflectance to characterize vegetation
# greenness while accounting for atmospheric and background effects.
evi <- function(r) {
  
  out <- 2.5 * ((r$B8 - r$B4) /
                  (r$B8 + 6 * r$B4 - 7.5 * r$B2 + 1))
  
  names(out) <- "evi"
  return(out)
}


# Green Normalized Difference Vegetation Index (GNDVI)
# Uses green and NIR reflectance to characterize vegetation greenness.
gndvi <- function(r) {
  
  out <- (r$B8 - r$B3) / (r$B8 + r$B3)
  
  names(out) <- "gndvi"
  return(out)
}


# Modified Chlorophyll Absorption Ratio Index (MCARI)
# Uses green, red, and red-edge reflectance to characterize chlorophyll
# absorption while reducing background effects.
mcari <- function(r) {
  
  out <- ((r$B5 - r$B4) - 0.2 * (r$B5 - r$B3)) *
    (r$B5 / r$B4)
  
  names(out) <- "mcari"
  return(out)
}


# Modified Normalized Difference Water Index (MNDWI)
# Uses green and SWIR reflectance to distinguish water from non-water.
mndwi <- function(r) {
  
  out <- (r$B3 - r$B11) / (r$B3 + r$B11)
  
  names(out) <- "mndwi"
  return(out)
}


# Normalized Difference Aquatic Vegetation Index (NDAVI)
# Uses blue and NIR reflectance to distinguish aquatic vegetation
# from surrounding water.
ndavi <- function(r) {
  
  out <- (r$B8 - r$B2) / (r$B8 + r$B2)
  
  names(out) <- "ndavi"
  return(out)
}


# Normalized Difference Moisture Index (NDMI)
# Uses narrow NIR and SWIR reflectance to characterize vegetation
# and surface moisture conditions.
ndmi <- function(r) {
  
  out <- (r$B8A - r$B11) / (r$B8A + r$B11)
  
  names(out) <- "ndmi"
  return(out)
}


# Normalized Difference Vegetation Index (NDVI)
# Uses red and NIR reflectance to characterize vegetation greenness.
ndvi <- function(r) {
  
  out <- (r$B8 - r$B4) / (r$B8 + r$B4)
  
  names(out) <- "ndvi"
  return(out)
}


# Soil Adjusted Vegetation Index (SAVI)
# Uses red and NIR reflectance with a soil adjustment factor (L = 0.5)
# to reduce background influences on vegetation greenness.
savi <- function(r) {
  
  out <- ((r$B8 - r$B4) / (r$B8 + r$B4 + 0.5)) * 1.5
  
  names(out) <- "savi"
  return(out)
}


# Visible Atmospherically Resistant Index (VARI)
# Uses visible bands to characterize vegetation greenness.
vari <- function(r) {
  
  out <- (r$B3 - r$B4) / (r$B3 + r$B4 - r$B2)
  
  names(out) <- "vari"
  return(out)
}


# Water Adjusted Vegetation Index (WAVI)
# Uses blue and NIR reflectance with a background adjustment factor
# to characterize aquatic vegetation.
wavi <- function(r) {
  
  out <- (r$B8 - r$B2) / (0.5 + r$B8 + r$B2)
  
  names(out) <- "wavi"
  return(out)
}