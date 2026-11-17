# =============================================================================
# Sentinel-2 Spectral Index Functions
# =============================================================================
# Description:
#   Functions to calculate spectral indices from Sentinel-2 surface
#   reflectance imagery (scaled by 10000 to reflectance units).
#
#   Valid reflectance values of 0 are preserved in the input imagery.
#   To prevent undefined index calculations (Inf, -Inf, NaN), denominators
#   equal to 0 are replaced with 0.0001 using the safe_denom() function.
#
#   Input:  terra SpatRaster with Sentinel-2 band names
#   Output: Single-layer SpatRaster named after the spectral index
# =============================================================================


# 1. Helper Function ----------------------------------------------------------

# nonnegative bands (e.g., NDVI), a zero denominator generally occurs when
# both contributing bands have zero reflectance.
#
# Dividing by 0 produces Inf, -Inf, or NaN, depending on the numerator.
# These nonfinite values can cause problems during random forest classification,
# including pixels being excluded from prediction and gaps in classification
# maps.
#
# The safe_denom() function replaces denominators equal to 0 with a small
# positive constant (0.0001) before division. All nonzero denominators are
# left unchanged, and existing NA values remain NA.
#
# This adjustment prevents division-by-zero errors while preserving the
# original spectral reflectance values. However, index values calculated
# using the replacement constant are approximations rather than mathematically
# defined index values and may be unusually large, particularly for VARI.

# Replace zero denominators with 0.0001
safe_denom <- function(x) {
  terra::ifel(x == 0, 0.0001, x)
}


# 2. Spectral Index Functions -------------------------------------------------

# 2.1. Enhanced Vegetation Index (EVI) ----------------------------------------

evi <- function(r) {

  denom <- safe_denom(r$B8 + 6 * r$B4 - 7.5 * r$B2 + 1)

  out <- 2.5 * ((r$B8 - r$B4) / denom)

  names(out) <- "evi"
  return(out)
}


# 2.2. Green Normalized Difference Vegetation Index (GNDVI) -------------------

gndvi <- function(r) {

  denom <- safe_denom(r$B8 + r$B3)

  out <- (r$B8 - r$B3) / denom

  names(out) <- "gndvi"
  return(out)
}


# 2.3. Modified Chlorophyll Absorption Ratio Index (MCARI) --------------------

mcari <- function(r) {
  
  denom <- safe_denom(r$B4)
  
  out <- ((r$B5 - r$B4) - 0.2 * (r$B5 - r$B3)) *
    (r$B5 / denom)
  
  names(out) <- "mcari"
  return(out)
}


# 2.4. Modified Normalized Difference Water Index (MNDWI) ---------------------

mndwi <- function(r) {
  
  denom <- safe_denom(r$B3 + r$B11)
  
  out <- (r$B3 - r$B11) / denom
  
  names(out) <- "mndwi"
  return(out)
}


# 2.5. Normalized Difference Aquatic Vegetation Index (NDAVI) -----------------

ndavi <- function(r) {

  denom <- safe_denom(r$B8 + r$B2)

  out <- (r$B8 - r$B2) / denom

  names(out) <- "ndavi"
  return(out)
}


# 2.6. Normalized Difference Moisture Index (NDMI) ----------------------------

ndmi <- function(r) {
  
  denom <- safe_denom(r$B8A + r$B11)
  
  out <- (r$B8A - r$B11) / denom
  
  names(out) <- "ndmi"
  return(out)
}


# 2.7. Normalized Difference Vegetation Index (NDVI) --------------------------

ndvi <- function(r) {

  denom <- safe_denom(r$B8 + r$B4)

  out <- (r$B8 - r$B4) / denom

  names(out) <- "ndvi"
  return(out)
}


# 2.8. Soil Adjusted Vegetation Index (SAVI) ----------------------------------

savi <- function(r) {

  denom <- safe_denom(r$B8 + r$B4 + 0.5)

  out <- 1.5 * ((r$B8 - r$B4) / denom)

  names(out) <- "savi"
  return(out)
}


# 2.9. Visible Atmospherically Resistant Index (VARI) -------------------------

vari <- function(r) {
  
  denom <- safe_denom(r$B3 + r$B4 - r$B2)
  
  out <- (r$B3 - r$B4) / denom
  
  names(out) <- "vari"
  return(out)
}


# 2.10. Water Adjusted Vegetation Index (WAVI) --------------------------------

wavi <- function(r) {

  denom <- safe_denom(0.5 + r$B8 + r$B2)

  out <- (r$B8 - r$B2) / denom

  names(out) <- "wavi"
  return(out)
}