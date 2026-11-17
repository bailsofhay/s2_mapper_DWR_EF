# =============================================================================
# Sentinel-2 Harmonized: Calculate Spectral Indices
# =============================================================================
#
# Description:
#   Calculates spectral indices from mosaicked and masked Sentinel-2
#   Harmonized surface reflectance imagery for aquatic vegetation mapping.
#
#   Spectral index equations are defined in the external index_functions.R
#   script and sourced before processing.
#
# Workflow:
#   1. Load spectral index functions.
#   2. Identify Sentinel-2 spectral mosaics without corresponding index files.
#   3. Calculate 10 spectral indices in alphabetical order.
#   4. Stack the indices into a multiband raster.
#   5. Export the index raster to the raw_mosaics/indices directory.
#
# Notes:
#   - Input imagery contains surface reflectance scaled between 0 and 1.
#   - Input rasters must have Sentinel-2 band names (e.g., B2, B3, B4).
#   - Previously processed mosaics are identified by output filename.
#   - Output rasters contain 10 spectral index layers.
#
# =============================================================================


# 1. Load Required Packages and Functions ------------------------------------

library(terra)
library(stringr)

# Load spectral index functions
source("C://Users/bdmor/Documents/s2_mapper_DWR_EF/index_functions.R")


# 2. Define Input and Output Directories --------------------------------------

# Directory containing mosaicked and masked Sentinel-2 spectral imagery
indir <- "E:/sentinel2_harmonized/raw_mosaics/spectra"

# Directory for calculated spectral indices
outdir <- "E:/sentinel2_harmonized/raw_mosaics/indices"

# Create output directory if it does not exist
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# 3. Identify Unprocessed Sentinel-2 Mosaics ----------------------------------

# List available spectral mosaics
mos_files <- list.files(
  path = indir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# Generate expected index filenames
expected_files <- basename(mos_files) |>
  str_replace("Bands", "Indices")

# Keep only mosaics without corresponding index files
mos_files <- mos_files[
  !file.exists(file.path(outdir, expected_files))
]

cat("Found", length(mos_files),
    "Sentinel-2 mosaics requiring spectral indices.\n")


# 4. Calculate Spectral Indices -----------------------------------------------

for (mos_file in mos_files) {
  
  cat("\nWorking on:", basename(mos_file), "\n")
  
  # 4.1. Load Spectral Mosaic -------------------------------------------------
  
  # Load mosaicked Sentinel-2 surface reflectance imagery
  r <- rast(mos_file, lyrs = c(1:4, 7:9))
  
  
  # 4.2. Calculate Spectral Indices -------------------------------------------
  
  # # Enhanced Vegetation Index
  # cat("\t", "Calculating EVI\n")
  # evi_rast <- evi(r)
  
  # # Green Normalized Difference Vegetation Index
  # cat("\t", "Calculating GNDVI\n")
  # gndvi_rast <- gndvi(r)
  
  # Modified Chlorophyll Absorption Ratio Index
  cat("\t", "Calculating MCARI\n")
  mcari_rast <- mcari(r)
  
  # Modified Normalized Difference Water Index
  cat("\t", "Calculating MNDWI\n")
  mndwi_rast <- mndwi(r)
  
  # # Normalized Difference Aquatic Vegetation Index
  # cat("\t", "Calculating NDAVI\n")
  # ndavi_rast <- ndavi(r)
  
  # Normalized Difference Moisture Index
  cat("\t", "Calculating NDMI\n")
  ndmi_rast <- ndmi(r)
  
  # # Normalized Difference Vegetation Index
  # cat("\t", "Calculating NDVI\n")
  # ndvi_rast <- ndvi(r)
  
  # # Soil Adjusted Vegetation Index
  # cat("\t", "Calculating SAVI\n")
  # savi_rast <- savi(r)
  
  # Visible Atmospherically Resistant Index
  cat("\t", "Calculating VARI\n")
  vari_rast <- vari(r)
  
  # # Water Adjusted Vegetation Index
  # cat("\t", "Calculating WAVI\n")
  # wavi_rast <- wavi(r)
  
  
  # 4.3. Stack Spectral Indices -----------------------------------------------
  
  # Combine the 10 spectral indices into a multiband raster
  ind <- c(
    mcari_rast,
    mndwi_rast,
    ndmi_rast,
    vari_rast
  )
  
  # # Combine the 10 spectral indices into a multiband raster
  # ind <- c(
  #   evi_rast,
  #   gndvi_rast,
  #   mcari_rast,
  #   mndwi_rast,
  #   ndavi_rast,
  #   ndmi_rast,
  #   ndvi_rast,
  #   savi_rast,
  #   vari_rast,
  #   wavi_rast
  # )
  
  # Assign band names in alphabetical order
  names(ind) <- c(
    "mcari", "mndwi", "ndmi", "vari"
  )
  
  # # Assign band names in alphabetical order
  # names(ind) <- c(
  #   "evi", "gndvi", "mcari", "mndwi", "ndavi",
  #   "ndmi", "ndvi", "savi", "vari", "wavi"
  # )
  
  
  # 4.4. Export Spectral Indices ----------------------------------------------
  
  # Generate output filename using the spectral mosaic filename
  outname <- basename(mos_file) |>
    str_replace("Bands", "Indices")
  
  outfile <- file.path(outdir, outname)
  
  # Write 10-band index raster to GeoTIFF
  writeRaster(
    ind,
    filename = outfile,
    overwrite = FALSE
  )
  
  cat("Completed:", basename(mos_file), "\n")
}

cat("\nSpectral index processing finished.\n")
