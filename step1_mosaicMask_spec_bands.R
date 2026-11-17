
# =============================================================================
# Sentinel-2 Harmonized: Mosaic and Mask Spectral Imagery
# =============================================================================
#
# Description:
#   Processes Sentinel-2 Harmonized imagery exported from Google Earth Engine
#   (GEE) for aquatic vegetation classification.
#
#   Each GEE export directory is expected to contain five GeoTIFF files:
#     - Four tiled, multiband spectral images
#     - One raster mask containing SCL and Cloud Score+ screening results
#
# Workflow:
#   1. Identify newly downloaded Sentinel-2 image directories.
#   2. Archive GeoTIFF files in the raw_spectra directory.
#   3. Load the spectral tiles and corresponding mask.
#   4. Mosaic the spectral tiles into a single multiband raster.
#   5. Apply the raster mask to the mosaicked imagery.
#   6. Replace zero reflectance values with 0.001.
#   7. Export the masked mosaic to the raw_mosaics directory.
#
# Notes:
#   - Export directories are removed only after successful file copying.
#   - Archived raw imagery is retained for reproducibility.
#   - Zero is a valid spectral reflectance value and is not treated as NoData.
#   - Zero reflectance values are replaced with 0.001 to support index
#     calculations and avoid division-by-zero errors.
#   - Masked pixels remain NA.
#   - All spectral tiles and masks must share a compatible raster grid.
#
# =============================================================================


# 1. Load Required Packages ---------------------------------------------------

library(terra)
library(stringr)


# 2. Define Input and Output Directories --------------------------------------

# Parent directory containing downloaded Sentinel-2 imagery
indir <- "E:/sentinel2_harmonized"

# Directory for archiving original GEE exports
raw_dir <- file.path(indir, "raw_spectra")

# Directory for processed spectral mosaics
outdir <- file.path(indir, "raw_mosaics", "spectra")

# Create output directories if they do not exist
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# 3. Identify Unprocessed Sentinel-2 Images -----------------------------------

# Identify GEE export directories that have not yet been archived
dirs <- list.dirs(indir, recursive = FALSE, full.names = FALSE)
dirs <- dirs[grepl("^S2[ABC]_SR_HARMONIZED_", dirs)]

cat("Found", length(dirs), "Sentinel-2 directories to process.\n")


# 4. Process Sentinel-2 Images ------------------------------------------------

for (drivedir in dirs) {
  
  cat("\nWorking on:", drivedir, "\n")
  
  
  # 4.1. Identify and Validate Input Files ------------------------------------
  
  # Each export should contain four spectral tiles and one mask
  input_dir <- file.path(indir, drivedir)
  
  files <- list.files(
    path = input_dir,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  band_files <- files[grepl("spectra", basename(files))]
  mask_files <- files[grepl("mask", basename(files))]
  
  # Skip incomplete exports --> GEE may have failed processing/exporting if 
  # there are not 5 files total.
  if (length(files) != 5 ||
      length(band_files) != 4 ||
      length(mask_files) != 1) {
    
    warning("Incomplete or unexpected files in: ", drivedir)
    next
  }
  
  
  # 4.2. Archive Original Sentinel-2 Files ------------------------------------
  
  # Move a copy of each export into raw_spectra for permanent storage
  archive_dir <- file.path(raw_dir, drivedir)
  
  dir.create(archive_dir, showWarnings = FALSE)
  
  archive_files <- file.path(archive_dir, basename(files))
  
  # Copy files and verify successful transfer
  copied <- file.copy(
    from = files,
    to = archive_files,
    overwrite = FALSE
  )
  
  if (!all(copied)) {
    warning("File copying failed for: ", drivedir)
    next
  }
  
  # Remove the original download directory after verifying the copies
  # Archived files remain available for future reprocessing
  unlink(input_dir, recursive = TRUE, force = TRUE)
  
  # Update file paths to archived copies
  band_files <- archive_files[grepl("spectra", basename(archive_files))]
  mask_files <- archive_files[grepl("mask", basename(archive_files))]
  
  
  # 4.3. Define Output Filename -----------------------------------------------
  
  # Preserve the existing GEE-based naming convention for mosaics
  specname <- basename(band_files[1]) |>
    substr(1, 47) |>
    paste0(".tif") |>
    str_remove("spectral")
  
  outfile <- file.path(outdir, specname)
  
  
  # 4.4. Load Raster Mask -----------------------------------------------------
  
  # Mask incorporates SCL and Cloud Score+ screening from GEE
  mask_rast <- rast(mask_files)
  
  # Interpret zero-valued mask pixels as NoData
  NAflag(mask_rast) <- 0
  
  
  # 4.5. Mosaic Spectral Band Tiles -------------------------------------------
  
  cat("\t", "Mosaicking spectral tiles...\n")
  
  # Load each multiband spectral tile
  # Zero reflectance values are retained during mosaicking
  r <- lapply(band_files, rast)
  
  # Combine spectral tiles into a SpatRasterCollection
  r <- sprc(r)
  
  # Mosaic overlapping tiles into one multiband raster
  out <- mosaic(r)
  
  # Note: Why the code does not set spectral band pixels with value of 0 to NA.
  # COPERNICUS/S2_SR_HARMONIZED, a value of 0 in the spectral reflectance 
  # bands represents an actual scaled reflectance value of 0 (or a negative 
  # value clipped to 0), rather than NA (NoData).
  
  # In Google Earth Engine (GEE), true NoData areas (such as pixels outside the 
  # satellite swath or granule footprint) are natively masked out via the 
  # image's internal validity mask rather than stored as unmasked 0s.
  
  # Spectral Reflectance Bands (B1–B12, B8A). These bands are stored as UINT16 
  # integers scaled by 10000. A value of 0 represents a valid pixel where 
  # surface reflectance is 0.0. Because the collection harmonizes post-January 
  # 2022 data (PROCESSING_BASELINE >= 04.00) by subtracting the 1000 DN offset 
  # and clamping negative results to the UINT16 floor, dark surfaces or 
  # over-corrected atmospheric pixels can legitimately hit 0
  
  # Scene Classification Map (SCL). Unlike the spectral bands, 0 in the raw 
  # Copernicus SCL specification stands for "No Data". According to the Earth 
  # Engine Data Catalog, GEE automatically masks out SCL == 0 upon ingestion, 
  # leaving valid SCL values ranging from 1 to 11.
  
  # https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED
  # https://clearsky.vision/knowledge/sentinel2-scaling-harmonization
  
  
  # 4.6. Apply Raster Mask ----------------------------------------------------
  
  cat("\t", "Applying SCL / Cloud Score+ mask...\n")
  
  # Exclude pixels flagged as NoData in the raster mask
  out <- mask(out, mask_rast)
  
  # Optional visual quality check (RGB composite)
  # plotRGB(out, r = 3, g = 2, b = 1, stretch = "lin")
  
  
  # # 4.7. Replace Zero Reflectance Values --------------------------------------
  # DO NOT RUN THIS!
  # cat("\t", "Replacing zero reflectance values...\n")
  # 
  # # Replace valid zero reflectance values with a small positive value
  # # to prevent division-by-zero errors when calculating spectral indices.
  # # Existing NA values are preserved.
  # out <- ifel(out == 0, 0.001, out)
  
  
  # 4.8. Export Processed Spectral Mosaic -------------------------------------
  
  cat("\t", "Writing masked spectral mosaic...\n")
  
  writeRaster(
    out,
    filename = outfile,
    overwrite = FALSE
  )
  
  cat("Completed:", drivedir, "\n")
}

cat("\nSentinel-2 mosaic processing finished.\n")
