library(terra)

# Configure GDAL curl options for S3 public bucket streaming
setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".vrt,.tif")

# Construct the full URL safely using paste0()
bucket <- "ctrees-tree-height-ca-2020.s3.us-west-2.amazonaws.com"
file   <- "ca-treeheight-utm-10.vrt"
vrt_url <- paste0("/vsicurl/https://", bucket, "/", file)

r <- rast(vrt_url)

# Reproject AOI and crop
aoi <- vect("E://sentinel2_harmonized/processing_boundary/boundary_w_ss_yolo_FINAL.shp")
aoi_utm10 <- project(aoi, r)

r_crop <- crop(r, aoi_utm10)

# Reference sentinel-2 mosaic image with correct res, crs, extent
ref = rast("E://sentinel2_harmonized/raw_mosaics/spectra/S2B_SR_HARMONIZED_Bands_20261005T185239.tif", lyrs = 3)

r_resamp = project(r_crop, ref, method = "bilinear", threads = T)

r_crop2 = crop(r_resamp, aoi)

r_mask = mask(r_crop2, aoi)

writeRaster(r_mask, "E://sentinel2_harmonized/riparian_mask/canopy_height_10m.tif")


# Note: 
# How to Cite: Sub-Meter Canopy Tree Height of California in 2020 by CTrees.org 
# was accessed on DATE from 
# https://registry.opendata.aws/ctrees-california-vhr-tree-height. 
# Fabien H. Wagner, Sophia Roberts, Alison L. Ritz, Griffin Carter, Ricardo 
# Dalagnol, Samuel Favrichon, Mayumi C.M. Hirye, Martin Brandt, Philippe Ciais, 
# Sassan Saatchi, Sub-meter tree height mapping of California using aerial 
# images and LiDAR-informed U-Net model, Remote Sensing of Environment, 2024, 
# https://doi.org/10.1016/j.rse.2024.114099.



