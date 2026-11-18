library(terra)

dem = rast("E://sentinel2_harmonized/riparian_mask/dem_bay_delta_10m_20250312_2mHoleFill.tif")

# raster of where elevation is >= 2.6m and can be reclassed (if veg) to riparian 
# class
rip = as.numeric(dem >= 2.6)
NAflag(rip) = 0

writeRaster(rip, "E://sentinel2_harmonized/riparian_mask/dem_over_threshold_binary.tif")
