library(terra)


# raw dem file
dem = rast("E://sentinel2_harmonized/riparian_mask/dem_bay_delta_10m_20250312/dem_bay_delta_10m_20250312.tif")

# boundary shapefile of area to make class maps
bound = vect("E://sentinel2_harmonized/processing_boundary/boundary_w_ss_yolo_FINAL.shp")

# reference sentinel-2 image 
ref = rast("E://sentinel2_harmonized/raw_mosaics/spectra/S2C_SR_HARMONIZED_Bands_20261010T185301.tif", lyrs = 1)


# resample/align dem to sentinel-2 mosaic
dem2 = project(dem, ref, method = 'bilinear', threads = T)

# crop dem to the same extent as mosaic
dem3 = crop(dem2, ref)

# mask dem to the boundary of the study area
dem_mask = mask(dem3, bound)
writeRaster(dem_mask, "E://sentinel2_harmonized/riparian_mask/dem_bay_delta_10m_20250312_aligned2S2.tif")

# there are holes throughout the dem. Fill in areas where 2m dem's occur
files = list.files(path = "E://sentinel2_harmonized/riparian_mask/2m_dem/", ".tif", include.dirs = T, full.names = T)

# make a raster collection to merge the different 2m dems
r = lapply(files, rast)
r = sprc(r)
r = merge(r, na.rm = T)

# resample and align to ref s2 mosaic
r2 = project(r, ref, method = 'bilinear', threads = T)
r2 = crop(r2, ref)

# identify pixels where dem_mask is NA but r2 has a value
x = which(is.na(dem_mask[]) & !is.na(r2[]))

dem_fill = dem_mask
dem_fill[x] = r2[x]

dem_fill = mask(dem_fill, bound)

# make reference raster to keep track of where there were no values --> for uncertainty calcs later
dummy = dem_fill
dummy[is.na(dummy[])] = -9999
dummy = mask(dummy, bound)
writeRaster(dummy, "E://sentinel2_harmonized/riparian_mask/dem_bay_delta_10m_20250312_withHoles.tif")


writeRaster(dem_fill, "E://sentinel2_harmonized/riparian_mask/dem_bay_delta_10m_20250312_2mHoleFill.tif")
