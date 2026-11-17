library(terra)


files = list.files(path = "C://Users/bdmor/Documents/persistence/data/raw/", "img|tif", include.dirs = T, full.names = T)

# only need 2017 onward
files = files[9:16]

years = 2017:2024

ref = rast("E://sentinel2_harmonized/raw_mosaics/spectra/S2A_SR_HARMONIZED_Bands_20260803T185831.tif", lyrs = 1)



rip_stack = rast()
for (i in 1:length(files)){
  cat("working on: ", years[i], "\n")
  
  r = rast(files[i])

  # 9 == riparian class
  out = as.numeric(r == 9)
  
  # make sure file matches sentinel-2 mosaic
  rip = resample(out, ref, method = 'near', threads = T)
  names(rip) = paste0("rip_", years[i])
  
  rip[is.na(rip[])] = 0
  
  # add to rip_stack
  rip_stack = c(rip_stack, rip)
}


# add data together to see how much rip is persistent throught the years
rip_sum = tapp(rip_stack, index = rep(1, 8), fun = sum)

# sums of 0 mean riparian was never classed at a pixel for all years
rip_sum2 = rip_sum
NAflag(rip_sum2) = 0


# need a raster to know where areas were always classed for all years (e.g., not missing flightlines)
common = ref*0
common[is.na(common[])] = 0

for (i in 1:length(files)){
  r = rast(files[i])*1
  r = project(r, ref, method = "near", threads = T)
  r[r[] > 0] = 1
  
  common = sum(common, r, na.rm = T)
}

# common values of 8 mean there was data collected for all 8 years
compare = as.numeric(common == 8)

# mask out pixels that weren't present all 8 years
rip_sum3 = mask(rip_sum2, compare, maskvalues = 0, updatevalue = NA)

# look into changes between rip
quantile(rip_sum3[], na.rm = T)

# 50% of rip pixels were only present at the same pixels 2 of 8 years
vals = values(rip_sum3, dataframe = T, na.rm = T)$X1
hist(vals)

rip_mask = as.numeric(rip_sum3 >= 4)
NAflag(rip_mask) = 0

writeRaster(rip_mask, "E://sentinel2_harmonized/riparian_mask/rip_mask_allYears.tif")
