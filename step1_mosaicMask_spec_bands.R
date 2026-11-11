library(terra)
library(dplyr)
library(stringr)
library(lubridate)

indir = "E://sentinel2_harmonized/"
setwd(indir)

outdir = "E://sentinel2_harmonized/"


dirs = list.files(path = "E://sentinel2_harmonized/",
                  "HARMONIZED_")

for (x in 1:length(dirs)){
  drivedir = dirs[x] # folders for specific Sentinel-2 image dates
  cat("Working on", drivedir, "\n")
  

  # should be 5 files: 4 spectral band tiles, 1 mask file
  files = list.files(paste(getwd(), drivedir, sep = "/"), 
                     ".tif$", include.dirs = T, full.names = T)
  
  # if there are < 5 files then something went wrong
  if (length(files) == 5){
    # make directory for s2 date in raw_spectra folder
    dir_out = paste(outdir, "raw_spectra", drivedir, sep = "/")
    dir.create(dir_out)
    
    # copy S2 date .tif files to the new folder
    copy_out = paste(dir_out, basename(files), sep = "/")
    file.copy(from = files, to = copy_out)
    
    # helps keep track of whats been mosaicked and what hasn't
    unlink(drivedir, recursive = TRUE, force = T) # remove the old folder and data

    # spectral tifs
    band_files = copy_out[grep("spectra", files)]
    # mask tif --> includes SCL + CloudScore+ mask
    mask_files = copy_out[grep("mask", files)]
    
    # Basename for output mosaics
    specname = basename(band_files[1]) %>%
      substr(1, 47) %>%
      paste0(".tif") %>%
      str_remove("spectral")
    
    # date of the s2 image
    date = basename(band_files[1]) %>%
      str_extract("[0-9]{8}T[0-9]{6}")%>%
      substr(1, 8) %>%
      as_date(tz = NULL)%>%
      as.character()
    
    print(date)
    
    # Load in the Mask raster
    mask = rast(mask_files)
    NAflag(mask) = 0 # set no data flag
    # plot(mask)
    
    
    # set up spectral band tiles for mosaic
    cat("Mosaicking spectral bands", "\n")
    r = lapply(band_files, rast)
    r = sprc(r)
    
    # mosaic spectral band tiles
    out = mosaic(r)
    NAflag(out) = 0

    # Apply mask to mosaicked spectral bands
    cat("Masking spectral bands", "\n")
    out = mask(out, mask) 
    
    
    # if needed, plot RGB of s2 image for sanity
    # plotRGB(out, 3, 2, 1, stretch = "lin")

    # write mosaicked/masked s2 image to geotiff file
    outfile = paste(outdir, "raw_mosaics", 'spectra', specname, sep = "/")
    writeRaster(out, outfile)
    
    cat("Completed mosaic and masking for", drivedir, "\n")
    
  } else {
    cat(drivedir, "is missing files. Please check.", "\n")
  }

}

