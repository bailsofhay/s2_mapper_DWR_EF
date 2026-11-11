library(terra)
library(stringr)
library(dplyr)

# mosaicked spetral bands
indir = "E://sentinel2_harmonized/raw_mosaics/spectra/"
setwd(indir)

# output directory for indicies to be calculated
outdir = "E://sentinel2_harmonized/raw_mosaics/indices/"

# list of input mosaic files
mos_files = list.files(indir, ".tif")

# dates of mosaic files
dates = str_extract(mos_files, "_[0-9]{8}") %>%
  str_remove("_")


# only want to calculate indices for mosaics that have not had them calculated yet



for (i in 1:length(mos_files)){
  r = rast(mos_files[i])
  
  outname = gsub("Bands", "Indicices", mos_files[i])
  
  # test = r
  # NAflag(test) = 0
  
  # Indices to calculate: mcari + mndwi + ndmi + vari
  mcari = ((r$B5-r$B4)-0.2*(r$B5-r$B3))*(r$B5/r$B4)
  names(mcari) = "mcari"
  cat("\t", "Completed: mcari", "\n")
  
  mndwi = (r$B3-r$B11)/(r$B3+r$B11)
  names(mndwi) = "mndwi"
  cat("\t", "Completed: mndwi", "\n")
  
  ndmi = (r$B8A-r$B11)/(r$B8A+r$B11) #normalized difference moisture index
  names(ndmi) = "ndmi"
  cat("\t", "Completed: ndmi", "\n")
  
  
  vari = (r$B3-r$B4)/(r$B3+r$B4-r$B2)
  names(vari) = "vari"
  cat("\t", "Completed: vari", "\n")
  
  ind = c(mcari, mndwi, ndmi, vari)
  
  writeRaster(ind,
              paste0(outdir, outname))
}