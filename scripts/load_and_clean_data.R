############################## Load required libraries #####################################
rm(list = ls())
suppressPackageStartupMessages({
  library(simpleSeg)
  library(cytomapper)
  library(dplyr)
  library(tidyverse)
})

############################## Set up parallell processing #####################################
nCores <- 40
BPPARAM <- simpleSeg:::generateBPParam(nCores)

############################## Read clinical information #####################################
clinicalData <- readxl::read_xlsx("../Data/oSCC_ICT_TMA1-2_2022.xlsx") |> as.data.frame()

clinicalData$`P16pos` <- clinicalData$`P16+ve`
clinicalData$TMA <- substr(clinicalData$Sample,1,2)
clinicalData$TMA <- substr(clinicalData$Sample,1,2)
clinicalData$Clean <- gsub("\\(A/B\\)", "", clinicalData$Sample)
clinicalData$Clean <- gsub(" ", "", clinicalData$Clean)
clinicalData$R1 <- unlist(lapply(strsplit(clinicalData$Clean,"/"), function(x)x[1]))
clinicalData$Patient <- paste0("Patient", seq_len(nrow(clinicalData)))
clinicalData$R2 <- paste0(clinicalData$TMA,
                          unlist(lapply(strsplit(clinicalData$Clean,"/"), function(x)x[2])))
clinicalData <- tidyr::pivot_longer(clinicalData, col = c("R1", "R2"), 
                                    names_to = "Replicate", values_to = "sampleID")
loc_dat <- data.table::fread('oSCC_ICT_TMA1-2_2022.csv')
colnames(loc_dat)[[4]] <- 'Location'
clinicalData <- merge(clinicalData, loc_dat)

############################## Read in imges #####################################
# The images are stored in the `images` folder within the `Data` folder. 
# Here we use `readImages()` from the `EBImage` package to read these into R. 
# If memory is a restricting factor, and the files are in a slightly different 
# format, you could use `loadImages()` from the `cytomapper` package to load all 
# of the tiff images into a `CytoImageList` object, which can store the images as h5 on-disk.
pathToImages <- "../Data/TMA1&2_tiffs"
# Get directories of images
imageDirs <- dir(pathToImages, full.names = TRUE)
names(imageDirs) <- dir(pathToImages, full.names = FALSE)
# Get files in each directory
files <- sapply(imageDirs, list.files, pattern = "tif", full.names = TRUE, simplify = FALSE)
# Read files with readImage from EBImage
images <- lapply(files, EBImage::readImage, as.is = TRUE)

images <- lapply(images, function(x){
  if("172Yb_172Yb_antiCy5.ome" %in% dimnames(x)[[3]]){
    dimnames(x)[[3]][dimnames(x)[[3]]=="172Yb_172Yb_antiCy5.ome"] <- "172Yb_172Yb_antiCy5_p40.ome"
  }
  x
})

# We will make use of the `on_disk` option to convert our images to a 
# `CytoImageList` with the images not held in memory.
# Store images in a CytoImageList with images on_disk as h5 files to save memory. 
dir.create("../Data/h5Files")
images <- cytomapper::CytoImageList(images, 
                                    on_disk = TRUE, 
                                    h5FilesPath = "../Data/h5Files", 
                                    BPPARAM = BPPARAM)
gc()

###################### Combine clinical data with images ######################
mcols(images) <- data.frame(imageID = names(images), sampleID = substr(names(images),1,4))
mcols(images) <-  dplyr::left_join(as.data.frame(mcols(images)), clinicalData)

###################### Get the marker information ######################
cytomapper::channelNames(images) <- gsub("\\.ome", "", cytomapper::channelNames(images))
cytomapper::channelNames(images) <- gsub("_antiCy5", "", cytomapper::channelNames(images))
n <- unlist(lapply(strsplit(cytomapper::channelNames(images), "_"), length))
cytomapper::channelNames(images) <- unlist(lapply(strsplit(cytomapper::channelNames(images), "_"), 
                                                  function(x)x[length(x)]))
markers <- cytomapper::channelNames(images)[n==3]

###################### Generate Segmentation masks ######################
masks <- simpleSeg(images,
                   nucleus = c("HH3", "DNA1", "DNA2"),
                   cellBody = "dilate",
                   transform = "sqrt",
                   sizeSelection = 40,
                   discSize = 3,
                   cores = nCores)

######################  Summarise cell features ######################
cells <- cytomapper::measureObjects(masks, 
                                    images, 
                                    img_id = "imageID", 
                                    BPPARAM = BPPARAM)

######################  Summarise cell features ######################
# Transform and normalise the marker expression of each cell type.
cells <- normalizeCells(cells,
                        transformation = "sqrt",
                        method = c("trim99","mean","PC1"),
                        assayIn = "counts",
                        markers = markers,
                        cores = nCores)

######################  Filter images based on responder ######################
useImages = unique(cells$sampleID[!is.na(cells$Responder)])
cells <- cells[,cells$imageID %in% useImages]

save(cells, markers, clinicalData, images, masks,
     file = "oSCC_cells_normalised.RData")


