if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBImage")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSMB")

library(BiocManager)
library(EBImage)
library(MSMB)


imagefile = system.file("images", "mosquito.png",
                        package = "MSMB")
mosq = readImage(imagefile)

display(mosq)

display(mosq, method = "raster")
text(x = 85, y = 800, label = "A mosquito",
     adj = 1.5, col = "orange", cex = 1.5)

imagefile = system.file("images", "hiv.png",
                        package = "MSMB")
hivc = readImage(imagefile)
display(hivc)

dim(hivc)
dim(mosq)

nuc = readImage(system.file("images", "nuclei.tif",
                            package = "EBImage"))
display(1- nuc, method = "raster", all = TRUE) #1-nuc exchanges the color of the white and black

display(1 - nuc, method = "raster", frame = 2)


class(mosq)
dim(mosq)

imageData(mosq)[1:3, 1:6]

class(nuc)
dim(imageData(nuc))

writeImage(hivc,'hivc.jpeg',quality=85)

mosqinv = normalize(-mosq)
display(mosqinv )

mosqcont=mosq*3
display(mosqcont)

mosqexp=mosq^(1/3)
display(mosqexp)

mosqcrop   = mosq[100:438, 112:550]
display(mosqcrop)

mosqthresh = mosq > 0.5
display(mosqthresh)
class(mosqthresh)
mosqthresh

mosqtransp = transpose(mosq)
display(mosqtransp)
display(t(mosq))

mosqrot   = EBImage::rotate(mosq, angle = 30)
display(mosqrot)

mosqshift = translate(mosq, v = c(40, 70))
display(mosqshift)

mosqflip  = flip(mosq)
mosqflop  = flop(mosq)

################

imagefiles = system.file("images", c("image-DAPI.tif", "image-FITC.tif", "image-Cy3.tif"), package="MSMB")
cells = readImage(imagefiles)
display(cells)
hist(cells[,,3])
display(cells[,,1]>0.1)
display(cells[,,1]>0.2)

apply(cells,3,range)

cells[,,1]   = 32 * cells[,,1]
cells[,,2:3] = 16 * cells[,,2:3]
apply(cells, 3, range)#after the transformation, the cells become visible
display(cells)

w = makeBrush(size = 51, shape = "gaussian", sigma = 7)
nucSmooth = filter2(getFrame(cells, 1), w)
display(w)


cellsSmooth = Image(dim = dim(cells))
sigma = c(1, 3, 3)
for(i in seq_along(sigma))
  cellsSmooth[,,i] = filter2( cells[,,i],
                              filter = makeBrush(size = 51, shape = "gaussian",
                                                 sigma = sigma[i]) )
display(cellsSmooth)

imageData(cellsSmooth[1,1:3,1:3])

py = seq(-1, +1, length.out = dim(cellsSmooth)[1])
px = seq(-1, +1, length.out = dim(cellsSmooth)[2])
illuminationGradient = Image(
  outer(py, px, function(x, y) exp(-(x^2+y^2))))
nucBadlyIlluminated = cellsSmooth[,,1] * illuminationGradient
display(nucBadlyIlluminated)

disc = makeBrush(21, "disc")# solution for the nucBadlyIlluminated
head(disc)
display(disc)
disc = disc / sum(disc)
display(disc)
localBackground = filter2(nucBadlyIlluminated, disc)#performed 2-d convolution on nBI and disc
display(localBackground)
offset = 0.02#how to decide off set
nucBadThresh = (nucBadlyIlluminated - localBackground > offset)
display(nucBadThresh)

display(nucBadlyIlluminated)

nucThresh =
  (cellsSmooth[,,1] - filter2(cellsSmooth[,,1], disc) > offset)

display(nucThresh)

display(filter2(cellsSmooth[,,1], disc))
#morphological operation
nucOpened = EBImage::opening(nucThresh,
                             kern = makeBrush(5, shape = "disc"))#why disc
display(nucOpened)

#segmentation
nucSeed = bwlabel(nucOpened)
table(nucSeed)

display(colorLabels(nucSeed))

nucMask = cellsSmooth[,,1] - filter2(cellsSmooth[,,1], disc) > 0
display(nucMask)
display(cellsSmooth[,,1])
nucMask = fillHull(nucMask)

display(nucMask)

nuclei = propagate(cellsSmooth[,,1], nucSeed, mask = nucMask)
display(nuclei)



#Voronoi tessellation
zeros        = Image(dim = dim(nuclei))#generate black image
display(zeros)
voronoiExamp = propagate(seeds = nuclei, x = zeros, lambda = 100)#lambda=0
voronoiPaint = paintObjects(voronoiExamp,  1-nucOpened)

display(voronoiPaint)

table(voronoiExamp)
ind = which(voronoiExamp == 13, arr.ind = TRUE)#?

#11.13
hist(log(cellsSmooth[,,3]) )
hist(log(cellsSmooth[,,3]), xlim = -c(3.6, 3.1), breaks = 300)

library("genefilter")
bgPars = function(x) {
  x    = log(x)#why take log
  loc  = half.range.mode( x )
  left = (x - loc)[ x < loc ]#return numbers that is smaller than loc in (X-loc)
  wid  = sqrt( mean(left^2) )#the standard deviation??
  c(loc = loc, wid = wid, thr = loc + 6*wid)
}
cellBg = apply(cellsSmooth, MARGIN = 3, FUN = bgPars) #margin=3 means the third graph?
cellBg

hist(log(cellsSmooth[,,3]), xlim = -c(3.6, 3.1), breaks = 300)
abline(v = cellBg[c("loc", "thr"), 3], col = c("brown", "red"))

cytoplasmMask = (cellsSmooth[,,2] > exp(cellBg["thr", 2])) |
  nuclei | (cellsSmooth[,,3] > exp(cellBg["thr", 3]))

cellbodies = propagate(x = cellsSmooth[,,3], seeds = nuclei,
                       lambda = 1.0e-2, mask = cytoplasmMask)

display(cytoplasmMask)

makeBrush(5,shape="disc")






hist(log(cellsSmooth[,,3]), xlim = -c(3.6, 3.1), breaks = 300)
abline(v = cellBg[c("loc", "thr"), 3], col = c("brown", "red"))

cytoplasmMask = (cellsSmooth[,,2] > exp(cellBg["thr", 2])) |
  nuclei | (cellsSmooth[,,3] > exp(cellBg["thr", 3]))
display(cytoplasmMask)

cellbodies = propagate(x = cellsSmooth[,,3], seeds = nuclei,
                       lambda = 1.0e-2, mask = cytoplasmMask)
display(cellbodies)


cellsColor = rgbImage(red   = cells[,,3],
                      green = cells[,,2],
                      blue  = cells[,,1])

nucSegOnNuc  = paintObjects(nuclei, tgt = toRGB(cells[,,1]),
                            col = "#ffff00")
display(nucSegOnNuc )
nucSegOnAll  = paintObjects(nuclei, tgt = cellsColor,
                            col = "#ffff00")
display(nucSegOnAll)
cellSegOnAll = paintObjects(cellbodies, tgt = nucSegOnAll,
                            col = "#ff0080")
display(cellSegOnAll)

#feature extraction
meanNucInt       = tapply(cells[,,1], nuclei, mean)
meanActIntInNuc  = tapply(cells[,,3], nuclei, mean)
meanActIntInCell = tapply(cells[,,3], cellbodies, mean)


F1 = computeFeatures(nuclei,     cells[,,1], xname = "nuc",
                     refnames = "nuc")
F2 = computeFeatures(cellbodies, cells[,,2], xname = "cell",
                     refnames = "tub")
F3 = computeFeatures(cellbodies, cells[,,3], xname = "cell",
                     refnames = "act")
dim(F1)

class(cellbodies)

F1






