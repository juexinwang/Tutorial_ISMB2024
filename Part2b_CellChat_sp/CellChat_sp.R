#  Single spatial transcriptomics Analysis using CellChat

###############################################################################
#### Part I: Data Input & Processing and Initialization of CellChat Object ####
###############################################################################

# Load in required libraries
ptm = Sys.time()
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Load data
load("data/visium_mouse_cortex_annotated.RData")

# Show the image and annotated spots
color.use <- scPalette(nlevels(visium.brain)); names(color.use) <- levels(visium.brain)
Seurat::SpatialDimPlot(visium.brain, label = T, label.size = 3, cols = color.use)

# Use normalized data matrix
data.input = Seurat::GetAssayData(visium.brain, slot = "data", assay = "SCT") 

# Define the meta data: 
# manually create a dataframe consisting of the cell labels
meta = data.frame(labels = Seurat::Idents(visium.brain), samples = "sample1", row.names = names(Seurat::Idents(visium.brain))) 

meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
unique(meta$samples) # check the sample labels


# load spatial transcriptomics information
spatial.locs = Seurat::GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 

# Spatial factors of spatial coordinates
scalefactors = jsonlite::fromJSON(txt =  'data/scalefactors_json.json')
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) # this value should approximately equal 100um for 10X Visium data




## Create CellChat Object using the loaded data ##
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat


## Set the Ligand-Receptor Interaction Database ##
# When analyzing mouse samples, we use the database 'CellChatDB.mouse'
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use



## Preprocessing the Expression Data for Cell-Cell Communication Analysis ##
# subset the expression data of signaling genes for saving computation cost
ptm = Sys.time()
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.mouse)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

###############################################################################
###### Part II: Inference of Cell-Cell Communication Network #########
###############################################################################
ptm = Sys.time()

# Compute the Communication Crobability and Infer Cellular Communication Network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                              contact.dependent = FALSE, contact.range = 100, nboot = 20)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#saveRDS(cellchat, file = "data/cellchat_inferred_object.rds")
#cellchat_alt <- readRDS("data/cellchat_inferred_object.rds")


# filter out the cell-cell communication if there are only few cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


# Extra: 'subsetCommunication' to easily access the inferred cell-cell communications of interest
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))


# Infer the Cell-Cell Communication at a Signaling Pathway Level
cellchat <- computeCommunProbPathway(cellchat)



# Calculate the Aggregated Cell-Cell Communication Network
cellchat <- aggregateNet(cellchat)


# Visualize aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions",  arrow.size = 0)
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",  arrow.size = 0)

# Do heatmap based on a single object, for both counts and weights
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")


###############################################################################
#### Part III: Visualization of Cell-Cell Communication Network ####
###############################################################################

# Visualization of Cell-Cell Communication at Different Levels
pathways.show <- c("IGF") 

# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)


# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Spatial Plot - Incoming signaling weight
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3.5)


# Visualize gene expression distribution - Take an input of a ligand-receptor pair
spatialFeaturePlot(cellchat, pairLR.use = "IGF1_IGF1R", point.size = 2, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

# Show same expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "IGF1_IGF1R", point.size = 2, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

# Save Cellchat object
saveRDS(cellchat, file = "data/cellchat_visium_mouse_cortex.rds")

runCellChatApp(cellchat)
