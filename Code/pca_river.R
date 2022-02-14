library(stpca)
library(SSN)
library(rgdal)
library(GWmodel)
library(geosphere)
library(ggplot2)
library(ggdendro)

# The structure of this code follows the structure of the paper. 

# ========================= #
#  load Geum River network  #
# ========================= #
shape_Geum <- readOGR("/home/kyu9510/pca_on_river/Data/ProcessedData/Geum.shp")

source("/home/kyu9510/pca_on_river/Code/source.R", chdir = TRUE)
binaryIDs_obj <- get_binaryIDs_stream(mouth_node=614 , shape_Geum) # defined in source.R, shape file needs RCH_ID, binaryID
# river mouth : node 614, called Geumgang Gapmun (금강갑문)

# get adjacency 
adjacency <- get_adjacency_stream(binaryIDs_obj) # defined in source.R 

shreve_order <- compute_shreve(adjacency) # defined in source.R 

# generate SSN obj
ssn_Geum <- importSSN_stream(shreve_obj = shreve_order, location="Full", multipleplaces=TRUE) 


# =============================== #
#  compute spatial weight matrix  #
# =============================== #
# afvcol="shreve" works well b/c shreve ratio = addfunccol ratio 
spatial.wt_Geum <- createWeightS_SC(ssndata=ssn_Geum, afvcol="shreve", adjacency=adjacency) 

colnames(spatial.wt_Geum) <-
  as.character(ssn_Geum@obspoints@SSNPoints[[1]]@point.data['X'][[1]])

rownames(spatial.wt_Geum) <-
  as.character(ssn_Geum@obspoints@SSNPoints[[1]]@point.data['X'][[1]])


# ========================== #
#  data load and preprocess  #
# ========================== #
data <- readRDS("/home/kyu9510/pca_on_river/Data/ProcessedData/Geum2020(extended).RDS")
# TN: Total nitrogen
# temp: Water temperature 
# ph: Hydrogen ion concentration
# elec: Electrical conductivity 
# O2: Dissolved oxygen
# BOD: Biochemical oxygen demand
# COD: Chemical oxygen demand
# susp: Suspended solid
# TP: Total phosphorus
# TOC: Total organic carbon
# flow: Flow volume
data$normaldata_TOC <- data$normaldata_TOC["2012-12/2020-11",] # extract from December, 2012 to November, 2020
head(data$normaldata_TOC)

# delete unused point (Samgacheon (삼가천), Seokcheon (석천), Mihocheon4 (미호천4))
# and Yudeungcheon A (유등천A) do not have observations)
remove.cand <- which(colnames(data$normaldata_TOC)%in%c("삼가천", "석천", "미호천4", "유등천A"))
data$normaldata_TOC <- data$normaldata_TOC[,-remove.cand ]

# unobserved value 0 -> NA
for(ii in 1:ncol(data$normaldata_TOC)){
  if(length(which(data$normaldata_TOC[,ii]==0))!=0){
    data$normaldata_TOC[which(data$normaldata_TOC[,ii]==0),ii] <- NA
  }
}

normaldata_TOC_log <- log(data$normaldata_TOC)

# check basics
data$eventplace$X[which(is.na(match(data$eventplace$X, colnames(data$normaldata_TOC))))]

dim(data$normaldata_TOC)

dim(data$eventplace)

# add day, week, month, year, season columns
normaldata_TOC_log <- as.data.frame(normaldata_TOC_log)
normaldata_TOC_log['day'] <- as.Date(rownames(as.data.frame(normaldata_TOC_log)), "%Y-%m-%d")
normaldata_TOC_log['month'] <- format(as.Date(rownames(as.data.frame(normaldata_TOC_log)), "%Y-%m-%d"), "%Y-%m")    

tmp <- sapply(strsplit(normaldata_TOC_log$month, "-"), "[[", 2)
tmp <- replace(tmp, tmp %in% c("03","04","05"), "spring")
tmp <- replace(tmp, tmp %in% c("06","07","08"), "summer")
tmp <- replace(tmp, tmp %in% c("09","10","11"), "fall")
tmp <- replace(tmp, tmp %in% c("12","01","02"), "winter")

normaldata_TOC_log['season'] <- paste(sapply(strsplit(normaldata_TOC_log$month, "-"), "[[", 1), tmp, sep='-')
normaldata_TOC_log['week'] <- strftime(normaldata_TOC_log$day, format = "%Y-%V")
normaldata_TOC_log['year'] <- format(as.Date(rownames(as.data.frame(normaldata_TOC_log)), "%Y-%m-%d"), "%Y")



# ============================== #
#  annual summer avg T-mode PCA  #
# ============================== #

# yearly average summer data
TOC_log_summer_year <- subset(aggregate(normaldata_TOC_log[endsWith(normaldata_TOC_log$season, "summer"),], 
                                       by = list(normaldata_TOC_log[endsWith(normaldata_TOC_log$season, "summer"),]$year), 
                                       FUN = mean, na.rm=T), select = -c(Group.1, day, month, season, week, year))


# first impute missing values
annually_avg_TOC_log_summer.Tmode.impute <- completeData(x = TOC_log_summer_year, pca.mode = "Tmode",
                                                         mean.plot = FALSE, estim.plot = FALSE,
                                                         mi.plot = FALSE, 
                                                         nbsim = 10, nboot = 10)


# transpose the completed data so that columns are monitoring sites and
# rows are time points
Y_avg_TOC_log_summer.Tmode.completedata <- 
  t(annually_avg_TOC_log_summer.Tmode.impute$res.MIPCA$res.imputePCA)

Y_avg_TOC_log_summer.Tmode.completedata <- 
  t(annually_avg_TOC_log_summer.Tmode.impute$res.MIPCA$res.imputePCA)

# make column order identical to the order of spatial.wt_Geum 
Y_avg_TOC_log_summer.Tmode.completedata <- 
  Y_avg_TOC_log_summer.Tmode.completedata[,colnames(spatial.wt_Geum)]

# check column names are in same format/order as spatial.wt_Geum
mean(colnames(Y_avg_TOC_log_summer.Tmode.completedata) == colnames(spatial.wt_Geum))

rownames(Y_avg_TOC_log_summer.Tmode.completedata) <- aggregate(normaldata_TOC_log[endsWith(normaldata_TOC_log$season, "summer"),], 
                                                              by = list(normaldata_TOC_log[endsWith(normaldata_TOC_log$season, "summer"),]$year), 
                                                              FUN = mean, na.rm=T)$Group.1

# rename from 1 to 127
colnames(Y_avg_TOC_log_summer.Tmode.completedata) <- 1:127
colnames(spatial.wt_Geum) <- 1:127
rownames(spatial.wt_Geum) <- 1:127

# colname : monitoring site
# 1 : 구량천 	2 : 가막 	3 : 진안천 	4 : 정자천 	5 : 용담 	6 : 용포 	7 : 무주남대천 	8 : 무주남대천-1 	9 : 무주남대천-2 	10 : 부리 	
# 11 : 제원 	12 : 봉황천 	13 : 제원A 	14 : 영동 	15 : 영동천1 	16 : 영동천2 	17 : 추풍령천 	18 : 초강1 	19 : 금계천 	20 : 초강2 	
# 21 : 이원 	22 : 옥천 	23 : 우산 	24 : 보청천1 	25 : 보청천2 	26 : 항건천 	27 : 보청천3 	28 : 보청천4 	29 : 상곡천 	30 : 추풍천 	
# 31 : 옥천천 	32 : 회인천 	33 : 주원천 	34 : 대청댐 	35 : 품곡천 	36 : 대청 	37 : 용호천 	38 : 현도 	39 : 두계천1 	40 : 봉곡2교 	
# 41 : 갑천1 	42 : 갑천2 	43 : 갑천3 	44 : 침산교 	45 : 유등천1 	46 : 옥계교 	47 : 대전천1 	48 : 대전천2 	49 : 대전천3 	50 : 유등천5 	
# 51 : 갑천4 	52 : 갑천5 	53 : 갑천5-1 	54 : 세종1 	55 : 청원-1 	56 : 백천 	57 : 미호천1 	58 : 칠장천 	59 : 미호천1-1 	60 : 한천 	
# 61 : 미호천2 	62 : 백곡천1 	63 : 백곡천2 	64 : 미호천7 	65 : 초평천 	66 : 미호천3 	67 : 보강천1 	68 : 보강천 	69 : 성암천 	70 : 무심천 	
# 71 : 무심천1 	72 : 무심천2 	73 : 무심천3 	74 : 석남천 	75 : 미호천8 	76 : 병천천A 	77 : 용두천 	78 : 병천천 	79 : 미호천5 	80 : 미호천5A 	
# 81 : 조천 	82 : 조천1 	83 : 월하천 	84 : 미호천6-1 	85 : 연기 	86 : 금남 	87 : 용수천-1 	88 : 용수천 	89 : 대교천 	90 : 대교천2 	
# 91 : 세종2 	92 : 공주1 	93 : 정안천 	94 : 곰나루 	95 : 금강 	96 : 유구천 	97 : 목면 	98 : 공주2 	99 : 부여 	100 : 지천 	
# 101 : 지천-1 	102 : 정동 	103 : 은산천 	104 : 부여1 	105 : 금천 	106 : 부여2 	107 : 석성천-1 	108 : 석성천 	109 : 석성천2 	110 : 성동 
# 111 : 논산천-1 	112 : 논산천1 	113 : 연산천 	114 : 노성천-1 	115 : 노성천 	116 : 논산천2 	117 : 방축천 	118 : 논산천4 	119 : 수철천  120 : 마산천 	
# 121 : 강경천 	122 : 강경 	123 : 산북천 	124 : 양화-1 	125 : 길산천 	126 : 길산천2 	127 : 금강갑문 	

# ================================ #
#  compute temporal weight matrix  #
# ================================ #
t_wt_rho <- vector(length = ncol(Y_avg_TOC_log_summer.Tmode.completedata))
for(i in 1:ncol(Y_avg_TOC_log_summer.Tmode.completedata)){
  t_wt_rho[i] <- abs(as.numeric(arima(Y_avg_TOC_log_summer.Tmode.completedata[,i], order = c(1,0,0), method = "ML")[[1]][1]))
}

summary(t_wt_rho)

annually_temporal_summer.wt_Geum <- createWeightT(n = nrow(Y_avg_TOC_log_summer.Tmode.completedata), rho = median(t_wt_rho))



# ============================ #
#  upstream distance d-matrix  #
# ============================ #
dmat <- matrix(0, nrow=nrow(spatial.wt_Geum), ncol=ncol(spatial.wt_Geum))
rownames(dmat) <- rownames(spatial.wt_Geum) 
colnames(dmat) <- colnames(spatial.wt_Geum)


updist_mat <- ssn_Geum@obspoints@SSNPoints[[1]]@network.point.coords
rownames(updist_mat) <- 1:127

for(i in 1 : nrow(dmat)){
  for(j in 1 : ncol(dmat)){
    if(spatial.wt_Geum[i,j]!=0){
      dmat[i,j] <- abs(updist_mat[colnames(dmat)[j],'DistanceUpstream'] -
                        updist_mat[colnames(dmat)[i],'DistanceUpstream'])
    }
    else{
      dmat[i,j] <- max(updist_mat['DistanceUpstream'])*100 # for calculation, set very large value instead of Inf
      # dmat[i,j] = Inf
    }
    
  }
}

# plot(ssn_Geum, "shreve", cex=1) # plot general plot
# plot(ssn_Geum, cex=0.8, col="red", xlab="lon", ylab="lat")


# =========================== #
#  flow-directed PCA results  #
# =========================== #
coords_Geum <- getSSNdata.frame(ssn_Geum, "Obs")[, c("경도.Degree.", "위도.Degree.")]
colnames(coords_Geum) <- c("lon","lat")

Y_avg_TOC_log_summer.Tmode.pca.all <- stpca(x = Y_avg_TOC_log_summer.Tmode.completedata,
                                            pca.mode = "Tmode", 
                                            spatial.wt = spatial.wt_Geum, 
                                            temporal.wt = annually_temporal_summer.wt_Geum,
                                            pca.wt = c("unweighted", "spatial", 
                                                       "temporal", "spatiotemporal"),
                                            scree = FALSE)

# visualize for flow-directed PCA
# score map
plot_stpca(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "map",
           river = shape_Geum, coords = coords_Geum)

plot_stpca(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "map",
           river = shape_Geum, coords = coords_Geum, pc.map=2)

plot_stpca(x = Y_avg_TOC_log_summer.Tmode.pca.uw, plots = "map",
           river = shape_Geum, coords = coords_Geum, pc.map=3)

# biplot
par(mfrow=c(1,1),mar=c(5,5,5,5))

# change sign of loadings
Y_avg_TOC_log_summer.Tmode.pca.all$spatial$loads[,2] <- -Y_avg_TOC_log_summer.Tmode.pca.all$spatial$loads[,2]
Y_avg_TOC_log_summer.Tmode.pca.all$spatial$scores[,2] <- -Y_avg_TOC_log_summer.Tmode.pca.all$spatial$scores[,2]

plot_river(Y_avg_TOC_log_summer.Tmode.pca.all, river = shape_Geum, plots="biplot", biplot.factor=30)

# glyph plot
par(mfrow=c(1,1),mar=c(4.5,4.5,4.5,4.5))
plot_stpca(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "glyph", 
           river = shape_Geum, coords = coords_Geum)


plot_stpca_zoom(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "glyph", 
                river = shape_Geum, coords = coords_Geum, lwd=5, 
                zoom=TRUE, zoomx=c(127.21,127.31), zoomy=c(36.21,36.5))



# ====================== #
#      PCA on river      #
# ====================== #
# spatial and temporal weights
sp.wt <- expm::sqrtm(solve(spatial.wt_Geum))

t.wt <- solve(expm::sqrtm(t(annually_temporal_summer.wt_Geum)))

# make spatial points dataframe
dls <- SpatialPointsDataFrame(coords_Geum, as.data.frame(t(Y_avg_TOC_log_summer.Tmode.completedata)))
colnames(dls@data) <- make.names(colnames(t(Y_avg_TOC_log_summer.Tmode.completedata)))
rownames(dls@data) <- rownames(spatial.wt_Geum)
n.var <- length(colnames(dls@data))

# evaluate the optimal bandwidth
bw <- bw.pca.river(dls, vars = colnames(dls@data), k=6, adaptive=TRUE, dMat = dmat, rwt = sp.wt, cwt = t.wt)

# implement the proposed PCA on river network  
pca.river.basic <- pca.river(dls, vars = colnames(dls@data),adaptive=TRUE,bw=bw,k=6,dMat=dmat, rwt = sp.wt, cwt = t.wt)

prop.var <- function(pca.obj, n.components) { # calculate var proportion
  if(n.components == 1){
    return((pca.obj$var[, 1:n.components] / 
              rowSums(pca.obj$var)) * 100)
  }else{
    return((rowSums(pca.obj$var[, 1:n.components]) / 
              rowSums(pca.obj$var)) * 100)
  }
}

var.pca.river.basic <- prop.var(pca.river.basic, 3) # 127 location, first 3 PCs'var proportion

dls$var.pca.river.basic <- var.pca.river.basic

loadings.pc1.pca.river.basic <- pca.river.basic$loadings[,,1] #127*8 matrix, 1st PC loadings on each location
win.item1.pca.river.basic <- max.col(abs(loadings.pc1.pca.river.basic)) # 127-vector, winning variable for pc1 on each location

dls$win.item1.pca.river.basic <- win.item1.pca.river.basic 


# summary of local PTV (proportion of total variance)
summary(pca.river.basic$local.PV[,1])
summary(pca.river.basic$local.PV[,2])
summary(pca.river.basic$local.PV[,3])
summary(pca.river.basic$local.PV[,1] + pca.river.basic$local.PV[,2])
summary(pca.river.basic$local.PV[,1] + pca.river.basic$local.PV[,2] + pca.river.basic$local.PV[,3])


loadings.pc2.pca.river.basic <- pca.river.basic$loadings[,,2]
loadings.pc12.pca.river.basic <- cbind(loadings.pc1.pca.river.basic, loadings.pc2.pca.river.basic)

# ward clustering result 
par(mfrow=c(2,2))

hc.ward <- hclust(dist(loadings.pc12.pca.river.basic), method='ward.D2')
plot(hc.ward, main=NULL, labels=FALSE, hang = -1, xlab="", ylab="", sub="")

plt <- c("#FF9933", "#3366FF", "#009966", "#FF6666", "#000000")
for(i in c(2,4,5)){
  n.cls <- i
  sub_grp <- cutree(hc.ward, k=n.cls) 
  scores <- sub_grp
  nums <- length(unique(scores))
  n.colors <- nums
  data.vec <- sort(scores)
  quants <- nums
  getPalette <- plt[1:n.colors]
  plot(coords_Geum[, 1], coords_Geum[, 2], type = "n", xlab = "Easting", 
       ylab = "Northing", main = paste("k =",i), cex.main = 1.5, 
       cex.axis = 1.5, cex.lab = 1.5, xlim = c(shape_Geum@bbox[1, 
       ]), ylim = c(shape_Geum@bbox[2, ]))
  lines(shape_Geum)
  points(coords_Geum[, 1], coords_Geum[, 2], pch = 19, cex = 1.2, 
         col = getPalette[cut(scores, breaks = quants)])
  
}

# k = 2 case 
n.cls <- 2
sub_grp <- cutree(hc.ward, k=n.cls) 
dls$cluster.pca.river.pc12 <- sub_grp


Y_avg_TOC_log_summer.Tmode.pca.all$PTV <- matrix(var.pca.river.basic, ncol=1) 
Y_avg_TOC_log_summer.Tmode.pca.all$group <- sub_grp
Y_avg_TOC_log_summer.Tmode.pca.all$winvar <- win.item1.pca.river.basic

Y_avg_TOC_log_summer.Tmode.pca.all$info <- c("PTV","group","winvar") 
Y_avg_TOC_log_summer.Tmode.pca.all$riverpca$loadings <- pca.river.basic$loadings
Y_avg_TOC_log_summer.Tmode.pca.all$riverpca$scores <- matrix(unlist(pca.river.basic$pca.river.scores), ncol=6, byrow=T)


# show local loading glyphs for PC 1 and 2, score glyphs for PC 1,2, and 3
plot_river(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "riverpca", river.glyph.loading=c(1,2),
           river.glyph.score = c(1,2,3), river = shape_Geum, coords = coords_Geum, r1=30, lwd=2.5)

# plot_river(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "riverpca", river.glyph=c(1,2),
#            river = shape_Geum, coords = coords_Geum, r1=30, lwd=2, zoom=TRUE, zoomx=c(127.21,127.31), zoomy=c(36.21,36.5))

plot_river(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "score", pc.map=1,
           river = shape_Geum, coords = coords_Geum)


plot_river(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "info",
           river = shape_Geum, coords = coords_Geum)


# ====================== #
#          CPCA          #
# ====================== #

# Flury's Hierarchy Model
S <- array(NA, c(n.var, n.var, n.cls), 
           dimnames = list(colnames(dls@data[,1:n.var]), colnames(dls@data[,1:n.var]), 
                           paste("Group",1:n.cls, sep = ""))) # sample covariance for each group 
nvec <- vector(length=n.cls) # vector of sample size

adjusted_data <- sp.wt %*% as.matrix(scale(dls@data[,1:n.var], center=TRUE)) %*% t.wt
rownames(adjusted_data) <- rownames(dls@data[,1:n.var])
colnames(adjusted_data) <- colnames(dls@data[,1:n.var])

for(i in 1:n.cls){
  target <- adjusted_data[dls@data$cluster.pca.river.pc12 == i,]
  nvec[i] <- nrow(target)
  S[,,i] <- cov(target)
}


# results with FG algorithm
pcpca.result.FG <- flury.hierarchy(data = as.matrix(dls@data[,1:n.var]), covmats = S, nvec = nvec, n.var = n.var, 
                                   cluster = dls$cluster.pca.river.pc12, alpha = 0.05, mode = "FG")


q <- 2

B <- array(NA, c(n.var, n.var, n.cls), 
          dimnames = list(colnames(dls@data[,1:n.var]), paste("PC",1:n.var, sep = ""), 
                          paste("Group",1:n.cls, sep = ""))) 

score <- array(NA, c(nrow(spatial.wt_Geum), n.var), 
              dimnames = list(colnames(spatial.wt_Geum), paste("PC",1:n.var, sep = "")))

for(i in 1:n.cls){
  B[,,i] <- t(solve(t.wt))%*%B.partial(covmats=S, nvec=nvec, B=pcpca.result.FG$B.cpc,
                                    commonvec.order = pcpca.result.FG$common_order, q=q)[,,i]
  score[dls$cluster.pca.river.pc12==i,] <- as.matrix(scale(dls@data[,1:n.var], 
                                                     center=TRUE)[dls$cluster.pca.river.pc12==i,])%*%t.wt%*%
    B.partial(covmats=S, nvec=nvec, B=pcpca.result.FG$B.cpc, commonvec.order = pcpca.result.FG$common_order, q=q)[,,i]
}



# apply CPC(2) model 
res <- cpc::cpcq.test(covmats = S, nvec = nvec, 
                     B = pcpca.result.FG$B.cpc[,pcpca.result.FG$common_order], q = q)

# eigenvalue & proportion of total variance 
evptv <- array(0, c(2, n.var, n.cls), 
              dimnames = list(c("Eigenvalues", "Percentage of total variation"),
                              paste("PC",c(1:n.var),sep=""),paste("Group",1:n.cls, sep = "")))

for(j in 1:n.cls){
  Fj <- t(B.partial(covmats=S, nvec=nvec, B=pcpca.result.FG$B.cpc,
                   commonvec.order = pcpca.result.FG$common_order, q=q)[,,j])%*%
    S[,,j]%*%B.partial(covmats=S, nvec=nvec, B=pcpca.result.FG$B.cpc,
                        commonvec.order = pcpca.result.FG$common_order, q=q)[,,j] # covariance of score 
  lambda.inv.sq.rt <- diag(1/sqrt(diag(Fj)))
  # Rj = lambda.inv.sq.rt%*%Fj%*%lambda.inv.sq.rt # correlation of score
  
  evptv[,,j][1,] <- diag(Fj)
  evptv[,,j][2,] <- diag(Fj)/sum(diag(Fj))*100
}

evptv

B[,2,] = - B[,2,]

B # loadings for each group including CPC

score[,2] = -score[,2]
Y_avg_TOC_log_summer.Tmode.pca.all$riverpca$scores <- score

# take a look at zoomed section 
plot_river(x = Y_avg_TOC_log_summer.Tmode.pca.all, plots = "riverpca", river.glyph.loading=c(1,2),
           river.glyph.score = c(1,2,3), river = shape_Geum, coords = coords_Geum, r1=10, lwd=6,
           zoom=TRUE, zoomx=c(127.21,127.31), zoomy=c(36.21,36.5))
