library(stringr)
library(SSN)
library(rgdal) # for readOGR function
library(riverdist) # for computing upstream distance
library(rivernet)
library(xts)
library(shp2graph) # for get_binaryIDs_stream function
library(spam)
library(scales)
library(stpca)
library(SSN)
library(rgdal)
library(GWmodel)

# from "get_binaryIDs_stream" (line 19) to "createWeightS_from_ssn" (line 646) was defined by Seoncheol Park 
# (https://github.com/SeoncheolPark/paper-StreamflowLifting)

# need RCH_ID in shape file
get_binaryIDs_stream <- function(mouth_node, shape, RCH_ID=NULL){
  if(is.null(RCH_ID)){
    RCH_ID <- shape$RCH_ID
  }
  
  rtNEL1 <-readshpnw(shape, ELComputed=TRUE, longlat=TRUE)
  igr1 <-nel2igraph(rtNEL1[[2]], rtNEL1[[3]])
  
  edgelist <- rtNEL1[[3]]
  mouth_node_index <- edgelist[mouth_node,c(2,3)]
  mouth_node_final <- c()
  if(sum(edgelist[,c(2,3)]== mouth_node_index[1])==1){
    mouth_node_final <- c(mouth_node_final,mouth_node_index[1])
    key_ind <- mouth_node_index[2]
  }else if(sum(edgelist[,c(2,3)]== mouth_node_index[2])==1){
    mouth_node_final <- c(mouth_node_final,mouth_node_index[2])
    key_ind <- mouth_node_index[1]
  }
  
  
  # result_binaryIDs <- data.frame(rid=c(1:length(shape$RCH_ID)), binaryID=rep("", length(shape$RCH_ID)))
  result_binaryIDs <- data.frame(rid=c(1:length(RCH_ID)), binaryID=rep("", length(RCH_ID)))
  result_binaryIDs$binaryID <- as.character(result_binaryIDs$binaryID)
  idx <- 1
  result_binaryIDs[mouth_node,2] <- as.character("1")
  key_RID_ind <- mouth_node
  pointsin <- c(1:nrow(result_binaryIDs))
  pointsin <- setdiff(pointsin, key_RID_ind)
  while(length(key_RID_ind)!=0){
    RID_candidate <- which(edgelist[,3]==key_ind[1])
    if(length(RID_candidate)==0){
      pointsin <- setdiff(pointsin, key_RID_ind[1])
    }else{
      if(result_binaryIDs[RID_candidate[1],2]==""){
        result_binaryIDs[RID_candidate[1],2] <- as.character(paste(result_binaryIDs[key_RID_ind[1],2],"0", sep=""))
        key_RID_ind <- c(key_RID_ind, RID_candidate[1])
        pointsin <- setdiff(pointsin, RID_candidate[1])
        if(length(RID_candidate)==2 & result_binaryIDs[RID_candidate[2],2]==""){
          result_binaryIDs[RID_candidate[2],2] <- as.character(paste(result_binaryIDs[key_RID_ind[1],2],"1", sep=""))
          key_RID_ind <- c(key_RID_ind, RID_candidate[2])
          pointsin <- setdiff(pointsin, RID_candidate[2])
        }
      }else if(length(RID_candidate)==2 & result_binaryIDs[RID_candidate[2],2]==""){
        result_binaryIDs[RID_candidate[2],2] <- as.character(paste(result_binaryIDs[key_RID_ind[1],2],"0", sep=""))
        key_RID_ind <- c(key_RID_ind, RID_candidate[2])
        pointsin <- setdiff(pointsin, RID_candidate[2])
      }
    }
    if(length(key_RID_ind)!=0){
      key_RID_ind <- key_RID_ind[-1]
      if(length(key_RID_ind)!=0){
        key_ind <- edgelist[key_RID_ind[1],2]
      }
    }
    # print(key_RID_ind)
  }
  return(result_binaryIDs)
}


# obtain binaryIDs_obj from get_binaryIDs_stream
get_adjacency_stream <- function (binaryIDs_obj, netID = 1) 
{
  binaryIDs <- binaryIDs_obj
  # binaryIDs <- get_binaryIDs(ssn_directory, net = netID)
  bid <- binaryIDs[, 2]
  rid <- binaryIDs[, 1]
  nch <- nchar(bid)
  bid.list <- split(bid, nch)
  rid.list <- split(rid, nch)
  n.segments <- nrow(binaryIDs)
  xy <- matrix(ncol = 2, nrow = (n.segments - 1))
  counter <- 0
  for (j in 1:(length(bid.list) - 1)) {
    bid.dn <- bid.list[[j]]
    bid.up <- bid.list[[j + 1]]
    rid.dn <- rid.list[[j]]
    rid.up <- rid.list[[j + 1]]
    bid.up.sub.vec <- bid.up
    rid.up.sub.vec <- rid.up
    for (i in 1:length(bid.dn)) {
      current.dn.bid <- bid.dn[i]
      current.dn.rid <- rid.dn[i]
      inner.count <- 1
      number.upstream <- 0
      n.bid.up <- length(bid.up.sub.vec)
      crit <- FALSE
      while (!crit) {
        if (n.bid.up > 0) {
          current.up.bid <- bid.up.sub.vec[inner.count]
          connected <- substr(current.up.bid, 1, nchar(current.dn.bid)) == 
            current.dn.bid
          if (connected) {
            counter <- counter + 1
            number.upstream <- number.upstream + 1
            xy[counter, ] <- c(current.dn.rid, rid.up.sub.vec[inner.count])
            rid.up.sub.vec <- rid.up.sub.vec[-inner.count]
            bid.up.sub.vec <- bid.up.sub.vec[-inner.count]
          }
          if (!connected) 
            inner.count <- inner.count + 1
          crit <- (number.upstream == 2) | ((number.upstream +  inner.count - 1) == n.bid.up)
        }
        if (n.bid.up == 0) 
          crit <- TRUE
      }
    }
  }
  xy[, 1] <- re_map_rid(rid_vector = xy[, 1], all_rid = rid)
  xy[, 2] <- re_map_rid(rid_vector = xy[, 2], all_rid = rid)
  add.one <- min(xy) == 0
  list(adjacency = spam(list(j = (xy[, 1] + add.one), i = (xy[, 2] + add.one), rep(1, nrow(xy))), nrow = n.segments, 
                        ncol = n.segments), rid_bid = cbind(rid, bid))
}


# compute shreve stream order
compute_shreve <- function(adjacency){
  # adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    # find next stream
    # case 1: Y자
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
    }
    # case 2: 1자
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    # print(pointsin)
  }
  
  return(shreve)
}



importSSN_stream_past <- function(shreve_obj, location="Miho", multipleplaces=FALSE, RID="original"){
  # location: "Miho", "Full" 
  # multipleplaces=FALSE -> fix and sum
  # multipleplaces=TRUE -> go ahead
  
  # netID = 1 (b/c there is one network)
  
  # (a) stream network data
  edges  <- readOGR("/home/kyu9510/pca_on_river/Data/KRF_3.0_Geumgang/KRF_ver3_LINE_금강수계.shp")
  # rownames(edges@data) <- edges@data[, "RCH_ID"]
  # reach ID (print rid)
  
  # (b) monitoring sites data
  sites <- read.csv("/home/kyu9510/pca_on_river/Data/ProcessedData/금강장소(엑셀편집)UTF8.csv",header=T)
  sites <- sites[which(sites$총량==1),]
  sites <- data.frame(name= sites$이름,x=sites$경도.Degree., y=sites$위도.Degree.)
  
  # rownames(sites@data) <- sites@data[,"pid"]
  # rownames(sites@coords) <- sites@data[, "pid"]
  
  # (c) water quality measurements data
  alldata1 <- read.table("/home/kyu9510/pca_on_river/Data/ProcessedData/수질20082013.txt",sep=',', header=T)
  alldata2 <- read.table("/home/kyu9510/pca_on_river/Data/ProcessedData/수질20132018.txt",sep=',', header=T)
  alldata <- rbind(alldata1, alldata2)
  
  # (d) preprocessed data 
  if(location=="Miho"){
    data <- readRDS("/home/kyu9510/pca_on_river/Data/ProcessedData/Geum(miho).RDS")
  }else if(location=="Full"){
    data <- readRDS("/home/kyu9510/pca_on_river/Data/ProcessedData/Geum.RDS")
  }else{
    stop("Wrong location")
  }
  shape <- data$shape
  
  data$streamdata <- data$streamdata["2012/2017",]
  # mid catchments
  region_index <- unique(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)))[c(1)]
  # the first 4 digits of ID = mid catchment code
  # 3011: 미호천
  # 3012: 금강공주
  # 3008: 대청댐
  # 3007: 보청천
  # 3005: 초강
  # 3009: 갑천
  # 3013: 논산천
  # 3006: 대청댐상류
  # 3014: 금강하구언
  # 3004: 영동천
  # 3003: 무주남대천
  # 3001: 용담댐
  # 3002: 용담댐하류
  # 3310: 대청댐하류
  # 3301: 만경강(x)
  # 3303: 직소천(x)
  # 3302: 동진강(x)
  # 3101: 삽교천(x)
  # 2005: 금강상류부분(필요없음)
  # 3202: 부남방조제(x)
  # 1101, 3201: 대호방조제(x)
  # 3203: 금강서해(x)
  # 5301, 5302: 주진천(x)
  
  # for(i in 1:length(region_index)){
  #   shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  #   plot(shape, main=region_index[i])
  #   plot(shape_sub, col="red", add=T, lwd=3)
  # }
  
  name_index <- c("미호천", "금강공주", "대청댐", "보청천", "초강", "갑천", "논산천", "대청댐상류", "금강하구언", "영동천", "무주남대천", "용담댐", "용담댐하류", "대청댐하류")
  name_subindex <- list()
  name_subindex[[1]] <- c("미호댐", "원남댐", "병천천상류", "병천천하류", "보강천", "미호천중류", "미호천상류", "작천보", "한천", "백곡천", "조천", "석화수위표", "미호천하류", "무심천", "백곡댐" ) #15개
  name_subindex[[2]] <- c("세종보", "지천상류", "지천합류후", "어천합류후", "용수천", "백제보", "석성천", "규암수위표", "금천", "논산천합류전", "유구천", "공주보", "대교천", "공주수위표") #14개
  name_subindex[[3]] <- c("대청댐", "대청댐조정지", "소옥천상류", "소옥천하류", "대청댐상류") # 5개
  name_subindex[[4]] <- c("항건천", "보청천중류", "삼가천합류후", "삼가천", "보청천상류", "보청천하류") # 6개
  name_subindex[[5]] <- c("석천", "초강상류", "초강하류") # 3개
  name_subindex[[6]] <- c("갑천하류", "유성수위표", "갑천상류", "대전천", "유등천상류", "유등천하류") # 6개
  name_subindex[[7]] <- c("노성천", "논산천하류", "강경천", "논산천상류", "탑정댐") # 5개
  name_subindex[[8]] <- c("보청천합류전") # 1개
  name_subindex[[9]] <- c("입포수위표", "길산천", "금강하구언") # 3개
  name_subindex[[10]] <- c("영동천", "봉황천하류", "초강합류후", "봉황천상류", "봉황천합류후", "호탄수위표") # 6개
  name_subindex[[11]] <- c("무주남대천상류", "무주남대천중류", "무주남대천하류") # 3개
  name_subindex[[12]] <- c("주자천", "용담댐", "정자천", "구량천", "진안천", "장계천합류후", "장계천", "진안천합류후") # 8개
  name_subindex[[13]] <- c("용담댐하류") # 1개
  name_subindex[[14]] <- c("미호천합류전", "매포수위표") # 2개
  region_subindex <- list()
  
  Geum_RCH_ID <- c()
  for(i in 1:length(region_index)){
    shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
    # plot(shape_sub)
    region_subindex[[i]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6)))
    for(j in 1:length(region_subindex[[i]])){
      shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
      # plot(shape_sub_sub, col="red", add=T)
      Geum_RCH_ID_Imsi <- cbind(shape_sub_sub@data, data.frame(중권역번호=region_index[i], 중권역이름=name_index[i], 소권역번호=region_subindex[[i]][j], 소권역이름=name_subindex[[i]][j]))
      Geum_RCH_ID <- rbind(Geum_RCH_ID, Geum_RCH_ID_Imsi)
    }
  }
  
  # when print subregion
  # i=8
  # shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  # plot(shape_sub)
  # region_subindex[[i]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6)))
  # region_subsubindex[[i]] <- list()
  #
  # j=1
  # shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
  # plot(shape_sub_sub, col="red", add=T)
  # region_subsubindex[[i]][[j]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=8)))
  #
  # for(k in 1:length(region_subsubindex[[i]][[j]])){
  #   shape_sub_sub_sub <- subset(shape_sub_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=8))==region_subsubindex[[i]][[j]][k])
  #   plot(shape_sub_sub_sub, col=k, add=T)
  # }
  
  ########################################
  ##Adjacency matrix
  ########################################
  # reference: https://cran.r-project.org/web/packages/riverdist/vignettes/riverdist_vignette.html
  
  # connections: Object of class "matrix", with "numeric" elements. Defined as a square matrix, with elements describing the type of connection detected between line segments.
  # • A value of 1 in element [i,j] indicates that the beginning of segment i is connected to the beginning of segment j.
  # • A value of 2 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j.
  # • A value of 3 in element [i,j] indicates that the end of segment i is connected to the beginning of segment j.
  # • A value of 4 in element [i,j] indicates that the end of segment i is connected to the end of segment j.
  # • A value of 5 in element [i,j] indicates that segments i and j are connected at both beginning and end.
  # • A value of 6 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j, and the end of segment i is connected to the beginning of segment j.
  # • A value of NA in element [i,j] indicates that segments i and j are not connected.
  if(location=="Miho"){
    flow_network_original <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="GeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    flow_network <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="thinGeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    # plot(flow_network)
  }else if(location=="Full"){
    flow_network_original <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="Geum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    flow_network <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="thinGeum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    # plot(flow_network)
  }else{
    stop("Wrong location")
  }
  
  adj_matrix <- flow_network$connections # adjacency of network
  
  # topologydots(rivers=flow_network) # topologydots: if connected -> green, if not -> red
  
  # Basic distance calculation in a non-messy river network
  # detectroute(start=4, end=2, rivers=flow_network)
  # riverdistance(startseg=4, startvert=1, endseg=2, endvert=1, rivers=flow_network, map=TRUE)
  
  # exmaple of subnetwork output 
  # i=1
  # flow_subnetwork <- trimriver(trimto = c(which(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=6)==region_subindex[[1]][i]) )), rivers=flow_network)
  # plot(flow_subnetwork)
  
  #하구 설정
  if(location=="Miho"){
    flow_network <- setmouth(seg=68,vert=2, rivers=flow_network)  # node 68 : river mouth 
    flow_network_original <- setmouth(seg=68,vert=20, rivers=flow_network_original) # node 68 : river mouth 
  }else if(location=="Full"){
    flow_network <- setmouth(seg=614,vert=2, rivers=flow_network) # node 68 : river mouth 
    flow_network_original <- setmouth(seg=614,vert=8, rivers=flow_network_original) # node 68 : river mouth 
  }else{
    stop("Invalid location")
  }
  
  #showends(seg=80, rivers=flow_network)
  
  # added stream lengths by myself 
  # reach_length <- shape@data$Shape_Leng[which(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==rregion_subindex[[1]][i]) )]
  
  # ind1 <- colnames(edges@data) == c("netID")
  # ind2 <- colnames(edges@data) == c("rid")
  # ind3 <- colnames(edges@data) == c("upDist")
  # if (is.factor(edges@data$netID)) {
  #   edges@data$netID <- as.character(edges@data$netID)
  # }
  
  if(location=="Miho"){
    
    if(multipleplaces==FALSE){
      RCH_ID_manual <- c(38, 30, 40, 44, 45, 63, 50, 46, 1, 88, 22, 87, 27, 84, 62,
                         105, 109, 83, 54, 81, 18, 16, 19, 80, 78, 90, 77, 56, 68)
    }else{
      RCH_ID_manual <- c(38, 30, 40, 44, 45, 63, 50, 46, 1, 88, 22, 87, 27, 84, 62,
                         105, 109, 83, 54, 81, 18, 16, 19, 81, 78, 90, 90, 56, 68)
    }
    
    # data$eventplace$RCH_ID <- shape@data$RCH_ID[c(38, 30, 40, 44, 45, 63, 50, 46, 1, 88, 22, 87, 27, 84, 62,
    #                                              105, 109, 83, 54, 81, 18, 16, 19, 81, 78, 90, 90, 56, 68)]
    
    data$eventplace$RCH_ID <- shape@data$RCH_ID[RCH_ID_manual]
    # data$eventplace$RCH_ID : 
    
    # match upstream of the nearest lines for each location
  }else if(location=="Full" & multipleplaces==FALSE){
    
    # which(data$eventplace$X=="미호천4" | data$eventplace$X=="유등천A" | data$eventplace$X=="삼가천" | data$eventplace$X=="석천")
    # included so that location index different to the written index 
    # 17: 추풍령천, 18: 초강1
    data$eventplace$RCH_ID[17] <- as.integer(30050108)
    # 23: 옥천, 24: 우산 (옥천으로 평균)
    # 36: 대청댐, 38: 대청
    data$eventplace$RCH_ID[36] <- as.integer(30080409); data$eventplace$RCH_ID[38] <- as.integer(30080419)
    # 43: 갑천1, 44: 갑천2
    data$eventplace$RCH_ID[44] <- as.integer(30090208)
    # 47: 침산교, 48: 유등천1
    data$eventplace$RCH_ID[47] <- as.integer(30090502)
    # 49: 옥계교, 50: 대전천1 (두개 평균을내어야 (대전천1로 평균))
    # 55: 갑천5, 56: 갑천5-1 (갑천5로 평균)
    # 58: 청원-1 59: 백천
    data$eventplace$RCH_ID[58] <- as.integer(30100203)
    # 79: 미호천8, 83: 미호천5
    data$eventplace$RCH_ID[79] <- as.integer(30111301)
    # 85: 조천, 86: 조천1
    data$eventplace$RCH_ID[86] <- as.integer(30111505)
    # 98: 곰나루, 99: 금강
    data$eventplace$RCH_ID[98] <- as.integer(30120532)
    
    
    # [1] 구량천       가막         진안천       정자천       용담         용포         무주남대천   무주남대천-1 무주남대천-2 부리        
    # [11] 제원         봉황천       제원A        영동         영동천1      영동천2      추풍령천     초강1        금계천       초강2       
    # [21] 이원         옥천         우산         보청천1      보청천2      항건천       보청천3      보청천4      상곡천       추풍천      
    # [31] 옥천천       회인천       주원천       대청댐       품곡천       대청         용호천       현도         두계천1      봉곡2교     
    # [41] 갑천1        갑천2        갑천3        침산교       유등천1      옥계교       대전천1      대전천2      대전천3      유등천5     
    # [51] 갑천4        갑천5        갑천5-1      세종1        청원-1       백천         미호천1      칠장천       미호천1-1    한천        
    # [61] 미호천2      백곡천1      백곡천2      미호천7      초평천       미호천3      보강천1      보강천       성암천       무심천      
    # [71] 무심천1      무심천2      무심천3      석남천       미호천8      병천천A      용두천       병천천       미호천5      미호천5A    
    # [81] 조천         조천1        월하천       미호천6-1    연기         금남         용수천-1     용수천       대교천       대교천2     
    # [91] 세종2        공주1        정안천       곰나루       금강         유구천       목면         공주2        부여         지천        
    # [101] 지천-1       정동         은산천       부여1        금천         부여2        석성천-1     석성천       석성천2      성동        
    # [111] 논산천-1     논산천1      연산천       노성천-1     노성천       논산천2      방축천       논산천4      수철천       마산천      
    # [121] 강경천       강경         산북천       양화-1       길산천       길산천2      금강갑문    
    # }else{
    #   
    #}
    
  }
  dist.pts <- rep(0, nrow(data$eventplace))
  dist.imsi <- c()
  for(jj in 1:nrow(data$eventplace)){
    pts.coords <- c(data$eventplace$경도.Degree.[jj], data$eventplace$위도.Degree.[jj])
    
    kk <- match(data$eventplace$RCH_ID[jj], shape@data$RCH_ID) # fixme !!
    
    # dist_upper <- apply(shape@data[,c(3,4)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
    # dist_lower <- apply(shape@data[,c(5,6)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
    #
    # kk <- sort(unique(union(which(dist_upper==min(c(dist_upper, dist_lower))), which(dist_lower==min(c(dist_upper, dist_lower))))))
    #
    # dist.imsiimsi <- c()
    # for(pp in 1:length(kk)){
    #   for(ll in 1:nrow(flow_network_original$lines[[kk[pp]]])){
    #     dist.compute <- sqrt(sum((flow_network_original$lines[[kk[pp]]][ll,]-pts.coords)^2))
    #     dist.imsiimsivec <- matrix(c(dist.compute, kk[pp], ll), nrow=1)
    #     dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
    #   }
    # }
    # dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
    # dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
    
    dist.imsiimsi <- c()
    for(ll in 1:nrow(flow_network_original$lines[[kk]])){
      dist.compute <- sqrt(sum((flow_network_original$lines[[kk]][ll,]-pts.coords)^2))
      dist.imsiimsivec <- matrix(c(dist.compute, kk, ll), nrow=1)
      dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
    }
    dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
    dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
  }
  colnames(dist.imsi) <- c("mindist", "rid", "segid")
  rownames(dist.imsi) <- data$eventplace$X
  
  if(location=="Miho"){
    # upDist should also be calculated
    DistanceUpstream_seg <- rep(0, nrow(shape@data))
    for(ii in 1:length(DistanceUpstream_seg)){
      DistanceUpstream_seg[ii] <- upstream(startseg = 68, endseg = ii, startvert = 20, endvert = 1, rivers=flow_network_original, net=TRUE)
    }
    # DistanceUpstream of upstream: use dist.imsi
    DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
    for(ii in 1:length(DistanceUpstream_obs)){
      DistanceUpstream_obs[ii] <- upstream(startseg = 68, endseg = dist.imsi[ii,2], startvert = 20, endvert = dist.imsi[ii,3], rivers=flow_network_original, net=TRUE)
    }
    shape_tinned <- readOGR("/home/kyu9510/pca_on_river/Data/ProcessedData/thinGeumMiho.shp")
  }else if(location=="Full"){
    # upDist should also be calculated
    DistanceUpstream_seg <- rep(0, nrow(shape@data))
    for(ii in 1:length(DistanceUpstream_seg)){
      DistanceUpstream_seg[ii] <- upstream(startseg = 614, endseg = ii, startvert = 8, endvert = 1, rivers=flow_network_original, net=TRUE)
    }
    # DistanceUpstream of upstream: use dist.imsi
    DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
    for(ii in 1:length(DistanceUpstream_obs)){
      DistanceUpstream_obs[ii] <- upstream(startseg = 614, endseg = dist.imsi[ii,2], startvert = 8, endvert = dist.imsi[ii,3], rivers=flow_network_original, net=TRUE)
    }
    shape_tinned <- readOGR("/home/kyu9510/pca_on_river/Data/ProcessedData/thinGeum.shp")
  }
  
  # (1)network.line.coords: network.line.coords
  network.line.coords <- data.frame(NetworkID=as.factor(rep(1,nrow(shape@data))), SegmentID=c(0:(nrow(shape@data)-1)), DistanceUpstream=DistanceUpstream_seg)
  # observeddata <- SSNPoints
  if(location=="Miho"){
    miho <- which(data$eventplace$X=="미호천4")
  }else if(location=="Full"){
    # 삼가천, 석천: no obs.points
    # 유등천A, 미호천4: data$eventplace$X exist but colnames(data$normaldata_TN) not
    miho <- which(data$eventplace$X=="미호천4" | data$eventplace$X=="유등천A" | data$eventplace$X=="삼가천" | data$eventplace$X=="석천")
    if(multipleplaces==FALSE){
      miho <- sort(c(miho, which(data$eventplace$X=="무주남대천" | data$eventplace$X=="무주남대천-1" | data$eventplace$X=="우산" | data$eventplace$X=="옥계교" | data$eventplace$X=="갑천5-1")))
    }else{
      miho <- sort(miho)
    }
  }
  # (2)obspoint: obspoints
  network.line.coords.obsdata <- data.frame(NetworkID=rep(1,(nrow(data$eventplace)-length(miho))), SegmentID=dist.imsi[-miho,2], DistanceUpstream=DistanceUpstream_obs[-miho])
  pointcoords <- cbind(data$eventplace$경도.Degree.[-miho], data$eventplace$위도.Degree.[-miho]); colnames(pointcoords) <- c("coords.x1", "coords.x2")
  # pointdata에 ratio needs to be added
  if(multipleplaces==TRUE & location=="Full"){
    if(RID=="original"){
      RCH_ID_manual <- c(403, 445, 792, 418, 807, 812, 869, 869, 869, 821,
                         826, 827, 830, 833, 788, 789, 228, 228, 861, 867,
                         662, 669, 669, 539, 513, 14, 516, 877, 222, 199,
                         507, 577, 921, 153, 530, 750, 90, 718, 174, 175,
                         579, 579, 584, 569, 569, 182, 182, 856, 855, 591,
                         589, 590, 590, 739, 508, 508, 235, 142, 237, 241,
                         242, 509, 247, 243, 1, 612, 93, 611, 128, 503, 
                         769, 773, 607, 325, 605, 76, 74, 77, 605, 602,
                         754, 754, 448, 592, 738, 737, 62, 49, 481, 476,
                         700, 720, 366, 752, 752, 901, 687, 677, 674, 890,
                         888, 928, 163, 663, 853, 660, 162, 841, 836, 647,
                         775, 554, 190, 560, 555, 551, 252, 646, 785, 919,
                         576, 645, 270, 630, 884, 883, 614) 
    }else if(RID=="corrected"){
      RCH_ID_manual <- c(403, 445, 792, 418, 807, 813, 869, 815, 815, 821,
                         826, 294, 830, 833, 789, 789, 858, 228, 861, 867,
                         662, 669, 669, 539, 520, 14, 516, 877, 222, 506,
                         507, 577, 921, 153, 530, 736, 90, 718, 174, 216,
                         579, 581, 584, 569, 569, 182, 182, 856, 572, 573,
                         589, 590, 590, 739, 739, 508, 235, 142, 245, 241,
                         242, 509, 247, 243, 1, 612, 611, 611, 128, 503, 
                         769, 773, 607, 325, 605, 76, 74, 77, 604, 602,
                         754, 754, 448, 592, 737, 737, 56, 728, 482, 476,
                         700, 723, 366, 733, 752, 904, 687, 677, 753, 37,
                         888, 672, 163, 663, 853, 658, 162, 839, 837, 647,
                         775, 554, 190, 560, 555, 551, 252, 646, 296, 575,
                         576, 645, 270, 628, 883, 883, 614)
    }
    ## (738, 737) -> (737, 751)
    dist.imsi <- c()
    lll <- 0
    for(jj in 1:nrow(data$eventplace)){
      pts.coords <- c(data$eventplace$경도.Degree.[jj], data$eventplace$위도.Degree.[jj])
      
      lll <- lll+1
      if((jj%in%miho)){
        lll <-lll - 1
      }
      kk <- RCH_ID_manual[lll]
      # kk <- match(RCH_ID_manual[lll], shape@data$RCH_ID) # to be fixed
      # kk <- match(data$eventplace$RCH_ID[jj], shape@data$RCH_ID) # to be fixed
      
      # dist_upper <- apply(shape@data[,c(3,4)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
      # dist_lower <- apply(shape@data[,c(5,6)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
      #
      # kk <- sort(unique(union(which(dist_upper==min(c(dist_upper, dist_lower))), which(dist_lower==min(c(dist_upper, dist_lower))))))
      #
      # dist.imsiimsi <- c()
      # for(pp in 1:length(kk)){
      #   for(ll in 1:nrow(flow_network_original$lines[[kk[pp]]])){
      #     dist.compute <- sqrt(sum((flow_network_original$lines[[kk[pp]]][ll,]-pts.coords)^2))
      #     dist.imsiimsivec <- matrix(c(dist.compute, kk[pp], ll), nrow=1)
      #     dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
      #   }
      # }
      # dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
      # dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
      
      dist.imsiimsi <- c()
      for(ll in 1:nrow(flow_network_original$lines[[kk]])){
        dist.compute <- sqrt(sum((flow_network_original$lines[[kk]][ll,]-pts.coords)^2))
        if(lll==9){
          dist.compute <- sqrt(sum((flow_network_original$lines[[kk]][ll,]-pts.coords-c(0.01,0))^2))
        }
        dist.imsiimsivec <- matrix(c(dist.compute, kk, ll), nrow=1)
        dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
      }
      dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
      dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
    }
    colnames(dist.imsi) <- c("mindist", "rid", "segid")
    rownames(dist.imsi) <- data$eventplace$X
    
    if(RID=="original"){
      dist.imsi[9,3] <- 126
    }else if(RID=="corrected"){
      dist.imsi[9,3] <- 41
    }
    
    DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
    for(ii in 1:length(DistanceUpstream_obs)){
      DistanceUpstream_obs[ii] <- upstream(startseg = 614, endseg = dist.imsi[ii,2], startvert = 8, endvert = dist.imsi[ii,3], rivers=flow_network_original, net=TRUE)
    }
    
    # data$eventplace$RCH_ID <- shape@data$RCH_ID[RCH_ID_manual]
    pointdata <- cbind(data$eventplace[-miho,], rid_old=dist.imsi[-miho,2], upDist=DistanceUpstream_obs[-miho], locID=c(1:(nrow(data$eventplace)-length(miho))), netID=rep(1,(nrow(data$eventplace)-length(miho))), pid=setdiff(c(1:nrow(data$eventplace)), miho), shreve=shreve_obj[dist.imsi[-miho,2]], mindist=dist.imsi[-miho,1], rid=RCH_ID_manual, segid=dist.imsi[-miho,3])
  }else{
    pointdata <- cbind(data$eventplace[-miho,], rid_old=dist.imsi[-miho,2], upDist=DistanceUpstream_obs[-miho], locID=c(1:(nrow(data$eventplace)-length(miho))), netID=rep(1,(nrow(data$eventplace)-length(miho))), pid=setdiff(c(1:nrow(data$eventplace)), miho), shreve=shreve_obj[dist.imsi[-miho,2]], mindist=dist.imsi[-miho,1], rid=dist.imsi[-miho,2], segid=dist.imsi[-miho,3])
  }
  pointbbox <- rbind(range(pointcoords[-miho,1]), range(pointcoords[-miho,2])); colnames(pointbbox) <- c("min", "max"); rownames(pointbbox) <- c("coords.x1", "coords.x2")
  network.SSNPoints <- list()
  network.SSNPoints[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.obsdata, point.coords = pointcoords, point.data=pointdata, points.bbox=pointbbox, proj4string=shape@proj4string)
  obspoints <- new("SSNPoints", SSNPoints=network.SSNPoints, ID="Obs")
  # (3)predpoint
  network.line.coords.preddata <- data.frame(NetworkID=rep(1,length(miho)), SegmentID=dist.imsi[miho,2], DistanceUpstream=DistanceUpstream_obs[miho])
  rownames(network.line.coords.preddata) <- rownames(dist.imsi)[miho]
  pointcoords.preddata <- cbind(data$eventplace$경도.Degree.[miho], data$eventplace$위도.Degree.[miho]); colnames(pointcoords.preddata) <- c("coords.x1", "coords.x2")
  # pointdata에 ratio 추가해야
  pointdata.preddata <- cbind(data$eventplace[miho,], rid_old=dist.imsi[miho,2], upDist=DistanceUpstream_obs[miho], locID=c(1:length(miho)), netID=rep(1,length(miho)), pid=miho, shreve=shreve_obj[dist.imsi[miho,2]], mindist=dist.imsi[miho,1], rid=dist.imsi[miho,2], segid=dist.imsi[miho,3])
  pointbbox.preddata <- rbind(range(pointcoords[miho,1]), range(pointcoords[miho,2])); colnames(pointbbox.preddata) <- c("min", "max"); rownames(pointbbox.preddata) <- c("coords.x1", "coords.x2")
  network.SSNPoints.preddata <- list()
  network.SSNPoints.preddata[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.preddata, point.coords = pointcoords.preddata, point.data=pointdata.preddata, points.bbox=pointbbox.preddata, proj4string=shape@proj4string)
  if(length(miho)==0){
    # (3-1) no predpt case
    predpoints <- new("SSNPoints", SSNPoints=list(), ID=character(0))
  }else{
    predpoints <- new("SSNPoints", SSNPoints=network.SSNPoints.preddata, ID="Preds")
  }
  # (4)SSNpath: path
  SSNpath <- "NULL"
  # (5)data: data.frame type ( rid,upDist,    Length, netID)
  # data:Object of class "data.frame". 
  # The number of rows in data should equal the number of lines in the lines object. Row names correspond to SegmentID values
  shapedata <- shape@data
  shapedata <- cbind(shapedata, rid=c(1:nrow(shape@data)), upDist=DistanceUpstream_seg, Length=shape@data$Shape_Leng, netID=rep(1,nrow(shape@data)), shreve=shreve_obj )
  # (6)lines: segment polygon
  shapelines <- shape@lines
  # (7)bbox
  shapebbox <- shape@bbox; colnames(shapebbox) = c("min", "max"); rownames(shapebbox) <- c("x", "y")
  # (8)roj4string
  shapeproj4string <- shape@proj4string
  
  Geum_original <- new("SpatialStreamNetwork", network.line.coords=network.line.coords, obspoints=obspoints, predpoints=predpoints, path=SSNpath, data=shapedata, lines=shapelines, bbox=shapebbox, proj4string=shapeproj4string)
  
  # Object of class "SSNPoints" with 2 slots
  #
  # @ SSNPoints: List of SSNPoint objects with 5 slots
  # @ network.point.coords: object of class "data.frame". Row names
  # represent point identifiers (pid) stored in the point.data
  # data.frame.
  # $ NetworkID: factor identifying the NetworkID of that point
  # $ SegmentID: factor identifying the unique stream segment of that point
  # $ DistanceUpstream: numeric value representing the cumulative
  # distance from the network outlet, the most downstream point
  # on a network, to that point
  # @ point.coords: numeric matrix or "data.frame" with x- and y-
  #   coordinates (each row is a point); row names represent point
  # identifiers (pid) stored in the point.data data.frame.
  # @ point.data: object of class "data.frame"; the number of rows in
  # data should equal the number of points in the
  # network.point.coords object; row names are set to the pid
  # attribute.
  # @ points.bbox: Object of class "matrix"; see Spatial-class
  # @ proj4string: Object of class "CRS"; see CRS-class
  # @ ID: character string representing the name of the observation points
  
  
  # network.line <- shape_tinned@lines
  # bbox, proj4string
  
  return(Geum_original)
}

# existing createWeightS function takes file path of ssn object as an input 
# createWeightS_from_ssn takes ssn object itself as an input
createWeightS_from_ssn <- function (ssndata, afvcol="shreve", adjacency) {
  # ssndata : should be SpatialStreamNetwork class
  # afvcol: character (which contains weight information)
  # adjacency: adjacency information (can be obtained from get_adjacency_stream)
  
  if(class(ssndata)[1]!="SpatialStreamNetwork"){
    stop("ssndata must be a SpatialStreamNetwork class")
  }else if( !(afvcol %in% names(ssndata@obspoints@SSNPoints[[1]]@point.data) )){
    stop("afvcol does not exist")
  }else if(!(all(c("upDist", "rid", "pid")%in%colnames(ssndata@obspoints@SSNPoints[[1]]@point.data)))){
    stop("you must check whether (upDist, rid, pid) are in your SpatialStreamNetwork object" )
  }
  
  n.data <- nrow(ssndata@obspoints@SSNPoints[[1]]@point.data)
  
  binary_rid <- adjacency$rid_bid[as.numeric(ssndata@obspoints@SSNPoints[[1]]@point.data$rid),2]
  nchar_vec <- nchar(binary_rid)
  
  
  distmat1 <- matrix(0, nrow = n.data, ncol = n.data)
  distmats <- list()
  pid.names <- list()
  
  # rule of distmat1
  # (downstream, upstream): 1
  # diagonal elt: 0
  # (in the same segment): (upstreamDist low, high): 1
  for(k in 1:n.data){
    rid_k <- ssndata@obspoints@SSNPoints[[1]]@point.data$rid[k]
    binary_rid_k <- binary_rid[k]
    
    
    upstream_candidate0 <- which(sapply(binary_rid, function(x) substr(x, start=1, stop=nchar(binary_rid_k)))==as.character(binary_rid_k))
    # (1) nchar_vec of candidate is long : set 1
    # which(nchar_vec[upstream_candidate0]>nchar_vec[k])
    distmat1[k,upstream_candidate0[which(nchar_vec[upstream_candidate0]>nchar_vec[k])]] <- 1
    # (2) equal -> look upstream dist
    # which(nchar_vec[upstream_candidate0]<=nchar_vec[k])
    upstream_candidate1 <- upstream_candidate0[which(nchar_vec[upstream_candidate0]<=nchar_vec[k])]
    
    upstream_candidate2 <- upstream_candidate1[which(ssndata@obspoints@SSNPoints[[1]]@point.data$upDist[upstream_candidate1]>ssndata@obspoints@SSNPoints[[1]]@point.data$upDist[k])]
    if(length(upstream_candidate2)!=0){
      distmat1[k,upstream_candidate2] <- 1
    }
  }
  distmats[[1]] <- distmat1
  pid.names[[1]] <- ssndata@obspoints@SSNPoints[[1]]@point.data$pid
  
  connectedness <- as.matrix(Matrix::bdiag(distmats))
  colnames(connectedness) <- unlist(pid.names)
  rownames(connectedness) <- unlist(pid.names)
  
  afvdata <- SSN::getSSNdata.frame(ssndata)
  afvdata <- afvdata[order(afvdata$pid), ]
  
  
  ordrow <- order(as.numeric(rownames(connectedness)))
  ordcol <- order(as.numeric(colnames(connectedness)))
  connectedness <- connectedness[ordrow, ordcol, drop = F]
  if (mean(afvdata$pid == colnames(connectedness)) != 1) {
    stop("data frame observations are not in same order as distance matrix")
  }
  x <- afvdata[, afvcol]
  y <- afvdata[, afvcol]
  weight <- matrix(0, nrow = length(x), ncol = length(x))
  for (i in 1:length(x)) {
    weight[i, ] <- sqrt(x/y[i])
  }
  weight.connected <- matrixcalc::hadamard.prod(connectedness, 
                                                weight)
  diag(weight.connected) <- rep(1, length(diag(weight.connected)))
  if (sum(range(weight.connected) %in% c(0, 1)) != 2) {
    stop("weights are not in range [0,1]")
  }
  
  # x2 <- SSN::importSSN(ssndata)
  x2 <- ssndata
  order.data <- SSN::getSSNdata.frame(x2)
  ordrow2 <- as.character(order.data$pid)
  ordcol2 <- as.character(order.data$pid)
  weight.connected <- weight.connected[ordrow2, ordcol2, drop = F]
  colnames(weight.connected) <- as.character(ssndata@obspoints@SSNPoints[[1]]@point.data['X'][[1]])
  rownames(weight.connected) <- as.character(ssndata@obspoints@SSNPoints[[1]]@point.data['X'][[1]])
  
  return(weight.connected)
}


## =============================================================================================== ##

# ======================================= #
#           functions by myself           #
# ======================================= #
# These functions are basically constructed by modifying existing functions of the "GWPCA" and "stpca" packages
# for the purpose of analyzing river network data. 


importSSN_stream <- function(data, flow_network_original, shreve_obj, mouth_node){
  # location: "Miho", "Full" 
  # multipleplaces=FALSE -> fix and sum
  # multipleplaces=TRUE -> go ahead
  
  # netID = 1 (b/c there is one network)
  
  shape <- data$shape
  
  
  ########################################
  ##Adjacency matrix
  ########################################
  # reference: https://cran.r-project.org/web/packages/riverdist/vignettes/riverdist_vignette.html
  
  # connections: Object of class "matrix", with "numeric" elements. Defined as a square matrix, with elements describing the type of connection detected between line segments.
  # • A value of 1 in element [i,j] indicates that the beginning of segment i is connected to the beginning of segment j.
  # • A value of 2 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j.
  # • A value of 3 in element [i,j] indicates that the end of segment i is connected to the beginning of segment j.
  # • A value of 4 in element [i,j] indicates that the end of segment i is connected to the end of segment j.
  # • A value of 5 in element [i,j] indicates that segments i and j are connected at both beginning and end.
  # • A value of 6 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j, and the end of segment i is connected to the beginning of segment j.
  # • A value of NA in element [i,j] indicates that segments i and j are not connected.
  
  # topologydots(rivers=flow_network) # topologydots: if connected -> green, if not -> red
  
  # Basic distance calculation in a non-messy river network
  # detectroute(start=4, end=2, rivers=flow_network)
  # riverdistance(startseg=4, startvert=1, endseg=2, endvert=1, rivers=flow_network, map=TRUE)
  
  # exmaple of subnetwork output 
  # i=1
  # flow_subnetwork <- trimriver(trimto = c(which(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=6)==region_subindex[[1]][i]) )), rivers=flow_network)
  # plot(flow_subnetwork)
  
  # set river mouth 
  flow_network_original <- setmouth(seg = mouth_node, vert=nrow(flow_network_original$lines[[mouth_node]]), rivers=flow_network_original)
  
  dist.pts <- rep(0, nrow(data$eventplace))
  dist.tmp <- c()
  for(i in 1:nrow(data$eventplace)){
    pts.coords <- c(data$eventplace$lon[i], data$eventplace$lat[i])
    
    k <- match(data$eventplace$RCH_ID[i], shape@data$RCH_ID) # fixme !! RCH_ID where the ith monitoring site exists
    
    dist.tmp2 <- c()
    for(j in 1:nrow(flow_network_original$lines[[k]])){ # flow_network_original$lines[[k]] : the coordinates consisting of the stream segment where the ith monitoring site exists 
      dist.compute <- sqrt(sum((flow_network_original$lines[[k]][j,] - pts.coords) ^ 2))
      dist.tmp2.vec <- matrix(c(dist.compute, k, j), nrow=1)
      dist.tmp2 <- rbind(dist.tmp2, dist.tmp2.vec)
    }
    dist.tmp2 <- dist.tmp2[which.min(dist.tmp2[,1]),] # choose the closest vertex 
    dist.tmp <- rbind(dist.tmp, dist.tmp2)
  }
  colnames(dist.tmp) <- c("mindist", "rid", "segid") 
  # mindist : minimum distance to the closest point among points on the segment where the observation exists
  # rid : river stream index where the observation exists
  # segid : closest vertex index of the kth river stream 
  rownames(dist.tmp) <- data$eventplace$X
  
  # upDist should also be calculated
  DistanceUpstream_seg <- rep(0, nrow(shape@data))
  for(i in 1:length(DistanceUpstream_seg)){
    DistanceUpstream_seg[i] <- upstream(startseg = mouth_node, endseg = i, startvert = nrow(flow_network_original$lines[[mouth_node]]), 
                                        endvert = 1, rivers = flow_network_original, net=TRUE) 
    # start from the river mouth, end to the most upstream point of the end segment
    # endvert = 1 means the first point of end segment (ith river stream)
  }
  # DistanceUpstream of upstream: use dist.tmp
  DistanceUpstream_obs <- rep(0, nrow(data$eventplace))
  for(i in 1:length(DistanceUpstream_obs)){
    DistanceUpstream_obs[i] <- upstream(startseg = mouth_node, endseg = dist.tmp[i,2], startvert = nrow(flow_network_original$lines[[mouth_node]]), 
                                        endvert = dist.tmp[i,3], rivers = flow_network_original, net = TRUE) 
    # measure distance until observation. 
    # actually, not exact distance to the observation but an approximate value
    # exact value is equal to the distance from the river mouth to the closest vertex on river stream 
  }
  
  # (1) network.line.coords: network.line.coords
  network.line.coords <- data.frame(NetworkID = as.factor(rep(1, nrow(shape@data))), SegmentID = c(1:nrow(shape@data)), 
                                    DistanceUpstream = DistanceUpstream_seg)
  # observeddata <- SSNPoints
  
  # (2) obspoint: obspoints
  network.line.coords.obsdata <- data.frame(NetworkID = rep(1, nrow(data$eventplace)), SegmentID = dist.tmp[, 2], DistanceUpstream = DistanceUpstream_obs)
  # SegmentID : segment ID where the observation exists
  # DistanceUpstream : approximate distance from the river mouth to the observation
  
  pointcoords <- cbind(data$eventplace$lon, data$eventplace$lat); colnames(pointcoords) <- c("coords.x1", "coords.x2")
  
  pointbbox <- rbind(range(pointcoords[, 1]), range(pointcoords[, 2]))
  colnames(pointbbox) <- c("min", "max"); rownames(pointbbox) <- c("coords.x1", "coords.x2")
  
  data$eventplace$RCH_ID <- shape@data$RCH_ID[dist.tmp[, 2]]
  pointdata <- cbind(data$eventplace, rid_old = dist.tmp[, 2], upDist = DistanceUpstream_obs, locID = c(1:nrow(data$eventplace)), 
                     netID = rep(1,nrow(data$eventplace)), pid = c(1:nrow(data$eventplace)), shreve = shreve_obj[dist.tmp[, 2]], 
                     mindist = dist.tmp[, 1], rid = dist.tmp[, 2], segid = dist.tmp[, 3])
  
  network.SSNPoints <- list()
  network.SSNPoints[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.obsdata, point.coords = pointcoords, point.data = pointdata, 
                                points.bbox = pointbbox, proj4string = shape@proj4string)
  obspoints <- new("SSNPoints", SSNPoints = network.SSNPoints, ID = "Obs")
  
  # (3) predpoint
  predpoints <- new("SSNPoints", SSNPoints = list(), ID = character(0))
  
  # (4) SSNpath: path
  SSNpath <- "NULL"
  # (5) data: data.frame type ( rid,upDist,    Length, netID)
  # data:Object of class "data.frame". 
  # The number of rows in data should equal the number of lines in the lines object. Row names correspond to SegmentID values
  shapedata <- shape@data
  shapedata <- cbind(shapedata, rid=c(1:nrow(shape@data)), upDist = DistanceUpstream_seg, Length = shape@data$Shape_Leng, 
                     netID = rep(1,nrow(shape@data)), shreve = shreve_obj)
  # (6) lines: segment polygon
  shapelines <- shape@lines
  # (7) bbox
  shapebbox <- shape@bbox; colnames(shapebbox) = c("min", "max"); rownames(shapebbox) <- c("x", "y")
  # (8) proj4string
  shapeproj4string <- shape@proj4string
  
  river_original <- new("SpatialStreamNetwork", network.line.coords = network.line.coords, obspoints = obspoints, predpoints = predpoints, 
                        path = SSNpath, data = shapedata, lines = shapelines, bbox = shapebbox, proj4string = shapeproj4string)
  
  # Object of class "SSNPoints" with 2 slots
  #
  # @ SSNPoints: List of SSNPoint objects with 5 slots
  # @ network.point.coords: object of class "data.frame". Row names
  # represent point identifiers (pid) stored in the point.data
  # data.frame.
  # $ NetworkID: factor identifying the NetworkID of that point
  # $ SegmentID: factor identifying the unique stream segment of that point
  # $ DistanceUpstream: numeric value representing the cumulative
  # distance from the network outlet, the most downstream point
  # on a network, to that point
  # @ point.coords: numeric matrix or "data.frame" with x- and y-
  #   coordinates (each row is a point); row names represent point
  # identifiers (pid) stored in the point.data data.frame.
  # @ point.data: object of class "data.frame"; the number of rows in
  # data should equal the number of points in the
  # network.point.coords object; row names are set to the pid
  # attribute.
  # @ points.bbox: Object of class "matrix"; see Spatial-class
  # @ proj4string: Object of class "CRS"; see CRS-class
  # @ ID: character string representing the name of the observation points
  
  
  # network.line <- shape_tinned@lines
  # bbox, proj4string
  
  return(river_original)
}


calculate_UpstreamDist <- function (ssndata, spatialwt)
{
  n <- nrow(spatialwt)
  updist_mat <- matrix(0, nrow = n, ncol = n)
  rownames(updist_mat) <- rownames(spatialwt) 
  colnames(updist_mat) <- colnames(spatialwt)
  
  updist_info <- ssndata@obspoints@SSNPoints[[1]]@network.point.coords
  rownames(updist_info) <- c(1:n)
  
  for (i in 1 : n) {
    for (j in 1 : n) {
      if (spatialwt[i,j]!=0) {
        updist_mat[i,j] <- abs(updist_info[colnames(updist_mat)[j],'DistanceUpstream'] -
                                 updist_info[colnames(updist_mat)[i],'DistanceUpstream'])
      } else {
        updist_mat[i,j] <- max(updist_info['DistanceUpstream'])*100 # for calculation, set very large value instead of Inf
      }
    }
  }
  
  return(updist_mat)
}


wpca.adjust <- function(x,wt,rwt,cwt,...) {
  local.center <- function(x,wt)  
    sweep(x,2,colSums(sweep(x,1,wt,'*'))/sum(wt))
  svd(rwt %*% sweep(local.center(x,wt),1,sqrt(wt),'*') %*% cwt,...)}


pca.cv.river <- function(bw,x,loc,k=2,kernel="bisquare",adaptive=FALSE,p=2, theta=0, longlat=F, dMat, rwt, cwt) 
{
  if (missing(dMat))
    DM.given<-F
  else
    DM.given<-T
  n <- nrow(loc) # number of locations
  m <- ncol(x) # number of variables
  w <- array(0,c(n,m,k)) # array for local loadings
  
  score <- 0
  pcafun <- wpca.adjust
  
  d <- matrix(0,n,m) # matrix for eigenvalues 
  
  miss <- c() # store locations with insufficient upstream connected sites
  for (i in 1:n) 
  {
    if (DM.given)
      dist.vi<-dMat[i,] # distance vector for ith location
    else
    {
      dist.vi<-gw.dist(loc, focus=i, p=p, theta=theta, longlat=longlat) 
    }
    wt <-gw.weight(dist.vi,bw,kernel,adaptive) # kernel weight vector 
    # wt[i] <- 0
    use <- wt>0
    # wt <-wt[use]
    if(length(wt[use])<m)
    {
      expr <- paste("Too small bandwidth at location: ", i)
      warning(paste(expr, "and the results can't be given there.", sep=", "))
      miss <- c(miss, i)
      next
    }
    temp <- pcafun(x,wt,rwt,cwt,nu=0,nv=k) # svd result of s^-1/2 Wi*^1/2 X T^-1/2
    w[i,,] <- t(solve(cwt))%*%temp$v # loadings after back transformation 
    d[i,] <- temp$d # eigenvalue
    
    # cv score
    score <- score + sum((x[i,] - x[i,] %*% cwt %*% temp$v %*% t(temp$v) %*% solve(cwt)))**2
  }
  
  # for insufficient locations, find the nearest downstream site
  for(j in miss){
    ord <- order(dMat[,j]) # index set in order of close point from j th point
    ind <- ord[which(!(ord %in% miss))[1]]
    d[j,] <- d[ind,]
    w[j,,] <- w[ind,,]
    score <- score + sum((x[j,] - x[j,] %*% cwt %*% t(cwt) %*% w[j,,] %*% t(w[j,,])))**2
  }
  
  if(adaptive)
    cat("Adaptive bandwidth(number of nearest neighbours):", bw, "CV score:", score, "\n")
  else
    cat("Fixed bandwidth:", bw, "CV score:", score, "\n")
  score
}


bw.pca.river <- function(data,vars,k=2,kernel="bisquare",adaptive=FALSE,p=2, theta=0, longlat=F,dMat,rwt,cwt){
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
  }
  else if (is(data, "data.frame")&&(!missing(dMat)))
    data<-data
  else
    stop("Given data must be a Spatial*DataFrame or data.frame object")
  data <- as(data, "data.frame")
  dp.n<-nrow(data)
  if (missing(dMat))
  {
    DM.given<-F
    if(dp.n <= 5000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop("Dimensions of dMat are not correct")
  }
  if(missing(vars))
    stop("Variables input error")
  ##########Extract data
  col.nm<-colnames(data)
  var.idx<-match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if(length(var.idx)==0) stop("Variables input doesn't match with data")
  x<-data[,var.idx]
  x<-as.matrix(x)
  var.nms<-colnames(x)
  var.n<-ncol(x)
  #########Find the range of the fixed bandwidth
  if(adaptive)
  {
    upper<-dp.n
    lower<-2
  }
  else
  {
    if(DM.given)
    {
      upper<-range(dMat)[2]
      lower<-upper/5000
    }
    ###!!!!!!!! Important note: if the distance matrix is not specified, the running time will be consuming very much by choosing the range of fixed bandwidth when the p is not 2; 
    ### because the range can't be decided by boundary box
    else
    {
      dMat<-NULL
      if (p==2)
      {
        b.box<-bbox(dp.locat)
        upper<-sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower<-upper/5000
      }
      else
      {
        upper<-0
        for (i in 1:dp.n)
        {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
          upper<-max(upper, range(dist.vi)[2])
        }
        lower<-upper/5000
      }
    }
    
  }
  ########################## Now the problem for the golden selection is too computationally heavy
  #Select the bandwidth by golden selection
  # Add this bit please ############################################################
  #  if(robust==FALSE)
  #    pcafun=wpca
  #  else
  #    pcafun=rwpca
  ##################################################################################
  bw<-NA
  bw <- gold(pca.cv.river, lower,upper, adapt.bw=adaptive, x, dp.locat, k, 
             kernel, adaptive, p, theta, longlat, dMat, rwt, cwt) 
  # bw<-NA
  #    if(approach=="cv"||approach=="CV")
  #       bw <- optimize(bw.cv,lower=lower,upper=upper,maximum=FALSE,X=x,Y=y,kernel=kernel,
  #       adaptive=adaptive, dp.locat=dp.locat, p=p, theta=theta, longlat=longlat,dMat=dMat,tol=.Machine$double.eps^0.25)
  #    else if(approach=="aic"||approach=="AIC"||approach=="AICc")
  #       bw<-optimize(bw.aic,lower=lower,upper=upper,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
  bw
}


pca.cv.contrib.river <- function(x,loc,bw, k=2,kernel="bisquare",adaptive=FALSE,p=2, theta=0, longlat=F,dMat,rwt,cwt)
{
  if (missing(dMat))
    DM.given<-F
  else
    DM.given<-T
  n <- nrow(loc) # number of locations
  m <- ncol(x) # number of variables
  w <- array(0,c(n,m,k)) # array for local loadings
  score <- numeric(n)
  pcafun <- wpca.adjust
  d <- matrix(0,n,m) # matrix for eigenvalues 
  miss <- c()
  
  for (i in 1:n) 
  {
    if (DM.given)
      dist.vi<-dMat[i,]
    else
    {
      dist.vi<-gw.dist(loc, focus=i, p=p, theta=theta, longlat=longlat) 
    }
    wt <-gw.weight(dist.vi,bw,kernel,adaptive)
    # wt[i] <- 0
    use <- wt>0
    wt <-wt[use]
    if(length(wt[use])<m)
    {
      expr <- paste("Too small bandwidth at location: ", i)
      warning(paste(expr, "and the results can't be given there.", sep=", "))
      miss <- c(miss, i)
      next
    }
    temp <- pcafun(x,wt,rwt,cwt,nu=0,nv=k) # svd result of s^-1/2 Wi*^1/2 X T^-1/2
    w[i,,] <- t(solve(cwt))%*%temp$v # loadings after back transformation 
    d[i,] <- temp$d # eigenvalue
    
    score[i] <- sum((x[i,] - x[i,] %*% cwt %*% temp$v %*% t(temp$v) %*% solve(cwt)))**2
  }
  score
}


# Main GWPCA function
# x is the data matrix,
# loc is a 2-column coordinate matrix
# bw is the bandwidth in fixed or adaptive form,
# k is the number of retained components
# eloc are locations to evaluate the loadings,  if different from loc - eloc is also a 2 column matrix
# pcafun is the weighted PCA function which can be in a basic or robust form (from before)
# This returns a list with <loadings> for each eloc; d for each eloc; and the bandwidth used	
#gwpca <- function(x,loc,bw,k=2,eloc=loc,pcafun=wpca,...) 

pca.river <- function (data, elocat, vars, k = 2, kernel = "bisquare",
                       adaptive = FALSE, bw, p = 2, theta = 0, longlat = F, cv = T, scores=T,
                       dMat, rwt, cwt)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    polygons <- NULL
    if(is(data, "SpatialPolygonsDataFrame"))
      polygons <- polygons(data)
  }
  else
    stop("Given data must be a Spatial*DataFrame or data.frame object")
  
  if (missing(elocat))
  {
    ep.given <- FALSE
    elocat<-coordinates(data)
  }
  else 
  {
    ep.given <- T
    if (is(elocat, "Spatial"))
    {
      espdf<-elocat
      elocat<-coordinates(espdf)
    }
    else if (is.numeric(elocat) && dim(elocat)[2] == 2)
      elocat<-elocat
    else
    {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      elocat<-dp.locat
    }
  }
  data <- as(data, "data.frame")
  dp.n<-nrow(data)
  ep.n<-nrow(elocat)
  if (missing(dMat))
  {
    DM.given<-F
    DM1.given<-F
    if(dp.n + ep.n <= 10000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=elocat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    DM1.given<-T 
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=ep.n)
      stop("Dimensions of dMat are not correct")
  }
  if(missing(vars))
    stop("Variables input error")
  if(missing(bw)||bw<=0)
    stop("Bandwidth is not specified incorrectly")
  len.var <- length(vars)
  ##########Extract data
  col.nm<-colnames(data)
  var.idx<-match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if(length(var.idx)==0) stop("Variables input doesn't match with data")
  x<-data[,var.idx]
  x<-as.matrix(x)
  var.nms<-colnames(x)
  var.n<-ncol(x)
  if(len.var > var.n)
    warning("Invalid variables have been specified, please check them again!")
  #####PCA
  pca.res <- princomp(x, cor=T, scores=scores)
  ############################# WPCA
  #The local loading for the principle components
  w <- array(0,c(ep.n,var.n,k))
  #Return the local scores
  if(scores)
    pca.river.scores <- list()
  else
    pca.river.scores <- NULL
  #
  d <- matrix(0,ep.n,var.n)
  
  pcafun = wpca.adjust
  miss = c()
  for (i in 1:ep.n) 
  {
    if (DM.given)
      dist.vi<-dMat[i,]
    else
    {
      if (ep.given)
        dist.vi<-gw.dist(dp.locat, elocat, focus=i, p, theta, longlat) 
      else
        dist.vi<-gw.dist(dp.locat, focus=i, p=p, theta=theta, longlat=longlat) 
    }
    wt <-gw.weight(dist.vi,bw,kernel,adaptive)
    use <- wt>0
    # wt <-wt[use]
    if(length(wt[use])<var.n)
    {
      expr <- paste("Too small bandwidth at location: ", i)
      warning(paste(expr, "and the results can't be given there.", sep=", "))
      miss <- c(miss, i)
      next
    }
    temp <- pcafun(x,wt,rwt,cwt,nu=0,nv=k)
    w[i,,] <- t(solve(cwt))%*%temp$v
    d[i,] <- temp$d
    
    #Calculate the local scores
    
  }
  if (!is.null(rownames(x))) dimnames(w)[[1]] <- rownames(x)
  if (!is.null(colnames(x))) dimnames(w)[[2]] <- colnames(x)
  dimnames(w)[[3]] <- paste("PC",1:k,sep='')
  CV<-numeric(dp.n)
  if(cv)
    CV<-pca.cv.contrib.river(x,dp.locat,bw,k,kernel,adaptive, p, theta, longlat,dMat,rwt,cwt)
  GW.arguments<-list(vars=vars,k=k, bw=bw, kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat, dp.n=dp.n, DM.given=DM.given,scores=scores)
  # And add and change this bit please #################################################
  
  for(j in miss){
    ord <- order(dMat[,j]) # index set in order of close point from j th point
    ind <- ord[which(!(ord %in% miss))[1]]
    d[j,] <- d[ind,]
    w[j,,] <- w[ind,,]
  }
  
  if(scores)
  {
    for(i in 1:ep.n){
      scores.i <- c()
      for(j in 1:k)
      {
        score <- as.numeric(x[i,]%*%cwt%*%t(cwt)%*%w[i,,j])
        scores.i <- cbind(scores.i, score)
      }
      pca.river.scores[[i]] <- scores.i
    }
  }
  
  d1 <- (d/(sum(wt)^0.5))^2
  
  local.PV <- d1[, 1:k]/rowSums(d1) * 100
  var.names <- c()
  for(i in 1:k)
    var.names <- c(var.names, paste(paste("Comp", i, sep="."), "PV", sep="_"))
  win.var.pc1 <- max.col(abs(w[,,1]))
  res.df <- data.frame(local.PV, rowSums(local.PV), vars[win.var.pc1])
  names(res.df) <- c(var.names, "local_CP", "win_var_PC1")
  
  if (!is.null(polygons))
  {
    rownames(res.df) <- sapply(slot(polygons, "polygons"),
                               function(i) slot(i, "ID"))
    SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=res.df,match.ID=F)
  }
  else
  {
    SDF <- SpatialPointsDataFrame(coords=elocat, data=res.df, proj4string=CRS(p4s), match.ID=FALSE)
  }
  ################################################
  timings[["stop"]] <- Sys.time()
  res <- list(pca=pca.res,loadings = w, SDF=SDF, pca.river.scores=pca.river.scores, var=d1,  
              local.PV=local.PV, GW.arguments = GW.arguments, CV = CV, timings=timings)
  class(res) <-"riverpca"
  invisible(res)
  ##################################################################################
}

print.pca.river<-function(x, ...)
{
  if(class(x) != "riverpca") stop("It's not a riverpca object")
  cat("   ***********************************************************************\n")
  cat("   *                           Basic Info                                *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  vars <- x$GW.arguments$vars
  cat("\n   Variables concerned: ",vars)
  cat("\n   The number of retained components: ",x$GW.arguments$k)
  dp.n<- dim(x$loadings)[1]
  cat("\n   Number of data points:",dp.n)
  ################################################################ Print Linear
  cat("\n   ***********************************************************************\n")
  cat("   *                Results of Principal Components Analysis               *\n")
  cat("   ***********************************************************************\n")
  print(summary(x$pca,  loadings = TRUE, cutoff=0))
  
  #########################################################################
  cat("\n   ***********************************************************************\n")
  cat("   *           Results of Principal Components Analysis on River         *\n")
  cat("   ***********************************************************************\n")
  cat("\n   *********************Model calibration information*********************\n")
  cat("   Kernel function for geographically weighting:", x$GW.arguments$kernel, "\n")
  if(x$GW.arguments$adaptive)
    cat("   Adaptive bandwidth for geographically and temporally  weighting: ", x$GW.arguments$bw, " (number of nearest neighbours)\n", sep="")
  else
    cat("   Fixed bandwidth for geographically and temporally weighting: ", x$GW.arguments$bw, "\n")
  if (x$GW.arguments$DM.given)
    cat("   Distance metric for geographically weighting: A distance matrix is specified for this model calibration.\n")
  else
  {
    if (x$GW.arguments$longlat)
      cat("   Distance metric for geographically weighting: Great Circle distance metric is used.\n")
    else if (x$GW.arguments$p==2)
      cat("   Distance metric for geographically weighting: Euclidean distance metric is used.\n")
    else if (x$GW.arguments$p==1)
      cat("   Distance metric for geographically weighting: Manhattan distance metric is used.\n") 
    else if (is.infinite(x$GW.arguments$p))
      cat("   Distance metric for geographically weighting: Chebyshev distance metric is used.\n")
    else 
      cat("   Distance metric for geographically weighting: A generalized Minkowski distance metric is used with p=",x$GW.arguments$p,".\n")
    if (x$GW.arguments$theta!=0&&x$GW.arguments$p!=2&&!x$GW.arguments$longlat)
      cat("   Coordinate rotation: The coordinate system is rotated by an angle", x$GW.arguments$theta, "in radian.\n")   
  } 
  
  cat("\n   ****************     Summary of river PCA information:    *****************\n")       
  var.names <- c()
  for(i in 1:x$GW.arguments$k)
    var.names <- c(var.names, paste("Comp", i, sep="."))
  cat("   Local variance: \n")
  local.SD <- t(apply(x$var[,1:x$GW.arguments$k], 2, summary))[,c(1:3,5,6)]
  rownames(local.SD) <- var.names
  if(x$GW.arguments$k==1) 
  { 
    local.SD <- matrix(local.SD, nrow=1)
    colnames(local.SD) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
    rownames(local.SD) <- var.names
  }
  rnames<-rownames(local.SD)
  for (i in 1:length(rnames))
    rnames[i]<-paste("   ",rnames[i],sep="")
  rownames(local.SD) <-rnames 
  printCoefmat(local.SD)
  cat("   Local Proportion of Variance: \n")
  local.PV <-  t(apply(as(x$SDF, "data.frame")[,1:(x$GW.arguments$k+1), drop=FALSE], 2, summary))[,c(1:3,5,6)]
  if(x$GW.arguments$k==1) 
  { 
    local.PV <- matrix(local.SD, nrow=1)
    colnames(local.PV) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
  }
  rownames(local.PV) <-c(rnames, paste("   ","Cumulative",sep=""))
  printCoefmat(local.PV)
  cat("\n   ***********************************************************************\n")
  cat("   Program stops at:", as.character(x$timings$stop), "\n")
  invisible(x)
}

plot_river <- function (x, plots = c("score", "biplot", "glyph", "info", "riverpca"),
                        pc.map = 1, probs = c(0, 0.25, 0.75, 1),
                        pc.biplot = c(1,2), biplot.factor = 100,
                        pc.glyph = c(1, 2, 3), river.glyph.loading = c(1,2,3), river.glyph.score = c(1,2,3), river, r1 = 10,
                        coords, pc.ts = c(1, 2, 3), pc.meanPM = c(1, 2,3),
                        meanPM.factor = 0.05, lwd=2, zoom=FALSE, zoomx=c(0,1), zoomy=c(0,1))
{
  pca.mode <- x$pca.mode
  if(zoom==FALSE)
  {
    zoomx <- c(river@bbox[1, ])
    zoomy <- c(river@bbox[2, ])
  }
  if (pca.mode == "Smode")
  {
    message("For river PCA, only support T-mode PCA")
    return()
  }
  else
  {
    if ("biplot" %in% plots) 
    {
      nplot <- length(x$pca.wt)
      if (nplot == 4) 
      {
        par(mfrow = c(2, 2),oma = c(2,2,2,2) + 0.1, mar = c(4,4,2,1) + 0.1)
      }
      else 
      {
        par(mfrow = c(1, nplot))
      }
      pca.type <- x$pca.wt
      xlims <- matrix(nrow = 2, ncol = length(pca.type))
      ylims <- matrix(nrow = 2, ncol = length(pca.type))
      
      for (i in 1:nplot) 
      {
        pca.type1 <- pca.type[i]
        scores <- x[[pca.type1]]$scores
        rownames(scores) <- x[["row.names"]]
        loadings <- x[[pca.type1]]$loads
        rownames(loadings) <- x[["col.names"]]
        var.val <- x[[pca.type1]]$var
        sd.val = sqrt(var.val)
        xlims[, i] <- range(scores[, pc.biplot[1]]/sd.val[pc.biplot[1]], 
                            loadings[, pc.biplot[1]] * sd.val[pc.biplot[1]]/biplot.factor * 
                              1.2)
        ylims[, i] <- range(scores[, pc.biplot[2]]/sd.val[pc.biplot[2]], 
                            loadings[, pc.biplot[2]] * sd.val[pc.biplot[2]]/biplot.factor * 
                              1.2)
      }
      for (i in 1:nplot) 
      {
        pca.type1 <- pca.type[i]
        scores <- x[[pca.type1]]$scores
        rownames(scores) <- x[["row.names"]]
        loadings <- x[[pca.type1]]$loads
        rownames(loadings) <- x[["col.names"]]
        var.val <- x[[pca.type1]]$var
        var.prop <- x[[pca.type1]]$var.prop
        pc1 = 100 * round(var.prop[pc.biplot[1]], digits = 2)
        pc2 = 100 * round(var.prop[pc.biplot[2]], digits = 2)
        xlab = paste("PC", pc.biplot[1], ": ", pc1, "%", 
                     sep = "")
        ylab = paste("PC", pc.biplot[2], ": ", pc2, "%", 
                     sep = "")
        sd.val = sqrt(var.val)
        scores.min = min(scores[, pc.biplot])
        scores.max = max(scores[, pc.biplot])
        xlim <- c(min(xlims[1, ]), max(xlims[2, ]))
        ylim <- c(min(ylims[1, ]), max(ylims[2, ]))
        plot(scores[, pc.biplot[1]]/sd.val[pc.biplot[1]], 
             scores[, pc.biplot[2]]/sd.val[pc.biplot[2]], 
             main = pca.type1, xlab = xlab, 
             ylab = ylab, type = "n", xlim = xlim, ylim = ylim, 
             col.axis = "black", cex.main = 1.5, cex.lab = 1.5, 
             cex.axis = 1.5)
        arrows(0, 0, loadings[, pc.biplot[1]] * sd.val[pc.biplot[1]]/biplot.factor, 
               loadings[, pc.biplot[2]] * sd.val[pc.biplot[2]]/biplot.factor, 
               length = 0.1, lwd = lwd, angle = 20, col = "red")
        text(loadings[, pc.biplot[1]] * sd.val[pc.biplot[1]]/biplot.factor * 
               1.2, loadings[, pc.biplot[2]] * sd.val[pc.biplot[2]]/biplot.factor * 
               1.2, labels = rownames(loadings), col = "red", 
             cex = 1.2, font=2)
        text(scores[, pc.biplot[1]]/sd.val[pc.biplot[1]], 
             scores[, pc.biplot[2]]/sd.val[pc.biplot[2]], 
             labels = rownames(scores), col = "blue", cex = 0.9, font=2)
        abline(0, 0, col = "black")
        abline(v = 0, col = "green")
      }
    }
    if ("glyph" %in% plots) 
    {
      if (length(pc.glyph) < 3) 
      {
        message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
      }
      nplot <- length(x$pca.wt)
      if (nplot == 4) 
      {
        par(mfrow = c(2, 2))
      }
      else 
      {
        par(mfrow = c(nplot, 1))
      }
      pca.type <- x$pca.wt
      
      for (i in 1:nplot) 
      {
        glyph.data <- x[[pca.type[i]]]
        scores <- glyph.data$scores[, pc.glyph]
        plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
             ylab = names(coords)[2], main = paste("Scores: weight =", 
                                                   pca.type[i]), cex.main = 1.5, cex.lab = 1.5, 
             cex.axis = 1.5, xlim = zoomx, 
             ylim = zoomy)
        lines(river, col = "lightgrey")
        GWmodel::glyph.plot(as.matrix(scores), coords, 
                            r1 = r1, add = TRUE)
        
      }
      message("For plots = glyph: did you check that coordinates for monitoring site ID's are in \nthe same order as data used for stpca()?  \nIf not then glyphs might not correspond to the correct monitoring site.")
    }
    
    if ("score" %in% plots)
    {
      nplot <- length(x$pca.wt)
      pca.type <- x$pca.wt
      if (nplot == 4) 
      {
        par(mfrow = c(2, 2))
      }
      for (i in 1:nplot) 
      {
        title <- paste("Scores PC", pc.map, "\n weight =", 
                       pca.type[i])
        map.data <- x[[pca.type[i]]]
        scores <- map.data$scores[, pc.map]
        n.colors <- length(probs) - 1
        data.vec <- sort(scores)
        quants <- quantile(data.vec, probs = probs)
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(n.colors)))(n.colors)
        par(mar = c(5, 4, 4, 2) + 0.1, mfrow=c(2,2))
        graphics::layout(matrix(1:2, ncol = 2), widths = c(0.85, 
                                                           0.15), heights = c(1, 1))
        plot(coords[, 1], coords[, 2], type = "n", xlab = "Easting", 
             ylab = "Northing", main = title, cex.main = 1.5, 
             cex.axis = 1.5, cex.lab = 1.5, xlim = c(river@bbox[1, 
             ]), ylim = c(river@bbox[2, ]))
        lines(river)
        points(coords[, 1], coords[, 2], pch = 19, cex = 2, 
               col = getPalette[cut(scores, breaks = quants)])
        par(mar = c(5.1, 0.4, 4.1, 0.5))
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(100)))(100)
        plot(NA, type = "n", ann = FALSE, xlim = c(0, 
                                                   2.1), ylim = c(0, 101), xaxt = "n", yaxt = "n", 
             bty = "n")
        rect(1, c(0:99), 1.5, c(1:100), col = c(getPalette), 
             border = NA)
        a <- n.colors
        b <- 100/n.colors
        text(x = 0.5, y = ((0:a) * b), labels = round(quants, 
                                                      2), cex = 1.5)
      }
    }
    if ("info" %in% plots) 
    {
      info.num <- length(x$info)
      info.type <- x$info
      for (i in 1:info.num) 
      {
        title <- paste("Plot for", info.type[i])
        map.data <- x[[info.type[i]]]
        if(info.type[i] == "PTV")
        {
          scores <- x[[info.type[i]]]
          n.colors <- length(probs) - 1
          data.vec <- sort(scores)
          quants <- quantile(data.vec, probs = probs)
        }
        else
        {
          scores <- x[[info.type[i]]]
          nums = length(unique(scores))
          if(nums==1) nums=2
          n.colors <- nums
          data.vec <- sort(scores)
          quants <- nums
        }
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(n.colors)))(n.colors)
        par(mar = c(5, 4, 4, 2) + 0.1)
        graphics::layout(matrix(1:2, ncol = 2), widths = c(0.85, 
                                                           0.15), heights = c(1, 1))
        plot(coords[, 1], coords[, 2], type = "n", xlab = "Easting", 
             ylab = "Northing", main = title, cex.main = 1.5, 
             cex.axis = 1.5, cex.lab = 1.5, xlim = c(river@bbox[1, 
             ]), ylim = c(river@bbox[2, ]))
        lines(river)
        points(coords[, 1], coords[, 2], pch = 19, cex = 2, 
               col = getPalette[cut(scores, breaks = quants)])
        par(mar = c(5.1, 0.4, 4.1, 0.5))
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(100)))(100)
        plot(NA, type = "n", ann = FALSE, xlim = c(0, 
                                                   2.1), ylim = c(0, 101), xaxt = "n", yaxt = "n", 
             bty = "n")
        rect(1, c(0:99), 1.5, c(1:100), col = c(getPalette), 
             border = NA)
        a <- n.colors
        b <- 100/n.colors
        text(x = 0.5, y = ((0:a) * b), labels = round(quants, 
                                                      2), cex = 1.5)
      }
    }
    if ("riverpca" %in% plots)
    {
      if (length(river.glyph.loading) < 3) 
      {
        message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
      }
      nplot <- length(river.glyph.loading)
      if (nplot == 4) 
      {
        par(mfrow = c(2, 2))
      }
      else 
      {
        par(mfrow = c(1,nplot), oma = c(2,2,2,2) + 0.1, mar = c(4,5,2,2) + 0.1)
      }
      pca.type <- "riverpca"
      for (i in 1:nplot) 
      {
        glyph.data <- x[[pca.type]]$loadings
        local.loadings <- glyph.data[,,river.glyph.loading[i]]
        # scores <- glyph.data$scores[, pc.glyph]
        plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
             ylab = names(coords)[2], main = paste("Loadings for PC", i), cex.main = 1.5, cex.lab = 1.5, 
             cex.axis = 1.5, xlim = zoomx, 
             ylim = zoomy)
        lines(river, col = "black")
        glyph.plot.river(as.matrix(local.loadings), coords, 
                         r1 = r1, add = TRUE, lwd=lwd)
        if(zoom==FALSE)
        {
          glyph.plot.river(rbind(t(x$spatiotemporal$loads[,i]),0), rbind(c(min(coords[,1])+0.1,max(coords[,2])-0.05),0), 
                           r1 = 40*r1, add = TRUE, lwd=lwd) # manually control
          text(min(coords[,1])+0.1,max(coords[,2])-0.2,expression(TPCA[ST]), cex=lwd, font=2) # manually control
        }
      }
      message("For plots = glyph: did you check that coordinates for monitoring site ID's are in \nthe same order as data used for stpca()?  \nIf not then glyphs might not correspond to the correct monitoring site.")
      
      if (length(river.glyph.score) < 3) 
      {
        message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
      }
      
      glyph.data <- x[[pca.type]]$scores
      local.scores <- glyph.data[,river.glyph.score]
      # scores <- glyph.data$scores[, pc.glyph]
      
      par(mfrow = c(1,2), oma = c(2,2,2,2) + 0.1, mar = c(4,5,2,2) + 0.1)
      
      plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
           ylab = names(coords)[2], main = paste("Scores: unweighted pca"), cex.main = 1.5, cex.lab = 1.5, 
           cex.axis = 1.5, xlim = zoomx, 
           ylim = zoomy)
      lines(river, col = "black")
      glyph.plot.river(as.matrix(x[["unweighted"]]$scores[,river.glyph.score]), coords, 
                       r1 = r1, add = TRUE, lwd=lwd)
      
      # plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
      #      ylab = names(coords)[2], main = NULL, cex.main = 1.5, cex.lab = 1.5, 
      #      cex.axis = 1.5, xlim = zoomx, 
      #      ylim = zoomy)
      # lines(river, col = "black")
      # glyph.plot.river(as.matrix(x[["spatiotemporal"]]$scores[,river.glyph.score]), coords, 
      #                  r1 = r1, add = TRUE, lwd=lwd)
      
      plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
           ylab = names(coords)[2], main = paste("Scores: river pca"), cex.main = 1.5, cex.lab = 1.5, 
           cex.axis = 1.5, xlim = zoomx, 
           ylim = zoomy)
      lines(river, col = "black")
      glyph.plot.river(as.matrix(local.scores), coords, 
                       r1 = r1, add = TRUE, lwd=lwd)
      
      message("For plots = glyph: did you check that coordinates for monitoring site ID's are in \nthe same order as data used for stpca()?  \nIf not then glyphs might not correspond to the correct monitoring site.")
      
    }
  }
}


glyph.plot.river <- function (ld, loc, r1 = 50, add = FALSE, alpha = 1, sep.contrasts = FALSE, lwd=2) 
{
  r <- max(max(loc[, 1]) - min(loc[, 1]), max(loc[, 2]) - min(loc[,2]))/r1
  glyph.plot1 <- function(ld, loc, r = max(max(loc[, 1]) - 
                                             min(loc[, 1]), max(loc[, 2]) - min(loc[, 2]))/r1, add = FALSE, 
                          alpha = 1, lwd) {
    rowmax <- function(z) z[cbind(1:nrow(z), max.col(abs(z)))]
    ld <- sweep(ld, 1, sign(rowmax(ld)), "*")
    ld.max <- max(ld)
    ld.min <- min(ld)
    n.col <- ncol(ld)
    n.row <- nrow(ld)
    angles <- (0:(n.col - 1)) * pi/n.col
    J <- 0 + (0+1i)
    disp <- exp((pi/2 - angles) * J) * r
    loc2 <- loc[, 1] + loc[, 2] * J
    ld.scaled <- 2 * ld/(max(ld.max, -ld.min))
    if (!add) 
      plot(loc, asp = 1, type = "n")
    points(loc, pch = 16, cex = 0.1, col = "black")
    for (i in 1:n.row) {
      for (j in 1:n.col) {
        l.from <- loc2[i]
        l.to <- loc2[i] + disp[j] * ld.scaled[i, j]
        col <- if (ld[i, j] > 0) {
          rgb(0, 0, 1, alpha)
        }
        else {
          rgb(1, 0, 0, alpha)
        }
        lines(Re(c(l.from, l.to)), Im(c(l.from, l.to)), 
              col = col, lwd=lwd)
      }
    }
  }
  glyph.plot2 <- function(ld, loc, r = max(max(loc[, 1]) - 
                                             min(loc[, 1]), max(loc[, 2]) - min(loc[, 2])) / r1, add = FALSE, 
                          alpha = 1, lwd) {
    rowmax <- function(z) z[cbind(1:nrow(z), max.col(abs(z)))]
    ld <- sweep(ld, 1, sign(rowmax(ld)), "*")
    ld.max <- max(ld)
    ld.min <- min(ld)
    n.col <- ncol(ld)
    n.row <- nrow(ld)
    angles <- (0:(n.col - 1)) * 2 * pi/n.col
    J <- 0 + (0+1i)
    disp <- exp((pi/2 - angles) * J) * r
    loc2 <- loc[, 1] + loc[, 2] * J
    ld.scaled <- abs(ld)/(max(ld.max))
    if (!add) 
      plot(loc, asp = 1, type = "n")
    points(loc, pch = 16, cex = 0.1, col = "black")
    for (i in 1:n.row) {
      for (j in 1:n.col) {
        l.from <- loc2[i]
        l.to <- loc2[i] + disp[j] * ld.scaled[i, j]
        col <- if (ld[i, j] > 0) {
          rgb(0, 0, 1, alpha)
        }
        else {
          rgb(1, 0, 0, alpha)
        }
        lines(Re(c(l.from, l.to)), Im(c(l.from, l.to)), 
              col = col, lwd=lwd)
      }
    }
  }
  if (sep.contrasts) {
    glyph.plot1(ld, loc, r, add, alpha, lwd)
  }
  else {
    glyph.plot2(ld, loc, r, add, alpha, lwd)
  }
}



plot_stpca_zoom <- function (x, plots = c("map", "biplot", "glyph", "ts", "meanPM"), 
                             pc.map = 1, probs = c(0, 0.25, 0.75, 1), pc.biplot = c(1, 
                                                                                    2), biplot.factor = 100, pc.glyph = c(1, 2, 3), river, 
                             r1 = 10, coords, pc.ts = c(1, 2, 3), pc.meanPM = c(1, 2, 
                                                                                3), meanPM.factor = 0.05, lwd, zoom=F, zoomx=c(0,1), zoomy=c(0,1)) 
{
  pca.mode <- x$pca.mode
  if(zoom==FALSE){
    zoomx <- c(river@bbox[1, ])
    zoomy <- c(river@bbox[2, ])
  }
  if ("biplot" %in% plots) {
    nplot <- length(x$pca.wt)
    if (nplot == 4) {
      par(mfrow = c(2, 2))
    }
    else {
      par(mfrow = c(1, nplot))
    }
    pca.type <- x$pca.wt
    xlims <- matrix(nrow = 2, ncol = length(pca.type))
    ylims <- matrix(nrow = 2, ncol = length(pca.type))
    for (i in 1:nplot) {
      pca.type1 <- pca.type[i]
      scores <- x[[pca.type1]]$scores
      rownames(scores) <- x[["row.names"]]
      loadings <- x[[pca.type1]]$loads
      rownames(loadings) <- x[["col.names"]]
      var.val <- x[[pca.type1]]$var
      sd.val = sqrt(var.val)
      xlims[, i] <- range(scores[, pc.biplot[1]]/sd.val[pc.biplot[1]], 
                          loadings[, pc.biplot[1]] * sd.val[pc.biplot[1]]/biplot.factor * 
                            1.2)
      ylims[, i] <- range(scores[, pc.biplot[2]]/sd.val[pc.biplot[2]], 
                          loadings[, pc.biplot[2]] * sd.val[pc.biplot[2]]/biplot.factor * 
                            1.2)
    }
    for (i in 1:nplot) {
      pca.type1 <- pca.type[i]
      scores <- x[[pca.type1]]$scores
      rownames(scores) <- x[["row.names"]]
      loadings <- x[[pca.type1]]$loads
      rownames(loadings) <- x[["col.names"]]
      var.val <- x[[pca.type1]]$var
      var.prop <- x[[pca.type1]]$var.prop
      pc1 = 100 * round(var.prop[pc.biplot[1]], digits = 2)
      pc2 = 100 * round(var.prop[pc.biplot[2]], digits = 2)
      xlab = paste("PC", pc.biplot[1], ": ", pc1, "%", 
                   sep = "")
      ylab = paste("PC", pc.biplot[2], ": ", pc2, "%", 
                   sep = "")
      sd.val = sqrt(var.val)
      scores.min = min(scores[, pc.biplot])
      scores.max = max(scores[, pc.biplot])
      xlim <- c(min(xlims[1, ]), max(xlims[2, ]))
      ylim <- c(min(ylims[1, ]), max(ylims[2, ]))
      plot(scores[, pc.biplot[1]]/sd.val[pc.biplot[1]], 
           scores[, pc.biplot[2]]/sd.val[pc.biplot[2]], 
           main = paste("Weighting =", pca.type1), xlab = xlab, 
           ylab = ylab, type = "n", xlim = xlim, ylim = ylim, 
           col.axis = "blue", cex.main = 1.5, cex.lab = 1.5, 
           cex.axis = 1.5)
      arrows(0, 0, loadings[, pc.biplot[1]] * sd.val[pc.biplot[1]]/biplot.factor, 
             loadings[, pc.biplot[2]] * sd.val[pc.biplot[2]]/biplot.factor, 
             length = 0.1, lwd = 2, angle = 20, col = "red")
      text(loadings[, pc.biplot[1]] * sd.val[pc.biplot[1]]/biplot.factor * 
             1.2, loadings[, pc.biplot[2]] * sd.val[pc.biplot[2]]/biplot.factor * 
             1.2, labels = rownames(loadings), col = "red", 
           cex = 1.2)
      text(scores[, pc.biplot[1]]/sd.val[pc.biplot[1]], 
           scores[, pc.biplot[2]]/sd.val[pc.biplot[2]], 
           labels = rownames(scores), col = "blue", cex = 0.7)
      abline(0, 0, col = "black")
      abline(v = 0, col = "green")
    }
  }
  if ("glyph" %in% plots) {
    if (length(pc.glyph) < 3) {
      message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
    }
    nplot <- length(x$pca.wt)
    # nplot=2
    if (nplot == 4) {
      par(mfrow = c(2, 2))
    }
    else {
      par(mfrow = c(1,nplot))
    }
    pca.type <- x$pca.wt
    # pca.type <- c("unweighted", "spatiotemporal")
    
    for (i in 1:nplot) {
      if (x$pca.mode == "Smode") {
        glyph.data <- x[[pca.type[i]]]
        loadings <- glyph.data$loads[, pc.glyph]
        plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
             ylab = names(coords)[2], main = paste("Loadings: weight =", 
                                                   pca.type[i]), cex.main = 1.5, cex.lab = 1.5, 
             cex.axis = 1.5, xlim = zoomx, 
             ylim = zoomy)
        lines(river, col = "lightgrey")
        glyph.plot.river(as.matrix(loadings), coords, 
                         r1 = r1, add = TRUE, lwd=lwd)
      }
      if (x$pca.mode == "Tmode") {
        glyph.data <- x[[pca.type[i]]]
        scores <- glyph.data$scores[, pc.glyph]
        plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
             ylab = names(coords)[2], main = paste("Scores: weight =", 
                                                   pca.type[i]), cex.main = 1.5, cex.lab = 1.5, 
             cex.axis = 1.5, xlim = zoomx, 
             ylim = zoomy)
        lines(river, col = "black")
        glyph.plot.river(as.matrix(scores), coords, 
                         r1 = r1, add = TRUE, lwd=lwd)
      }
    }
    message("For plots = glyph: did you check that coordinates for monitoring site ID's are in \nthe same order as data used for stpca()?  \nIf not then glyphs might not correspond to the correct monitoring site.")
  }
  if ("ts" %in% plots) {
    if (x$pca.mode == "Tmode") {
      message("For plots = ts: ts plot should only be used when pca.mode = Smode")
    }
    if (x$pca.mode == "Smode") {
      nplot <- length(pc.ts)
      pca.type <- x$pca.wt
      par(mfcol = c(nplot, length(pca.type)))
      min.y <- numeric(length(pca.type))
      max.y <- numeric(length(pca.type))
      for (i in 1:length(pca.type)) {
        min.y[i] <- min(x[[pca.type[i]]]$scores[, pc.ts])
        max.y[i] <- max(x[[pca.type[i]]]$scores[, pc.ts])
      }
      for (i in 1:length(pca.type)) {
        ts.data <- x[[pca.type[i]]]
        for (j in pc.ts) {
          plot(1:nrow(ts.data$scores), ts.data$scores[, 
                                                      j], type = "l", xaxt = "n", ylim = c(min(min.y), 
                                                                                           max(max.y)), xlab = "Date", ylab = "score", 
               main = paste("weight = ", pca.type[i], ", PC = ", 
                            j, sep = ""))
          axis(side = 1, at = round(quantile(1:nrow(ts.data$scores), 
                                             probs = c(0, 0.33, 0.66, 1))), labels = rownames(ts.data$scores)[round(quantile(1:nrow(ts.data$scores), 
                                                                                                                             probs = c(0, 0.33, 0.66, 1)))])
        }
      }
    }
  }
  if ("meanPM" %in% plots) {
    if (x$pca.mode == "Tmode") {
      message("For plots = meanPM: meanPM plot should only be used when pca.mode = Smode")
    }
    if (x$pca.mode == "Smode") {
      nplot <- length(pc.meanPM)
      pca.type <- x$pca.wt
      par(mfcol = c(nplot, length(pca.type)))
      min.y <- numeric(length(pca.type))
      max.y <- numeric(length(pca.type))
      for (i in 1:length(pca.type)) {
        meanPM.data <- x[[pca.type[i]]]
        min.y[i] <- min((x[["mean.timeseries"]] + (meanPM.factor * 
                                                     meanPM.data$scores[, i])), (x[["mean.timeseries"]] - 
                                                                                   (meanPM.factor * meanPM.data$scores[, i])))
        max.y[i] <- max((x[["mean.timeseries"]] + (meanPM.factor * 
                                                     meanPM.data$scores[, i])), (x[["mean.timeseries"]] - 
                                                                                   (meanPM.factor * meanPM.data$scores[, i])))
      }
      for (i in 1:length(pca.type)) {
        meanPM.data <- x[[pca.type[i]]]
        for (j in pc.meanPM) {
          plot(1:nrow(meanPM.data$scores), x[["mean.timeseries"]], 
               type = "l", xaxt = "n", ylim = c(min(min.y), 
                                                max(max.y)), xlab = "Date", ylab = "score", 
               main = paste("weight = ", pca.type[i], ", PC = ", 
                            j, sep = ""), lwd = 2)
          lines(1:nrow(meanPM.data$scores), x[["mean.timeseries"]])
          points(1:nrow(meanPM.data$scores), x[["mean.timeseries"]] + 
                   (meanPM.factor * meanPM.data$scores[, j]), 
                 pch = "+", col = "red")
          points(1:nrow(meanPM.data$scores), x[["mean.timeseries"]] - 
                   (meanPM.factor * meanPM.data$scores[, j]), 
                 pch = "-", col = "blue")
          axis(side = 1, at = round(quantile(1:nrow(meanPM.data$scores), 
                                             probs = c(0, 0.33, 0.66, 1))), labels = rownames(meanPM.data$scores)[round(quantile(1:nrow(meanPM.data$scores), 
                                                                                                                                 probs = c(0, 0.33, 0.66, 1)))])
        }
      }
    }
  }
  if ("map" %in% plots) {
    nplot <- length(x$pca.wt)
    pca.type <- x$pca.wt
    for (i in 1:nplot) {
      if (x$pca.mode == "Smode") {
        title <- paste("Loadings PC", pc.map, "\n weight =", 
                       pca.type[i])
        map.data <- x[[pca.type[i]]]
        loadings <- map.data$loads[, pc.map]
        n.colors <- length(probs) - 1
        data.vec <- sort(loadings)
        quants <- quantile(data.vec, probs = probs)
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(n.colors)))(n.colors)
        par(mar = c(5, 4, 4, 2) + 0.1)
        graphics::layout(matrix(1:2, ncol = 2), widths = c(0.85, 
                                                           0.15), heights = c(1, 1))
        plot(coords[, 1], coords[, 2], type = "n", xlab = "Easting", 
             ylab = "Northing", main = title, cex.main = 1.5, 
             cex.axis = 1.5, cex.lab = 1.5, xlim = zoomx, ylim = zoomy)
        lines(river)
        points(coords[, 1], coords[, 2], pch = 19, cex = 2, 
               col = getPalette[cut(loadings, breaks = quants)])
        par(mar = c(5.1, 0.4, 4.1, 0.5))
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(100)))(100)
        plot(NA, type = "n", ann = FALSE, xlim = c(0, 
                                                   2.1), ylim = c(0, 101), xaxt = "n", yaxt = "n", 
             bty = "n")
        rect(1, c(0:99), 1.5, c(1:100), col = c(getPalette), 
             border = NA)
        text(x = 0.5, y = ((0:3) * 33), labels = round(quants, 
                                                       2), cex = 1.5)
      }
      if (x$pca.mode == "Tmode") {
        title <- paste("Scores PC", pc.map, "\n weight =", 
                       pca.type[i])
        map.data <- x[[pca.type[i]]]
        scores <- map.data$scores[, pc.map]
        n.colors <- length(probs) - 1
        data.vec <- sort(scores)
        quants <- quantile(data.vec, probs = probs)
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(n.colors)))(n.colors)
        par(mar = c(5, 4, 4, 2) + 0.1)
        graphics::layout(matrix(1:2, ncol = 2), widths = c(0.85, 
                                                           0.15), heights = c(1, 1))
        plot(coords[, 1], coords[, 2], type = "n", xlab = "Easting", 
             ylab = "Northing", main = title, cex.main = 1.5, 
             cex.axis = 1.5, cex.lab = 1.5, xlim = zoomx, ylim = zoomy)
        lines(river)
        points(coords[, 1], coords[, 2], pch = 19, cex = 2, 
               col = getPalette[cut(scores, breaks = quants)])
        par(mar = c(5.1, 0.4, 4.1, 0.5))
        getPalette <- (grDevices::colorRampPalette(sp::bpy.colors(100)))(100)
        plot(NA, type = "n", ann = FALSE, xlim = c(0, 
                                                   2.1), ylim = c(0, 101), xaxt = "n", yaxt = "n", 
             bty = "n")
        rect(1, c(0:99), 1.5, c(1:100), col = c(getPalette), 
             border = NA)
        a <- n.colors
        b <- 100/n.colors
        text(x = 0.5, y = ((0:a) * b), labels = round(quants, 
                                                      2), cex = 1.5)
      }
    }
  }
}


flury.hierarchy <- function(data, covmats, nvec, n.var, cluster, alpha=0.05, mode="FG"){
  # flury's AIC
  if(mode == "FG"){
    g <- cpc::FG
  }else if(mode == "stepwise"){
    g <- cpc::stepwisecpc
  }
  
  varnames <- colnames(covmats[,,1]) # variable names
  no.groups <- dim(covmats)[3] # number of groups
  
  B.cpc = g(covmats <- covmats, nvec = nvec)$B # common eigenvectors under CPC model
  rownames(B.cpc) <- varnames
  colnames(B.cpc) <- paste("PC",c(1:length(varnames)),sep="")
  common_order <- cpc::findcpc(covmats=covmats, B=B.cpc, # common eigenvectors order for cpcq model 
                               plotting=FALSE)$commonvec.order
  
  flury_res <- cpc::flury.test(covmats=covmats, nvec=nvec, B=B.cpc, # flury's hierarchy model res 
                               p = n.var, qmax = n.var - 2, 
                               commonvec.order = common_order)
  flury_res <- cbind(flury_res, rep(0, n.var+2), rep(0, n.var+2), rep(0, n.var+2))
  colnames(flury_res)[7] <- "total.chi.square"
  colnames(flury_res)[8] <- "total.df"
  colnames(flury_res)[9] <- "criterion"
  
  res.equal <- cpc::equal.test(covmats=covmats, nvec=nvec) # equality test
  
  equal <- 0 # if equality model holds, becomes 1
  prop <- 0 # if proportionality model holds, becomes 1
  cpc <- 0 # if cpc model holds, becomes 1
  
  res <- list()
  res[[1]] <- equal ; res[[2]] <- prop ; res[[3]] <- cpc
  res[[4]] <- B.cpc ; res[[5]] <- common_order ; res[[6]] <- flury_res  
  names(res) <- c("equal", "prop", "cpc", "B.cpc", "common_order", "flury_res")
  flury_res[1,7] <- res.equal$chi.square
  flury_res[1,8] <- res.equal$df
  flury_res[1,9] <- qchisq(1-alpha, res.equal$df)
  
  cat("####### Hypothesis Testing #######", end= "\n")
  if(res.equal$chi.square > 
     qchisq(1-alpha, res.equal$df)){
    cat("1 Reject Equality Model", end="\n")
    res.prop <- cpc::prop.test(covmats=covmats, nvec=nvec) # proportionality test
    
    flury_res[2,7] <- res.prop$chi.square
    flury_res[2,8] <- res.prop$df
    flury_res[2,9] <- qchisq(1-alpha, res.prop$df)
    
    if(res.prop$chi.square >
       qchisq(1-alpha, res.prop$df)){
      cat("2 Reject Proportionality Model", end="\n")
      res.cpc <- cpc::cpc.test(covmats=covmats, nvec=nvec) # cpc test
      
      flury_res[3,7] <- res.cpc$chi.square
      flury_res[3,8] <- res.cpc$df
      flury_res[3,9] <- qchisq(1-alpha, res.cpc$df)
      
      if(res.cpc$chi.square >
         qchisq(1-alpha, res.cpc$df)){
        cat("3 Reject CPC Model", end="\n")
        
        q.common <- 0
        if(n.var > 2){
          for(q in (n.var-2):1){
            res.cpcq <- cpc::cpcq.test(covmats = covmats, nvec = nvec, 
                                       B = B.cpc[,common_order], q = q) # cpcq test
            
            flury_res[(n.var+2-q),7] <- res.cpcq$chi.square
            flury_res[(n.var+2-q),8] <- res.cpcq$df
            flury_res[(n.var+2-q),9] <- qchisq(1-alpha, res.cpcq$df)
            
            if(res.cpcq$chi.square > qchisq(1-alpha, res.cpcq$df)){
              cat(paste(2+n.var-q, " Reject CPC(",q,") Model", sep=""), end="\n")
              q.common <- 0
            }else{
              cat(paste(2+n.var-q, " Accept CPC(",q,") Model", sep=""), end="\n")
              q.common <- q # determine cpc(q) model 
              res[[7]] <- q.common ; res[[8]] = res.cpcq$covmats.cpcq
              res[[9]] <- B.partial(covmats=covmats, nvec=nvec, B=B.cpc,
                                    commonvec.order = common_order, q=q.common)
              
              score <- matrix(0, nrow = sum(nvec), ncol = n.var)
              R.F.partial <- array(0, c(n.var, n.var, no.groups), dimnames = list(c(),c(),paste("Group",1:no.groups, sep = "")))
              evptv <- array(0, c(2, n.var, no.groups), 
                             dimnames = list(c("Eigenvalues", "Percentage of total variation"),
                                             paste("PC",c(1:n.var),sep=""),paste("Group",1:no.groups, sep = "")))
              for(j in 1:no.groups){
                Fj <- t(res[[9]][,,j])%*%covmats[,,j]%*%res[[9]][,,j] # covariance of score 
                lambda.inv.sq.rt <- diag(1/sqrt(diag(Fj)))
                Rj <- lambda.inv.sq.rt%*%Fj%*%lambda.inv.sq.rt # correlation of score
                
                R.F.partial[,,j] <- Fj
                R.F.partial[,,j][lower.tri(R.F.partial[,,j])] <- Rj[lower.tri(Rj)]
                
                score[cluster==j,] <- data[cluster==j,]%*%res[[9]][,,j]
                
                evptv[,,j][1,] <- diag(Fj)
                evptv[,,j][2,] <- diag(Fj)/sum(diag(Fj))*100
              }
              
              score <- cbind(score, cluster)
              colnames(score) <- c(paste("PC",c(1:n.var),sep=""), "cluster")
              res[[10]] <- R.F.partial
              res[[11]] <- score
              res[[12]] <- evptv
              
              names(res)[c(7:12)] <- c("q.common", "covmats.cpcq", "B.partial", "R.F.partial", "scores","EVPTV")
              dimnames(res[[8]]) <- list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
              dimnames(res[[9]]) <- list(varnames,paste("PC",c(1:n.var),sep=""),paste("Group",1:no.groups, sep = ""))
              
              break
            }
          }
        }
        
      }else{
        cat("3 Accept CPC Model", end="\n")
        cpc <- 1
        res[[3]] <- cpc
        res[[7]] <- res.cpc$covmats.cpc
        names(res)[7] <- "covmats.cpc"
        dimnames(res[[7]]) <- list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
        
        score <- data%*%B.cpc
        R.F.cpc <- array(0, c(n.var, n.var, no.groups), dimnames = list(c(),c(),paste("Group",1:no.groups, sep = "")))
        evptv <- array(0, c(2, n.var, no.groups), 
                       dimnames = list(c("Eigenvalues", "Percentage of total variation"),
                                       paste("PC",c(1:n.var),sep=""),paste("Group",1:no.groups, sep = "")))
        for(j in 1:no.groups){
          Fj <- t(B.cpc)%*%covmats[,,j]%*%B.cpc # covariance of score 
          lambda.inv.sq.rt <- diag(1/sqrt(diag(Fj)))
          Rj <- lambda.inv.sq.rt%*%Fj%*%lambda.inv.sq.rt # correlation of score
          
          R.F.cpc[,,j] <- Fj
          R.F.cpc[,,j][lower.tri(R.F.cpc[,,j])] <- Rj[lower.tri(Rj)]
          
          evptv[,,j][1,] <- diag(Fj)
          evptv[,,j][2,] <- diag(Fj)/sum(diag(Fj))*100
        }
        score <- cbind(score, cluster)
        colnames(score) <- c(paste("PC",c(1:n.var),sep=""), "cluster")
        res[[8]] <- R.F.cpc
        res[[9]] <- score
        res[[10]] <- evptv
        
        names(res)[c(7:10)] <- c("covmats.cpc", "R.F.cpc", "scores", "EVPTV")
      }
    }else{
      cat("2 Accept Proportionality Model", end="\n")
      prop <- 1
      res[[2]] <- prop
      res[[7]] <- res.prop$covmats.prop
      names(res)[7] <- "covmats.prop"
      dimnames(res[[7]]) <- list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
    }
  }else{
    cat("1 Accept Equality Model", end="\n")
    equal <- 1
    res[[1]] <- equal
    res[[7]] <- res.equal$covmats.equal
    names(res)[7] <- "covmats.equal"
    dimnames(res[[7]]) <- list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
  }
  
  cat("\n")
  print(flury_res)
  cat("\n")
  
  cat("####### Model selection #######", end="\n")
  if(equal == 1){
    cat(paste("1 Hypothesis testing : ", "Equality Model", sep=""), end="\n")
  }else if(prop == 1){
    cat(paste("1 Hypothesis testing : ", "Proportionality Model", sep=""), end="\n")
  }else if(cpc == 1){
    cat(paste("1 Hypothesis testing : ", "CPC Model", sep=""), end="\n")
  }else if(q.common > 0){
    cat(paste("1 Hypothesis testing : ", "CPC(",q.common,") Model", sep=""), end="\n")
  }else if(q.common==0){
    cat(paste("1 Hypothesis testing : ", "Heterogeneity Model", sep=""), end="\n")
  }
  
  cat(paste("2 AIC criterion :", as.character(flury_res$Model[order(flury_res$AIC)[1]]), "Model"), end="\n")
  cat(paste("3 Partial chi square :", as.character(flury_res$Model[order(abs(1-flury_res$Chi2.div.df))[1]]), "Model","\n"), end="\n")
  
  return(res)
}

