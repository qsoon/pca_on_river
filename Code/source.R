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


# from "createWeight_SC" (line 21) to "compute_flow_vol" (line 2334) was defined by Seoncheol Park 
# (https://github.com/SeoncheolPark/paper-StreamflowLifting)


# existing createWeightS function takes file path of ssn object as an input 
# createWeightS_SC takes ssn object itself as an input
createWeightS_SC <- function (ssndata, afvcol="addfunccol", adjacency) {
  # ssndata : should be SpatialStreamNetwork class
  # afvcol: character (which contains weight information)
  # adjacency: adjacency information (can be obtained from get_adjacency_stream)
  
  adjacency_old <- adjacency
  
  # x@obspoints@SSNPoints[[1]]@point.data$addfunccol
  
  if(class(ssndata)[1]!="SpatialStreamNetwork"){
    stop("ssndata must be a SpatialStreamNetwork class")
  }else if( !(afvcol %in% names(ssndata@obspoints@SSNPoints[[1]]@point.data) )){
    stop("afvcol does not exist")
  }else if(!(all(c("upDist", "rid", "pid")%in%colnames(ssndata@obspoints@SSNPoints[[1]]@point.data)))){
    stop("you must check whether (upDist, rid, pid) are in your SpatialStreamNetwork object" )
  }
  
  #for(i in 1:length(unique(ssndata@network.line.coords$NetworkID))){
  #  
  #}
  n.data <- nrow(ssndata@obspoints@SSNPoints[[1]]@point.data)
  
  binary_rid <- adjacency_old$rid_bid[as.numeric(ssndata@obspoints@SSNPoints[[1]]@point.data$rid),2]
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
  return(weight.connected)
}


########################################
##CODES related to the SSN and smnet package
########################################
initint2_stream <- function(example_network, adjacency, linear=TRUE, MyRivernetwork=NULL){
  #example_network: should be SpatialStreamNetwork class from SSN package
  if(linear==FALSE & is.null(MyRivernetwork)==TRUE){
    stop("For nonlinear stream network, riverdist object is needed.")
  }
  # example_network: "SpatialStreamNetwork" obj
  # adjacency: list type
  # linear: indicate whether given stream network is linear or not 
  setClass("ILines", slots=list(Lines="list", ID="character", weight="numeric", range="numeric"))
  
  I_lines <- list()
  I <- matrix(rep(0, length(unique(example_network@obspoints@SSNPoints[[1]]@point.data$pid))), nrow=1)
  colnames(I) <- unique(example_network@obspoints@SSNPoints[[1]]@point.data$pid)
  for(j in 1:length(unique(example_network@obspoints@SSNPoints[[1]]@point.data$pid))){
    I_lines[[j]] <- list()
  }
  
  network_SegmentID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID # Segment ID of the observation station 
  for(i in 1:length(unique(example_network@network.line.coords$SegmentID))){
    # function matching each segment to the nearest point according to stream network structure
    n_inthe_Segment <- sum(i==network_SegmentID) # to check whether observation exists or not on the ith segment
    if(n_inthe_Segment==0){
      # case when there is no point on the segment
      # find a stream segement where there exists obs.point among the nearest streams along downstream direction 
      ind_imsi <- i
      # pt_downstream_dist <- 0
      while(TRUE){
        # find downstream for the stream
        mouth_index_imsi <- sum(as.matrix(adjacency$adjacency)[ind_imsi,]==1) # number of tributaries downstream.
        if(mouth_index_imsi==0){
          # means that we can not find obs.points along the downstream direction
          break
        }else if(mouth_index_imsi==1){
          # assume that the number of tributaries downstream != 2
          mouth_seg_imsi <- which(as.matrix(adjacency$adjacency)[ind_imsi,]==1)
          if(sum(mouth_seg_imsi==network_SegmentID)==0){
            # have to move on to the next segment
            ind_imsi <- mouth_seg_imsi
            # pt_downstream_dist
          }else{
            # case when there exist a few obs.points for the segment
            # find the nearest index on the segment
            if(sum(mouth_seg_imsi==network_SegmentID)==1){
              # when there exists only one obs.point on the segment
              ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi)
              I[ind] <- I[ind] + example_network@data$Length[i]
              I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range=c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
              #I_lines[[ind]] <- c(I_lines[[ind]], example_network@lines[[i]])
              I_lines[[ind]] <- c(I_lines[[ind]], I_lines_new)
            }else{
              # when there exist several obs.points on the segment
              ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi) # obs coordinates on the segment
              # we have to identify the nearest place from the segment
              # example_network@lines : lines are ordered upstream -> downstream
              pt_upstream <- example_network@lines[[mouth_seg_imsi]]@Lines[[1]]@coords[1,]
              nearest_pt_index <- ind[which.min(as.matrix(dist(rbind(pt_upstream,example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]) ))[-1,1])]
              I[nearest_pt_index] <- I[nearest_pt_index] + example_network@data$Length[i]
              I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range=c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
              #I_lines[[nearest_pt_index]] <- c(I_lines[[nearest_pt_index]], example_network@lines[[i]])
              I_lines[[nearest_pt_index]] <- c(I_lines[[nearest_pt_index]], I_lines_new)
            }
            break
          }
        }
      }
    }else if(n_inthe_Segment==1){
      # when there is only one point on the segment
      ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == i)
      I[ind] <- I[ind] + example_network@data$Length[i]
      I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range= c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
      #I_lines[[ind]] <- c(I_lines[[ind]], example_network@lines[[i]])
      I_lines[[ind]] <- c(I_lines[[ind]], I_lines_new)
    }else{
      # when there are several points on the segment
      ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == i)
      
      # example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]
      node_order_index <- order(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[ind], decreasing = TRUE)
      
      if(linear==TRUE){
        # if linear, use upDist
        ordered_coords <- example_network@obspoints@SSNPoints[[1]]@point.coords[ind[node_order_index],]
        
        middle_pts <- c() # vector for centroids
        for(j in 1: (length(ind)-1)){
          # calculate centroids
          middle_pts_imsi <- c((ordered_coords[j,1]+ordered_coords[(j+1),1])/2, (ordered_coords[j,2]+ordered_coords[(j+1),2])/2)
          middle_pts <- rbind(middle_pts, middle_pts_imsi)
        }
        #end_pts <- example_network@lines[[i]]@Lines[[1]]@coords[order(example_network@lines[[i]]@Lines[[1]]@coords[,1]),]
        end_pts <- example_network@lines[[i]]@Lines[[1]]@coords[,]
        middle_pts <- rbind(end_pts[1,], middle_pts, end_pts[2,])
        
        # use end points and mid point to calculate lengths
        for(k in 1:length(ind)){
          I[ind[node_order_index[k]]] <- I[ind[node_order_index[k]]] + sqrt(sum((middle_pts[k,]-middle_pts[(k+1),])^2))
          example_network_lines_imsi <- example_network@lines[[i]]
          example_network_lines_imsi@ID <- paste(example_network_lines_imsi@ID, "-", k, sep="")
          example_network_lines_imsi@Lines[[1]]@coords <- rbind(middle_pts[k,], middle_pts[(k+1),])
          I_lines_new <- new("ILines", Lines=example_network_lines_imsi@Lines, ID=example_network_lines_imsi@ID, weight=1, range=c(1, nrow(example_network@lines[[i]]@Lines[[1]]@coords)))
          #I_lines[[ind[node_order_index[k]]]] <- c(I_lines[[ind[node_order_index[k]]]], example_network_lines_imsi)
          I_lines[[ind[node_order_index[k]]]] <- c(I_lines[[ind[node_order_index[k]]]], I_lines_new)
        }
      }else if(linear==FALSE){
        # when given stream set is not linear
        # assume that each partition has the same distance within each segment and divide 
        
        
        #example_network@obspoints@SSNPoints[[1]]@point.data$rid[ind]
        #example_network@obspoints@SSNPoints[[1]]@point.data$segid[ind]
        
        dist_combined <- rep(0, nrow(example_network@lines[[i]]@Lines[[1]]@coords))
        for(j in 1:nrow(example_network@lines[[i]]@Lines[[1]]@coords)){
          dist_imsi_vec <- rep(0, length(ind))
          for(k in 1:length(ind)){
            dist_imsi_vec[k] <- riverdistance(startseg = i, endseg = i, startvert = min(example_network@obspoints@SSNPoints[[1]]@point.data$segid[ind[k]], j), endvert = max(example_network@obspoints@SSNPoints[[1]]@point.data$segid[ind[k]], j), rivers= MyRivernetwork)
          }
          dist_combined[j] <- ind[which.min(dist_imsi_vec)]
        }
        
        
        
        for(l in 1:length(ind)){
          range_imsi <- range(which(dist_combined==ind[l]))
          I[ind[l]] <- I[ind[l]] + riverdistance(startseg = i, endseg = i, startvert = range_imsi[1], endvert = range_imsi[2], rivers= MyRivernetwork)
          if(range_imsi[1]!=1){
            I[ind[l]] <- I[ind[l]] + riverdistance(startseg = i, endseg = i, startvert = (range_imsi[1]-1), endvert = range_imsi[1], rivers= MyRivernetwork)/2
          }
          if(range_imsi[2]!=nrow(example_network@lines[[i]]@Lines[[1]]@coords)){
            I[ind[l]] <- I[ind[l]] + riverdistance(startseg = i, endseg = i, startvert = range_imsi[2], endvert = (range_imsi[2]+1), rivers= MyRivernetwork)/2
          }
          I_lines_new <- new("ILines", Lines=example_network@lines[[i]]@Lines, ID=example_network@lines[[i]]@ID, weight=1, range=range(which(dist_combined==ind[l])))
          #I_lines[[ind[node_order_index[k]]]] <- c(I_lines[[ind[node_order_index[k]]]], example_network_lines_imsi)
          I_lines[[ind[l]]] <- c(I_lines[[ind[l]]], I_lines_new)
        }
      }
      
      
    }
  }
  return(list(I=I, I_lines=I_lines))
  #I: area of node i
  #I_lines: shows that each node i has which segement
}

getnbrs_stream3 <- function(remove, pointsin, example_network, I_lines, adjacency, upperExtra=FALSE){
  # prototype 3
  r = which(remove==pointsin)
  
  nbrs_upstream <- c()
  nbrs_downstream <- c()
  
  network_SegmentID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[r]
  n_inthe_Segment <- sum(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID==network_SegmentID)
  r_length <- nchar(adjacency$rid_bid[network_SegmentID,2])
  
  #index <- setdiff()
  
  network_indexID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[-r]
  
  #samestream_index: have to check whether upstream or downstream after looking at upDist
  samestream_index <- setdiff(which(adjacency$rid_bid[network_indexID,2]==adjacency$rid_bid[network_SegmentID,2]), remove)
  #upstream_index: upstream
  upstream_index <- setdiff(str_which( substr(adjacency$rid_bid[network_indexID,2], start=1, stop=r_length), adjacency$rid_bid[network_SegmentID,2]), remove)
  upstream_index <- setdiff(upstream_index, samestream_index)
  
  if(length(samestream_index)!=0){
    samestream_index_imsi <- match(samestream_index, pointsin)
    for(i in 1:length(samestream_index_imsi)){
      if(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[r] < example_network@obspoints@SSNPoints[[1]]@point.data$upDist[samestream_index_imsi[i]]){
        # add neighbors to upstream
        nbrs_upstream <- c(nbrs_upstream, samestream_index[i])
      }else if(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[r] > example_network@obspoints@SSNPoints[[1]]@point.data$upDist[samestream_index_imsi[i]]){
        # add neighbors to downstream
        nbrs_downstream <- c(nbrs_downstream, samestream_index[i])
      }
    }
  }
  
  if(length(nbrs_upstream)>=2){
    
  }
  if(length(nbrs_downstream)>=2){
    
  }
  
  if(length(nbrs_upstream)==0 & length(upstream_index)!=0){
    # if there is no upstream neighbor, find upstream neighbor
    cand_length <- nchar(adjacency$rid_bid[network_indexID,2])[-remove]
    nbrs_cand <- pointsin[nchar(adjacency$rid_bid[network_indexID,2]) > nchar(adjacency$rid_bid[network_SegmentID,2])]
    nbrs_cand <- setdiff(nbrs_cand,remove)
    nbrs_cand <- intersect(nbrs_cand, upstream_index)
    while(TRUE){
      diff_vec <- nchar(adjacency$rid_bid[network_indexID[nbrs_cand],2])- nchar(adjacency$rid_bid[network_SegmentID,2])
      upstream_candidate <- nbrs_cand[which(diff_vec==min(diff_vec))]
      if(length(upstream_candidate)>1){
        upstream_candidate <- upstream_candidate[which.min(example_network@obspoints@SSNPoints[[1]]@point.data$upDist[match(upstream_candidate, pointsin)])]
      }
      
      nbrs_upstream <- c(nbrs_upstream, upstream_candidate)
      nchar_upstream_cand <- nchar(adjacency$rid_bid[network_indexID[upstream_candidate],2])
      up_upstream <- str_which( substr(adjacency$rid_bid[network_indexID[nbrs_cand],2], start=1, stop=nchar_upstream_cand ), adjacency$rid_bid[network_indexID[upstream_candidate],2])
      if(length(up_upstream)!=0){
        nbrs_cand <- nbrs_cand[-up_upstream]
      }
      if(length(nbrs_cand)==0){
        break
      }
    }
  }
  if(length(nbrs_downstream)==0){
    # finding downstream neighbour
    # when there is no point on the segment
    # find a stream segement where there exists obs.point among the nearest streams along downstream direction 
    ind_imsi <- network_SegmentID 
    #pt_downstream_dist <- 0
    while(TRUE){
      # find downstream
      mouth_index_imsi <- sum(as.matrix(adjacency$adjacency)[ind_imsi,]==1) # number of tributaries downstream
      if(mouth_index_imsi==0){
        break
      }else if(mouth_index_imsi==1){
        # assume that the number of tributaries downstream != 2
        mouth_seg_imsi <- which(as.matrix(adjacency$adjacency)[ind_imsi,]==1)
        if(sum(mouth_seg_imsi==network_indexID)==0){
          # have to move on to the next segment
          ind_imsi <- mouth_seg_imsi
          #pt_downstream_dist
        }else{
          # when there are several monitoring sites on the segment
          # find the nearest index on the segment
          if(sum(mouth_seg_imsi==network_indexID)==1){
            # when there is only one monitoring site on the segment
            ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi)
            nbrs_downstream <- c(nbrs_downstream, pointsin[ind])
          }else{
            # when there are multiple monitoring site on the segment
            ind <- which(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID == mouth_seg_imsi) #해당 segment들에 들어있는 obs들의 좌표들
            # have to identify the place nearest from the segment
            # example_network@lines : lines are ordered upstream -> downstream
            pt_upstream <- example_network@lines[[mouth_seg_imsi]]@Lines[[1]]@coords[1,]
            nearest_pt_index <- ind[which.min(as.matrix(dist(rbind(pt_upstream,example_network@obspoints@SSNPoints[[1]]@point.coords[ind,]) ))[-1,1])]
            nbrs_downstream <- c(nbrs_downstream, pointsin[nearest_pt_index])
          }
          break
        }
      }
    }
  }
  return(list(nbrs=sort(c(nbrs_upstream, nbrs_downstream)), nbrs_upstream=sort(nbrs_upstream), nbrs_downstream=sort(nbrs_downstream)))
}


# to be fixed
getnbrs_stream <- function(remove, pointsin, example_network, I_lines, adjacency, upperExtra=FALSE){
  r = which(remove==pointsin)

  network_SegmentID <- example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[r]
  n_inthe_Segment <- sum(example_network@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID==network_SegmentID)
  
  coords_stream <- c()
  for(ii in 1:length(I_lines[[remove]])){
    # ii: all stream indices in r
    range_up <- I_lines[[remove]][[ii]]@range[1]
    range_down <- I_lines[[remove]][[ii]]@range[2]
    if(range_up!=1){
      range_up <- range_up - 1
    }
    if(range_down!=nrow(I_lines[[remove]][[ii]]@Lines[[1]]@coords)){
      range_down <- range_down + 1
    }
    coords_stream <- rbind(coords_stream, I_lines[[remove]][[ii]]@Lines[[1]]@coords[range_up,], I_lines[[remove]][[ii]]@Lines[[1]]@coords[range_down,])
  }
  coords_stream <- coords_stream[!duplicated(coords_stream),]
  
  n_test <- setdiff(example_network@obspoints@SSNPoints[[1]]@point.data$pid, remove)
  # remains after checking r
  
  stream_candidate <- c()
  for(jj in 1:length(n_test)){
    # jj: remains after removing r among nodes 
    lengths_segment <- length(I_lines[[n_test[jj]]])
    match_result <- c()
    for(kk in 1:length(I_lines[[n_test[[jj]]]])){
      # kk: all stream indices on the node jj
      match_result_imsi <- apply(I_lines[[n_test[[jj]]]][[kk]]@Lines[[1]]@coords, 1, function(a) sum(which(a%in%coords_stream)))
      match_result <- c(match_result, match_result_imsi)
    }
    if(sum(match_result)!=0){
      if(upperExtra==TRUE){
        # upperExtra==TRUE: point on the tributary of Y network : neighbor
        stream_candidate <- c(stream_candidate, n_test[jj])
      }else{
        # upperExtra==FALSE: point on the tributary of Y network : NOT neighbor
        segment_index <- as.numeric(example_network@obspoints@SSNPoints[[1]]@point.data$rid[[n_test[jj]]])
        # segment_index <- as.numeric(example_network@obspoints@SSNPoints[[1]]@point.data$rid[jj])
        segment_index_r <- as.numeric(example_network@obspoints@SSNPoints[[1]]@point.data$rid[r])
        if(segment_index==segment_index_r){
          # if two points exist on the same segment -> neighbors
          stream_candidate <- c(stream_candidate, n_test[jj])
        }else{
          distUpstream_index <- example_network@network.line.coords$DistanceUpstream[n_test[jj]]
          #distUpstream_index <- example_network@network.line.coords$DistanceUpstream[segment_index]
          distUpstream_index_r <- example_network@network.line.coords$DistanceUpstream[segment_index_r]
          
          bid_index <- adjacency$rid_bid[segment_index,2]
          bid_index_r <- adjacency$rid_bid[segment_index_r,2]
          if(distUpstream_index_r>distUpstream_index){
            # point to be removed exists upstream than neighbors
            if(str_detect(bid_index_r, bid_index) | str_detect(bid_index, bid_index_r)  ){
              adj_imsi <- as.matrix(adjacency$adjacency)
              if(sum(as.matrix(adjacency$adjacency)[,segment_index_r])>=2 | sum(adj_imsi[,as.numeric(as.character(network_SegmentID))])==0){
                # if upstream of removed point >=2 -> add neighbors to downstream
                # or removed point exist the most upstream -> add neighbors
                stream_candidate <- c(stream_candidate, n_test[jj])
              }
            }
          }else{
            # point to be removed exists downstream than neighbors
            if(str_detect(bid_index, bid_index_r) | str_detect(bid_index_r, bid_index)){
              stream_candidate <- c(stream_candidate, n_test[jj])
            }
          }
        }
      }
    }
  }

  rid_stream_candidate <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[match(stream_candidate,pointsin)]
  rid_stream_r <- example_network@obspoints@SSNPoints[[1]]@point.data$rid[r]
  if(!all(rid_stream_candidate!=rid_stream_r)){
    stream_candidate <- stream_candidate[which(rid_stream_candidate==rid_stream_r)]
  }
  
  return(stream_candidate)

}


# construct streamPred function
# basic rule: inverse distance weight on segment, else : regression along the stream
streamPred <- function(pointsin, X=example_network, coeff=sims@obspoints@SSNPoints[[1]]@point.data$Sim_Values, nbrs, remove, adj=adjacency, method=c("IDW","LC"), intercept=FALSE, neighbours=NULL, adaptive=FALSE){
  #prediction weight를 반환하는 함수다
  #method: IDW: inverse distance weight, LC: local constant
  
  #X=example_network; coeff=sims@obspoints@SSNPoints[[1]]@point.data$Sim_Values; adj=adjacency
  
  Xcoords <- X@obspoints@SSNPoints[[1]]@point.coords
  index <- match(nbrs, pointsin)
  r <- which(pointsin == remove)
  Xneighbours <- Xcoords[index,]
  Xremove <- Xcoords[r,]
  #Xneighbours <- Xcoords[nbrs,]
  #Xremove <- Xcoords[remove,]
  
  if (intercept) {
    #Xneighbours <- cbind(1, Xneighbours)
    #Xremove <- as.row(c(1, Xremove))
  }
  
  
  if(sum(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)!=0){
    updist_remove <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
    reach_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
    #updist_remove <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==r)]
    #reach_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[which(X@obspoints@SSNPoints[[1]]@point.data$pid==r)]
    addfunccol_remove <- X@obspoints@SSNPoints[[1]]@point.data$shreve[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
  }else{
    updist_remove <- X@predpoints@SSNPoints[[1]]@point.data$upDist[length(X@predpoints@SSNPoints[[1]]@point.data$upDist)]
    reach_remove <- X@predpoints@SSNPoints[[1]]@point.data$rid[length(X@predpoints@SSNPoints[[1]]@point.data$rid)]
    addfunccol_remove <- X@predpoints@SSNPoints[[1]]@point.data$shreve[length(X@predpoints@SSNPoints[[1]]@point.data$addfunccol)]
  }
  rid_remove <- reach_remove
  rid_nbrs <- X@obspoints@SSNPoints[[1]]@point.data$rid[index]
  
  #length_vec <- rep(0, length(nbrs)) # length vec의 inverse weight?
  weight_vec <- rep(0, length(nbrs)) # calculate weight from the beginning -> normalize
  
  if(length(nbrs)>=2){
    
    for(i in 1:length(nbrs)){
      nb_index <- nbrs[i]
      #nb_index <- index[i]
      if(method=="IDW" | sum(rid_nbrs==rid_remove)==2){
        # have to calculate distances between each nbr point and removeed point
        
        updist_nb_index <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==nb_index)]
        #updist_remove <- X@obspoints@SSNPoints[[1]]@point.data$upDist[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
        
        weight_vec[i] <- 1/abs(updist_nb_index-updist_remove)
        
      }else if(method=="LC"){
        # weight based on segment?
        
        addfunccol_nb_index <- X@obspoints@SSNPoints[[1]]@point.data$shreve[which(X@obspoints@SSNPoints[[1]]@point.data$pid==nb_index)]
        #addfunccol_remove <- X@obspoints@SSNPoints[[1]]@point.data$addfunccol[which(X@obspoints@SSNPoints[[1]]@point.data$pid==remove)]
        
        weight_vec[i] <- min(addfunccol_nb_index, addfunccol_remove)/max(addfunccol_nb_index, addfunccol_remove)
      }
    }
    #normalizing weight vector
    weights <- weight_vec/sum(weight_vec)
  }else{
    # length ==1 -> don't need to fix
    mm <- 0
    bhat <- 1
    weights <- 1
    pred <- coeff[nbrs]
  }
  
  pred = sum(weights*coeff[nbrs])
  coeff = coeff[remove]
  return(list(weights = weights, pred = pred, coeff = coeff))
}



PointsUpdate_stream <- function (X=example_network, coeff, nbrs, index=NULL, remove, pointsin, weights, lengths=I, lengths_lines=I_lines, method="LC", stream=TRUE, linear=TRUE) 
{
  #X=example_network; lengths=I; lengths_lines=I_lines; method="LC"; index=match(nbrs, pointsin)
  r <- which(pointsin == remove)
  N <- length(pointsin)
  weights_imsi <- (1/weights)/sum(1/weights)
  if(stream==FALSE){
    # if not stream network
    if ((r >= 2) & (r <= (N - 1))) {
      lengths[index] <- as.row(lengths[index])
      weights <- as.row(weights)
      lengths[index] <- lengths[index] + lengths[r] * weights
    }else{
      if (r == 1) {
        lengths[2] <- lengths[2] + lengths[1]
      }
      if (r == N) {
        lengths[N - 1] <- lengths[N - 1] + lengths[N]
      }
    }
    alpha <- matrix(0, 1, length(nbrs))
    if (length(nbrs) >= 2) {
      alpha <- lengths[r] * lengths[index]/(sum(lengths[index]^2))
      coeff[pointsin[index]] <- coeff[pointsin[index]] + alpha * coeff[remove]
    }
    else {
      # q <- which(pointsin == nbrs)
      q <- match(nbrs, pointsin)
      alpha <- lengths[r]/lengths[q]
      coeff[pointsin[q]] <- coeff[pointsin[q]] + alpha * coeff[remove]
    }
    return(list(coeff = coeff, lengths = lengths, r = r, N = N, 
                weights = weights, alpha = alpha))
  }else{
    # if point to be removed has one stream
    # else if(lengths_lines[[remove]][[1]]@Lines==1)
    # stream==TRUE
    rid_nbrs <- X@obspoints@SSNPoints[[1]]@point.data$rid[index]
    rid_remove <- X@obspoints@SSNPoints[[1]]@point.data$rid[r]
    
    # I update
    # I_new <- I
    # I_new[nbrs] <- I_new[nbrs] + weights*I[remove]
    lengths[index] <- as.row(lengths[index])
    weights <- as.row(weights)
    lengths[index] <- lengths[index] + lengths[r] * weights
    
    if(linear==TRUE){
      remove_upper <- lengths_lines[[remove]][[1]]@Lines[[1]]@coords[1,]
      remove_lower <- lengths_lines[[remove]][[1]]@Lines[[1]]@coords[2,]
      remove_id <- strsplit(lengths_lines[[remove]][[1]]@ID, "-")[[1]][1]
    }
    
    
    # if stream network
    if(length(nbrs)>=2){
      if(linear==TRUE){
        streamcoord_new <- t(weights_imsi[order(X@network.line.coords$DistanceUpstream[rid_nbrs], decreasing=T)])%*%rbind(remove_upper, remove_lower)
      }
      if(sum(rid_nbrs==rid_remove)==2){
        # removed point and its neighbors are located in the same segment
        if(method=="IDW"){
          
        }else if(method=="LC"){
          # (1) update lengths_lines
          
          # have to classify which neighbors touch upper stream and which neighbors touch lower stream
          for(i in 1:length(nbrs)){
            nbr_index <- nbrs[i] # caution: lengths_lines did not be removed
            for(j in 1:length(lengths_lines[[nbr_index]])){
              if(strsplit(lengths_lines[[nbr_index]][[j]]@ID,"-")[[1]][1] == remove_id){
                if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_upper)){
                  lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,] <- streamcoord_new
                }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_lower)){
                  lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,] <- streamcoord_new
                }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_lower)){
     
                }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_upper)){
                 
                }
              }
            }
          }
          

        }
      }else if(length(unique(c(rid_nbrs, rid_remove)))>=3){
        # point to be removed exists downstream
        if(length(lengths_lines[[remove]])==1){
          # point to be removed has only one segment
          
          for(i in 1:length(nbrs)){
            nbr_index <- nbrs[i]
            org_index <- length(lengths_lines[[nbr_index]])
            # two reaches is located in different stream

            for(j in 1:length(lengths_lines[[remove]])){
              lengths_lines[[nbr_index]][[org_index+1]] <- lengths_lines[[remove]][[j]]
              lengths_lines[[nbr_index]][[org_index+1]]@weight <- weights[i]
            }
          }
          
        }
      }else{
        # reach index of neighboring points and reach index of remove point totally different
      }
      
    }else{
      if(length(nbrs)==1){
        # special trtment
        # divide two cases: reach index equal (on the same stream segment), not equal
        nbr_index <- nbrs
        # don't need to manage streamcoord_new 
        if(rid_nbrs==rid_remove){
          # two reaches are located in the same stream
          for(j in 1:length(lengths_lines[[nbr_index]])){
            if(linear==TRUE){
              if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_upper)){
                lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,] <- remove_lower
              }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_lower)){
                lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,] <- remove_upper
              }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[2,]==remove_lower)){
  
              }else if(all(lengths_lines[[nbr_index]][[j]]@Lines[[1]]@coords[1,]==remove_upper)){

              }
            }else{
              if(lengths_lines[[nbr_index]][[j]]@ID==I_lines[[remove]][[1]]@ID){
                lengths_lines[[nbr_index]][[j]]@range <- range(c(lengths_lines[[nbr_index]][[j]]@range, I_lines[[remove]][[1]]@range))
              }
            }
          }
        }else{
          org_index <- length(lengths_lines[[nbr_index]])
          # two reaches are located in different streams
        
          for(j in 1:length(lengths_lines[[remove]])){
            lengths_lines[[nbr_index]][[(org_index+j)]] <- lengths_lines[[remove]][[j]]
          }
        }
      }
    }

    alpha <- matrix(0, 1, length(nbrs))
    if (length(nbrs) >= 2) {
      alpha <- lengths[r] * lengths[index]/(sum(lengths[index]^2))
      coeff[pointsin[index]] <- coeff[pointsin[index]] + alpha * coeff[remove]
    }
    else {
      q <- which(pointsin == nbrs)
      alpha <- lengths[r]/lengths[q]
      coeff[pointsin[q]] <- coeff[pointsin[q]] + alpha * coeff[remove]
    }
    
    return(list(coeff = coeff, lengths = lengths, lengths_lines = lengths_lines, r = r, N = N, weights = weights, alpha = alpha))
  }
  
}


#lifting forward step
fwtnp_stream <- function(example_network, I, I_lines, adjacency, nkeep=2, linear=FALSE, do.W = TRUE, varonly = FALSE){
  #do.W = TRUE; varonly = FALSE; linear=FALSE; nkeep=2
  
  removelist <- NULL
  lengthsremove <- NULL
  neighbrs <- list()
  gamlist <- list()
  alphalist <- list()
  schemehist <- NULL
  interhist <- NULL
  clolist <- NULL
  W <- v <- NULL
  
  n <- length(I_lines)
  
  if ((do.W == 1) & (varonly == 1)) {
    varonly <- FALSE
  }
  ex <- do.W + varonly
  if (ex == 1) {
    W <- diag(n)
  }
  if (varonly) {
    v <- rep(1, times = n)
  }
  
  
  pointsin <- matrix(1:n, 1, n)
  
  #nkeep=2
  j<- 1
  
  while(TRUE){
    print(paste("Level ", j))
    if(j==1){
      remove <- order(abs(I))[1] # use initial I
    }else{
      remove <- order(abs(lengths))[1]
    }
    remove <- pointsin[remove] # return original index
    removelist[j] <- remove
    print(paste("remove point is ", remove))
    
    #points(example_network@obspoints@SSNPoints[[1]]@point.coords[which(pointsin==remove),1] , example_network@obspoints@SSNPoints[[1]]@point.coords[which(pointsin==remove),2], cex=2, pch=5, col="red" )
    
    #X_dataframe <- as.data.frame(t(X))
    
    #have to make getnbrs function
    nbrs_obj <- getnbrs_stream3(remove, pointsin, example_network, I_lines, adjacency)
    nbrs <- nbrs_obj$nbrs
    nbrs_upstream <- nbrs_obj$nbrs_upstream
    nbrs_downstream <- nbrs_obj$nbrs_downstream
    print(paste("neighbours of ", remove, " are ", nbrs))
    
    # res <- LocalPred(pointsin, X, coeff, nbrs, remove, intercept, neighbours)
    if(j==1){
      coeff=example_network@obspoints@SSNPoints[[1]]@point.data$Sim_Values
      res <- streamPred(pointsin,  X=example_network, coeff=example_network@obspoints@SSNPoints[[1]]@point.data$Sim_Values, nbrs, remove, adj=adjacency, method="LC", intercept=FALSE, neighbours=NULL)
    }else{
      res <- streamPred(pointsin,  X=example_network, coeff=coeff, nbrs, remove, adj=adjacency, method="LC", intercept=FALSE, neighbours=NULL)
    }
    
    if (length(res) == 2) {
      l <- res[[1]]
      clolist[j] <- res[[2]][[1]]
      nbrs <- res[[2]][[2]]
      index <- res[[2]][[3]]
    }else {
      l <- res
    }
    neighbrs[[j]] <- nbrs
    weights <- l[[1]]
    pred <- l[[2]]
    if (length(l) == 3) {
      scheme <- NULL
      int <- NULL
      details <- NULL
    }else {
      scheme <- l[[5]]
      int <- l[[4]]
      details <- l[[6]]
    }
    coeff[remove] <- coeff[remove] - pred
    
    # construct PointsUpdate function: 
    # l1 <- PointsUpdate(X, coeff, nbrs, index, remove, pointsin, weights, lengths)
    # res <- streamPred(pointsin,  X=example_network, coeff=sims@obspoints@SSNPoints[[1]]@point.data$Sim_Values, nbrs, remove, adj=adjacency, method="LC", intercept=FALSE, neighbours=NULL)
    
    l1 <- PointsUpdate_stream(X=example_network, coeff, nbrs, index=match(nbrs, pointsin), remove, pointsin, weights, lengths=I, lengths_lines=I_lines, method="LC", stream=TRUE, linear=linear)
    
    
    coeff <- l1$coeff
    lengths <- l1$lengths
    r <- l1$r
    weights <- l1$weights
    N <- l1$N
    alpha <- l1$alpha
    if (ex){
      if (varonly) {
        W[r, ] <- W[r, ] - colSums(as.vector(weights) *  matrix(W[index, ], nrow = length(nbrs)))
        W[index, ] <- W[index, ] + matrix(alpha) %*% W[r, ]
        v[remove] <- sum(W[r, ]^2)
        np <- setdiff(1:length(pointsin), r)
        W <- W[np, ]
      }else {
        W[remove, ] <- W[remove, ] - colSums(as.vector(weights) * matrix(W[nbrs, ], nrow = length(nbrs)))
        W[nbrs, ] <- W[nbrs, ] + matrix(alpha) %*% W[remove,  ]
      }
    }
    
    lengthsremove[j] <- lengths[r]
    gamlist[[j]] <- weights
    alphalist[[j]] <- alpha
    schemehist[j] <- scheme
    interhist[j] <- int
    lengths <- lengths[setdiff(1:length(pointsin), r)]
    pointsin <- setdiff(pointsin, remove)
    
    ## example network itself needs to be updated: refer subsetSSN function
    #if(j==1){
    #  #example_network@predpoints <- example_network@obspoints
    #  #example_network@predpoints@SSNPoints[[1]]@network.point.coords <- example_network@predpoints@SSNPoints[[1]]@network.point.coords[remove,]
    #  #example_network@predpoints@SSNPoints[[1]]@point.coords <- example_network@predpoints@SSNPoints[[1]]@point.coords[remove,]
    #  #example_network@predpoints@SSNPoints[[1]]@point.data <- example_network@predpoints@SSNPoints[[1]]@point.data[remove,]
    #}else{
    #  
    #}
    example_network@predpoints@SSNPoints[[1]]@network.point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@network.point.coords, example_network@obspoints@SSNPoints[[1]]@network.point.coords[r,])
    example_network@predpoints@SSNPoints[[1]]@point.coords <- rbind(example_network@predpoints@SSNPoints[[1]]@point.coords, matrix(example_network@obspoints@SSNPoints[[1]]@point.coords[r,], nrow=1, ncol=2))
    example_network@predpoints@SSNPoints[[1]]@point.data <- rbind(example_network@predpoints@SSNPoints[[1]]@point.data, example_network@obspoints@SSNPoints[[1]]@point.data[r,c(1:ncol(example_network@predpoints@SSNPoints[[1]]@point.data))])
    #if(nrow(example_network@predpoints@SSNPoints[[1]]@network.point.coords)!=1){
    #example_network@predpoints@SSNPoints[[1]]@network.point.coords <- example_network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(example_network@predpoints@SSNPoints[[1]]@network.point.coords),]
    #example_network@predpoints@SSNPoints[[1]]@point.coords <- matrix(example_network@predpoints@SSNPoints[[1]]@point.coords[nrow(example_network@predpoints@SSNPoints[[1]]@point.coords),], nrow=ncol(example_network@predpoints@SSNPoints[[1]]@point.coords), ncol=1)
    #example_network@predpoints@SSNPoints[[1]]@point.data <- example_network@predpoints@SSNPoints[[1]]@point.data[nrow(example_network@predpoints@SSNPoints[[1]]@point.data),]
    #}
    
    example_network@obspoints@SSNPoints[[1]]@network.point.coords <- example_network@obspoints@SSNPoints[[1]]@network.point.coords[-r,]
    example_network@obspoints@SSNPoints[[1]]@point.coords <- example_network@obspoints@SSNPoints[[1]]@point.coords[-r,]
    example_network@obspoints@SSNPoints[[1]]@point.data <- example_network@obspoints@SSNPoints[[1]]@point.data[-r,]
    
    #I update
    I <- l1$lengths[-r]
    #I_lines update (it is impossible to subset list typelist)
    I_lines <- l1$lengths_lines
    
    
    ##plotting new data object
    SSN::plot.SpatialStreamNetwork(x=example_network, main=paste("Level", j), VariableName = "Sim_Values", color.palette=rainbow(12), nclasses=10, breaktype = "quantile", brks = NULL, PredPointsID = NULL, add = FALSE, addWithLegend = FALSE, lwdLineCol = "addfunccol", lwdLineEx = 3, lineCol = "black")
    #need to be corrected
    
    #if(length(pointsin)<=nkeep | sum(as.matrix(adjacency$adjacency)[pointsin, pointsin])==0){
    if(length(pointsin)<=nkeep){
      break
    }else{
      #next level
      
      j = j+1
    }
  }
  if (varonly) {
    v[pointsin] <- rowSums(W^2)
    W <- NULL
  }
  return(list(network=example_network, I=I, I_lines=I_lines, adjacency=adjacency, coeff=coeff, lengths=lengths, lengthsremove = lengthsremove, pointsin = pointsin, removelist = removelist, 
              neighbrs = neighbrs, schemehist = schemehist, interhist = interhist, 
              clolist = clolist, gamlist = gamlist, alphalist = alphalist, 
              W = W, v=v))
}

#lifting backward step


UndoPointsUpdate_stream <- function (X, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, I_lines, stream=TRUE) 
{
  #I_lines=X$I_lines; X=X$network;  stream=TRUE
  
  
  # have to construct a function reversing stream network
  gamweights_imsi <- (1/gamweights)/sum(1/gamweights)
  if(stream==TRUE){
    # X$network@obspoints needs to be reconstructed
    #locID(=number of itself), upDist, pid(=number of itself), rid, ratio, shreve, addfunccol, X, X2, Sim_Values
    #point.data_reconstruct <- c()
    
    rid_remove <- X@predpoints@SSNPoints[[1]]@point.data$rid[length(X@predpoints@SSNPoints[[1]]@point.data$rid)]
    rid_nbrs <- X@obspoints@SSNPoints[[1]]@point.data$rid[index]
    
    remove_upper <- I_lines[[remove]][[1]]@Lines[[1]]@coords[1,]
    remove_lower <- I_lines[[remove]][[1]]@Lines[[1]]@coords[2,]
    remove_id <- strsplit(I_lines[[remove]][[1]]@ID, "-")[[1]][1]
    
    remove_coords <- X@predpoints@SSNPoints[[1]]@point.coords[nrow(X@predpoints@SSNPoints[[1]]@point.coords),]
    
    if(length(nbrs)>=2){
      if(sum(rid_nbrs==rid_remove)==2){
        # removed points and neighbors are located on the same segment
        coords_nbrs <- X@obspoints@SSNPoints[[1]]@point.coords[index,]
        #midpt_nbrs <- c(mean(coords_nbrs[,1]), mean(coords_nbrs[,2]))
        
        nbrs_order_index <- order(X@obspoints@SSNPoints[[1]]@point.data$upDist[index], decreasing=T)
        ordered_nbrs_coords <- X@obspoints@SSNPoints[[1]]@point.coords[index[nbrs_order_index],]
        
        #middle_pts <- c() # vector for centroids
        middle_pts_imsi <- c((ordered_nbrs_coords[1,1]+remove_coords[1])/2,  (ordered_nbrs_coords[1,2]+remove_coords[2])/2)
        middle_pts_imsi2 <- c((ordered_nbrs_coords[2,1]+remove_coords[1])/2,  (ordered_nbrs_coords[2,2]+remove_coords[2])/2)
        middle_pts <- rbind(middle_pts_imsi, middle_pts_imsi2)
        
        for(ii in 1:length(nbrs)){
          ID_remove <- c()
          for(kk in 1:length(I_lines[[remove]])){
            ID_remove <- c(ID_remove, strsplit(I_lines[[remove]][[1]]@ID, "-")[[kk]][1])
          }
          ll <- length(I_lines[[nbrs[ii]]])
          while(TRUE){
            if(strsplit(I_lines[[nbrs[ii]]][[ll]]@ID, "-")[[1]][1] %in% ID_remove){
              if(length(strsplit(I_lines[[nbrs[ii]]][[ll]]@ID, "-")[[1]])==1){
                I_lines[[nbrs[ii]]] <- I_lines[[nbrs[ii]]][-ll]
              }else{
                if(X@obspoints@SSNPoints[[1]]@point.data$upDist[index[ii]]>X@predpoints@SSNPoints[[1]]@point.data$upDist[length(X@predpoints@SSNPoints[[1]]@point.data$upDist)]){
                  I_lines[[nbrs[ii]]][[ll]]@Lines[[1]]@coords[2,] <- middle_pts[1,]
                }else{
                  I_lines[[nbrs[ii]]][[ll]]@Lines[[1]]@coords[1,] <- middle_pts[2,]
                }
              }
              ll <- ll - 1
            }else{
              ll <- ll - 1
            }
            if(ll==0){
              break
            }
          }
        }
        
      }else if(length(unique(c(rid_nbrs, rid_remove)))>=3){
        # point to be removed exists downstream
        for(ii in 1:length(nbrs)){
          ID_remove <- c()
          for(kk in 1:length(I_lines[[remove]])){
            ID_remove <- c(ID_remove, strsplit(I_lines[[remove]][[1]]@ID, "-")[[kk]][1])
          }
          ll <- length(I_lines[[nbrs[ii]]])
          while(TRUE){
            if(strsplit(I_lines[[nbrs[ii]]][[ll]]@ID, "-")[[1]][1] %in% ID_remove){
              I_lines[[nbrs[ii]]] <- I_lines[[nbrs[ii]]][-ll]
              ll <- ll - 1
            }else{
              break
            }
          }
        }
      }
      
    }else if(length(nbrs)==1){
      if(rid_remove==rid_nbrs){
        # two stream edges are equal
      }else{
        # two stream edges are different 
        ID_remove <- c()
        for(kk in 1:length(I_lines[[remove]])){
          ID_remove <- c(ID_remove, strsplit(I_lines[[remove]][[1]]@ID, "-")[[kk]][1])
        }
        ll <- length(I_lines[[nbrs]])
        while(TRUE){
          if(strsplit(I_lines[[nbrs]][[ll]]@ID, "-")[[1]][1] %in% ID_remove){
            I_lines[[nbrs]] <- I_lines[[nbrs]][-ll]
            ll <- ll - 1
          }else{
            break
          }
        }
      }
    }
  }
  
  alpha <- matrix(0, 1, length(nbrs))
  #if ((r > 1) & (r <= N)) {
  alpha <- lengths[index] * lengthrem/(sum(lengths[index]^2))
  coeff[nbrs] <- coeff[nbrs] - alpha * coeff[remove]
  lengths[index] <- as.row(lengths[index])
  prod <- gamweights * lengthrem
  prod <- as.row(prod)
  lengths[index] <- lengths[index] - prod
  #}
  #if ((r == 1) | (r == (N + 1))) {
  #  q <- which(pointsin == nbrs)
  #  alpha <- lengthrem/lengths[q]
  #  coeff[pointsin[q]] <- coeff[pointsin[q]] - alpha * coeff[remove]
  #  lengths[q] <- lengths[q] - lengthrem
  #}
  
  return(list(coeff = coeff, lengths = lengths, alpha = alpha, I_lines=I_lines))
}

invtnp_stream <- function (X, coeff, lengths, lengthsremove, pointsin, removelist, 
                           neighbrs, schemehist, interhist, nadd = length(X) - 2, intercept = TRUE, 
                           neighbours = 1,  closest = FALSE, LocalPred = LinearPred, length_diff) {
  
  # inverse lifting step
  if (is.list(X)) {
    coeff <- X$coeff
    lengths <- X$lengths
    lengthsremove <- X$lengthsremove
    removelist <- X$removelist
    neighbrs <- X$neighbrs
    pointsin <- X$pointsin
    schemehist <- X$schemehist
    interhist <- X$interhist
    #X <- X$network
  }
  
  #X <- as.row(X)
  coeff <- as.row(coeff)
  #n <- length(X)
  n <- length(coeff)
  N <- length(pointsin)
  m <- length(removelist) #n = N+m
  #d <- neighbours
  if (nadd > 0) {
    for (j in 1:nadd) {
      N <- length(pointsin)
      remove <- removelist[m - j + 1]
      lengthrem <- lengthsremove[m - j + 1]
      nbrs <- neighbrs[[m - j + 1]]
      index <- NULL
      index <- match(nbrs, pointsin)
      
      #decide r index
      #B <- (X[remove] > X[nbrs])
      B <- (remove > nbrs)
      nt <- sum(B)
      if (nt == 0) {
        r <- which(pointsin == nbrs[1])
      }
      if (nt == length(nbrs)) {
        r <- which(pointsin == nbrs[length(nbrs)]) +  1
      }
      if ((nt > 0) & (nt < length(nbrs))) {
        r <- which(pointsin == nbrs[nt + 1])
      }
      if (is.null(schemehist) == FALSE) {
        if (schemehist[m - j + 1] == "Linear") {
          res <- LinearPred(pointsin, X, coeff, nbrs, 
                            remove, intercept = interhist[m - j + 1], 
                            neighbours)
        }
        if (schemehist[m - j + 1] == "Quad") {
          res <- QuadPred(pointsin, X, coeff, nbrs, remove, 
                          intercept = interhist[m - j + 1], neighbours)
        }
        if (schemehist[m - j + 1] == "Cubic") {
          res <- CubicPred(pointsin, X, coeff, nbrs, 
                           remove, intercept = interhist[m - j + 1], 
                           neighbours)
        }
      }else {
        #schemehist가 null인 경우
        
        res <- streamPred(pointsin,  X=X$network, coeff=coeff, nbrs, remove, adj=X$adjacency, method="LC", intercept=FALSE, neighbours=NULL)
        
        #res <- LocalPred(pointsin, X, coeff, nbrs, remove,intercept, neighbours)
      }
      if (length(res) == 2) {
        l <- res[[1]]
      }else {
        l <- res
      }
      gamweights <- l$weights
      #UndoPointsUpdate함수의 stream version을 만들어야
      #l1 <- UndoPointsUpdate(X, coeff, nbrs, index, remove,  r, N, pointsin, gamweights, lengths, lengthrem)
      l1 <- UndoPointsUpdate_stream(X=X$network, coeff, nbrs, index, remove, r, N, pointsin, gamweights, lengths, lengthrem, I_lines=X$I_lines, stream=TRUE)
      X$I_lines <- l1$I_lines
      coeff <- l1$coeff
      lengths <- l1$lengths
      pred <- sum(as.column(gamweights) * coeff[nbrs])
      coeff[remove] <- coeff[remove] + pred
      removelist <- setdiff(removelist, remove)
      
      if (r == 1) {
        lengths <- c(lengthrem, lengths)
        pointsin <- c(remove, pointsin)
        #update example network
        X$network@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(X$network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),], X$network@obspoints@SSNPoints[[1]]@network.point.coords)
        X$network@obspoints@SSNPoints[[1]]@point.coords <- rbind(X$network@predpoints@SSNPoints[[1]]@point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),], X$network@obspoints@SSNPoints[[1]]@point.coords)
        if(length_diff<=0){
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X$network@obspoints@SSNPoints[[1]]@point.data)
        }else{
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(cbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X=NA, X2=NA, Sim_Values=NA) , X$network@obspoints@SSNPoints[[1]]@point.data)
        }
      }
      if (r == (N + 1)) {
        lengths <- c(lengths, lengthrem)
        pointsin <- c(pointsin, remove)
        #update example network
        X$network@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@network.point.coords, X$network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),])
        X$network@obspoints@SSNPoints[[1]]@point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@point.coords, X$network@predpoints@SSNPoints[[1]]@point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),])
        if(length_diff<=0){
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data, X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),])
        }else{
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data, cbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X=NA, X2=NA, Sim_Values=NA))
        }
      }
      if ((r > 1) & (r < (N + 1))) {
        lengths <- c(lengths[1:(r - 1)], lengthrem, lengths[r:N])
        pointsin <- c(pointsin[1:(r - 1)], remove, pointsin[r:N])
        #update example network
        X$network@obspoints@SSNPoints[[1]]@network.point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@network.point.coords[c(1:(r - 1)),], X$network@predpoints@SSNPoints[[1]]@network.point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),], X$network@obspoints@SSNPoints[[1]]@network.point.coords[c(r:N),])
        X$network@obspoints@SSNPoints[[1]]@point.coords <- rbind(X$network@obspoints@SSNPoints[[1]]@point.coords[c(1:(r - 1)),], X$network@predpoints@SSNPoints[[1]]@point.coords[nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),], X$network@obspoints@SSNPoints[[1]]@point.coords[c(r:N),])
        if(length_diff<=0){
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data[c(1:(r - 1)),], X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X$network@obspoints@SSNPoints[[1]]@point.data[c(r:N),])
        }else{
          X$network@obspoints@SSNPoints[[1]]@point.data <- rbind(X$network@obspoints@SSNPoints[[1]]@point.data[c(1:(r - 1)),], cbind(X$network@predpoints@SSNPoints[[1]]@point.data[nrow(X$network@predpoints@SSNPoints[[1]]@point.data),], X=NA, X2=NA, Sim_Values=NA), X$network@obspoints@SSNPoints[[1]]@point.data[c(r:N),])
        }
      }
      #delete points from predpoint
      X$network@predpoints@SSNPoints[[1]]@network.point.coords <- X$network@predpoints@SSNPoints[[1]]@network.point.coords[-nrow(X$network@predpoints@SSNPoints[[1]]@network.point.coords),]
      X$network@predpoints@SSNPoints[[1]]@point.coords <- X$network@predpoints@SSNPoints[[1]]@point.coords[-nrow(X$network@predpoints@SSNPoints[[1]]@point.coords),]
      X$network@predpoints@SSNPoints[[1]]@point.data <- X$network@predpoints@SSNPoints[[1]]@point.data[-nrow(X$network@predpoints@SSNPoints[[1]]@point.data),]
    }
  }
  
  
  
  return(list(coeff = coeff, lengths = lengths, lengthsremove = lengthsremove, 
              pointsin = pointsin, removelist = removelist))
}


########################################
##Data-adaptive stream flow weight selection
########################################  



########################################
##Network poltting function in the future
########################################




########################################
##SpatialStreamNetwork object making function in SSN
########################################
importSSN_stream <- function(shreve_obj, location="Miho", multipleplaces=FALSE, RID="original"){
  # location: "Miho", "Full" 
  # multipleplaces=FALSE -> fix and sum
  # multipleplaces=TRUE -> go ahead
  
  #netID = 1 (b/c there is one network)
  
  #(a) stream network data
  edges  <- readOGR("/home/kyu9510/pca_on_river/Data/KRF_3.0_Geumgang/KRF_ver3_LINE_금강수계.shp")
  #rownames(edges@data) <- edges@data[, "RCH_ID"]
  #reach ID (print rid)
  
  #(b) 관측 장소 자료
  sites <- read.csv("/home/kyu9510/pca_on_river/Data/ProcessedData/금강장소(엑셀편집)UTF8.csv",header=T)
  sites <- sites[which(sites$총량==1),]
  sites <- data.frame(name= sites$이름,x=sites$경도.Degree., y=sites$위도.Degree.)
  
  #rownames(sites@data) <- sites@data[,"pid"]
  #rownames(sites@coords) <- sites@data[, "pid"]
  
  #(c) 수질 데이터 자료
  alldata1 <- read.table("/home/kyu9510/pca_on_river/Data/ProcessedData/수질20082013.txt",sep=',', header=T)
  alldata2 <- read.table("/home/kyu9510/pca_on_river/Data/ProcessedData/수질20132018.txt",sep=',', header=T)
  alldata <- rbind(alldata1, alldata2)
  
  #(d) 예전에 편집한 자료
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
  #3011: 미호천
  #3012: 금강공주
  #3008: 대청댐
  #3007: 보청천
  #3005: 초강
  #3009: 갑천
  #3013: 논산천
  #3006: 대청댐상류
  #3014: 금강하구언
  #3004: 영동천
  #3003: 무주남대천
  #3001: 용담댐
  #3002: 용담댐하류
  #3310: 대청댐하류
  #3301: 만경강(x)
  #3303: 직소천(x)
  #3302: 동진강(x)
  #3101: 삽교천(x)
  #2005: 금강상류부분(필요없음)
  #3202: 부남방조제(x)
  #1101, 3201: 대호방조제(x)
  #3203: 금강서해(x)
  #5301, 5302: 주진천(x)
  
  #for(i in 1:length(region_index)){
  #  shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  #  plot(shape, main=region_index[i])
  #  plot(shape_sub, col="red", add=T, lwd=3)
  #}
  
  name_index <- c("미호천", "금강공주", "대청댐", "보청천", "초강", "갑천", "논산천", "대청댐상류", "금강하구언", "영동천", "무주남대천", "용담댐", "용담댐하류", "대청댐하류")
  name_subindex <- list()
  name_subindex[[1]] <- c("미호댐", "원남댐", "병천천상류", "병천천하류", "보강천", "미호천중류", "미호천상류", "작천보", "한천", "백곡천", "조천", "석화수위표", "미호천하류", "무심천", "백곡댐" ) #15개
  name_subindex[[2]] <- c("세종보", "지천상류", "지천합류후", "어천합류후", "용수천", "백제보", "석성천", "규암수위표", "금천", "논산천합류전", "유구천", "공주보", "대교천", "공주수위표") #14개
  name_subindex[[3]] <- c("대청댐", "대청댐조정지", "소옥천상류", "소옥천하류", "대청댐상류") #5개
  name_subindex[[4]] <- c("항건천", "보청천중류", "삼가천합류후", "삼가천", "보청천상류", "보청천하류") #6개
  name_subindex[[5]] <- c("석천", "초강상류", "초강하류") #3개
  name_subindex[[6]] <- c("갑천하류", "유성수위표", "갑천상류", "대전천", "유등천상류", "유등천하류") #6개
  name_subindex[[7]] <- c("노성천", "논산천하류", "강경천", "논산천상류", "탑정댐") #5개
  name_subindex[[8]] <- c("보청천합류전") #1개
  name_subindex[[9]] <- c("입포수위표", "길산천", "금강하구언") #3개
  name_subindex[[10]] <- c("영동천", "봉황천하류", "초강합류후", "봉황천상류", "봉황천합류후", "호탄수위표") #6개
  name_subindex[[11]] <- c("무주남대천상류", "무주남대천중류", "무주남대천하류") #3개
  name_subindex[[12]] <- c("주자천", "용담댐", "정자천", "구량천", "진안천", "장계천합류후", "장계천", "진안천합류후") #8개
  name_subindex[[13]] <- c("용담댐하류") #1개
  name_subindex[[14]] <- c("미호천합류전", "매포수위표") #2개
  region_subindex <- list()
  
  Geum_RCH_ID <- c()
  for(i in 1:length(region_index)){
    shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
    #plot(shape_sub)
    region_subindex[[i]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6)))
    for(j in 1:length(region_subindex[[i]])){
      shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
      #plot(shape_sub_sub, col="red", add=T)
      Geum_RCH_ID_Imsi <- cbind(shape_sub_sub@data, data.frame(중권역번호=region_index[i], 중권역이름=name_index[i], 소권역번호=region_subindex[[i]][j], 소권역이름=name_subindex[[i]][j]))
      Geum_RCH_ID <- rbind(Geum_RCH_ID, Geum_RCH_ID_Imsi)
    }
  }
  
  # when print subregion
  #i=8
  #shape_sub <- subset(shape, sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==region_index[i]) )
  #plot(shape_sub)
  #region_subindex[[i]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6)))
  #region_subsubindex[[i]] <- list()
  #
  #j=1
  #shape_sub_sub <- subset(shape_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=6))==region_subindex[[i]][j])
  #plot(shape_sub_sub, col="red", add=T)
  #region_subsubindex[[i]][[j]] <-unique(sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=8)))
  #
  #for(k in 1:length(region_subsubindex[[i]][[j]])){
  #  shape_sub_sub_sub <- subset(shape_sub_sub, sapply(as.character(shape_sub@data$RCH_ID), function(x) substr(x, start=1, stop=8))==region_subsubindex[[i]][[j]][k])
  #  plot(shape_sub_sub_sub, col=k, add=T)
  #}
  
  ########################################
  ##Adjacency matrix
  ########################################
  #reference: https://cran.r-project.org/web/packages/riverdist/vignettes/riverdist_vignette.html
  
  #connections: Object of class "matrix", with "numeric" elements. Defined as a square matrix, with elements describing the type of connection detected between line segments.
  #• A value of 1 in element [i,j] indicates that the beginning of segment i is connected to the beginning of segment j.
  #• A value of 2 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j.
  #• A value of 3 in element [i,j] indicates that the end of segment i is connected to the beginning of segment j.
  #• A value of 4 in element [i,j] indicates that the end of segment i is connected to the end of segment j.
  #• A value of 5 in element [i,j] indicates that segments i and j are connected at both beginning and end.
  #• A value of 6 in element [i,j] indicates that the beginning of segment i is connected to the end of segment j, and the end of segment i is connected to the beginning of segment j.
  #• A value of NA in element [i,j] indicates that segments i and j are not connected.
  if(location=="Miho"){
    flow_network_original <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="GeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    flow_network <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="thinGeumMiho", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    #plot(flow_network)
  }else if(location=="Full"){
    flow_network_original <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="Geum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    flow_network <- line2network(path = "/home/kyu9510/pca_on_river/Data/ProcessedData/", layer="thinGeum", tolerance =  0.000001, reproject ="+proj=longlat +datum=WGS84", supplyprojection = NULL)
    #plot(flow_network)
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
  #reach_length <- shape@data$Shape_Leng[which(sapply(as.character(shape@data$RCH_ID), function(x) substr(x, start=1, stop=4)==rregion_subindex[[1]][i]) )]
  

  
  #ind1 <- colnames(edges@data) == c("netID")
  #ind2 <- colnames(edges@data) == c("rid")
  #ind3 <- colnames(edges@data) == c("upDist")
  #if (is.factor(edges@data$netID)) {
  #  edges@data$netID <- as.character(edges@data$netID)
  #}
  
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
    #17: 추풍령천, 18: 초강1
    data$eventplace$RCH_ID[17] <- as.integer(30050108)
    #23: 옥천, 24: 우산 (옥천으로 평균)
    #36: 대청댐, 38: 대청
    data$eventplace$RCH_ID[36] <- as.integer(30080409); data$eventplace$RCH_ID[38] <- as.integer(30080419)
    #43: 갑천1, 44: 갑천2
    data$eventplace$RCH_ID[44] <- as.integer(30090208)
    #47: 침산교, 48: 유등천1
    data$eventplace$RCH_ID[47] <- as.integer(30090502)
    #49: 옥계교, 50: 대전천1 (두개 평균을내어야 (대전천1로 평균))
    #55: 갑천5, 56: 갑천5-1 (갑천5로 평균)
    #58: 청원-1 59: 백천
    data$eventplace$RCH_ID[58] <- as.integer(30100203)
    #79: 미호천8, 83: 미호천5
    data$eventplace$RCH_ID[79] <- as.integer(30111301)
    #85: 조천, 86: 조천1
    data$eventplace$RCH_ID[86] <- as.integer(30111505)
    #98: 곰나루, 99: 금강
    data$eventplace$RCH_ID[98] <- as.integer(30120532)
    
    
    #[1] 구량천       가막         진안천       정자천       용담         용포         무주남대천   무주남대천-1 무주남대천-2 부리        
    #[11] 제원         봉황천       제원A        영동         영동천1      영동천2      추풍령천     초강1        금계천       초강2       
    #[21] 이원         옥천         우산         보청천1      보청천2      항건천       보청천3      보청천4      상곡천       추풍천      
    #[31] 옥천천       회인천       주원천       대청댐       품곡천       대청         용호천       현도         두계천1      봉곡2교     
    #[41] 갑천1        갑천2        갑천3        침산교       유등천1      옥계교       대전천1      대전천2      대전천3      유등천5     
    #[51] 갑천4        갑천5        갑천5-1      세종1        청원-1       백천         미호천1      칠장천       미호천1-1    한천        
    #[61] 미호천2      백곡천1      백곡천2      미호천7      초평천       미호천3      보강천1      보강천       성암천       무심천      
    #[71] 무심천1      무심천2      무심천3      석남천       미호천8      병천천A      용두천       병천천       미호천5      미호천5A    
    #[81] 조천         조천1        월하천       미호천6-1    연기         금남         용수천-1     용수천       대교천       대교천2     
    #[91] 세종2        공주1        정안천       곰나루       금강         유구천       목면         공주2        부여         지천        
    #[101] 지천-1       정동         은산천       부여1        금천         부여2        석성천-1     석성천       석성천2      성동        
    #[111] 논산천-1     논산천1      연산천       노성천-1     노성천       논산천2      방축천       논산천4      수철천       마산천      
    #[121] 강경천       강경         산북천       양화-1       길산천       길산천2      금강갑문    
    #}else{
    #  
    #}
    
  }
  dist.pts <- rep(0, nrow(data$eventplace))
  dist.imsi <- c()
  for(jj in 1:nrow(data$eventplace)){
    pts.coords <- c(data$eventplace$경도.Degree.[jj], data$eventplace$위도.Degree.[jj])
    
    kk <- match(data$eventplace$RCH_ID[jj], shape@data$RCH_ID) # fixme !!

    #dist_upper <- apply(shape@data[,c(3,4)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
    #dist_lower <- apply(shape@data[,c(5,6)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
    #
    #kk <- sort(unique(union(which(dist_upper==min(c(dist_upper, dist_lower))), which(dist_lower==min(c(dist_upper, dist_lower))))))
    #
    #dist.imsiimsi <- c()
    #for(pp in 1:length(kk)){
    #  for(ll in 1:nrow(flow_network_original$lines[[kk[pp]]])){
    #    dist.compute <- sqrt(sum((flow_network_original$lines[[kk[pp]]][ll,]-pts.coords)^2))
    #    dist.imsiimsivec <- matrix(c(dist.compute, kk[pp], ll), nrow=1)
    #    dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
    #  }
    #}
    #dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
    #dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
    
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
    #삼가천, 석천: no obs.points
    #유등천A, 미호천4: data$eventplace$X exist but colnames(data$normaldata_TN) not
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
    ##(738, 737) -> (737, 751)
    dist.imsi <- c()
    lll <- 0
    for(jj in 1:nrow(data$eventplace)){
      pts.coords <- c(data$eventplace$경도.Degree.[jj], data$eventplace$위도.Degree.[jj])
      
      lll <- lll+1
      if((jj%in%miho)){
        lll <-lll - 1
      }
      kk <- RCH_ID_manual[lll]
      #kk <- match(RCH_ID_manual[lll], shape@data$RCH_ID) # to be fixed
      #kk <- match(data$eventplace$RCH_ID[jj], shape@data$RCH_ID) # to be fixed
  
      #dist_upper <- apply(shape@data[,c(3,4)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
      #dist_lower <- apply(shape@data[,c(5,6)], 1, function(x) sqrt(sum((x-pts.coords)^2)) )
      #
      #kk <- sort(unique(union(which(dist_upper==min(c(dist_upper, dist_lower))), which(dist_lower==min(c(dist_upper, dist_lower))))))
      #
      #dist.imsiimsi <- c()
      #for(pp in 1:length(kk)){
      #  for(ll in 1:nrow(flow_network_original$lines[[kk[pp]]])){
      #    dist.compute <- sqrt(sum((flow_network_original$lines[[kk[pp]]][ll,]-pts.coords)^2))
      #    dist.imsiimsivec <- matrix(c(dist.compute, kk[pp], ll), nrow=1)
      #    dist.imsiimsi <- rbind(dist.imsiimsi,dist.imsiimsivec)
      #  }
      #}
      #dist.imsiimsi <- dist.imsiimsi[which.min(dist.imsiimsi[,1]),]
      #dist.imsi <- rbind(dist.imsi, dist.imsiimsi)
      
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
    
    #data$eventplace$RCH_ID <- shape@data$RCH_ID[RCH_ID_manual]
    pointdata <- cbind(data$eventplace[-miho,], rid_old=dist.imsi[-miho,2], upDist=DistanceUpstream_obs[-miho], locID=c(1:(nrow(data$eventplace)-length(miho))), netID=rep(1,(nrow(data$eventplace)-length(miho))), pid=setdiff(c(1:nrow(data$eventplace)), miho), shreve=shreve_obj[dist.imsi[-miho,2]], mindist=dist.imsi[-miho,1], rid=RCH_ID_manual, segid=dist.imsi[-miho,3])
  }else{
    pointdata <- cbind(data$eventplace[-miho,], rid_old=dist.imsi[-miho,2], upDist=DistanceUpstream_obs[-miho], locID=c(1:(nrow(data$eventplace)-length(miho))), netID=rep(1,(nrow(data$eventplace)-length(miho))), pid=setdiff(c(1:nrow(data$eventplace)), miho), shreve=shreve_obj[dist.imsi[-miho,2]], mindist=dist.imsi[-miho,1], rid=dist.imsi[-miho,2], segid=dist.imsi[-miho,3])
  }
  pointbbox <- rbind(range(pointcoords[-miho,1]), range(pointcoords[-miho,2])); colnames(pointbbox) <- c("min", "max"); rownames(pointbbox) <- c("coords.x1", "coords.x2")
  network.SSNPoints <- list()
  network.SSNPoints[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.obsdata, point.coords = pointcoords, point.data=pointdata, points.bbox=pointbbox, proj4string=shape@proj4string)
  obspoints <- new("SSNPoints", SSNPoints=network.SSNPoints, ID="Obs")
  #(3)predpoint
  network.line.coords.preddata <- data.frame(NetworkID=rep(1,length(miho)), SegmentID=dist.imsi[miho,2], DistanceUpstream=DistanceUpstream_obs[miho])
  rownames(network.line.coords.preddata) <- rownames(dist.imsi)[miho]
  pointcoords.preddata <- cbind(data$eventplace$경도.Degree.[miho], data$eventplace$위도.Degree.[miho]); colnames(pointcoords.preddata) <- c("coords.x1", "coords.x2")
  #pointdata에 ratio 추가해야
  pointdata.preddata <- cbind(data$eventplace[miho,], rid_old=dist.imsi[miho,2], upDist=DistanceUpstream_obs[miho], locID=c(1:length(miho)), netID=rep(1,length(miho)), pid=miho, shreve=shreve_obj[dist.imsi[miho,2]], mindist=dist.imsi[miho,1], rid=dist.imsi[miho,2], segid=dist.imsi[miho,3])
  pointbbox.preddata <- rbind(range(pointcoords[miho,1]), range(pointcoords[miho,2])); colnames(pointbbox.preddata) <- c("min", "max"); rownames(pointbbox.preddata) <- c("coords.x1", "coords.x2")
  network.SSNPoints.preddata <- list()
  network.SSNPoints.preddata[[1]] <- new("SSNPoint", network.point.coords = network.line.coords.preddata, point.coords = pointcoords.preddata, point.data=pointdata.preddata, points.bbox=pointbbox.preddata, proj4string=shape@proj4string)
  if(length(miho)==0){
    #(3-1) no predpt case
    predpoints <- new("SSNPoints", SSNPoints=list(), ID=character(0))
  }else{
    predpoints <- new("SSNPoints", SSNPoints=network.SSNPoints.preddata, ID="Preds")
  }
  #(4)SSNpath: path
  SSNpath <- "NULL"
  #(5)data: data.frame type ( rid,upDist,    Length, netID)
  #data:Object of class "data.frame". 
  #The number of rows in data should equal the number of lines in the lines object. Row names correspond to SegmentID values
  shapedata <- shape@data
  shapedata <- cbind(shapedata, rid=c(1:nrow(shape@data)), upDist=DistanceUpstream_seg, Length=shape@data$Shape_Leng, netID=rep(1,nrow(shape@data)), shreve=shreve_obj )
  #(6)lines: segment polygon
  shapelines <- shape@lines
  #(7)bbox
  shapebbox <- shape@bbox; colnames(shapebbox) = c("min", "max"); rownames(shapebbox) <- c("x", "y")
  #(8)roj4string
  shapeproj4string <- shape@proj4string
  
  Geum_original <- new("SpatialStreamNetwork", network.line.coords=network.line.coords, obspoints=obspoints, predpoints=predpoints, path=SSNpath, data=shapedata, lines=shapelines, bbox=shapebbox, proj4string=shapeproj4string)
  
  #Object of class "SSNPoints" with 2 slots
  #
  #@ SSNPoints: List of SSNPoint objects with 5 slots
  #@ network.point.coords: object of class "data.frame". Row names
  #represent point identifiers (pid) stored in the point.data
  #data.frame.
  #$ NetworkID: factor identifying the NetworkID of that point
  #$ SegmentID: factor identifying the unique stream segment of that point
  #$ DistanceUpstream: numeric value representing the cumulative
  #distance from the network outlet, the most downstream point
  #on a network, to that point
  #@ point.coords: numeric matrix or "data.frame" with x- and y-
  #  coordinates (each row is a point); row names represent point
  #identifiers (pid) stored in the point.data data.frame.
  #@ point.data: object of class "data.frame"; the number of rows in
  #data should equal the number of points in the
  #network.point.coords object; row names are set to the pid
  #attribute.
  #@ points.bbox: Object of class "matrix"; see Spatial-class
  #@ proj4string: Object of class "CRS"; see CRS-class
  #@ ID: character string representing the name of the observation points
  
  
  #network.line <- shape_tinned@lines
  #bbox, proj4string
  
  return(Geum_original)
}


get_binaryIDs_usual <- function(mouth_node, shape){
  # "get_binaryIDs_stream" needs RCH_ID, binaryID in shape file 
  # new version (generalized): for shape files in stpca package
  rtNEL1 <-readshpnw(shape, ELComputed=TRUE, longlat=TRUE) # length  7
  igr1 <-nel2igraph(rtNEL1[[2]], rtNEL1[[3]]) #shape file to graph
  
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
  
}

# mouth_node=68; shape=shape_miho
# need RCH_ID,binaryID in shape file
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
  
  
  #result_binaryIDs <- data.frame(rid=c(1:length(shape$RCH_ID)), binaryID=rep("", length(shape$RCH_ID)))
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
    #print(key_RID_ind)
  }
  return(result_binaryIDs)
}

# functions for get_adjacency_stream 
re_map_rid <- function(rid_vector, all_rid){
  # in cases where multiple networks are present, gaps could occur in the rid sequence
  # this function maps the rid vector in a simple way so that it agrees with the adjacency matrix
  mapped_rid <- vector("numeric", length = length(rid_vector))
  rng_rid    <- range(all_rid)
  new_key    <- 1:length(unique(all_rid))
  old_rid    <- sort(all_rid)
  for(i in 1:length(rid_vector)){
    mapped_rid[i] <- new_key[which(old_rid == rid_vector[i])]
  }
  mapped_rid
}

# identical to get_adjacency but obtain binaryIDs_obj from get_binaryIDs_stream
get_adjacency_stream <- function (binaryIDs_obj, netID = 1) 
{
  binaryIDs <- binaryIDs_obj
  #binaryIDs <- get_binaryIDs(ssn_directory, net = netID)
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

#compute shreve stream order
compute_shreve <- function(adjacency){
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y자
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
    }
    #case 2: 1자
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  return(shreve)
}


#compute shreve stream order
compute_shreve_and_dist <- function(adjacency, example_network, type="equal", scalevec=c(0.2, 1.5), logdata=FALSE){
  #type="lengthprop": shreve_and_dist[first_stream] is defined by proportional to their lengths
  #type="equal": shreve_and_dist[first_stream] is defined by equal weight (== shreve dist style)
  #scalevec: final weight를 scaling (default = (0.2,1.5), (refer to O'Donnell's method))
  
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  shreve_and_dist <- rep(0, nrow(adj_mat))
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  if(type=="lengthprop"){
    shreve_and_dist[first_stream] <- example_network@data$Length[first_stream]
  }else if(type=="equal"){
    shreve_and_dist[first_stream] <- 1
  }else{
    stop("wrong types")
  }
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y shape
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
      shreve_and_dist[next_stream1[i]] <- sum(shreve_and_dist[which(adj_mat[,next_stream1[i]]!=0)]) + example_network@data$Length[next_stream1[i]]
    }
    #case 2: 1 shape
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
      shreve_and_dist[next_stream2[j]] <- sum(shreve_and_dist[which(adj_mat[,next_stream2[j]]!=0)]) + example_network@data$Length[next_stream2[j]]
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  #scaling
  #weight_vec_candidate2 <- scales::rescale(log(weight_vec_candidate2$distweight), to=c(0.2, 1.5))
  if(logdata==TRUE){
    shreve_and_dist <- scales::rescale(log(sqrt(shreve_and_dist)), to=c(scalevec[1], scalevec[2]))
    #shreve_and_dist <- scales::rescale(log(shreve_and_dist), to=c(scalevec[1], scalevec[2]))
  }else{
    shreve_and_dist <- scales::rescale((sqrt(shreve_and_dist)), to=c(scalevec[1], scalevec[2]))
  }
  return(list(shreve=shreve, distweight=shreve_and_dist))
}


compute_shreve_weighted <- function(adjacency, weight){
  # compute_shreve_and_dist takes length of example_network as an input
  # instead, take weight object as an input
  # extend to weighted version
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  shreve_and_dist <- rep(0, nrow(adj_mat))
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  shreve_and_dist[first_stream] <- weight[first_stream]
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y shape
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
      shreve_and_dist[next_stream1[i]] <- sum(shreve_and_dist[which(adj_mat[,next_stream1[i]]!=0)]) + weight[next_stream1[i]]
    }
    #case 2: 1 shape
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
      shreve_and_dist[next_stream2[j]] <- sum(shreve_and_dist[which(adj_mat[,next_stream2[j]]!=0)]) + weight[next_stream2[j]]
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  return(list(shreve=shreve, distweight=shreve_and_dist))
}


########################################
##Graphics
########################################

source("/home/kyu9510/pca_on_river/Code/iwanthue.R", chdir = TRUE)

plot.SpatialStreamNetwork.SC <- 
  function (x, VariableName = NULL, color.palette = NULL, nclasses = NULL, 
            breaktype = "quantile", brks = NULL, PredPointsID = NULL, 
            add = FALSE, addWithLegend = FALSE, lwdLineCol = NULL, lwdLineEx = 1, 
            lineCol = "black", I_lines=NULL, pointsin=NULL, ...) 
  {
    #my addition
    if(is.null(pointsin)){
      pointsin <- c(1:length(I_lines))
    }
    if (missing(lwdLineEx)) 
      lwdLineEx <- 1
    if (missing(lwdLineCol)) {
      x@data$lineWidth <- rep(1, nrow(x@data))
      lwdLineCol <- "lineWidth"
    }
    if (is.null(as.list(match.call()[-1])$pch)) {
      plch = 19
    }
    else plch <- as.list(match.call()[-1])$pch
    if (is.null(as.list(match.call()[-1])$cex)) {
      chex = 1
    }
    else chex <- as.list(match.call()[-1])$cex
    if (is.null(as.list(match.call()[-1])$col)) {
      colr = "black"
    }
    else colr <- as.list(match.call()[-1])$col
    par.orig <- par(no.readonly = TRUE)
    if (!is.null(PredPointsID)) {
      # when prediction Point ID : input
      for (i in 1:length(x@predpoints@ID)) {
        if (x@predpoints@ID[i] == PredPointsID) {
          if (add == FALSE & addWithLegend == FALSE) {
            plot(x@bbox[1, ], x@bbox[2, ], type = "n",   ...)
            for (j in 1:length(x@lines)) for (k in 1:length(x@lines[[j]])) if (is.null(lwdLineCol)) 
              lines((x@lines[[j]]@Lines[[k]]@coords), col = lineCol,  ...)
            else lines(x@lines[[j]]@Lines[[k]]@coords, 
                       lwd = lwdLineEx * x@data[i, lwdLineCol], 
                       col = lineCol, ...)
          }
          if (add == TRUE) {
            par(new = TRUE)
            plot(x@bbox[1, ], x@bbox[2, ], type = "n", 
                 bty = "n", xlab = "", ylab = "", ...)
          }
          if (addWithLegend == TRUE) {
            par(new = TRUE)
            layout(matrix(1:2, nrow = 1), widths = c(4, 
                                                     1))
            par(mar = c(5, 5, 3, 0))
            par(mfg = c(1, 1))
            plot(x@bbox[1, ], x@bbox[2, ], type = "n", 
                 bty = "n", xlab = "", ylab = "", ...)
          }
          points(x@predpoints@SSNPoints[[i]]@point.coords, 
                 pch = plch, cex = chex, col = colr)
        }
      }
      par(par.orig)
    }else if (is.null(VariableName)) {
      #variable name = null case
      plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
      for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
        lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol, 
              ...)
      else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * 
                   x@data[i, lwdLineCol], col = lineCol, ...)
      points(x@obspoints@SSNPoints[[1]]@point.coords, pch = plch, 
             cex = chex, col = colr)
      par(par.orig)
    }else {
      # most cases 
      layout(matrix(1:2, nrow = 1), widths = c(4, 1))
      par(mar = c(5, 5, 3, 0))
      plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
      #for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
      #  lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol,  ...)
      #else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * 
      #             x@data[i, lwdLineCol], col = lineCol, ...)
      
      #define color palette
      #color.palette = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      #set.seed(1)
      #col.stream <- sample(color.palette, length(I_lines))
      #pie(rep(1,length(I_lines)), col=col.stream)
      col.stream <- iwanthue(n=length(I_lines), hmin=0, hmax=360, cmin=30, cmax=80, lmin=35, lmax=80, plot=FALSE, random=1) 
      
      for (i in 1:length(pointsin)){
        for (j in 1:length(I_lines[[pointsin[i]]])){
          lines(I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[c(I_lines[[pointsin[i]]][[j]]@range[1]:I_lines[[pointsin[i]]][[j]]@range[2]),], col=col.stream[i], lwd=lwdLineEx*x@obspoints@SSNPoints[[1]]@point.data[i, lwdLineCol])
          if(I_lines[[pointsin[i]]][[j]]@range[1]!=1){
            mid_pts <- (I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[i]][[j]]@range[1],] + I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[(I_lines[[pointsin[i]]][[j]]@range[1]-1),])/2
            lines(rbind(mid_pts, I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[pointsin[i]]][[j]]@range[1],]), col=col.stream[i], lwd=lwdLineEx*x@obspoints@SSNPoints[[1]]@point.data[i, lwdLineCol])
          }
          if(I_lines[[pointsin[i]]][[j]]@range[2]!=nrow(I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords)){
            mid_pts <- (I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[pointsin[i]]][[j]]@range[2],] + I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[(I_lines[[pointsin[i]]][[j]]@range[2]+1),])/2
            lines(rbind(mid_pts, I_lines[[pointsin[i]]][[j]]@Lines[[1]]@coords[I_lines[[pointsin[i]]][[j]]@range[2],]), col=col.stream[i], lwd=lwdLineEx*x@obspoints@SSNPoints[[1]]@point.data[i, lwdLineCol])
          }
        }
      }
      
      data <- x@obspoints@SSNPoints[[1]]@point.data
      if (is.null(nclasses)) 
        nclasses <- 10
      lower.breaks <- matrix(0, nrow = nclasses, ncol = 1)
      upper.breaks <- matrix(0, nrow = nclasses, ncol = 1)
      if (breaktype == "quantile") {
        brks <- quantile(data[, VariableName], probs = (1:(nclasses -    1))/nclasses, na.rm = T)
        lower.breaks <- c(min(data[, VariableName], na.rm = T),   brks)
        upper.breaks <- c(brks, max(data[, VariableName],     na.rm = T))
      }
      if (breaktype == "even") {
        brks <- min(data[, VariableName]) + (max(data[, VariableName]) - min(data[, VariableName])) * (1:(nclasses - 1))/nclasses
        lower.breaks <- c(min(data[, VariableName], na.rm = T),   brks)
        upper.breaks <- c(brks, max(data[, VariableName],  na.rm = T))
      }
      if (breaktype == "user") {
        if (is.null(brks)) 
          return("Must specify brks if breaktype = user")
        minD <- min(data[, VariableName], na.rm = TRUE)
        maxD <- max(data[, VariableName], na.rm = TRUE)
        brks <- as.vector(unlist(brks))
        if (minD < min(brks)) 
          brks <- c(brks, minD)
        if (maxD > max(brks)) 
          brks <- c(brks, maxD)
        brks <- sort(unique(unlist(brks)))
        nclasses <- length(brks) - 1
        lower.breaks <- brks[1:nclasses]
        upper.breaks <- brks[2:(nclasses + 1)]
      }
      if (length(color.palette) == 0) 
        color.palette <- rainbow(nclasses, start = 0.66, 
                                 end = 0.99)
      for (j in 1:nclasses) {
        jmax <- upper.breaks[j]
        jmin <- lower.breaks[j]
        indj <- data[, VariableName] >= jmin & data[, VariableName] <=   jmax
        points(x@obspoints@SSNPoints[[1]]@point.coords[indj, , drop = F], col = color.palette[j], pch = plch,  cex = chex)
        #my addition
        text(x@obspoints@SSNPoints[[1]]@point.coords[indj, 1 , drop = F]+0.005, x@obspoints@SSNPoints[[1]]@point.coords[indj, 2 , drop = F]+0.005, pointsin[which(indj)])
      }
      dec.dig <- 2
      left <- as.character(as.numeric(as.integer(lower.breaks * 
                                                   10^dec.dig))/10^dec.dig)
      rght <- as.character(as.numeric(as.integer(upper.breaks * 
                                                   10^dec.dig))/10^dec.dig)
      leglabs <- paste(left, "to", rght)
      par(mar = c(0, 0, 0, 0))
      plot(c(0, 0), c(1, 1), type = "n", xaxt = "n", yaxt = "n", 
           xlab = "", ylab = "", bty = "n")
      legend(x = -1, y = 1.1, legend = leglabs, bty = "n", 
             pch = rep(plch, times = length(leglabs)), col = color.palette, 
             cex = 0.8)
      par(par.orig)
      return(invisible(data.frame(lower.breaks = lower.breaks, 
                                  upper.breaks = upper.breaks)))
    }
  }



########################################
##CODES not directly related to the SSN and smnet package
########################################

find_orthogonal <- function(x,y){
  # x: vector to be projected
  # y: point to be projected
  cc = crossprod(x,y)/crossprod(x)
  return(c(x[1]*cc,x[2]*cc))
}


forward_line <- function(line_remove, line_nbrs, lengths_remove, lengths_nbrs, adj_matrix, remove, nbrs, type="b2"){
  if(type=="a1" | type=="b1"){
    # case a1
    # upstream: line_remove, downstream: line_nbrs
    newpt <- (line_remove[1,] + find_orthogonal(x=(line_nbrs[2,]- line_remove[1,]), y=(line_remove[2,] - line_remove[1,]) ) )
    
    coarse_line <- rbind(line_remove[1,], line_nbrs[2,])
    detail_line <- rbind(line_remove[2,], newpt)
    
    coarse_lengths <- lengths_nbrs + lengths_remove
    detail_lengths <- lengths_remove - coarse_lengths
    
    #adjacency matrix change
    adj_matrix_imsi <- adj_matrix
    
    if(type=="b1"){
      
      adj_matrix_imsi[nbrs,which(adj_matrix[remove,]%in%c(2))] <- adj_matrix[remove,which(adj_matrix[remove,]%in%c(2))]
      adj_matrix_imsi[which(adj_matrix[,remove]%in%c(2)),nbrs] <- adj_matrix[remove,which(adj_matrix[,remove]%in%c(2))]
    }
    # change detail part negatives
    adj_matrix_imsi[remove,] <- -abs(adj_matrix_imsi[remove,])
    adj_matrix_imsi[,remove] <- -abs(adj_matrix_imsi[,remove])
    
  }else if(type=="a2" | type=="b2"){
    # case b2
    # upstream: line_nbrs, downstream: line_remove
    #test_a <- rbind(line_nbrs[1,], line_remove[2,])
    #test_b <- coeff_lines[[remove]][1,]
    
    newpt <- (line_nbrs[1,] + find_orthogonal(x=(line_remove[2,] - line_nbrs[1,]) , y= (line_nbrs[2,] - line_nbrs[1,]) ) )
    
    coarse_line <- rbind(line_nbrs[1,], line_remove[2,])
    detail_line <- rbind(line_remove[1,], newpt)
    
    coarse_lengths <- lengths_nbrs + lengths_remove
    detail_lengths <- lengths_remove - coarse_lengths
    
    
    #adjacency matrix change
    adj_matrix_imsi <- adj_matrix
    
    if(type=="b2"){
      
      adj_matrix_imsi[nbrs,which(adj_matrix[remove,]%in%c(3,4))] <- adj_matrix[remove,which(adj_matrix[remove,]%in%c(3,4))]
      adj_matrix_imsi[which(adj_matrix[,remove]%in%c(2,4)),nbrs] <- adj_matrix[remove,which(adj_matrix[,remove]%in%c(2,4))]
    }
    # change detail part negatives
    adj_matrix_imsi[remove,] <- -abs(adj_matrix_imsi[remove,])
    adj_matrix_imsi[,remove] <- -abs(adj_matrix_imsi[,remove])
  }
  return(list(coarse_line=coarse_line, detail_line=detail_line, coarse_lengths=coarse_lengths, detail_lengths=detail_lengths, adj_matrix=adj_matrix_imsi, type=type))
}




forward_line_Y <- function(line_remove, line_coremove, line_nbrs, lengths_remove, lengths_coremove, lengths_nbrs, adj_matrix, remove, coremove, nbrs){
  # constructed for the Y shape
  
  
  #plot(rbind(line_remove, line_coremove, line_nbrs), ylim=c(36.57288-(0.04621/2), 36.60227+(0.04621/2)))
  #segments(x0=line_remove[1,1], y0=line_remove[1,2], x1=line_remove[2,1], y1=line_remove[2,2], col="red")
  #segments(x0=line_coremove[1,1], y0=line_coremove[1,2], x1=line_coremove[2,1], y1=line_coremove[2,2], col="green")
  #segments(x0=line_nbrs[1,1], y0=line_nbrs[1,2], x1=line_nbrs[2,1], y1=line_nbrs[2,2], col="blue")
  
  # 1st removal
  newpt_imsi <- -(-(line_remove[1,]-line_remove[2,])-(line_coremove[1,]-line_coremove[2,]))+line_remove[2,]
  
  result_line_remove <- -line_remove
  result_line_coremove <- rbind(newpt_imsi, line_coremove[2,])
  
  #adjacency matrix update
  adj_matrix[remove, ] <- -abs(adj_matrix[remove,])
  adj_matrix[,remove] <- -abs(adj_matrix[,remove])
  
  #update length
  lengths_remove <- lengths_remove - lengths_coremove
  lengths_coremove <- lengths_remove + lengths_coremove
  
  #2nd removal
  result_new <- forward_line(line_remove=result_line_coremove, line_nbrs=line_nbrs, lengths_remove=lengths_coremove, lengths_nbrs=lengths_nbrs, adj_matrix, remove=coremove, nbrs=nbrs, type="b1")
  
  type <- c("c1")
  
  #points(newpt_imsi[1], newpt_imsi[2], col=5)
  
  #test_a <- rbind(coeff_lines[[nbrs]][1,], coeff_lines[[remove]][2,])
  #test_b <- coeff_lines[[remove]][1,]
  #project(test_b , (test_a[2,] -test_a[1,]))
  #
  #coeff_lines[[nbrs]][1,] + find_orthogonal(x=(test_a[2,] -test_a[1,]) , y= (coeff_lines[[nbrs]][2,]-coeff_lines[[nbrs]][1,]) )
  
  #plot(rbind(test_a,test_b), ylim=c(36.581, 36.609), xlim=c(127.292, 127.310))
  #segments(x0=test_a[1,1], y0=test_a[1,2], x1=test_a[2,1], y1=test_a[2,2])
  #segments(x0=coeff_lines[[nbrs]][1,1], y0=coeff_lines[[nbrs]][1,2], x1=coeff_lines[[nbrs]][2,1], y1=coeff_lines[[nbrs]][2,2], col="blue")
  #segments(x0=coeff_lines[[remove]][1,1], y0=coeff_lines[[remove]][1,2], x1=coeff_lines[[remove]][2,1], y1=coeff_lines[[remove]][2,2], col="red")
  #
  #points(127.30154,  36.58159, col="green", pch=16)
  #points(127.3016, 36.60801, col="green", pch=16)
  
  
  return(list(coarse_line=result_new$coarse_line, detail_line=result_new$detail_line, coarse_lengths=result_new$coarse_lengths, detail_lengths=result_new$detail_lengths, adj_matrix=result_new$adj_matrix, type=type))
}


########################################
##CODES related to graphics
########################################
#reference: https://stackoverflow.com/questions/20127282/r-color-scatterplot-points-by-z-value-with-legend
scatter_fill <- function (x, y, z,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),zlim=c(min(z),max(z)),
                          nlevels = 121, plot.title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", 
                          yaxs = "i", las = 1, y.axes.label=TRUE,
                          axes = TRUE, y.axes=TRUE, frame.plot = axes, smallplot=c(0.85,0.9,0.25,0.75), coldefault=colorRampPalette(c("cyan", "green", "yellow", "red", "black")), plot.legend=TRUE, ...) 
{
  old.par <- par(no.readonly = TRUE)
  bigpar <- old.par$plt
  
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #on.exit(par(par.orig))
  #w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  w <- (3 + mar.orig[2L]) * par("csi") * 2.24
  #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  
  #zlim=c(min(z, na.rm=T), max(z, na.rm=T))
  # choose colors to interpolate
  levels <- seq(zlim[1],zlim[2],length.out = nlevels)
  #col <- colorRampPalette(c("red","yellow","dark green"))(nlevels) 
  #following ODonnell's Approach
  col<-coldefault(nlevels)
  colz <- col[cut(z,nlevels)]  
  
  par(bigpar)
  mar <- mar.orig
  #mar[4L] <- 1
  #mar[4L] <- 0.1
  par(mar = mar)
  
  par(las = 0)
  # points
  plot(x,y,type = "n",xaxt='n',yaxt='n',xlab="",ylab="",xlim=xlim,ylim=ylim,bty="n")
  points(x,y,col = colz,xaxt='n',yaxt='n',xlab="",ylab="",bty="n",...)
  
  ## options to make mapping more customizable
  
  if (missing(plot.axes)) {
    if (axes & y.axes==TRUE) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1, labels = FALSE)
      Axis(y, side = 2, labels = FALSE)
    }else if (axes & y.axes==FALSE & y.axes.label==TRUE) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      #Axis(y, side = 2)
    }else if(axes & y.axes==FALSE & y.axes.label==FALSE){
      title(main = "", xlab = "", ylab = "")
      Axis(x, side=1, labels=FALSE, tick=FALSE)
      Axis(y, side=2, labels=FALSE, tick=FALSE)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  
  if(plot.legend==TRUE){
    #par(las = las)
    par(las = 0)
    mar <- mar.orig
    #mar[4L] <- mar[2L]
    #mar[2L] <- 1
    #mar[2L] <- 0.1
    par(mar = mar)
    #   
    #plot.new()
    if(!is.null(smallplot)){
      par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    }
    
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    
    par(las = 2)
    rect(0, levels[-length(levels)], 1, levels[-1L],col=col,border=col) 
    if (missing(key.axes)) {if (axes){axis(4)}}
    else key.axes
    box()
    if (!missing(key.title)) 
      key.title
  }
  
  par(new = TRUE, pty = "m", plt = bigpar, err = -1)
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  par(bigpar)
  mar <- mar.orig
  #mar[4L] <- 1
  #mar[4L] <- 0.1
  par(mar = mar)
  
  par(las = 0)
  plot(x,y,type = "n",xaxt='n',yaxt='n',xlab="",ylab="",xlim=xlim,ylim=ylim,bty="n")
  
  invisible()
}


#compute shreve stream order
compute_flow_vol <- function(adjacency, example_network, scalevec=c(0.2, 1.5), logdata=FALSE){
  #type="lengthprop": shreve_and_dist[first_stream] is defined by proportional to their lengths
  #type="equal": shreve_and_dist[first_stream] is defined by equal weight (== shreve dist style)
  #scalevec: weight scaling (default = (0.2,1.5)  (refer to O'Donnell method))
  
  #adjacency: adjacency object obtained from get_adjacency_stream
  adj_mat <- as.matrix(adjacency$adjacency)
  
  shreve <- rep(0, nrow(adj_mat))
  
  shreve_and_dist <- rep(0, nrow(adj_mat))
  
  first_stream <- which(colSums(adj_mat)==0)
  shreve[first_stream] <- 1
  shreve_and_dist[first_stream] <- example_network@data$Length[first_stream]
  
  pointsin <- c(1:nrow(adj_mat))
  pointsin <- setdiff(pointsin, first_stream)
  
  while(length(pointsin)!=0){
    #find next stream
    #case 1: Y shape
    next_stream1 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==2 ))]
    for(i in 1:length(next_stream1)){
      shreve[next_stream1[i]] <- sum(shreve[which(adj_mat[,next_stream1[i]]!=0)])
      shreve_and_dist[next_stream1[i]] <- sum(shreve_and_dist[which(adj_mat[,next_stream1[i]]!=0)])
    }
    #case 2: 1 shape
    next_stream2 <- pointsin[which(sapply(pointsin, function(x) sum(adj_mat[first_stream,x])==1 ) & sapply(pointsin, function(x) sum(adj_mat[,x])==1 ))]
    for(j in 1:length(next_stream2)){
      shreve[next_stream2[j]] <- sum(shreve[which(adj_mat[,next_stream2[j]]!=0)])
      shreve_and_dist[next_stream2[j]] <- sum(shreve_and_dist[which(adj_mat[,next_stream2[j]]!=0)])
    }
    next_stream <- c(next_stream1, next_stream2)
    pointsin <- setdiff(pointsin, next_stream)
    first_stream <- c(first_stream, next_stream)
    #print(pointsin)
  }
  
  #scaling
  #weight_vec_candidate2 <- scales::rescale(log(weight_vec_candidate2$distweight), to=c(0.2, 1.5))
  # if(logdata==TRUE){
  #   shreve_and_dist <- scales::rescale(log(sqrt(shreve_and_dist)), to=c(scalevec[1], scalevec[2]))
  #   #shreve_and_dist <- scales::rescale(log(shreve_and_dist), to=c(scalevec[1], scalevec[2]))
  # }else{
  #   shreve_and_dist <- scales::rescale((sqrt(shreve_and_dist)), to=c(scalevec[1], scalevec[2]))
  # }
  return(list(shreve=shreve, flow=shreve_and_dist))
}

## =============================================================================================== ##

####### functions by myself ###########
# These functions are basically constructed by modifying existing functions of the "GWPCA" and "stpca" packages
# for the purpose of analyzing river network data. 

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

plot_river <- function(x, plots = c("score", "biplot", "glyph", "info", "riverpca"), 
                       pc.map = 1, probs = c(0, 0.25, 0.75, 1), 
                       pc.biplot = c(1,2), biplot.factor = 100, 
                       pc.glyph = c(1, 2, 3), river.glyph.loading = c(1,2,3), river.glyph.score = c(1,2,3), river, r1 = 10, 
                       coords, pc.ts = c(1, 2, 3), pc.meanPM = c(1, 2,3), 
                       meanPM.factor = 0.05, lwd=2, zoom=FALSE, zoomx=c(0,1), zoomy=c(0,1)){
  pca.mode <- x$pca.mode
  if(zoom==FALSE){
    zoomx <- c(river@bbox[1, ])
    zoomy <- c(river@bbox[2, ])
  }
  if (pca.mode == "Smode"){
    message("For river PCA, only support T-mode PCA")
    return()
  }
  else{
  if ("biplot" %in% plots) {
    nplot <- length(x$pca.wt)
    if (nplot == 4) {
      par(mfrow = c(2, 2),oma = c(2,2,2,2) + 0.1, mar = c(4,4,2,1) + 0.1)
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
  if ("glyph" %in% plots) {
    if (length(pc.glyph) < 3) {
      message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
    }
    nplot <- length(x$pca.wt)
    if (nplot == 4) {
      par(mfrow = c(2, 2))
    }
    else {
      par(mfrow = c(nplot, 1))
    }
    pca.type <- x$pca.wt
    for (i in 1:nplot) {
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
  if ("score" %in% plots) {
    nplot <- length(x$pca.wt)
    pca.type <- x$pca.wt
    if (nplot == 4) {
      par(mfrow = c(2, 2))
    }
    for (i in 1:nplot) {
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
  if ("info" %in% plots) {
    info.num <- length(x$info)
    info.type <- x$info
    for (i in 1:info.num) {
      title <- paste("Plot for", info.type[i])
      map.data <- x[[info.type[i]]]
      if(info.type[i] == "PTV"){
        scores <- x[[info.type[i]]]
        n.colors <- length(probs) - 1
        data.vec <- sort(scores)
        quants <- quantile(data.vec, probs = probs)
      }
      else{
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
  if ("riverpca" %in% plots) {
    if (length(river.glyph.loading) < 3) {
      message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
    }
    nplot <- length(river.glyph.loading)
    if (nplot == 4) {
      par(mfrow = c(2, 2))
    }
    else {
      par(mfrow = c(1,nplot), oma = c(2,2,2,2) + 0.1, mar = c(4,5,2,2) + 0.1)
    }
    pca.type <- "riverpca"
    for (i in 1:nplot) {
      glyph.data <- x[[pca.type]]$loadings
      local.loadings <- glyph.data[,,river.glyph.loading[i]]
      # scores <- glyph.data$scores[, pc.glyph]
      plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
           ylab = names(coords)[2], main = NULL, cex.main = 1.5, cex.lab = 1.5, 
           cex.axis = 1.5, xlim = zoomx, 
           ylim = zoomy)
      lines(river, col = "black")
      glyph.plot.river(as.matrix(local.loadings), coords, 
                       r1 = r1, add = TRUE, lwd=lwd)
      if(zoom==FALSE){
      glyph.plot.river(rbind(t(x$spatiotemporal$loads[,i]),0), rbind(c(min(coords[,1])+0.1,max(coords[,2])-0.05),0), 
                       r1 = 40*r1, add = TRUE, lwd=lwd) # manually control
      text(min(coords[,1])+0.1,max(coords[,2])-0.2,expression(TPCA[ST]), cex=lwd, font=2) # manually control
      }
    }
    message("For plots = glyph: did you check that coordinates for monitoring site ID's are in \nthe same order as data used for stpca()?  \nIf not then glyphs might not correspond to the correct monitoring site.")
    
    if (length(river.glyph.score) < 3) {
      message("For plots = glyph: glyph should be used to display 3 or \n more principal components simultaneously.")
    }
    
    glyph.data <- x[[pca.type]]$scores
    local.scores <- glyph.data[,river.glyph.score]
    # scores <- glyph.data$scores[, pc.glyph]
    
    par(mfrow = c(1,2), oma = c(2,2,2,2) + 0.1, mar = c(4,5,2,2) + 0.1)
    
    plot(coords[, 1], coords[, 2], type = "n", xlab = names(coords)[1], 
         ylab = names(coords)[2], main = NULL, cex.main = 1.5, cex.lab = 1.5, 
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
         ylab = names(coords)[2], main = NULL, cex.main = 1.5, cex.lab = 1.5, 
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
                                             min(loc[, 1]), max(loc[, 2]) - min(loc[, 2]))/r1, add = FALSE, 
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


flury.hierarchy = function(data, covmats, nvec, n.var, cluster, alpha=0.05, mode="FG"){
  # flury's AIC
  if(mode == "FG"){
    g = cpc::FG
  }else if(mode == "stepwise"){
    g = cpc::stepwisecpc
  }
  
  varnames = colnames(covmats[,,1]) # variable names
  no.groups = dim(covmats)[3] # number of groups
  
  B.cpc = g(covmats = covmats, nvec = nvec)$B # common eigenvectors under CPC model
  rownames(B.cpc) = varnames
  colnames(B.cpc) = paste("PC",c(1:length(varnames)),sep="")
  common_order = cpc::findcpc(covmats=covmats, B=B.cpc, # common eigenvectors order for cpcq model 
                              plotting=FALSE)$commonvec.order
  
  flury_res = cpc::flury.test(covmats=covmats, nvec=nvec, B=B.cpc, # flury's hierarchy model res 
                              p = n.var, qmax = n.var - 2, 
                              commonvec.order = common_order)
  flury_res = cbind(flury_res, rep(0, n.var+2), rep(0, n.var+2), rep(0, n.var+2))
  colnames(flury_res)[7] = "total.chi.square"
  colnames(flury_res)[8] = "total.df"
  colnames(flury_res)[9] = "criterion"
  
  res.equal = cpc::equal.test(covmats=covmats, nvec=nvec) # equality test
  
  equal = 0 # if equality model holds, becomes 1
  prop = 0 # if proportionality model holds, becomes 1
  cpc = 0 # if cpc model holds, becomes 1
  
  res = list()
  res[[1]] = equal ; res[[2]] = prop ; res[[3]] = cpc
  res[[4]] = B.cpc ; res[[5]] = common_order ; res[[6]] = flury_res  
  names(res) = c("equal", "prop", "cpc", "B.cpc", "common_order", "flury_res")
  flury_res[1,7] = res.equal$chi.square
  flury_res[1,8] = res.equal$df
  flury_res[1,9] = qchisq(1-alpha, res.equal$df)
  
  cat("####### Hypothesis Testing #######", end= "\n")
  if(res.equal$chi.square > 
     qchisq(1-alpha, res.equal$df)){
    cat("1 Reject Equality Model", end="\n")
    res.prop = cpc::prop.test(covmats=covmats, nvec=nvec) # proportionality test
    
    flury_res[2,7] = res.prop$chi.square
    flury_res[2,8] = res.prop$df
    flury_res[2,9] = qchisq(1-alpha, res.prop$df)
    
    if(res.prop$chi.square >
       qchisq(1-alpha, res.prop$df)){
      cat("2 Reject Proportionality Model", end="\n")
      res.cpc = cpc::cpc.test(covmats=covmats, nvec=nvec) # cpc test
      
      flury_res[3,7] = res.cpc$chi.square
      flury_res[3,8] = res.cpc$df
      flury_res[3,9] = qchisq(1-alpha, res.cpc$df)
      
      if(res.cpc$chi.square >
         qchisq(1-alpha, res.cpc$df)){
        cat("3 Reject CPC Model", end="\n")
        
        q.common = 0
        if(n.var > 2){
          for(q in (n.var-2):1){
            res.cpcq = cpc::cpcq.test(covmats = covmats, nvec = nvec, 
                                      B = B.cpc[,common_order], q = q) # cpcq test
            
            flury_res[(n.var+2-q),7] = res.cpcq$chi.square
            flury_res[(n.var+2-q),8] = res.cpcq$df
            flury_res[(n.var+2-q),9] = qchisq(1-alpha, res.cpcq$df)
            
            if(res.cpcq$chi.square > qchisq(1-alpha, res.cpcq$df)){
              cat(paste(2+n.var-q, " Reject CPC(",q,") Model", sep=""), end="\n")
              q.common = 0
            }else{
              cat(paste(2+n.var-q, " Accept CPC(",q,") Model", sep=""), end="\n")
              q.common = q # determine cpc(q) model 
              res[[7]] = q.common ; res[[8]] = res.cpcq$covmats.cpcq
              res[[9]] = B.partial(covmats=covmats, nvec=nvec, B=B.cpc,
                                   commonvec.order = common_order, q=q.common)
              
              score = matrix(0, nrow = sum(nvec), ncol = n.var)
              R.F.partial = array(0, c(n.var, n.var, no.groups), dimnames = list(c(),c(),paste("Group",1:no.groups, sep = "")))
              evptv = array(0, c(2, n.var, no.groups), 
                            dimnames = list(c("Eigenvalues", "Percentage of total variation"),
                                            paste("PC",c(1:n.var),sep=""),paste("Group",1:no.groups, sep = "")))
              for(j in 1:no.groups){
                Fj = t(res[[9]][,,j])%*%covmats[,,j]%*%res[[9]][,,j] # covariance of score 
                lambda.inv.sq.rt = diag(1/sqrt(diag(Fj)))
                Rj = lambda.inv.sq.rt%*%Fj%*%lambda.inv.sq.rt # correlation of score
                
                R.F.partial[,,j] = Fj
                R.F.partial[,,j][lower.tri(R.F.partial[,,j])] = Rj[lower.tri(Rj)]
                
                score[cluster==j,] = data[cluster==j,]%*%res[[9]][,,j]
                
                evptv[,,j][1,] = diag(Fj)
                evptv[,,j][2,] = diag(Fj)/sum(diag(Fj))*100
              }
              
              score = cbind(score, cluster)
              colnames(score) = c(paste("PC",c(1:n.var),sep=""), "cluster")
              res[[10]] = R.F.partial
              res[[11]] = score
              res[[12]] = evptv
              
              names(res)[c(7:12)] = c("q.common", "covmats.cpcq", "B.partial", "R.F.partial", "scores","EVPTV")
              dimnames(res[[8]]) = list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
              dimnames(res[[9]]) = list(varnames,paste("PC",c(1:n.var),sep=""),paste("Group",1:no.groups, sep = ""))
              
              break
            }
          }
        }
        
      }else{
        cat("3 Accept CPC Model", end="\n")
        cpc = 1
        res[[3]] = cpc
        res[[7]] = res.cpc$covmats.cpc
        names(res)[7] = "covmats.cpc"
        dimnames(res[[7]]) = list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
        
        score = data%*%B.cpc
        R.F.cpc = array(0, c(n.var, n.var, no.groups), dimnames = list(c(),c(),paste("Group",1:no.groups, sep = "")))
        evptv = array(0, c(2, n.var, no.groups), 
                      dimnames = list(c("Eigenvalues", "Percentage of total variation"),
                                      paste("PC",c(1:n.var),sep=""),paste("Group",1:no.groups, sep = "")))
        for(j in 1:no.groups){
          Fj = t(B.cpc)%*%covmats[,,j]%*%B.cpc # covariance of score 
          lambda.inv.sq.rt = diag(1/sqrt(diag(Fj)))
          Rj = lambda.inv.sq.rt%*%Fj%*%lambda.inv.sq.rt # correlation of score
          
          R.F.cpc[,,j] = Fj
          R.F.cpc[,,j][lower.tri(R.F.cpc[,,j])] = Rj[lower.tri(Rj)]
          
          evptv[,,j][1,] = diag(Fj)
          evptv[,,j][2,] = diag(Fj)/sum(diag(Fj))*100
        }
        score = cbind(score, cluster)
        colnames(score) = c(paste("PC",c(1:n.var),sep=""), "cluster")
        res[[8]] = R.F.cpc
        res[[9]] = score
        res[[10]] = evptv
        
        
        names(res)[c(7:10)] = c("covmats.cpc", "R.F.cpc", "scores", "EVPTV")
      }
    }else{
      cat("2 Accept Proportionality Model", end="\n")
      prop = 1
      res[[2]] = prop
      res[[7]] = res.prop$covmats.prop
      names(res)[7] = "covmats.prop"
      dimnames(res[[7]]) = list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
    }
  }else{
    cat("1 Accept Equality Model", end="\n")
    equal = 1
    res[[1]] = equal
    res[[7]] = res.equal$covmats.equal
    names(res)[7] = "covmats.equal"
    dimnames(res[[7]]) = list(varnames,varnames,paste("Group",1:no.groups, sep = ""))
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

