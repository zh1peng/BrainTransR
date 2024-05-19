# Function to generate a permutation map from a set of cortical regions of interest to itself, 
# while (approximately) preserving contiguity and hemispheric symmetry.
# The function is based on a rotation of the FreeSurfer projection of coordinates
# of a set of regions of interest on the sphere.
#
# Inputs:
# coord.l       coordinates of left hemisphere regions on the sphere        array of size [n(LH regions) x 3]
# coord.r       coordinates of right hemisphere regions on the sphere       array of size [n(RH regions) x 3]
# nrot          number of rotations (default = 10000)                       scalar
#
# Output:
# perm.id      array of permutations, from set of regions to itself        array of size [n(total regions) x nrot]
#
# required library: matrixStats (for rowMins function)
#
# Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018
#       Updated on 16/10/2019 with permutation scheme that uniformly samples the space of permutations on the sphere
#       See github repo (@frantisekvasa) and references within for details

library(matrixStats)
library(dplyr)
rotate.parcellation.lh = function(coord.l,nrot=10000) {
  # check that coordinate dimensions are correct
  if (!all(dim(coord.l)[2]==3)) {
    if (all(dim(coord.l)[1]==3)) {
      print('transposing coordinates to be of dimension nROI x 3')
      coord.l = t(coord.l)
    }
  }
  
  nroi.l = dim(coord.l)[1]   # n(regions) in the left hemisphere
  nroi = nroi.l      # total n(regions)
  
  perm.id = array(0,dim=c(nroi,nrot)); # initialise output array
  r = 0; c = 0; # count successful (r) and unsuccessful (c) iterations
  
  # UPDATED 16/10/2019 - set up updated permutation scheme 
  I1 = diag(3); I1[1,1] = -1;
  # main loop -  use of "while" is to ensure any rotation that maps to itself is excluded (this is rare, but can happen)
  while (r < nrot) {
    
    # UPDATED 16/10/2019
    A = matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
    qrdec = qr(A)       # QR decomposition
    TL = qr.Q(qrdec)    # Q matrix
    temp = qr.R(qrdec)  # R matrix
    TL = TL%*%diag(sign(diag(temp)))
    if (det(TL)<0) {
      TL[,1] = -TL[,1]
    }
    # reflect across the Y-Z plane for right hemisphere

    coord.l.rot = coord.l %*% TL; # transformed (rotated) left coordinates
    dist.l = array(0,dim=c(nroi.l,nroi.l));
    for (i in 1:nroi.l) { # left
      for (j in 1:nroi.l) {
        dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
      }
    }
   
    
    # LEFT
    # calculate distances, proceed in order of "most distant minimum"
    # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
    # as this region is the hardest to match and would only become harder as other regions are assigned
    temp.dist.l = dist.l
    rot.l = c(); ref.l = c();
    #tba.r = tba.c = 1:nroi.l # rows and columns that are yet "to be assigned"
    for (i in 1:nroi.l) {
      # max(min) (described above)
      ref.ix = which( rowMins(temp.dist.l,na.rm=T) == max(rowMins(temp.dist.l,na.rm=T),na.rm=T) )   # "furthest" row
      rot.ix = which( temp.dist.l[ref.ix,] == min(temp.dist.l[ref.ix,],na.rm=T) ) # closest region
      
      # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
      # ref.ix = which(nanmean(temp.dist.l,2)==nanmax(nanmean(temp.dist.l,2)))    # "furthest" row
      # rot.ix = which(temp.dist.l(ref.ix,:)==nanmin(temp.dist.l(ref.ix,:)))      # closest region    
      ref.l = c(ref.l,ref.ix) # store reference and rotated indices
      rot.l = c(rot.l,rot.ix)
      temp.dist.l[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
      temp.dist.l[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
      #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
      #temp.dist.l[ref.ix,] = NA
    }
    

    # mapping is x->y
    # collate vectors from both hemispheres + sort mapping according to "reference" vector
    b = sort(ref.l,index.return=T); 
    ref.l.sort = ref.l[b$ix]; rot.l.sort = rot.l[b$ix];
    
    # verify that permutation worked (output should be vector with values 1:nroi = 1:(nroi_l+nroi_r))
    if (!all(sort(rot.l.sort,decreasing=F)==c(1:nroi))) {
      #save.image('~/Desktop/perm_error.RData')
      browser("permutation error")
    }
    
    # verify that permutation does not map to itself
    if (!all(rot.l.sort==c(1:nroi))) {
      r = r+1
      perm.id[,r] = rot.l.sort # if it doesn't, store it
    } else {
      c = c+1
      print(paste('map to itself n. ',toString(c),sep=''))
    }
    
    # track progress
    if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
  }
  return(perm.id)
}


create_lh_perm.id <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/perm_id_weights',  
                              atlas=c('desikan',
                                      'schaefer100',
                                      'schaefer200',
                                      'schaefer300'),
                                        nreps=10000,
                              null_type=c('random_brain','spin_brain')){
  atlas=match.arg(atlas)
  null_type=match.arg(null_type)
  n.region <- switch(atlas,
    'desikan'       = 34,
    'schaefer100'   = 50,
    'schaefer200'   = 100,
    'schaefer300'   = 150,
    NA_integer_
  )


  if (null_type=='spin_brain'){
  centroid_csv=sprintf('%s/%s_centroid.csv',data_path,atlas)

  # check if centroid_csv exist
  if (!file.exists(centroid_csv)){
   stop(sprintf('centroid_csv %s does not exist',centroid_csv))
  }

  centroid_data=read.csv(centroid_csv)
  centroid.l=centroid_data %>%  select(starts_with('centroid'))
  perm.id=rotate.parcellation.lh(as.matrix(centroid.l),
                                 nrot=nreps+100)
  }else if (null_type=='random_brain') {
  id.mat=matrix(rep(c(1:n.region),times=nreps+100),nrow = n.region,ncol = nreps+100)
  perm.id=apply(id.mat, MARGIN = 2, function(x) sample(x, replace = FALSE, size = length(x)))
  }

  print(paste0('find duplciated columns: ', sum(duplicated(t(perm.id)))))
  perm.id.nodup=perm.id[, !duplicated(t(perm.id))]
  perm.id.final=perm.id.nodup[,1:nreps]
  return(perm.id.final)
}




# ===================== generate perm.id =====================
data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/perm_id_weights'

for (atlas_i in c(
        'schaefer300')){

for (null_type in c('random_brain','spin_brain')){
perm.id=create_lh_perm.id(data_path=data_path, 
                         atlas=atlas_i,
                         nreps=10000,
                         null_type=null_type)
savename=sprintf('%s/%s_%s_perm_id.rds',data_path,atlas_i,null_type)
saveRDS(perm.id, savename)
 }
}


# ===================== generate Moran_w for each atlas =====================
for (atlas_i in c('desikan',
        'schaefer100',
        'schaefer200',
        'schaefer300')){

 centroid_csv=sprintf('%s/%s_centroid.csv',data_path,atlas_i)

  # check if centroid_csv exist
  if (!file.exists(centroid_csv)){
   stop(sprintf('centroid_csv %s does not exist',centroid_csv))
  }

 centroid_data=read.csv(centroid_csv)
 centroid.l=centroid_data %>%select(starts_with('centroid'))
 rownames(centroid.l)=centroid_data$Row
 dist.mat=as.matrix(dist(centroid.l))
 w=1/dist.mat
 diag(w)= 0
 savename=sprintf('%s/%s_Moran_w.rds',data_path,atlas_i)
 saveRDS(w, savename)
}