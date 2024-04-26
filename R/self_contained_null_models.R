load_permid <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/perm_id_weights',
                             atlas=c('desikan',
                                    'schaefer100',
                                    'schaefer200',
                                    'schaefer300'),
                             type=c('spin_brain','random_brain'),
                             perm.n=5000){
# load pre-generated perm.id
atlas=match.arg(atlas)
type=match.arg(type)

if (perm.n>10000){stop('>10k perm.id is supported')} 

 perm.id.file=sprintf('%s/%s_%s_perm_id.rds',data_path, atlas,type)
 
 # check if perm.id.file exist
  if (!file.exists(perm.id.file)){
    stop(sprintf('perm.id file %s does not exist',perm.id.file))
  }

 perm.id=readRDS(perm.id.file)   
  return(perm.id[,1:perm.n])
}

generate_null_brain_data <- function(brain_data,perm.id){
  # before put perm.id in the function 
  # need to make sure no duplicate in perm.id.
  region.n=dim(perm.id)[1]
  perm.n=dim(perm.id)[2]
  if (!dim(brain_data)[1]==region.n){stop('The number of regions in brain data nd perm.id is not matched ')}
  # null_brain_data=matrix(NA,nrow=region.n, ncol=perm.n)
  # for (idx in c(1:perm.n)){
  #   tmp_brain_data=brain_data
  #   null_brain_data[,idx]=brain_data[perm.id[,idx]]
  # }
  null_brain_data=sapply(c(1:perm.n),function(idx){
    brain_data[perm.id[,idx]]
  })
  rownames(null_brain_data)=rownames(brain_data)
  colnames(null_brain_data)=paste0('null_',c(1:perm.n))
  return(null_brain_data)
}
