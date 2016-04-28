naindex = which(is.na(geno),arr.ind = TRUE)
naindex = naindex[order(naindex[,2]),]
for(i in 1:dim(naindex)[1] ){
  row = naindex[i,1]
  col = naindex[i,2]
  if(col==2){
    geno[row,col] = geno[row,col+1]
  }else if(col==2058){
    geno[row,col] = geno[row,col-1]
  }else{
    left = abs(pmap$pos[col-1]-pmap$pos[col-2])
    right = abs(pmap$pos[col-1]-pmap$pos[col])
    if(left > right){
      geno[row,col] = geno[row, col-1]
    }else{
      geno[row,col] = geno[row, col+1]
    }
  }
}