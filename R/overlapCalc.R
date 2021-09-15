#' @title Calculates nucleotide overlap only
#' @description  The function calculates the nucleotide overlap between two bed files (query and hit).
#' @param bed1 The file path of a query bed file to be compared to the database files.
#' @param bed2 The file path of a bed file in the database files to be compared to the query file.
#' @export
#' @return A value of the nucleotide overlap. To be used as an internal function for PValue_nt and GPSMatch_nt.
#' @examples
#' overlapCalc("/dir/bed1.txt","/dir/bed2.txt")


### Calculate the overlap length between two bed files
overlapCalc = function(bed1,pwd){
  overlap = fread(cmd = paste("bedtools intersect -wo -a", bed1, "-b",pwd))
  overlap = as.data.frame(overlap)
  if (dim(overlap)[1] == 0) {
    num_overlap =0
  } else{
    num_overlap = sum(overlap[,ncol(overlap)])
  }

  jaccard_list = c(num_overlap)
  return(jaccard_list)
}
