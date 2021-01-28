#' @title Calculates nucleotide overlap and total length
#' @description  The function calculates the nucleotide overlap between two bed files (query and hit) and their respective nucleotide lengths.
#' @param bed1 The file path of a query bed file to be compared to the database files.
#' @param bed2 The file path of a bed file in the database files to be compared to the query file.
#' @export
#' @return A vector of the nucleotide overlap and their respective nucleotide lengths. To be used as an internal function for PValue_nt and GPSMatch_nt.
#' @examples
#' jaccardCalc("/dir/bed1.txt","/dir/folder_dir")

jaccardCalc = function(bed1,pwd){
  overlap = fread(cmd = paste("bedtools intersect -wo -a", bed1, "-b",pwd))
  overlap = as.data.frame(overlap)
  if (dim(overlap)[1] == 0) {
    num_overlap =0
  } else{
    num_overlap = sum(overlap[,ncol(overlap)])
  }

  #calculate total length
  overlapA = fread(cmd = paste("bedtools intersect -a", bed1, "-b",bed1, "-wo"))
  overlapA = as.data.frame(overlapA)
  num1 = sum(overlapA[,ncol(overlapA)])

  #overlapA[order(overlapA$V1,overlapA$V2),]
  #write.table(overlapA, "/Users/xianggao/Desktop/AmyHighschool/BaoNorthwestern/overlapA_0.95_merge_BroadHistoneH1hescH3k27ac.csv", sep = "\t")

  overlapB = fread(cmd = paste("bedtools intersect -a", pwd, "-b",pwd, "-wo"))
  overlapB = as.data.frame(overlapB)
  num2 = sum(overlapB[,ncol(overlapB)])

  jaccard_list = c(num_overlap,num1,num2)
  return(jaccard_list)
}
