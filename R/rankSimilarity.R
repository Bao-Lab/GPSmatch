#' @title Ranks the similarities between a given query bed file with multiple bed files (database files) based on the nucleotides length
#' @description The function compares a given query bed file with multiple bed files (database files), and ranks the relative similarity between each file pairing, computed using jaccard indexes, where 0 has no similarities and 1 has an identical file. The function can also provide a significance value (p-value) of the similarity index based on user selection.
#' @param bed1 The file path of a query bed file to be compared to the database files.
#' @param database_dir The directory of a folder containing database files to be used for comparison with the query file.
#' @param output_path The file path where the output .csv file will be generated.
#' @export
#' @return A dataframe that shows the similarities of the query file to the database files ranked from greatest to least.
#' @examples
#' rankSimilarity("/dir/bed1.txt","/dir/database_dir","/dir/output_path")

rankSimilarity = function(bed1,database_dir,output_path){
  print("rankSimilarity start")

  levels = list.files(database_dir)

  indexes = list()
  for (i in 1:length(levels)){
    ##jaccard
    pwd = paste(database_dir,levels[i],sep="/")

    jaccard_list = jaccardCalc(bed1, pwd)
    num_overlap = jaccard_list[1]
    num1=jaccard_list[2]
    num2=jaccard_list[3]

    #jaccard index
    num_total = num1 + num2 - num_overlap
    jaccard_id = num_overlap/num_total
    #jaccard_index

    #percentage --- ANB/A and ANB/B
    #ANB/A
    A_similarity = (num_overlap)/num1
    B_similarity = (num_overlap)/num2

    ##calc jaccard
    indexes[[i]] = c(bedfile = levels[i],jaccard_index = jaccard_id,percentage_A = A_similarity,percentage_B=B_similarity)
  }
  jaccard_df = as.data.frame(do.call(rbind,indexes))
  jaccard_df = jaccard_df[order(as.numeric(as.character(jaccard_df$jaccard_index)), decreasing = T),]

  write.csv(jaccard_df,paste(output_path,"/",gsub("^.*/", "", bed1),"_jaccard_nt.csv",sep=""))
  print(paste(output_path,"/",gsub("^.*/", "", bed1),"_jaccard_nt.csv",sep=""))
  print("file created!")
}
