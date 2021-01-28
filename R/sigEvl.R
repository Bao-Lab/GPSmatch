#' @title Computes the p-value and pi-scores for top ten hits from GPSmatch_nt()
#' @description The function evaluates the significance of the jaccard scores by calculating p-values and pi-scores.
#' @param n The number of background files generated in order to compute the p-value. As the n increases, the p-value will become more reliable, but the user should be aware that this will significantly increase the computing time. We have set a default n of 100.
#' @param genome The file path of a genome file, which should be tab delimited and structured as follows: <chromName><TAB><chromSize>. A pre-formatted hg19 genome file can be found on the Github.
#' @param bed1 The file path of a query bed file to be compared to the database files.
#' @param bed1_jaccard_output The file path of the output file generated from the GPSmatch_nt() function
#' @param folder_dir The directory of a folder containing database files to be used for comparison with the query file.
#' @param output_path The file path where the output .csv file will be generated.
#' @export
#' @return A dataframe that shows the similarities of the query file to the database files ranked from greatest to least, tailored to user method specification.
#' @examples
#' rankBedSimilarity(100, "/dir/genome-file.txt", "/dir/bed1.txt","/dir/bed1_jaccard_output.txt","/dir/database_dir","/dir/output_path")

sigEvl = function(n,genome,bed1,bed1_jaccard_output,database_dir, output_path){
  final_values = list()
  indexes = list()
  file = read.csv(bed1_jaccard_output)
  for (i in 1:10){
    #for each of the top hit files in folder_dir calculate original input file Jaccard index
    hitfile = paste(database_dir,file[i,2],sep = "/")

    ##generates, A and B... then A' B
    real_jaccard = jaccardCalc(bed1, hitfile)
    num_overlap = real_jaccard[1]
    num1=real_jaccard[2]
    num2=real_jaccard[3]

    #jaccard index
    num_total = num1 + num2 - num_overlap
    jaccard_id = num_overlap/num_total

    #percentage --- ANB/A and ANB/B
    #ANB/A
    A_similarity = (num_overlap)/num1
    B_similarity = (num_overlap)/num2
    #Summary
    indexes[[i]] = c(bedfile = as.character(file[i,2]),jaccard_index = jaccard_id,percentage_A = A_similarity,percentage_B=B_similarity)


    #generate background files from query:
    background_dir = paste(database_dir,"background",sep="_")
    dir.create(background_dir)
    print(paste("generating background files for", file[i,2]))
    background_list = list()

    for (j in 1:(n)){
      background <- fread(cmd = paste("bedtools shuffle -i", bed1, "-g", genome))
      write.table(background, paste(background_dir,"/background_file",j,".txt",sep = ""),sep="\t",row.names=F, col.names = F, quote = F)
    }
    print("finished generating background files")

    ######## Calculate background jaccard
    background_n = list.files(background_dir)
    background_jaccard = list()
    for (j in 1:length(background_n)){
      query_background = paste(background_dir,background_n[j],sep="/")
      jaccard_list = jaccardCalc(query_background,hitfile)
      num_overlap = jaccard_list[1]
      num1=jaccard_list[2]
      num2=jaccard_list[3]

      #jaccard index
      num_total = num1 + num2 - num_overlap
      jaccard_id = num_overlap/num_total
      background_jaccard[[j]] = c(bedfile = as.character(file[i,2]),jaccard_index = jaccard_id)
    }

    jaccard_df = as.data.frame(do.call(rbind,background_jaccard))
    bed = as.numeric(as.character(jaccard_df[,2]))   ## jacard scores of shuffled background

    #find the p-value
    num = as.numeric(as.character(match(as.numeric(as.character(indexes[[i]][2])),sort(c(bed,as.numeric(as.character(indexes[[i]][2]))),decreasing=T))))
    value = num/(n+1)

    #finding the pi-score = mean ratio * (-log10(p-value))
    if (mean(bed) == 0){
      piscore = (as.numeric(as.character(indexes[[i]][2])))/(0.000001) * (-log(value))
    } else{
      piscore = (as.numeric(as.character(indexes[[i]][2])))/(mean(bed)) * (-log(value))
    }
    final_values[[i]] = c(bedfile = as.character(file[i,2]),jaccard_index = as.numeric(as.character(indexes[[i]][2])),pi_score = piscore, p_value = value, percentage_A = as.numeric(as.character(indexes[[i]][3])), percentage_B = as.numeric(as.character(indexes[[i]][4])),summary(bed))

    unlink(background_dir, recursive = TRUE)
  }
  values_df = as.data.frame(do.call(rbind,final_values))
  values_df = values_df[order(as.numeric(as.character(values_df$jaccard_index)),decreasing = T),]
  write.csv(values_df,paste(output_path,"/",gsub("^.*/", "", bed1),"_jaccard_pval.csv",sep=""))
}

