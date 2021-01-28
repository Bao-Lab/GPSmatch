#' @title Calculate significance of bed file similarity
#' @description The function compare two bed files (query and hit) to check if the observed degree of similarity between the query and the hit files may occur purely by chance.
#' @param n The number of background files generated in order to compute the p-value. As the n increases, the p-value will become more reliable, but the user should be aware that this will significantly increase the computing time. We have set a default n of 100.
#' @param bed1 The file path of a query bed file to be compared to the database files.
#' @param genome The file path of a genome file, which should be tab delimited and structured as follows: <chromName><TAB><chromSize>. A pre-formatted hg19 genome file can be found on the Github.
#' @param bed2 The file path of a bedfile in the database, or any bed file.
#' @param output_path The file path where the output .csv file will be generated
#' @export
#' @return A .csv file showing the p-value, pi-score, and a five number summary of the background file distribution compared to the original jaccard index.
#' @examples
#' indivPVal(100,"/dir/bed1.txt","/dir/genome.txt","/dir/folder_dir","/dir/output_folder")

indivPVal = function(n=100, bed1, genome, bed2,output_path){
  #3levels = list.files(folder_dir)
  final_values = list()

  #generate background files from query:
  background_dir = paste(bed1,"background",sep="_")
  dir.create(background_dir)
  print("generating background files")
  background_list = list()
  for (i in 1:(n)){
    background <- fread(cmd = paste("bedtools shuffle -i", bed1, "-g", genome))
    write.table(background, paste(background_dir,"/background_file",i,".txt",sep = ""),sep="\t",row.names=F, col.names = F, quote = F)
  }
  print("finished generating background files")


  ##generates real jaccard for bed1 and bed2
  jaccard_list = jaccardCalc(bed1, bed2)
  num_overlap = jaccard_list[1]
  num1=jaccard_list[2]
  num2=jaccard_list[3]

  num_total = num1 + num2 - num_overlap
  real_jaccard_id = num_overlap/num_total

  ######## shuffled files and bed2
  background_n = list.files(background_dir)
  background_jaccard = list()
  for (j in 1:length(background_n)){
      #################################################watch
      bed1shuffle= paste(background_dir,background_n[j],sep="/")

      jaccard_list = jaccardCalc(bed1shuffle, bed2)
      num_overlap = jaccard_list[1]
      num1=jaccard_list[2]
      num2=jaccard_list[3]

      num_total = num1 + num2 - num_overlap
      shuffle_jaccard_id = num_overlap/num_total

      #jaccard_id = jaccard(paste(background_dir,background_n[j],sep="/"),pwd)
      background_jaccard[[j]] = c(bed1shuffle,jaccard_index = shuffle_jaccard_id)

      # print(j)
    }

    jaccard_df = as.data.frame(do.call(rbind,background_jaccard))

    bed = as.numeric(as.character(jaccard_df[,2]))

    #find the p-value
    num = as.numeric(as.character(match(real_jaccard_id,sort(c(bed,real_jaccard_id),decreasing=T))))
    value = num/(n+1)

    #finding the pi-score = mean ratio * (-log10(p-value))
    if (mean(bed) == 0){
      piscore = (real_jaccard_id)/(0.000001) * (-log(value))
    }
    else{
      piscore = (real_jaccard_id)/(mean(bed)) * (-log(value))
    }



    #final_values[[i]] = c(bed1,jaccard_index = real_jaccard_id,pi_score = piscore, p_value = value,summary(bed))
    final_values = c(bed1,bed2,jaccard_index = real_jaccard_id,pi_score = piscore, p_value = value,summary(bed))
    names(final_values) <- c("query", "hit", "Jaccard index","pi_score", "p-value", "Shuffle_Min jaccard", "Shuffle_1st Quarter jaccard", "Shuffle_Median jaccard", "Shuffle_Mean jaccard", "Shuffle_3rd Quarter jaccard","Shuffle_Max jaccard")
  unlink(background_dir, recursive = TRUE)


  write.csv(final_values,paste(output_path,"/",gsub("^.*/", "", bed1), "_",gsub("^.*/", "", bed2),"_pvalue_nt.csv",sep=""))
  print(paste(output_path,"/",gsub("^.*/", "", bed1), "_",gsub("^.*/", "", bed2),"_pvalue.csv",sep=""))
  print("file created!")
}

