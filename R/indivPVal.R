#' @title Calculate significance of bed file similarity
#' @description The function compare two bed files (query and hit) to check if the observed degree of similarity between the query and the hit files may occur purely by chance.
#' @param n The number of background files generated in order to compute the p-value. As the n increases, the p-value will become more reliable, but the user should be aware that this will significantly increase the computing time. We have set a default n of 2000.
#' @param bed1 The file path of a query bed file to be compared to the database files.
#' @param genome The file path of a genome file, which should be tab delimited and structured as follows: <chromName><TAB><chromSize>. A pre-formatted hg19 genome file can be found on the Github.
#' @param bed2 The file path of a bedfile in the database, or any bed file.
#' @param output_path The file path where the output .csv file will be generated
#' @export
#' @return A .csv file showing the p-value, pi-score, and a five number summary of the background file distribution compared to the original jaccard index.
#' @examples
#' indivPVal(n=n,bed1="/dir/bed1.txt",genome ="/dir/genome.txt",bed2 ="/dir/bed2.txt", output_path = "/dir/output_folder")

indivPVal = function(n=2000, bed1= bed1, genome =genome, bed2 = bed2, output_path = output_path){
  if(is_missing(n)) {
    n = 2000
  }
  if(is_missing(bed1)) {
    stop("The bed1 information is missing. Please use indivPVal(n=2000,bed1='/dir/bed1.txt', genome ='/dir/genome.txt',bed2 ='/dir/bed2.txt', output_path = '/dir/output_folder')")
  }
  if(is_missing(genome)) {
    stop("The genome information is missing. If you are working with human You can use the human genome BED file within this package './Example/hg19_formatted_genomebedfile.txt'")
  }
  if(is_missing(bed2)) {
    stop("The bed2 file is missing. Please use the following format indivPVal(n=2000,bed1='/dir/bed1.txt', genome ='/dir/genome.txt',bed2 ='/dir/bed2.txt', output_path = '/dir/output_folder')")
  }
  if(is_missing(output_path)) {
    stop("The output_path is missing. Please use the following format indivPVal(n=2000,bed1='/dir/bed1.txt', genome ='/dir/genome.txt',bed2 ='/dir/bed2.txt', output_path = '/dir/output_folder')")
  }


  #final_values = list()
  #calculate size of the total binding region in query bed file
  query_length = overlapCalc(bed1,bed1)


  #generate background files from query:
  #step1: create background file directory
  background_dir = paste(bed1,"background",sep="_")
  if (dir.exists(background_dir)){
    unlink(background_dir, recursive = TRUE)
  }
  dir.create(background_dir)
  print("generating background files")

  #step2: write n randomly generated files.
  for (j in 1:(n)){
    background <- fread(cmd = paste("bedtools shuffle -i", bed1, "-g", genome))
    write.table(background, paste(background_dir,"/background_file",j,".txt",sep = ""),sep="\t",row.names=F, col.names = F, quote = F)
  }
  print("finished generating background files")

  ##calculate the jaccard index for original file bed1 and bed2
  num_overlap = overlapCalc(bed1, bed2)
  num1=query_length
  hit_length = overlapCalc(bed2, bed2)
  num2=hit_length

  #jaccard index
  num_total = num1 + num2 - num_overlap
  jaccard_id = num_overlap/num_total

  #percentage --- ANB/A and ANB/B
  #ANB/A
  A_similarity = (num_overlap)/num1
  B_similarity = (num_overlap)/num2
  #Summary
  indexes = c(queryfile = bed1,hitfile = bed2, jaccard_index = jaccard_id,percentage_A = A_similarity,percentage_B=B_similarity)


  ######## Calculate background jaccard
  background_n = list.files(background_dir)
  background_jaccard = list()

  background_jaccard = lapply(background_n, FUN = Background_jaccardfun, hitfile=bed2, background_dir=background_dir,  query_length=query_length,  hit_length=hit_length)

  jaccard_df = as.data.frame(do.call(rbind,background_jaccard))
  bed = as.numeric(as.character(jaccard_df[,3]))   ## jaccard scores of shuffled background

  #find the p-value
  num = as.numeric(as.character(match(as.numeric(as.character(indexes[3])),sort(c(bed,as.numeric(as.character(indexes[3]))),decreasing=T))))
  value = num/(n+1)

  #finding the pi-score = mean ratio * (-log10(p-value))
  if (mean(bed) == 0){
    piscore = (as.numeric(as.character(indexes[3])))/(0.000001) * (-log(value))
  } else{
    piscore = (as.numeric(as.character(indexes[3])))/(mean(bed)) * (-log(value))
  }
  final_values = c(query = bed1, hit = bed2,jaccard_index = as.numeric(as.character(indexes[3])),pi_score = piscore, p_value = value, percentage_A = as.numeric(as.character(indexes[4])), percentage_B = as.numeric(as.character(indexes[5])),summary(bed))

unlink(background_dir, recursive = TRUE)
write.csv(final_values,paste(output_path,"/",gsub("^.*/", "", bed1),"_pair_jaccard_pval.csv",sep=""))
print("finished output")
}


