#' @title Computes the p-value and pi-scores for top hits from rankSimilarity()
#' @description The function evaluates the significance of the jaccard index scores by calculating p-values and pi-scores.
#' @param n The number of background files generated in order to compute the p-value. As the n increases, the p-value will become more reliable, but the user should be aware that this will significantly increase the computing time. We have set a default n of 2000.
#' @param tophit The number of top hits that the user are interested to compute p-value and pi-scores. Those top hits were generated and ordered by function rankSimilarity() at step 2 and were written to file queryBed_jaccard_output.
#' @param genome The file path of a genome file, which should be tab delimited and structured as follows: <chromName><TAB><chromSize>. A pre-formatted hg19 genome file can be found on the Github.
#' @param queryBed The file path of a query bed file to be compared to the database files.
#' @param queryBed_jaccard_output The file path of the output file generated from the rankSimilarity() function
#' @param database_dir The directory of a folder containing database files to be used for comparison with the query file.
#' @param output_path The file path where the output .csv file will be generated.
#' @export
#' @return A dataframe that shows the similarities of the query file to the database files ranked from greatest to least, tailored to user method specification.
#' @examples
#' sigEvl(n=2000, tophit=5,genome="/dir/genome-file.txt", queryBed="/dir/bed1.txt",queryBed_jaccard_output="/dir/bed1_jaccard_output.txt",database_dir="/dir/database_dir",output_path="/dir/output_path")


sigEvl = function(n=n,tophit=tophit,genome=genome,queryBed=queryBed,queryBed_jaccard_output= queryBed_jaccard_output,database_dir=database_dir, output_path=output_path){
  if(is_missing(n)) {
    n = 2000
  }
  if(is_missing(tophit)) {
    tophit = 1
  }
  if(is_missing(genome)) {
    stop("The genome information is missing. If you are working with human You can use the human genome BED file within this package './Example/hg19_formatted_genomebedfile.txt'")
  }
  if(is_missing(queryBed)) {
    stop("The queryBed is missing. Please use the following format sigEvl(n=2000, tophit=1, genome='/dir/genome-file.txt', queryBed='/dir/bed1.txt',queryBed_jaccard_output='/dir/bed1_jaccard_output.txt',database_dir='/dir/database_dir',output_path='/dir/output_path')")
  }
  if(is_missing(queryBed_jaccard_output)) {
    stop("The queryBed_jaccard_output from rankSimilarity() is missing. Please use the following format sigEvl(n=2000, tophit=1, genome='/dir/genome-file.txt', queryBed='/dir/bed1.txt',queryBed_jaccard_output='/dir/bed1_jaccard_output.txt',database_dir='/dir/database_dir',output_path='/dir/output_path')")
  }
  if(is_missing(database_dir)) {
    stop("The database_dir  is missing. Please use the following format sigEvl(n=2000, tophit=1, genome='/dir/genome-file.txt', queryBed='/dir/bed1.txt',queryBed_jaccard_output='/dir/bed1_jaccard_output.txt',database_dir='/dir/database_dir',output_path='/dir/output_path')")
  }
  if(is_missing(output_path)) {
    stop("The output_path is missing. Please use the following format sigEvl(n=2000, tophit=1, genome='/dir/genome-file.txt', queryBed='/dir/bed1.txt',queryBed_jaccard_output='/dir/bed1_jaccard_output.txt',database_dir='/dir/database_dir',output_path='/dir/output_path')")
  }

  final_values = list()
  indexes = list()
  file = read.csv(queryBed_jaccard_output)

  #calculate size of the total binding region in query bed file
  query_length = overlapCalc(queryBed,queryBed)

  #generate background files from query:
  #step1: create background file directory
  background_dir = paste(database_dir,"background",sep="_")
  if (dir.exists(background_dir)){
    unlink(background_dir, recursive = TRUE)
  }

  dir.create(background_dir)
  print("generating background files")
  #background_list = list()

  #step2: write n randomly generated files.
  for (j in 1:(n)){
    background <- fread(cmd = paste("bedtools shuffle -i", queryBed, "-g", genome))
    write.table(background, paste(background_dir,"/background_file",j,".txt",sep = ""),sep="\t",row.names=F, col.names = F, quote = F)
  }
  print("finished generating background files")

  #define how many top hits will be used to evaluate p-value and Pi-score
  if (nrow(file)< tophit){
    top = nrow(file)
  }else{
    top = tophit
  }

  for (i in 1:top){
    #for each of the top hit files in folder_dir calculate original input file Jaccard index
    hitfile = paste(database_dir,as.character(file[i,2]),sep = "/")

    hit_length = overlapCalc(hitfile,hitfile)

    num_overlap = overlapCalc(queryBed, hitfile)
    num1=query_length
    num2=hit_length

    #jaccard index
    num_total = num1 + num2 - num_overlap
    jaccard_id = num_overlap/num_total

    #percentage --- ANB/A and ANB/B
    #ANB/A
    A_similarity = (num_overlap)/num1
    B_similarity = (num_overlap)/num2
    #Summary
    indexes[[i]] = c(queryfile = queryBed,hitfile = as.character(file[i,2]), jaccard_index = jaccard_id,percentage_A = A_similarity,percentage_B=B_similarity)


    ######## Calculate background jaccard
    background_n = list.files(background_dir)
    background_jaccard = list()

    background_jaccard = lapply(background_n, FUN = Background_jaccardfun, hitfile=hitfile, background_dir=background_dir,  query_length=query_length,  hit_length=hit_length)

    jaccard_df = as.data.frame(do.call(rbind,background_jaccard))
    bed = as.numeric(as.character(jaccard_df[,3]))   ## jaccard scores of shuffled background

    #find the p-value
    num = as.numeric(as.character(match(as.numeric(as.character(indexes[[i]][3])),sort(c(bed,as.numeric(as.character(indexes[[i]][3]))),decreasing=T))))
    value = num/(n+1)

    #finding the pi-score = mean ratio * (-log10(p-value))
    if (mean(bed) == 0){
      piscore = (as.numeric(as.character(indexes[[i]][3])))/(0.000001) * (-log(value))
    } else{
      piscore = (as.numeric(as.character(indexes[[i]][3])))/(mean(bed)) * (-log(value))
    }
    final_values[[i]] = c(query = queryBed, hit = as.character(file[i,2]),jaccard_index = as.numeric(as.character(indexes[[i]][3])),pi_score = piscore, p_value = value, percentage_A = as.numeric(as.character(indexes[[i]][4])), percentage_B = as.numeric(as.character(indexes[[i]][5])),summary(bed))

  }
  unlink(background_dir, recursive = TRUE)
  values_df = as.data.frame(do.call(rbind,final_values))
  values_df = values_df[order(as.numeric(as.character(values_df$jaccard_index)),decreasing = T),]

  write.csv(values_df,paste(output_path,"/",gsub("^.*/", "", queryBed),"_jaccard_pval.csv",sep=""))
  print("finished output")
}

#Calculate the jaccard index between background file and query file
Background_jaccardfun = function(backgroundfile =backgroundfile, hitfile=hitfile, background_dir =background_dir, query_length = query_length, hit_length = hit_length){
  query_background = paste(background_dir,backgroundfile,sep="/")

  num_overlap = overlapCalc(query_background, hitfile)
  # We use the size of the query as an approximation for the size of the similated query file to spped up the computation.
  num1=query_length
  num2=hit_length

  #jaccard index
  num_total = num1 + num2 - num_overlap
  jaccard_id = num_overlap/num_total
  #line = c(bedfile = as.character(file[i,2]),jaccard_index = jaccard_id)

  ##calc jaccard
  line = c(bedfile = backgroundfile, hitfile = hitfile, jaccard_index = jaccard_id)
  return (line)
}
