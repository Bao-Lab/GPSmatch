#' @title Pre-compute the size of binding sites of each bed file in database
#' @description The function compute the size of binding sites of each bed file in database and save the output(dataframe containing two columns:filename and binding size) into a ../database_dir_bindingsize.csv file
#' @param database_dir directory path including the name of the database folder
#' @param output_path directory path of the output file, save the output dataframe to .csv file(naming: "database_dir"_bindingsize.csv) into the user provided "output_path"
#' @export
#' @examples
#' databaseFileSize(database_dir="/dir/database_dir", output_path="output_path" )



databaseFileSize = function(database_dir=database_dir, output_path=output_path){
  if(is_missing(database_dir)) {
    stop("database directory is missing. Please use the following format DatabaseFileSize('/dir/database_dir','output_path' )")
  }
  if(is_missing(output_path)) {
    stop("output_path is missing. Please use the following format DatabaseFileSize('/dir/database_dir','output_path' )")
  }


  print("calculating the length of binding region for each file in the database ...")

  # Extract the names of all Bed Files In Database
#  levels = list.files(database_dir, recursive = TRUE, pattern=".csv")
  levels = list.files(database_dir, recursive = TRUE)

#for each Bed file in the database, use databaseJaccard function to find out its length of binding sites.
  bindingsize = lapply(levels, FUN = databaseJaccard, database_dir = database_dir)
  #Construct a dataframe containing file name and the size of binding
  file_size = data.frame(levels, unlist(bindingsize))
  colnames(file_size) = c("filename", "size")
  #find the database name from its full path, and will use it as part of the name of output
  Split = str_split(database_dir, "/")
  databasename = Split[[1]][length(Split[[1]])]

  Output = paste(output_path, "/",databasename, "_BindingSize.csv", sep = "")
  write.csv(file_size, Output, row.names = FALSE)
}


#calculate size of the total binding region in each database bed file
databaseJaccard = function(databasefile=databasefile, database_dir=database_dir){
  pwd = paste(database_dir,databasefile,sep="/")
  databaefile_jaccard = overlapCalc(pwd,pwd)
}


