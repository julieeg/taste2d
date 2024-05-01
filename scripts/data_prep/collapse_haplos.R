
# load packages
library(tidyverse) ; library(data.table) ; library(vcfR)


file_names <- paste0("chunk",10:39,"_haplos_id.csv")
file_path <- "../data/processed/genos/supertasters"

files <- do.call(rbind.data.frame, lapply(as.list(paste0(file_path, "/", file_names)), function(file) fread(file)))


# write csv
fwrite(files, "../data/processed/genos/supertasters/supertasters_haplos_id.csv", col.names=T, row.names=F, quote=F)

#EOF