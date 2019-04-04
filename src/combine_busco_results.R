library(data.table)



FindBuscoResult <- function(x){
    list.files(paste("output/050_busco", x, sep = "/"),
               recursive = FALSE,
               pattern = "full_table_",
               full.names = TRUE)
}

all_busco_dirs <- list.dirs("output/050_busco", recursive = FALSE)
busco_dirs <- grep("^run_", basename(all_busco_dirs), value = TRUE)

busco_files <- sapply(busco_dirs, FindBuscoResult)
names(busco_files) <- gsub("full_table_([^[:space:]]+).tsv",
                           "\\1",
                           basename(busco_files))

busco_results_list <- lapply(busco_files, fread, skip = 4, fill = TRUE)
busco_results <- rbindlist(busco_results_list, idcol = "assembly")

fwrite(busco_results, "test/busco_results_combined.csv")
