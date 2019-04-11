library(data.table)
library(ggplot2)

new_dir <- "/cifs/ro_deardenlabarchive/tomharrop/data/nanopore_runs/asw_47/basecall_guppy2.3.7/basecalled"
old_dir <- "/cifs/ro_deardenlabarchive/tomharrop/data/nanopore_runs/asw_47/basecall/basecalled"

ss_lists <- lapply(list(new = new_dir,
            old = old_dir),
       list.files,
       pattern = "sequencing_summary.txt",
       recursive = TRUE,
       full.names = TRUE)

ss_table_list <- lapply(ss_lists, function(x)
    lapply(x, function(y)
        cbind(fread(y), fc = basename(dirname(y)))))

# release memory
bc_results <- rbindlist(lapply(ss_table_list, rbindlist),
                        idcol = "basecall_run",
                        fill = TRUE)

# release memory
rm(ss_table_list)

ggplot(cbind(bc_results, z = 1), 
       aes(x = sequence_length_template,
           y = mean_qscore_template,
           z = z)) +
    theme_minimal() +
    facet_grid(. ~ basecall_run) +
    scale_fill_viridis_c(guide = guide_colorbar(title = "Log counts")) +
    scale_x_log10() + 
    stat_summary_hex(fun=function(z){log(sum(z))},
                     bins = 50)
