#!/usr/bin/env Rscript
# Juliane Weller, 01 Jan 2024, Analyzing the NovaSeq data for Dipper, sequencing started with dark cycles in the primer

# Required libraries
library(ShortRead)
library(stringr)
library(tidyverse)


# A function to generate reverse complement nucleotides
rc <- function(x) {toupper(spgs::reverseComplement(x))}

# MiSeq
count_edits <- function(file, outdir, library) {
    #Read in the sequencing reads
    print(paste0(Sys.time(), " - Read FastQ file..."))
    filename = sub("_R.*$", "", basename(file))
    print(filename)
    reads <- readFastq(file)
    reads <- as.character(reads@sread)
#     reads <- head(reads, n = 1000)
    
    # extract the sequences
    print(paste0(Sys.time(), " - Extract oligo components..."))
    count_table <- tibble(raw_sequence = reads, 
                        spacer = paste0("G", substr(raw_sequence, start = 1, stop = 19)), # different to MiSeq due to dark cycles
                        extension = str_extract(raw_sequence, paste0("GAGTCGGTGC","(.*?)","(CG.{3}TTCTA|.{5}TT.TAT.T.GTT)"))%>% # I loosened the requirement here a bit
                          substr(start = 11, stop = nchar(.)-10),
                        scaffold = str_extract(raw_sequence, paste0("(.*?)","GAGTCGGTGC"))%>%
                          substr(start = 20, stop = nchar(.)),) # due to the dark cycle, this needed to be at 20 not 21
    
    # only keep sequences that have all components
    print(paste0(Sys.time(), " - Generate readcounts..."))
    count_table_complete = filter(count_table, !is.na(spacer) & !is.na(extension) & !is.na(scaffold)) %>%
        mutate(extension_rc = rc(extension), pbs_short = substr(extension_rc, start = 1, stop = 9)) # Take the shortest potential PBS, even tough some might be longer, but this is for checking if this substring can be found in the spacer
    
    
    print(paste0("Spacer is na: ", sum(is.na(count_table$spacer))))
    print(paste0("Extension is na: ", sum(is.na(count_table$extension))))
    print(paste0("Scaffold is na: ", sum(is.na(count_table$scaffold))))
    
    perc_fulloligo <- nrow(count_table_complete)/nrow(count_table)*100
    print(paste0("Percent of reads with PBS, extension and target: ", perc_fulloligo))
    
    # identify the types of recombinations that might have occured
    recombination_df <- count_table_complete %>% 
        mutate(spacer_match_pbs = ifelse(str_detect(spacer, pbs_short), T, F), 
                recombination = ifelse(spacer_match_pbs == T, "good", "no_match"))
    
    # some stats
    perc_norecomb <- recombination_df %>% group_by(recombination) %>% summarise(n = n()) %>% mutate(total = sum(n), fraction = n/total)
    print(perc_norecomb)

    # For first pass, filter for no recombination - I'm using the df with recombination, as I'll filter for it later
    # no_recombination <- filter(recombination_df, recombination == 'good')
    
    # Categorize the scaffold, don't do exact matching but identify type
    no_recombination_scfld <- recombination_df %>% # I'm using the one with recombination!
        mutate(scaffold_type = case_when(
            str_detect(scaffold, "GTTTAAGAGC(.*)GAAAAAGT") ~ "improved",
            str_detect(scaffold, "GTTTGAGAGC(.*)GAAAAAGT") ~ "paste",
            TRUE ~ "unknown"))
    
    # Generate editing outcomes
    print(paste0(Sys.time(), " - Generate count table..."))
    
    editing_table <- no_recombination_scfld %>% 
        group_by(spacer, extension_rc, scaffold_type) %>% 
        summarise(n = n())
    
    print(paste0(Sys.time(), " - Save count table..."))
    write.csv(editing_table, paste0(outdir, filename, "_counts_newpipe.csv"))
    print(paste0(Sys.time(), " - Count table saved"))
    
    total_roi <- sum(editing_table$n)
    perc_oligos <- total_roi/nrow(count_table)*100
    perc_library_covered <- nrow(editing_table)/nrow(library)*100 # this is not accurate, as some things might be mutations and not actually part of the library   
    
    write_tsv(tibble(sample = filename, total_roi, perc_fulloligo, 
                     recomb_good = perc_norecomb$fraction[perc_norecomb$recombination == "good"],
                     recomb_bad = perc_norecomb$fraction[perc_norecomb$recombination == "no_match"], 
                     perc_oligos, perc_library_covered), paste0(outdir, filename, "_metrics.tsv"))
    
    print(paste0(Sys.time(), " - File saved to: ", outdir, filename, "_counts_newpipe.csv"))
}


### RUN 
args <-commandArgs(trailingOnly=TRUE)
dipper_library <- read.csv("/library_features.csv")

count_edits(args[1], args[2], dipper_library)