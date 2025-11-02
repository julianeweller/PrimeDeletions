#!/usr/bin/env Rscript
# Juliane Weller, 04 Dec 2024, Analyzing the Avanti data for Dunnock

# Required libraries
library(ShortRead)
library(stringr)
library(tidyverse)

# A function to generate reverse complement nucleotides
rc <- function(x) {toupper(spgs::reverseComplement(x))}

# Avanti
count_edits <- function(file, outdir) {
    
    # Read in the sequencing reads from the file that has been stitched together
    print(paste0(Sys.time(), " - Read FastQ file..."))
    filename = sub("_stitched.fastq.*$", "", basename(file))
    print(filename)

    reads <- readFastq(file)
    reads <- as.character(reads@sread)

    # Define matching patterns
    extension_pattern <- "GAGTCGGTGC(.*?)CGCGGTTCTA"
    spacer2_pattern <- "(.*?)GTTTAAGAGC"
    
    # extract the sequences
    print(paste0(Sys.time(), " - Extract oligo components..."))

    count_table <- tibble(raw_sequence = reads) %>%
        mutate(
            raw_seq_peg2 = substr(raw_sequence, start = 150, stop = nchar(raw_sequence)),
            spacer1 = substr(raw_sequence, start = 1, stop = 20),
            spacer2 = str_extract(raw_seq_peg2, spacer2_pattern)%>%
                          substr(start = nchar(.)-29, stop = nchar(.)-10),
            extension1 = str_extract(raw_sequence, extension_pattern)%>%
                          substr(start = 11, stop = nchar(.)-10),
            extension2 = str_extract(raw_seq_peg2, extension_pattern)%>%
                          substr(start = 11, stop = nchar(.)-10))
            
    
    # only keep sequences that have all components
    print(paste0(Sys.time(), " - Filter for components..."))
    count_table_complete = filter(count_table, !is.na(spacer1) & !is.na(extension1) & !is.na(spacer2) & !is.na(extension2)) %>%
        mutate(
            extension1_rc = rc(extension1), pbs1 = substr(extension1_rc, start = 1, stop = 13),
            extension2_rc = rc(extension2), pbs2 = substr(extension2_rc, start = 1, stop = 13)) %>%
        filter(!is.na(pbs1) & !is.na(pbs2))  # Additional filter step to only keep rows where we have long enough extensions
    

    print(paste0("Spacer1 is na: ", sum(is.na(count_table$spacer1))))
    print(paste0("Spacer2 is na: ", sum(is.na(count_table$spacer2))))
    print(paste0("Extension1 is na: ", sum(is.na(count_table$extension1))))
    print(paste0("Extension2 is na: ", sum(is.na(count_table$extension2))))
    
    perc_fulloligo <- nrow(count_table_complete)/nrow(count_table)*100
    print(paste0("Percent of reads with all components: ", perc_fulloligo))
    
    # identify the types of recombinations that might have occured
    recombination_df <- count_table_complete %>%
        mutate(spacer1_match_pbs1 = str_detect(spacer1, pbs1),
            spacer1_match_pbs2 = str_detect(spacer2, pbs2), 
            recombination = case_when(
                spacer1_match_pbs1 & spacer1_match_pbs2 ~ "good",
                spacer1_match_pbs1 ~ "mismatch_pegRNA1",
                spacer1_match_pbs2 ~ "mismatch_pegRNA2",
                TRUE ~ "mismatch_both"))
               
    # some stats
    perc_norecomb <- recombination_df %>%
        count(recombination) %>%
        mutate(total = sum(n), fraction = n / total)           
    print(perc_norecomb)
    
    # Generate editing outcomes
    print(paste0(Sys.time(), " - Generate count table..."))
    
    editing_table <- recombination_df %>% 
        group_by(spacer1, extension1, spacer2, extension2) %>% 
        summarise(n = n())
    
    print(paste0(Sys.time(), " - Save count table..."))
    write.csv(editing_table, paste0(outdir, filename, "_counts.csv"))
    print(paste0(Sys.time(), " - Count table saved"))
    
    total_roi <- sum(editing_table$n)
    perc_oligos <- total_roi/nrow(count_table)*100
    print(total_roi)
       
    
    recomb_stats <- tibble(
        sample = filename, total_roi = total_roi, 
        perc_fulloligo = perc_fulloligo, 
        recomb_good = ifelse(any(perc_norecomb$recombination == "good"), perc_norecomb$fraction[perc_norecomb$recombination == "good"], 0),
        recomb_bad_pegRNA1 = ifelse(any(perc_norecomb$recombination == "mismatch_pegRNA1"), perc_norecomb$fraction[perc_norecomb$recombination == "mismatch_pegRNA1"], 0),
        recomb_bad_pegRNA2 = ifelse(any(perc_norecomb$recombination == "mismatch_pegRNA2"), perc_norecomb$fraction[perc_norecomb$recombination == "mismatch_pegRNA2"], 0),
        recomb_bad_both = ifelse(any(perc_norecomb$recombination == "mismatch_both"), perc_norecomb$fraction[perc_norecomb$recombination == "mismatch_both"], 0),
        perc_oligos = perc_oligos)
    
    write_tsv(recomb_stats, paste0(outdir, filename, "_metrics.tsv"))
    
    print(paste0(Sys.time(), " - File saved to: ", outdir, filename, "_counts.csv"))
}


### RUN 
args <-commandArgs(trailingOnly=TRUE)
print("Start the script...")
count_edits(args[1], args[2])