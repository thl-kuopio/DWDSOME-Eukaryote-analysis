######################### BIOINFORMATICS#############################################################
###########################18S analysis #############################################################
## removing the adaptors and bad quality reads using cutadapt and trimmomatic using default parameters
## using FLASH software to merge the forward and reverse reads
flash2 read.1.fastq read.2.fastq
##uses the split_libraries_fastq command to convert from fastq to fasta with quality filtering while definding sample names for each sample: resequenced individuals were given the suffix “R”. Because the command was done per sample, a dummy mapping file was used. ##
##combine the sequences##
cat Eu_slout_q20/slout_single_sample_q20_18s*/seqs.fna > final_demultiplexed_seqs_18S.fasta 
##identify (using usearch wrapped within QIIME) then exlucde chimeras##
identify_chimeric_seqs.py -i final_demultiplexed_seqs_18S.fasta -m usearch61 -o Eu_usearch_checked_chimeras/ -r greengenes/gg_13_8_otus/rep_set/97_otus.fasta
filter_fasta.py -f final_demultiplexed_seqs_18S.fasta -o final_chimeric_rmv_seqs_18S.fasta -s bact_usearch_checked_chimeras/chimeras.txt -n
##cluster sequences into OTUs, first by closed-reference against the greengenes database then by de novo clustering at 100% subsampling ##
pick_open_reference_otus.py -i final_chimeric_rmv_seqs_18S.fasta  -r Silva_database/128.fasta -o Euk_open_ref_picked_otus -a -O 4 --percent_subsample 1
## filter spurious OTUs from total sequences
filter_otus_from_otu_table.py -i Euk_open_ref_picked_otus/otu_table_mc2_w_tax_no_pynast_failures.biom -o Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001.biom --min_count_fraction 0.00001
## remove all chloroplasts/mitochondria then all empty OTUs
filter_taxa_from_otu_table.py -i Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001.biom -o Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom -n  c__Chloroplast,f__mitochondria
filter_otus_from_otu_table.py -i Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom -o Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom -n 1

########## community analyses###############
## alpha richness across samples ##
alpha_rarefaction.py -i Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -m map.18S.txt -o alpha_rare_e2873_DWDSOME -a -e 2873 -O 6 -p alpha_params.txt
## beta diversity across samples ##
beta_diversity_through_plots.py -i Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -m map.18S.txt -o bdiv_e2873_DWDSOME -a -t rep_set.tre -e 2873 -O 6 -p beta_params.txt
#############################################
summarize_taxa.py -i Euk_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -o Euk_tax/

############################################################################################################################
