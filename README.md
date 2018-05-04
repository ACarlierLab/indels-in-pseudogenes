# indels-in-pseudogenes
Scripts and data for the investigation of indel biases in bacterial genomes

This is a collection of Python scripts used for the analysis of indel profiles in pseudogenes in: 
"Patterns of nucleotide deletion and insertion inferred from bacterial pseudogenes", 2018, by Bram Danneels, Marta Pinto-Carbó and Aurélien Carlier

Curated alignments discussed in the paper are saved under the Burkholderia_curated-final.tar.gz archive.

The steps to replicate the analyses described in the manuscript are as follows:

1. Download sets of 3 genomes in genbank format and save in a directory called "genomes". Pseudogenes must be already annotated in the genomes using the /pseudo qualifier for CDS features
2. Run Pseudo_genome_align.py as follows:
    $ python Pseudo_genome_align.py genomes
    
    The script extracts sequences of annotated pseudogenes and performs a blastx search on intergenic regions to extract unannotated pseudogenes. Candidate pseudogenes are aligned to the best blastx hit using tfasty v3.6 (Pearson, 2000) to verify the presence of at least one inactivating mutation (frameshift, early stop codon or truncation) affecting at least 20% of the reading frame as described in Carlier et al. 2014.
    Predicted proteins are extracted for each gene, including pseudogenes, and orthologs are computed using OrthoMCL v1.4 (Li et al. 2003). Pseudogenes as well as upstream and downstream flanking genes are extracted and aligned to syntenic functional orthologs from the two other reference genomes using MAFFT v7.2 with the E-INS-i algorithm with default parameters. Alignments are further sliced to remove upstream and downstream features, leaving only the alignment region containg the pseudogene and orthologous functional references. 
    
    Output: - alignment_whole: Directory containing MAFFT alignments including upstream and downstream anchor sequences. Alignments are saved in the Clustal format.
            - alignments : directory containing the sliced alignments including solely the gene references and pseudogene regions. Alignments are saved in the Clustal format.
    
    Various temporary files useful for diagnostic:
            - orthomcl.out: output of orthomcl
            - one_to_one.txt: set of single copy orthologs from the 3 query genomes, extracted from the orthomcl output.
            - syntenic_blocks.txt: File containg the coordinates of syntenic blocks of genes in the 3 target genomes as inferred from the relative positions of single copy orthologs in one_to_one.txt
            
3. Manually curate the files in the alignments directory to retain only pseudogenes which contain at least 3 unambiguous indels. Save the curated alignments in a separate directory, for example called "curated_alignments"

4. Extract information about nucleotide insertions and deletions in pseudogenes by running the del_align_analysis.py Python script:
        $ python del_align_analysis.py [curated_alignments]
        
        Output: - ins_statistics.txt : text file containing the sizes of predicted insertions (one value per line)
                - del_statistics.txt : text file containing the sizes of predicted insertions (one value per line)
                - GCstatistics.txt: text file containing average %GC of pseudogene (with alignment gaps ignored), %GC of closest functional orthologs, number of indels in pseudogene and % identity between the functional reference and the pseudogene.  
                
5. Extract information about the presence of direct repeats at indel sites using Python script del_align_analysis_DR.py
        $ python del_align_analysis_DR.py
        
        Output: - del_statistics_DR.txt : text file containing three tab-separated columns:
                        1. size of indel (negative for insertion, positive value for deletion)
                        2. number of indels in pseudogene
                        3. Alignment file name (derived from pseudogene name)
                        
                - DR_statistics.txt : text file containing 4 tab-separated columns:
                        1. Name of pseudogene
                        2. motif in insertion
                        3. repeated motif
                        4. size of motif

6. To extract the data used for Figure S3 (ratio of 3n deletions per number of indels), run the Python script indel_per_pseudo.py on the output of del_align_analysis_DR.txt:
    $ python indel_per_pseudo del_statistics_DR.txt
    
        Output: - indel_by_gene.txt: tab-delimited text file with 4 columns

7. To perform bootstrapping analysis, first create separate files for deletions and insertions, containing overview tables of number of deletions of certain length. Each row should represent the size distribution for a given dataset (e.g. species, pooled data). Row- and columnames should be present. Values should be separated by tabs.
	    $ python bootstrap_analysis.py [input table] [indel type (either 'Deletion' or 'Insertion')]
	
	    Output: - bootstrap.out: text file containing results of the bootstrapping analysis for each set (row) of data in 3 tab-separated columns:
				    1. Indel size
				    2. Mean (n=10000) of the differences between boostrapped values and power-law estimates
				    3. 95 % confidence interval on the means
				- [dataset]_[indeltype].png: png file with the size distribution of deletion/insertion of certain dataset, power-law fit, and bootstrapped values + 95% confidence interval

8. To perform motif detection, execute the motif_search.py script in the folder containing the curated alignments in CLUSTAL format. Output is printed to STDOUT by default.
	    $ python motif_search.py
	
	    Output: For every pseudogene (= alignemnt) the following will be displayed:
				- Name of alignment file (flanked by 5 '=' characters)
				- Part of the alignment containing the gap and 10 nt flanking sequence (if present), each deletion seperated by a line of '-' characters
				- Forward motifs search results (Forward motif = motif where the start of the deleted sequence matches with he sequence right after the gap)
				- Reverse motifs search results (Reverse motif = motif where the end of the deleted sequence matches with the sequence right before the gap)
				The file is ended with overview statistics of detected deletions and motifs (# deletions, # small deletions, # forward motifs, # reverse motifs, # deletions with detected motif)
				
	
