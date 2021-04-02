

library(tidyverse)
library(GenomicRanges)

#Specify tissue-type(s) here
tissue.list <- c("Uterus")
smtsd.list <- c("Uterus")


###########################################################################################
############################  NON-TISSUE SPECIFIC CODE ####################################
############################## NEEDS TO RUN JUST ONCE #####################################
###########################################################################################


#################################################################################################
#################################### Reading in the full file ###################################
#################################################################################################

gtex.samp.attrs <- read_tsv("./ERR-datafiles/v7_eqtl_data/gtex7_sampleAttributes_groomed.txt")

# Separate out the DNA and RNA samples for genotype and RNASeq respectively
#gtex.dna.samples <- gtex.samp.attrs %>% filter(ANALYTE=="DNA")
gtex.dna.samples <- gtex.samp.attrs %>% filter(SMAFRZE=="WGS")

# Add a donor_id column (first two parts of SAMPID) to match with 
# genotype data with RNASeq
gtex.dna.samples <- gtex.dna.samples %>% 
  separate(SAMPID,sep="-",into=c("donor1","donor2","idp3","idp4","idp5","idp6"),remove=F) %>% 
  unite(donor_id, donor1:donor2, sep="-") %>% dplyr::select(-(idp3:idp6))


#gtex.rna.samples <- gtex.samp.attrs %>% filter(ANALYTE=="RNA")
gtex.rna.samples <- gtex.samp.attrs %>% filter(SMAFRZE=="RNASEQ")

# Add a donor_id column (first two parts of SAMPID) to match with 
# RNASeq data with genotype
gtex.rna.samples <- gtex.rna.samples %>% 
  separate(SAMPID,sep="-",into=c("donor1","donor2","idp3","idp4","idp5","idp6"),remove=F) %>% 
  unite(donor_id, donor1:donor2, sep="-") %>% dplyr::select(-(idp3:idp6))

###########################################################################################
###########################################################################################
#################################################################################################
###
rnaseq.file <- "./ERR-datafiles/v7_eqtl_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"

rnaseq.data <- read_tsv(rnaseq.file, skip=2)
gene.names <- rnaseq.data$Description
names(gene.names) <- rnaseq.data$Name




##########################################################################################################
###################################### -ENHANCER AND GENCODE FILE PATH SETTINGS -#########################
##########################################################################################################
enhancer.filePath <- "./ERR-datafiles/combinedEnhancerAnnotations.txt"
gencode.filepath <- "./ERR-datafiles/gencode_info.txt"

#gtex.ase.filepath <- "./gtex_ASE_uterus_head.txt"
##########################################################################################################
###################################### - ENCODE ANNOTATIONS -#############################################
##########################################################################################################

# Read in the gencode gene annotations (previously processed by script 'gencode.r')
gencode.info <- read_tsv(gencode.filepath)

# Get the required info (chromosome name, start and end coordinates, gene ids, and gene names)
# The rows are labelled by gene IDs for easy selection downstream.
gencode.gr <- GRanges(seqnames = gencode.info$Chr,
                      ranges = IRanges(start = gencode.info$Start,
                                       end = gencode.info$End,
                                       names = gencode.info$ID),
                      strand = gencode.info$Strand,
                      gene_id = gencode.info$ID,
                      gene_name = gencode.info$Name,
                      annot_type=gencode.info$AnnotType,
                      evidence=gencode.info$Evidence,
                      type=gencode.info$Type
)


##########################################################################################################
###################################### - ENHANCER ANNOTATIONS -#############################################
##########################################################################################################


########Enhancer Data
#enhancer.filePath <- "./enhancers_50rows.txt"
#enhancer.filePath <- "./3_Ehnacer.merge.method.txt"

# Read in enhancer coordinate files
enhancers <- read_csv(enhancer.filePath, col_names = F, col_types = "iciic")
names(enhancers) <- c("enhID", "Chr", "Start", "End", "source")
enhancers$enhID <- paste0("enh",enhancers$enhID)
enhancers$Chr <- paste0("chr",enhancers$Chr)

enhancers.gr <- GRanges(seqnames = enhancers$Chr,
                        ranges = IRanges(start = enhancers$Start,
                                         end = enhancers$End,
                                         names = paste("EN",1:nrow(enhancers),sep="_")),
                        #strand = "*",
                        source = enhancers$source
)

#cmnt
## COSMIC gene annotations...read in List of COSMIC genes for checking for genes in cosmic
cosmic.info <- read_tsv("./ERR-datafiles/cosmic_gene_info.txt")
## Also get DO brca and prad annotations
brca.prad.info <- read_tsv("./ERR-datafiles/DO/brca_prad_gene_info.txt")


#unique.enhancers.gr <- unique(enhancers.gr)
######################### END OF NON-TISSUE SPECIFIC CODE #################################
###########################################################################################
###########################################################################################



############################ TISSUE SPECIFIC CODE BEGINS ##################################
###########################################################################################
###########################################################################################




counts.table <- data.frame(Tissue=character(), enhancer_eqtl_counter=numeric(), multiTarget_eqtl_count=numeric(),
                           all_eqtls_in_promoter_count=numeric(), p_eqtls_count=numeric(), multiGene_pEqtls_count=numeric(),
                           common_genes_count=numeric(), peqtls_in_commonGenes_count=numeric(), genePairs_oppSlopes_count=numeric(),
                           geneA_enhancer_count=numeric(), bothGenes_enhancer_count=numeric(), stringsAsFactors = FALSE)

#Can be made to loop over multiple tissues
for(k in 1:length(tissue.list)){

  current.tissue <- tissue.list[k]
  current.smtsd <- smtsd.list[k]
  
  
  cat(paste0("--> Current tissue is ",current.tissue,"\n"))
  # Extract the sample attributes for RNASeq data
  rna.samps.tissue <- filter(gtex.rna.samples, SMTSD == current.smtsd) %>% dplyr::select(SAMPID, donor_id)
  # Save the sample IDs 
  rna.sampids.tissue <- rna.samps.tissue$SAMPID
  
  
  
  eqtl.filePath <- paste0("./ERR-datafiles/v7_eqtl_data/",current.tissue,".v7.signif_variant_gene_pairs.txt")

  
  eqtls <- read_tsv(eqtl.filePath,col_types="cci__ddd____")
  
  eqtls <- eqtls %>% 
    separate(variant_id, c("Chr", "Pos", "Allele1", "Allele2", "del"), convert = T, remove = F) %>%
    dplyr::select(-c(del))
  
  
  
  ###
  eqtls.gr <- GRanges(seqnames = paste0("chr",eqtls$Chr),
                      ranges = IRanges(start = eqtls$Pos, 
                                       width = 1)
  )
  
  values(eqtls.gr) <- dplyr::select(eqtls, 4:10)
  # eQTL id is kept same as the variant id for ease of retrieval later
  eqtls.gr$eqtl_id <- eqtls$variant_id
  
  #Limit eQTLs to 200kb to keep interactions within contact domains
  dom.eqtls.gr <- eqtls.gr[abs(eqtls.gr$tss_distance)<=200000]
  
  
  # Make a GrangesList of all eqtls based on the unique eqtl id
  eqtls.grList <- split(eqtls.gr, eqtls.gr$eqtl_id)
  # Make a shorter list of eqtls that have at least 2 targets
  multiTarget.eqtls.grList <- eqtls.grList[lapply(eqtls.grList, length) > 1] #keepCount
  multiTarget.eqtl.ids <- names(multiTarget.eqtls.grList)
  #counter
  multitarget.eqtl.counter <- length(multiTarget.eqtls.grList)
  #############################################
  
  ###############################################################################################
  ############################ eQTLs occuring within Enhancers ##################################
  ###############################################################################################
  #countThis
  enhancer.hits <- findOverlaps(eqtls.gr, enhancers.gr)
  eqtls.in.enhancers <- eqtls.gr[queryHits(enhancer.hits)]
  eqtls.in.enhancers$enhancer_id <- names(enhancers.gr)[subjectHits(enhancer.hits)]
  eqtls.in.enhancers$enhancer_start <- start(enhancers.gr)[subjectHits(enhancer.hits)]
  eqtls.in.enhancers$enhancer_end <- end(enhancers.gr)[subjectHits(enhancer.hits)] #keepCount
  #counter
  enhancer.eqtl.counter <- length(eqtls.in.enhancers)
  
  # Also make a GrangesList grouped by enhancer ID #### restart from here
  enhancer.eqtls.grList <- split(eqtls.in.enhancers, eqtls.in.enhancers$enhancer_id)
  #multiTarget.enhancers.grList <- enhancer.eqtls.grList[lapply(enhancer.eqtls.grList, length) > 1]
  multiTarget.enhancers.grList <- enhancer.eqtls.grList[lapply(enhancer.eqtls.grList, length) > 1]
  ###############################################################################################
  
  
  ###############################################################################################
  ############################ eQTLs occuring within Promoters ##################################
  ###############################################################################################
  
  # Narrow down the search by looking only at genes that are targetted by eQTLs in 
  # both Enhancer and promoter regions.
  
  # Get a genomic ranges object for promoter ranges from gencode
  # This returns eqtls that don't necessarily target genes targeted by promoter
  
  promoters.gr <- promoters(gencode.gr, downstream = 1000)
  promoter.hits <- findOverlaps(eqtls.gr, promoters.gr)
  all.eqtls.in.promoters <- eqtls.gr[queryHits(promoter.hits)] #keepCount
  #counter
  all.eqtls.in.promoter.counter <- length(all.eqtls.in.promoters)
  
  # Subset promoter eQTLs that target the gene itself
  eqtls.in.promoters <- GRanges()
  for(i in 1:length(all.eqtls.in.promoters)){
    current.eqtl <- all.eqtls.in.promoters[i]
    location <- start(current.eqtl)
    target.promoter <- promoters(gencode.gr[current.eqtl$gene_id], downstream=1000)
    if((location >= start(target.promoter)) & (location <= end(target.promoter))){
      eqtls.in.promoters <- c(eqtls.in.promoters, current.eqtl) #keepCount
    }
  }
  #counter
  p.eqtls.counter <- length(eqtls.in.promoters)
  #counter
  multiGene.p.eqtls.counter <-intersect(all.eqtls.in.promoters$eqtl_id,multiTarget.eqtl.ids)
  
  
  # Get a list of genes targeted by both enhancer and promoter eQTLs
  common.genes <- intersect(eqtls.in.enhancers$gene_id, eqtls.in.promoters$gene_id) #keepCount
  #counter
  common.genes.counter <- length(common.genes)
  
  common.promoter.eqtls <- eqtls.in.promoters[eqtls.in.promoters$gene_id %in% common.genes] #keepCount
  #counter
  peqtls.in.common.genes.counter <- length(common.promoter.eqtls)
  ###############################################################################################
  
    
  promoter.eqtl.targetPairs <- GRanges()
  #common.enhancer.eqtls <- GRangesList()
  common.enhancer.eqtls <- list()
  pair_counter <- 0
  #counter
  genepairs.oppSlopes.counter <- 0 #302
  #counter
  geneA.enhancer.counter <- 0
  #counter
  bothGenes.enhancer.counter <- 0
  
  for(i in 1:length(common.promoter.eqtls)){
    #for(i in 275:300){
    #Find other targets of each eQTL
    current.eqtl <- common.promoter.eqtls[i]
    current.target <- current.eqtl$gene_id
    current.slope <- current.eqtl$slope
    if(current.eqtl$eqtl_id %in% multiTarget.eqtl.ids){
      all.targets <- multiTarget.eqtls.grList[[current.eqtl$eqtl_id]]
      all.targets <- all.targets[-which(all.targets$gene_id %in% current.target)]
      for(j in 1:length(all.targets)){
        current.alt <- all.targets[j]  # <-The eqtl for the possible alternate target
        current.genedist <- abs(start(gencode.gr[current.target]) - start(gencode.gr[current.alt$gene_id]))
        #Check if other gene has opposite slope
        if((current.slope/current.alt$slope < 0) & (current.genedist >=5000)){
          #counter
          genepairs.oppSlopes.counter <- genepairs.oppSlopes.counter + 1
          # ## Check if current gene is in COSMIC
          # if(current.target %in% cosmic.info$gene_id) {
          #   current.eqtl$geneA_inCOSMIC <- 1
          # } else {current.eqtl$geneA_inCOSMIC <- 0}
          # ####
          # ## Check if current gene is in DO
          # if(current.target %in% brca.prad.info$gene_id) {
          #   cancer_type <- brca.prad.info$cancerType[brca.prad.info$gene_id %in% current.target]
          #   if(length(cancer_type)>1){current.eqtl$geneA_inDO <- "both"}
          #   else{current.eqtl$geneA_inDO <- cancer_type}
          # } else {current.eqtl$geneA_inDO <- NA}
          ####
          
          
          # Save eqtl details for alternate target gene
          current.eqtl$alt_target_id <- current.alt$gene_id
          current.eqtl$alt_tss_distance <- current.alt$tss_distance
          current.eqtl$alt_pvalue <- current.alt$pval_nominal
          current.eqtl$alt_slope <- current.alt$slope
          
          # ## Check if alternate gene is in COSMIC
          # if(current.alt$gene_id %in% cosmic.info$gene_id) {
          #   current.eqtl$geneB_inCOSMIC <- 1
          # } else {current.eqtl$geneB_inCOSMIC <- 0}
          # ###
          # ## Check if alternate gene is in DO
          # if(current.alt$gene_id %in% brca.prad.info$gene_id) {
          #   cancer_type <- brca.prad.info$cancerType[brca.prad.info$gene_id %in% current.alt$gene_id]
          #   if(length(cancer_type)>1){current.eqtl$geneB_inDO <- "both"}
          #   else{current.eqtl$geneB_inDO <- cancer_type}
          #   
          # } else {current.eqtl$geneB_inDO <- NA}
          ####
          
          # Assign a pair_id for cross referencing with enhancers
          pair_counter <- pair_counter + 1
          current.eqtl$pair_id <- paste("pair", pair_counter, sep="_")
          current.eqtl$common_enhancer <- "no"
          #promoter.eqtl.targetPairs <- c(promoter.eqtl.targetPairs, current.eqtl)
          
          # Check if any enhancers regulating target gene also regulates alternate target
          target.enhancer.list <- unique(eqtls.in.enhancers[eqtls.in.enhancers$gene_id %in% current.target]$enhancer_id)
          # For each enhancer, see if both genes are targeted 
          target.pair <- c(current.target, current.alt$gene_id) #Changed this vector to include both genes
          shared.enhancer.counter <- 0
          geneA.enhancer.counter <- 0
          
          common.enhancer.sublist <- GRangesList()
          for(k in 1:length(target.enhancer.list)){
            
            current.enhancer <- target.enhancer.list[k]
            enhancer.target.eqtls <- multiTarget.enhancers.grList[[current.enhancer]]
            enhancer.target.list <- unique(enhancer.target.eqtls$gene_id)
            
            
            
            ##################################################
            ## Added to annotate enhancer link to geneA, both or none
            
            #To begin, set common_enhancer variable to "none"
            current.eqtl$common_enhancer <- "none"
            
            if(target.pair[1] %in% enhancer.target.list){
              geneA.enhancer.counter <- geneA.enhancer.counter + 1
              current.eqtl$common_enhancer <- "geneA_only"
              geneA.enhancer.counter <- geneA.enhancer.counter + 1
              common.enhancer.sublist[[current.enhancer]] <- 
                enhancer.target.eqtls[which(enhancer.target.eqtls$gene_id %in% target.pair[1])]
            }
            
            ##################################################
            
            
            if(sum(target.pair %in% enhancer.target.list)==2){
              shared.enhancer.counter <- shared.enhancer.counter + 1
              # Change common_enhancer value to "yes"
              #current.eqtl$common_enhancer <- "yes"
              current.eqtl$common_enhancer <- "both_genes"
              ##############
              #counter
              geneA.enhancer.counter <- geneA.enhancer.counter - 1
              if(geneA.enhancer.counter <= 0) {geneA.enhancer.counter <- 0}
              bothGenes.enhancer.counter <- bothGenes.enhancer.counter + 1
              ##############
              #enhancer.pair.id <- paste(current.eqtl$pair_id, paste0("enh",shared.enhancer.counter),sep="_")
              common.enhancer.sublist[[current.enhancer]] <- 
                enhancer.target.eqtls[which(enhancer.target.eqtls$gene_id %in% target.pair)]
            }#common enhancer if loop end
            
          }#alt target if block end
          current.eqtl$no_shared_enhancers <- shared.enhancer.counter
          current.eqtl$no_geneA_enhancers <- geneA.enhancer.counter
          promoter.eqtl.targetPairs <- c(promoter.eqtl.targetPairs, current.eqtl)
          common.enhancer.eqtls[[current.eqtl$pair_id]] <- common.enhancer.sublist
        }#Check for opposite slopes if block end
      }#For all targets of promoter eqtl end
      
    }#Multi-target eqtl if block end
    
  }#Promoter eqtl loop end
  
  
  #### Adding Information about the Slope, etc
  
   #Get only the genes we are interested in for this tissue type
   gene.subset.list <- unique(c(promoter.eqtl.targetPairs$gene_id, promoter.eqtl.targetPairs$alt_target_id))
   gene.subset.inds <- which(rnaseq.data$Name %in% gene.subset.list)
  
  
   tissue.rnaseq <- rnaseq.data %>% dplyr::select(one_of(rna.sampids.tissue))
   tissue.rnaseq <- tissue.rnaseq[gene.subset.inds,]
   tissue.rnaseq <- as_tibble(t(tissue.rnaseq))
  
   names(tissue.rnaseq) <- as.character(names(gene.names[gene.subset.inds]))
  
  ####################################################################################

  #### Adding Information about the Slope
  expr.corr <- numeric(length(promoter.eqtl.targetPairs))
   for(i in 1:length(promoter.eqtl.targetPairs)){
      paired.eqtls <- promoter.eqtl.targetPairs[i]
  
      pro.gene <- paired.eqtls$gene_id
      alt.gene <- paired.eqtls$alt_target_id
  
      expr.corr[i] <- cor(tissue.rnaseq[,pro.gene],tissue.rnaseq[,alt.gene])

      promoter.eqtl.targetPairs$rnaseq_corr <- expr.corr
   }



  ####################################################################################
  ################# ADDING INDIVIDUAL E-EQTLs TO EACH ROW OF PAIR ####################
  ####################################################################################

   results.table <- data.frame(pair_id=character(),geneA_id=character(),geneA_name=character(),
                               geneB_id=character(), geneB_name=character(),
                               rnaseq_corr=numeric(),
                               geneA_COSMIC=numeric(),geneA_DO=character(),
                               geneB_COSMIC=numeric(),geneB_DO=character(),
                               peqtl_id=character(), gtex_maf=numeric(),#peqtl_rsid=character(),peqtl_CADD_score=numeric(),
                               peqtl_slopeA=numeric(), peqtl_pvalueA=numeric(),peqtl_tss_dist_A=numeric(),
                               peqtl_slopeB=numeric(), peqtl_pvalueB=numeric(),peqtl_tss_dist_B=numeric(),
                               no_shared_enhancers=numeric(),no_geneA_enhancers=numeric(),enhancer_id=character(),e_eqtl_geneA=character(),stringsAsFactors = FALSE)
                               #e_eqtl_geneB=character(),stringsAsFactors = FALSE)



   for(i in 1:length(promoter.eqtl.targetPairs)){
     current.pair <- promoter.eqtl.targetPairs[i]
  
  
  
     current.geneA <- current.pair$gene_id
     current.geneB <- current.pair$alt_target_id
  
  
     current.common.enhancers <- common.enhancer.eqtls[[current.pair$pair_id]]
     for(each.enhancer in names(current.common.enhancers)){
       which.geneA <- which(current.common.enhancers[[each.enhancer]]$gene_id %in% current.geneA)
       geneA.e_eqtls <- current.common.enhancers[[each.enhancer]][which.geneA]
       selected.geneA.e_eqtl <- geneA.e_eqtls[which.min(geneA.e_eqtls$pval_nominal)]
  
       # which.geneB <- which(current.common.enhancers[[each.enhancer]]$gene_id %in% current.geneB)
       # geneB.e_eqtls <- current.common.enhancers[[each.enhancer]][which.geneB]
       # selected.geneB.e_eqtl <- geneB.e_eqtls[which.min(geneB.e_eqtls$pval_nominal)]
  
       # Arrange the information for pEQTL in the preferred order for table (see later for column names)
       results.table <- results.table %>%
         add_row(pair_id=current.pair$pair_id, geneA_id=current.geneA, geneA_name=gencode.info$Name[gencode.info$ID %in% current.geneA],
                 geneB_id=current.geneB, geneB_name=gencode.info$Name[gencode.info$ID %in% current.geneB],
                 rnaseq_corr=current.pair$rnaseq_corr,
                 geneA_COSMIC=current.pair$geneA_inCOSMIC, geneA_DO=current.pair$geneA_inDO,
                 geneB_COSMIC=current.pair$geneB_inCOSMIC, geneB_DO=current.pair$geneB_inDO,
                 peqtl_id=current.pair$eqtl_id, gtex_maf=current.pair$maf,#peqtl_rsid=current.pair$pEQTL_rsid, peqtl_CADD_score=current.pair$cadd_scores,
                 peqtl_slopeA=current.pair$slope, peqtl_pvalueA=current.pair$pval_nominal, peqtl_tss_dist_A=current.pair$tss_distance,
                 peqtl_slopeB=current.pair$alt_slope, peqtl_pvalueB=current.pair$alt_pvalue, peqtl_tss_dist_B=current.pair$alt_tss_distance,
                 no_shared_enhancers=current.pair$no_shared_enhancers,no_geneA_enhancers=current.pair$no_geneA_enhancers,enhancer_id=each.enhancer, e_eqtl_geneA=selected.geneA.e_eqtl$eqtl_id)#, e_eqtl_geneB=selected.geneB.e_eqtl$eqtl_id)
  
  
     }
   }

   results.table <- tbl_df(results.table)
  
  ####################################################################################
  ####################################################################################
  
  ## Create a output directory for the tissue in "outputFiles" with the tissue name
  dir.create(paste0("./output/",current.tissue))
  output.filename <- paste0("./output/",current.tissue,"/",current.tissue,
                            "_allPE_eqtlTargetPairs_commonEnhancers.Rdat")

  result.filename <- paste0("./output/",current.tissue,"/",current.tissue,
                            "_allPE_eqtlTargetPairs_commonEnhancers_Table.txt")

  save(common.enhancer.eqtls, promoter.eqtl.targetPairs, file=output.filename)
  write_tsv(results.table, result.filename)
  
  ####################################################################################
  ###########################  COUNTS TABLE FOR OUTPUT ###############################
  ####################################################################################

  counts.table <- counts.table %>%
    add_row(Tissue=current.tissue, enhancer_eqtl_counter=enhancer.eqtl.counter, multiTarget_eqtl_count=multitarget.eqtl.counter,
            all_eqtls_in_promoter_count=all.eqtls.in.promoter.counter, p_eqtls_count=p.eqtls.counter, multiGene_pEqtls_count=multiGene.p.eqtls.counter,
            common_genes_count=common.genes.counter, peqtls_in_commonGenes_count=peqtls.in.common.genes.counter, genePairs_oppSlopes_count=genepairs.oppSlopes.counter,
            geneA_enhancer_count=geneA.enhancer.counter, bothGenes_enhancer_count=bothGenes.enhancer.counter)

  output.filename <- paste0("./",current.tissue,"_countsTable_results.txt")
  ####################################################################################
  ####################################################################################
  

  write_tsv(counts.table, output.filename)   
}

