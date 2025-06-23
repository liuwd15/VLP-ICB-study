# For some rare cases, a few CD4+T cells are in the CD8+ clusters or vice verse. 
# This step make sure that cells with same clonotype have same CD4 or CD8 annotation.
mode <- function(x) names(sort(-table(x)))[1]
correct_cd_subset <- function(tcr_obj)
{
    expanded_tcr <- tcr_obj@meta.data %>% group_by(CTaa, orig.ident) %>% summarize(count=n()) %>% 
        filter(count>1 & !is.na(CTaa)) %>% select(CTaa, orig.ident) 
    for(i in 1:nrow(expanded_tcr))
    {
        cells_tcr_i <- which(tcr_obj$CTaa == expanded_tcr$CTaa[i] & tcr_obj$orig.ident == expanded_tcr$orig.ident[i])
        cells_tcr_i_cd8 <- which(tcr_obj$CTaa == expanded_tcr$CTaa[i] & tcr_obj$orig.ident == expanded_tcr$orig.ident[i] & tcr_obj$CDsubset=='CD8')
        cd_exp <- FetchData(tcr_obj, c('Cd8a','Cd8b1','Cd4','Foxp3'), cells=rownames(tcr_obj@meta.data)[cells_tcr_i])
        if(sum(cd_exp[,'Cd8a']+cd_exp[,'Cd8b1']>0) > sum(cd_exp[,'Cd4']+cd_exp[,'Foxp3']>0) | length(cells_tcr_i_cd8)/length(cells_tcr_i)>0.5)
        {
            row_matched_tcr <- tcr_obj@meta.data[cells_tcr_i,] 
            row_matched_tcr$CDsubset <- 'CD8'
            all_annotations <- row_matched_tcr$annotation
            if(sum(str_detect(all_annotations, 'CD4\\+')>0))
                row_matched_tcr$annotation[str_detect(all_annotations, 'CD4\\+')] <- mode(all_annotations[str_detect(all_annotations, '(CD8\\+|CD4/8\\+)')])
            tcr_obj@meta.data[cells_tcr_i,] <- row_matched_tcr 
        }
        else if(sum(cd_exp[,'Cd8a']+cd_exp[,'Cd8b1']>0) < sum(cd_exp[,'Cd4']+cd_exp[,'Foxp3']>0) | length(cells_tcr_i_cd8)/length(cells_tcr_i)<0.5 |
               sum(cd_exp[,'Cd8a']+cd_exp[,'Cd8b1']) <= sum(cd_exp[,'Cd4']+cd_exp[,'Foxp3']))
        {
            row_matched_tcr <- tcr_obj@meta.data[cells_tcr_i,] 
            row_matched_tcr$CDsubset <- 'CD4'
            all_annotations <- row_matched_tcr$annotation
            if(sum(str_detect(all_annotations, 'CD8\\+')>0))
                row_matched_tcr$annotation[str_detect(all_annotations, 'CD8\\+')] <- mode(all_annotations[str_detect(all_annotations, '(CD4\\+|CD4/8\\+)')])
            tcr_obj@meta.data[cells_tcr_i,] <- row_matched_tcr 
        }
    }

    return(tcr_obj)
}

# Split single-cell duel TCRs into two TCRs
split_duel_tcr <- function(df)
{
    new_df <- df
    for(i in 1:nrow(df))
    {
        if(str_detect(df$cdr3a[i], ';'))
        {
            a1 <- str_split(df$cdr3a[i], ';')[[1]]
            new_df[i, 'cdr3a'] <- a1[1]
            newrow <- nrow(new_df)
            new_df[newrow+1, ] <- df[i, ]
            new_df[newrow+1, 'cdr3a'] <- a1[2]
        }
        if(str_detect(df$cdr3b[i], ';'))
        {
            b1 <- str_split(df$cdr3b[i], ';')[[1]]
            new_df[i, 'cdr3b'] <- b1[1]
            newrow <- nrow(new_df)
            new_df[newrow+1, ] <- df[i, ]
            new_df[newrow+1, 'cdr3b'] <- b1[2]
        }
    }
    new_df <- new_df %>% arrange(peptide, cdr3a, cdr3b) %>% distinct(peptide, cdr3a, cdr3b, mhc, .keep_all = TRUE)
    return(new_df)
}

# mapping TCR genes to IMGT alleles for B16F10 mice
tcr_gene2allele <- function(df)
{
    #https://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=Mus_musculus&group=TRAV
    trav_alleles <- c('TRAV1*02','TRAV10*02','TRAV10D*03','TRAV10N*01','TRAV11*01','TRAV12-1*06','TRAV12-2*02','TRAV12-3*05',
                      'TRAV12D-3*04','TRAV12N-1*01','TRAV12N-2*01','TRAV13-1*02','TRAV13-2*03',
                      'TRAV13-4/DV7*02','TRAV13D-1*04','TRAV13D-3*02','TRAV13D-4*03','TRAV13N-1*01','TRAV13N-2*01','TRAV13N-3*01','TRAV13N-4*01',
                      'TRAV14-1*01','TRAV14-2*02','TRAV14-3*04','TRAV14D-3/DV8*08','TRAV14N-1*01','TRAV14N-2*01','TRAV14N-3*01',
                      'TRAV15-1/DV6-1*01','TRAV15-2/DV6-2*02','TRAV15D-1/DV6D-1*03','TRAV15N-2*01','TRAV16*06','TRAV16D/DV11*02','TRAV16N*01','TRAV17*02',
                      'TRAV18*02','TRAV19*03',
                      'TRAV2*01','TRAV21/DV12*03',
                      'TRAV3-1*03','TRAV3-3*01','TRAV3-4*02','TRAV3N-3*01',
                      'TRAV4-2*02','TRAV4-3*03','TRAV4-4/DV10*02','TRAV4N-3*01','TRAV4N-4*01',
                      'TRAV5-1*02','TRAV5-4*02','TRAV5N-4*01',
                      'TRAV6-1*03','TRAV6-2*02','TRAV6-3*02','TRAV6-4*03','TRAV6-5*04','TRAV6-7/DV9*04','TRAV6D-3*03','TRAV6D-4*01','TRAV6D-5*02',
                      'TRAV6D-6*02','TRAV6N-5*01','TRAV6N-6*01','TRAV6N-7*01',
                      'TRAV7-1*02','TRAV7-3*04','TRAV7-4*02','TRAV7-5*03','TRAV7-6*02','TRAV7D-2*02','TRAV7D-3*03','TRAV7D-4*01',
                      'TRAV7N-4*01','TRAV7N-5*01','TRAV7N-6*01',
                      'TRAV8-1*03','TRAV8D-1*01','TRAV8D-2*03','TRAV8N-2*01',
                      'TRAV9-1*02','TRAV9-2*02','TRAV9-4*02','TRAV9D-1*03','TRAV9D-4*03','TRAV9N-2*01','TRAV9N-3*01','TRAV9N-4*01')
    names(trav_alleles) <- trav_alleles %>% str_replace('\\*[0-9]+','') %>% str_replace('/','-')
    #https://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=Mus_musculus&group=TRBV
    trbv_alleles <- c('TRBV1*01','TRBV10*01','TRBV12-1*01','TRBV12-2*01','TRBV13-1*02','TRBV13-2*01','TRBV13-3*01','TRBV14*01','TRBV15*01',
                      'TRBV16*01','TRBV17*01','TRBV19*01','TRBV2*01','TRBV20*01','TRBV23*01','TRBV26*01','TRBV29*01',
                      'TRBV3*01','TRBV30*01','TRBV31*01','TRBV5*01')
    names(trbv_alleles) <- str_replace(trbv_alleles,'\\*[0-9]+','')
    df$va <- trav_alleles[df$va]
    df$vb <- trbv_alleles[df$vb]
    return(df)
}

# assigning the pMTnet prediction precentile rank and logit to single cells
assign_percentile_logit <- function(tcr_obj, pmtnet_out_cd8=NULL, pmtnet_out_cd4=NULL, antigen='Ag')
{
    col_percentile <- paste0('percentile_', antigen)
    col_logit <- paste0('logit_', antigen)
    tcr_obj@meta.data[[col_percentile]] <- Inf
    tcr_obj@meta.data[[col_logit]] <- -Inf
    tcr_obj@meta.data[['MHC']] <- NA

    if(!is.null(pmtnet_out_cd8))
    {
        for(i in 1:nrow(pmtnet_out_cd8))
        {
            matched_tcr_cd8 <- which(str_detect(tcr_obj$CTaa, pmtnet_out_cd8$cdr3a[i]) & str_detect(tcr_obj$CTaa, pmtnet_out_cd8$cdr3b[i]))
            if( length(matched_tcr_cd8)>0 & mean(tcr_obj$CDsubset[matched_tcr_cd8] == 'CD8')>0.5)
            {
                if(max(tcr_obj@meta.data[[col_percentile]][matched_tcr_cd8]) > pmtnet_out_cd8$percentile_rank_0[i])
                {
                    tcr_obj@meta.data[[col_percentile]][matched_tcr_cd8] <- pmtnet_out_cd8$percentile_rank_0[i]
                    tcr_obj@meta.data[[col_logit]][matched_tcr_cd8] <- pmtnet_out_cd8$logit[i]
                    tcr_obj@meta.data[['MHC']][matched_tcr_cd8] <- 'MHC-I'
                }
            }
        }
    }
    if(!is.null(pmtnet_out_cd4))
    {
        for(i in 1:nrow(pmtnet_out_cd4))
        {
            matched_tcr_cd4 <- which(str_detect(tcr_obj$CTaa, pmtnet_out_cd4$cdr3a[i]) & str_detect(tcr_obj$CTaa, pmtnet_out_cd4$cdr3b[i]))
            if( length(matched_tcr_cd4)>0 & mean(tcr_obj$CDsubset[matched_tcr_cd4] == 'CD4')>0.5)
            {
                if(max(tcr_obj@meta.data[[col_percentile]][matched_tcr_cd4]) > pmtnet_out_cd4$percentile_rank_0[i])
                {
                    tcr_obj@meta.data[[col_percentile]][matched_tcr_cd4] <- pmtnet_out_cd4$percentile_rank_0[i]
                    tcr_obj@meta.data[[col_logit]][matched_tcr_cd4] <- pmtnet_out_cd4$logit[i]
                    tcr_obj@meta.data[['MHC']][matched_tcr_cd4] <- 'MHC-II'
                }
            }
        }
    }
    tcr_obj@meta.data[[col_percentile]][is.infinite(tcr_obj@meta.data[[col_percentile]])] <- NA
    tcr_obj@meta.data[[col_logit]][is.infinite(tcr_obj@meta.data[[col_logit]])] <- NA

    return(tcr_obj)
}