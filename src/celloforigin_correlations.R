#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
samplelist=args[1]
samples=fread(samplelist, header=F)

#tabula sapiens + vento tormo (batch corrected values) - broad cell types
#proteinAtlas <- read.table('ref/tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_broad_cell_types_vtannotated_normbatchcorrected_expression_data.txt', header=T,as.is=T,sep="\t")

#tabula sapiens + vento tormo (batch corrected values) - broad cell types stratified by tissue
#proteinAtlas <- read.table('ref/tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_tissue_general_cell_types_vtannotated_normbatchcorrected_expression_data.txt', header=T,as.is=T,sep="\t")

proteinAtlasfile = args[2]
proteinAtlas <- read.table(proteinAtlasfile, header=T,as.is=T,sep="\t")

#batch corrected values are not normally distributed (perform log normalization)
rownames(proteinAtlas) <- proteinAtlas$feature_id
ndata <- proteinAtlas[,-1]
#logndata <- log2(ndata) #originally what snyder did as bulk rna-seq was not log transformed
#dim(logndata)

tLabels <- as.data.table(colnames(ndata))
colnames(tLabels) <- c("RName")
fftColumns <- 29:52 # 160-222
selFreq <- c("193","196","199")

for (i in 1:nrow(samples)){
        tryCatch({
sample=samples[i]$V1
print(sample)

file_path <- Sys.glob(paste0("body/", sample, "/fft_summaries/fft_", sample, "*_WPS.tsv.gz"))
fdata <- read.table(file_path, as.is=TRUE, sep="\t", header=TRUE, comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    ndata2 <- ndata[fdata[,1],]

res <- cor(rowMeans(fdata[,selFreq]),ndata2[,order(names(ndata2))],use="pairwise.complete.obs")
    res.df <- as.data.frame(t(res))
    colnames(res.df) <- c("correlation")
    res.df$RName <- rownames(res.df)
    res.df <- merge(res.df, tLabels, by="RName")
    res.df <- res.df[order(res.df$correlation),]
    res.df$rank <- rep(1:nrow(res.df))

    sampleshort=substr(sample, 0, 14)
    write.table(res.df, paste0("body/", sample, "/fft_summaries/", sampleshort, "_tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbatch_broad_cell_types_vtannotated_normbatchcorrected_expression_Ave193-199bp_correlation.csv"), sep=",", quote=F, row.names=F, col.names=T)
}, error=function(e){})
        }
