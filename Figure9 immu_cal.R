#setwd('/pub1/data/mg_projects/projects/web_script/R/')
library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/37082fba10fc297522cb3156cfe1113f/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/37082fba10fc297522cb3156cfe1113f',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()
tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "GEO Data press"))
  logs=c(logs,paste0('run immu_cal.R-',basename(opt$outfile)))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 100)
  exp_path=unlist(data$exp_path)
  cancerCode=unlist(data$cancerCode)
  array=as.numeric(unlist(data$array))#0,1=array
  method=unlist(data$method)
  #method='timer'
  #array=1
  #cancerCode='STAD'
  dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                        ,na.strings="NA",data.table = F)
  
  row.names(dat)=dat[,1]
  dat=dat[,-1]
  library(hgu133plus2.db)
  if(sum(row.names(dat)%in%mappedRkeys(hgu133plus2SYMBOL))<2000){
    logs=c(logs,paste0('error:Too few genes'))
  }else{
    library(IOBR)
    deconvolute_timer.default<-function (args) 
    {
      cancers = check_cancer_types(args)
      TimerINFO("Loading immune gene expression")
      immune <- immuneCuratedData
      immune.geneExpression <- immune$genes
      immune.cellTypes <- immune$celltypes
      outlier.genes <- sort(GetOutlierGenes(cancers))
      print(paste("Outlier genes:", paste(outlier.genes, collapse = " ")))
      dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
      if (!dir.exists(paste(args$outdir, "/results", sep = ""))) {
        dir.create(paste(args$outdir, "/results", sep = ""))
      }
      abundance.score.matrix <- c()
      #pdf(paste(args$outdir, "/results/output.pdf", sep = ""))
      for (i in 1:nrow(cancers)) {
        cancer.expFile <- cancers[i, 1]
        cancer.category <- cancers[i, 2]
        cancer.expression <- ParseInputExpression(cancer.expFile)
        index <- !(row.names(cancer.expression) %in% outlier.genes)
        cancer.expression <- cancer.expression[index, , drop = FALSE]
        cancer.colnames <- colnames(cancer.expression)
        TimerINFO(paste("Removing the batch effect of", cancer.expFile))
        #for (j in 1:length(cancer.colnames)) {
        #  DrawQQPlot(cancer.expression[, j], immune.geneExpression[, 
        #                                                           1], name = cancer.colnames[j])
        #}
        tmp <- RemoveBatchEffect(cancer.expression, immune.geneExpression, 
                                 immune.cellTypes)
        cancer.expNorm <- tmp[[1]]
        immune.expNormMedian <- tmp[[3]]
        #for (j in 1:length(cancer.colnames)) {
        #  DrawQQPlot(cancer.expNorm[, j], immune.expNormMedian[, 
        #                                                       1], name = paste("After batch removing and aggregating for", 
        #                                                                        cancer.colnames[j]))
        #}
        gene.selected.marker <- cancer_type_genes[[which(names(cancer_type_genes) == 
                                                           cancer.category)]]
        gene.selected.marker <- intersect(gene.selected.marker, 
                                          row.names(cancer.expNorm))
        XX = immune.expNormMedian[gene.selected.marker, c(-4)]
        YY = cancer.expNorm[gene.selected.marker, , drop = FALSE]
        for (j in 1:length(cancer.colnames)) {
          fractions <- GetFractions.Abbas(XX, YY[, j])
          #barplot(fractions, cex.names = 0.8, names.arg = names(fractions), 
          #        xlab = "cell type", ylab = "abundance", main = paste("Abundance estimation for", 
          #                                                             cancer.colnames[j]))
          #box()
          abundance.score.matrix <- cbind(abundance.score.matrix, 
                                          fractions)
          colnames(abundance.score.matrix)[ncol(abundance.score.matrix)] <- cancer.colnames[j]
        }
      }
      #dev.off()
      write.table(abundance.score.matrix, paste(args$outdir, "/results/score_matrix.txt", 
                                                sep = ""), sep = "\t", quote = FALSE, row.names = TRUE, 
                  col.names = NA)
      return(abundance.score.matrix)
    }
    deconvo_timer<-function (eset, project = NULL, indications = NULL) 
    {
      indications = tolower(indications)
      checkmate::assert("indications fit to mixture matrix", length(indications) == 
                          ncol(eset))
      args = new.env()
      args$outdir = tempdir()
      args$batch = tempfile()
      lapply(unique(indications), function(ind) {
        tmp_file = tempfile()
        tmp_mat = eset[, indications == ind, drop = FALSE] %>% 
          as_tibble(rownames = "gene_symbol")
        readr::write_tsv(tmp_mat, tmp_file)
        cat(paste0(tmp_file, ",", ind, "\n"), file = args$batch, append = TRUE)
      })
      results <- deconvolute_timer.default(args)[, make.names(colnames(eset))]
      colnames(results) <- colnames(eset)
      results <- as.data.frame(t(results))
      colnames(results) <- paste(colnames(results), "_TIMER", sep = "")
      colnames(results) <- gsub(colnames(results), pattern = "\\.", replacement = "\\_")
      colnames(results) <- gsub(colnames(results), pattern = "\\ ", replacement = "\\_")
      if (!is.null(project)) {
        results$project <- project
        results <- results[, c(ncol(results), 1:ncol(results) - 1)]
      }
      results <- tibble::rownames_to_column(results, var = "ID")
      return(results)
    }
    array=(array==1)
    tp=cancerCode
    tp=tolower(tp)
    fst_exp=as.matrix(dat)
    logs=c(logs,paste0('start run immu,method=',method))
    if(method=='timer'|method=='T'){
      logs=c(logs,paste0('run timer'))
      imu<-deconvo_timer(eset=fst_exp,indications = rep(tp,ncol(fst_exp)))
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_TIMER$'),'',colnames(imu))
    }else if(method=='quantiseq'|method=='Q'){
      logs=c(logs,paste0('run quantiseq'))
      imu<-deconvo_quantiseq(fst_exp,tumor=TRUE,arrays = array,scale_mrna=T)
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_quantiseq$'),'',colnames(imu))
    }else if(method=='mcpcounter'|method=='M'){
      logs=c(logs,paste0('run mcpcounter'))
      imu<-deconvo_mcpcounter(fst_exp)
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_MCPcounter$'),'',colnames(imu))
    }else if(method=='estimate'|method=='ES'){
      logs=c(logs,paste0('run estimate'))
      ptf='illumina'
      if(array) ptf='affymetrix'
      imu<-deconvo_estimate(fst_exp,platform=ptf)
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_estimate$'),'',colnames(imu))
    }else if(method=='ips'|method=='I'){
      logs=c(logs,paste0('run ips'))
      imu<-deconvo_ips(fst_exp,plot = F)
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_IPS$'),'',colnames(imu))
    }else if(method=='epic'|method=='EP'){
      logs=c(logs,paste0('run epic'))
      imu<-deconvo_epic(eset = fst_exp,tumor = TRUE)
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_EPIC$'),'',colnames(imu))
    }else if(method=='xcell'|method=='X'){
      logs=c(logs,paste0('run xcell'))
      imu<-deconvo_xcell(fst_exp,arrays = array)
      logs=c(logs,paste0('run end'))
      colnames(imu)=gsub(paste0('_xCell$'),'',colnames(imu))
    }else if(method=='cibersort'|method=='C'){
      logs=c(logs,paste0('run cibersort'))
      fst_exp1=fst_exp[,which(apply(fst_exp[row.names(fst_exp)%in%row.names(lm22),], 2,sd)>0)]
      cibersort_result<-deconvo_cibersort(eset = fst_exp1,arrays = FALSE,absolute = FALSE, perm = 1000)
      logs=c(logs,paste0('run end'))
      narn=setdiff(colnames(fst_exp),colnames(fst_exp1))
      if(length(narn)>0){
        nam<-cbind(narn,matrix(rep(NA,25*length(narn)),ncol = 25))
        colnames(nam)=colnames(cibersort_result)
        cibersort_result=rbind(cibersort_result,nam)
      }
      imu=cibersort_result
    }
    logs=c(logs,paste0('outputing immu result'))
    colnames(imu)=gsub(paste0('_',toupper(method),'$'),'',colnames(imu))
    write.table(imu,file = paste0(opt$outfile,'/immueScore.txt')
                ,row.names = F,col.names = T,quote = F,sep = '\t')
    logs=c(logs,paste0('outputed'))
  }
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})


