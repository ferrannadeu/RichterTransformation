snvsClassification <- function(mat, context=3, contextColumn=NULL, mutTypeColumn=NULL, title="Classification of SNVs", subtitle=TRUE, yScale="percentage"){
  if(is.null(contextColumn)){
    mat <- mat[,1:4]
    colnames(mat) <- c("CHR", "POSITION", "REF", "ALT")
    if(!startsWith(as.character(mat$CHR[1]), "chr")){ mat$CHR <- paste0("chr", mat$CHR) }
  }
  else{
    mat <- as.data.frame(mat)
    mat1 <- mat[,1:4]
    mat1 <- cbind(mat1, mat[,  c(contextColumn, mutTypeColumn)])
    mat <- mat1
    colnames(mat) <- c("CHR", "POSITION", "REF", "ALT", "context", "mutType")
  }
  
  if(context==3){ s <- 0 }
  else if(context==5){ s <- 1 }
  else{ stop("Context must be 3 or 5!") }
  
  mat <- mat[nchar(mat$REF) == 1 & nchar(mat$ALT)==1, ] # remove potential non SNVs in matrix
  if(nrow(mat)==0){ p <- ggplot(NULL) }
  
  else{
    if(is.null(contextColumn)){
      gr1 <- GRanges(paste0(mat$CHR),IRanges(start=as.numeric(mat$POSITION)-1-s, end=as.numeric(mat$POSITION)+1+s), strand = "+")
  
      refbase <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr1)
      refbase <- as.data.frame(refbase)$x
      mat$context <- refbase
      
      # Convert them to the 6 types of base substitution types that are distinguished by convention:
      mat$mutType <- NA
      
      mat$mutType[mat$REF == "C" & mat$ALT == "A"] <- "C>A"
      mat$mutType[mat$REF == "C" & mat$ALT == "G"] <- "C>G"
      mat$mutType[mat$REF == "C" & mat$ALT == "T"] <- "C>T"
      mat$mutType[mat$REF == "T" & mat$ALT == "A"] <- "T>A"
      mat$mutType[mat$REF == "T" & mat$ALT == "C"] <- "T>C"
      mat$mutType[mat$REF == "T" & mat$ALT == "G"] <- "T>G"
      
      mat$mutType[mat$REF == "G" & mat$ALT == "T"] <- "C>A"
      mat$context[mat$REF == "G" & mat$ALT == "T"] <- reverse(chartr("ATGC", "TACG", mat$context[mat$REF == "G" & mat$ALT == "T"]))
      mat$mutType[mat$REF == "G" & mat$ALT == "C"] <- "C>G"
      mat$context[mat$REF == "G" & mat$ALT == "C"] <- reverse(chartr("ATGC", "TACG", mat$context[mat$REF == "G" & mat$ALT == "C"]))
      mat$mutType[mat$REF == "G" & mat$ALT == "A"] <- "C>T"
      mat$context[mat$REF == "G" & mat$ALT == "A"] <- reverse(chartr("ATGC", "TACG", mat$context[mat$REF == "G" & mat$ALT == "A"]))
      mat$mutType[mat$REF == "A" & mat$ALT == "T"] <- "T>A"
      mat$context[mat$REF == "A" & mat$ALT == "T"] <- reverse(chartr("ATGC", "TACG", mat$context[mat$REF == "A" & mat$ALT == "T"]))
      mat$mutType[mat$REF == "A" & mat$ALT == "G"] <- "T>C"
      mat$context[mat$REF == "A" & mat$ALT == "G"] <- reverse(chartr("ATGC", "TACG", mat$context[mat$REF == "A" & mat$ALT == "G"]))
      mat$mutType[mat$REF == "A" & mat$ALT == "C"] <- "T>G"
      mat$context[mat$REF == "A" & mat$ALT == "C"] <- reverse(chartr("ATGC", "TACG", mat$context[mat$REF == "A" & mat$ALT == "C"]))
      
    }
    
    nuc <- c("A", "T", "C", "G")
    if(context == 3){ 
      possibilities <- sort(do.call(paste0, expand.grid(nuc, nuc))) 
      possibilities <- paste0(substr(possibilities, start = 1, stop = 1+s), c(rep("C",16*3), rep("T",16*3)), substr(possibilities, start = 2+s, stop = 2+s*2))
    }
    else if(context == 5){ 
      possibilities <- sort(do.call(paste0, expand.grid(nuc, nuc, nuc, nuc))) 
      possibilities <- paste0(substr(possibilities, start = 1, stop = 1+s), c(rep("C",256*3), rep("T",256*3)), substr(possibilities, start = 2+s, stop = 2+s*2))
    }
    
    mmat <- data.frame(matrix(c(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each=ifelse(context==3, 16, 256)), possibilities), ncol = 2))
    colnames(mmat) <- c("mutType", "context")
    mmat <- merge(mmat, melt(table(mat[, c("mutType", "context")])), all.x=T)
    colnames(mmat) <- c("mutType", "context", "count")
    mmat$count <- as.numeric(as.character(mmat$count))
    mmat$count[is.na(mmat$count)] <- 0
    
    
    if(subtitle){ subtitle <-  paste0("SNVs input = ", nrow(mat), " | SNVs annotated = ", sum(mmat$count)) }
    else{ subtitle <- NULL }
    if(yScale== "count"){ yLab <- "Count" }
    else{ 
      yLab <- "Percentage (%)"
      mmat$count <- mmat$count/sum(mmat$count)*100
    }
    
    if(context==3){
      p<-ggplot(data=mmat, aes(x=context, y=count, fill=mutType, width=0.66)) +
        geom_bar(stat="identity", color="black") + facet_grid(.~mutType, scales = "free_x", space="free") + theme_classic() +
        scale_y_continuous(expand = c(0,0), limits = c(0, max(mmat$count)*1.05)) +
        theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 8)) +
        scale_fill_manual(values=c("#03BCEE", "#010101", "#E32926", "#999999", "#A1CE63", "#EBC6C4")) +
        labs(title=title, subtitle=subtitle, x=NULL, y=yLab)
    }
    
    else if(context==5){
      mmat$context <- paste0(substr(mmat$context,1,2), "[-]", substr(mmat$context,4,5))
      
      p<-ggplot(data=mmat, aes(x=context, y=count, fill=mutType, width=0.66)) +
        geom_bar(stat="identity", color="black") + facet_grid(mutType~., scales = "free_x", space="free") + theme_classic() +
        scale_y_continuous(expand = c(0,0), limits = c(0, max(mmat$count)*1.05)) +
        theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 5)) +
        scale_fill_manual(values=c("#03BCEE", "#010101", "#E32926", "#999999", "#A1CE63", "#EBC6C4")) +
        labs(title=title, subtitle=subtitle, x=NULL, y=yLab)
    }
  }
  return(list(mat, p))
}

fitting_sigs3 <- function(cancer_signatures, threshold, mut_mat, all_sig, patient_list=NULL, threshold_add=NULL, mutated_cases=NULL, addsigs=TRUE) { 

  cancer_signaturesPresent <- cancer_signatures[,colnames(cancer_signatures) %in% all_sig]
  cancer_signaturesPresent_bk <- cancer_signaturesPresent
  num <- ncol(cancer_signaturesPresent_bk) + 2
  
  contributionMatrix <- data.frame(matrix(0, ncol = num, nrow = ncol(mut_mat)))
  colnames(contributionMatrix) <- c("Case", "Accuracy", colnames(cancer_signaturesPresent_bk))
  
  for (i in 1:ncol(mut_mat)) {
    
    if(is.null(patient_list)){ 
      signaturesToStart <- cancer_signaturesPresent_bk
      signaturesNotPresent <- NULL
    } else{
      # For clonal fitting, start with detected signatures within current patient
      case <- strsplit(colnames(mut_mat[,i,drop=FALSE]), "_")[[1]][1]
      sigPatient <- patient_list[[case]]
      signaturesToStart <- cancer_signaturesPresent_bk[,colnames(cancer_signaturesPresent_bk) %in% sigPatient]
      signaturesNotPresent <- cancer_signaturesPresent_bk[,!colnames(cancer_signaturesPresent_bk) %in% sigPatient]
    }
    cancer_signaturesPresent <- signaturesToStart
    
    # Fit signatures to current sample
    fit_j <- fit_to_signatures(as.matrix(mut_mat[,i],1), as.matrix(cancer_signaturesPresent))
    v <- fit_j$contribution[fit_j$contribution[,1] > 0]
    names(v) <- rownames(fit_j$contribution)[fit_j$contribution[,1] > 0]
    
    cancer_signaturesPresent <- as.matrix(cancer_signaturesPresent[, names(v)])
    colnames(cancer_signaturesPresent) <- names(v) 
    
    c <- cosineSimilarity(mut_mat[,i], fit_j$reconstructed[,1])
    
    # While there are signatures to remove...
    if(length(v)>1){
    for(z in 1:(length(v)-1)){
      cosReduction <- NULL
      # Analyse cosine similarities for signatures removing one in every round
      for(j in colnames(cancer_signaturesPresent)){
        # Calculate cos sim reduction if we remove signature j...
        signs <- matrix(cancer_signaturesPresent[, colnames(cancer_signaturesPresent)!=j], ncol=ncol(cancer_signaturesPresent)-1)
        colnames(signs) <- colnames(cancer_signaturesPresent)[colnames(cancer_signaturesPresent)!=j]
        fit_res <- fit_to_signatures(matrix(mut_mat[,i], ncol=1), signs)
        temp2Cos <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
        cosReduction <- c(cosReduction, c-temp2Cos)
      }
      names(cosReduction) <- colnames(cancer_signaturesPresent)
      
      # After removing all possible signatures, if the minimum cosine similarity reduction is below the input threshold, 
      # we remove the corresponding mutational signature
      if(min(cosReduction) < threshold){ 
        clnames <- colnames(cancer_signaturesPresent)
        cancer_signaturesPresent <- matrix(cancer_signaturesPresent[, colnames(cancer_signaturesPresent) != names(cosReduction)[which.min(cosReduction)]], ncol=ncol(cancer_signaturesPresent)-1)
        colnames(cancer_signaturesPresent) <- clnames[clnames != names(cosReduction)[which.min(cosReduction)]]
        if (is.null(signaturesNotPresent)){
          signaturesNotPresent <- cancer_signaturesPresent_bk[,colnames(cancer_signaturesPresent_bk)==names(cosReduction)[which.min(cosReduction)],drop=FALSE]
        } else {
          signaturesNotPresent <- cbind(signaturesNotPresent,
                                      cancer_signaturesPresent_bk[,colnames(cancer_signaturesPresent_bk)==names(cosReduction)[which.min(cosReduction)],drop=FALSE])
        }
        fit_res <- fit_to_signatures(matrix(mut_mat[,i], ncol=1), cancer_signaturesPresent)
        c <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
      } else{  
        break 
      }
    }
    }
    
    fit_res <- fit_to_signatures(as.matrix(mut_mat[,i],1), as.matrix(cancer_signaturesPresent))
    c <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
    
    if(addsigs){
    # Append Signature.1 and Signature.5 if they have been removed, as we know all samples have them
    if(! "SBS1" %in% colnames(cancer_signaturesPresent)){
      signs <- cbind(cancer_signaturesPresent, "SBS1"=cancer_signaturesPresent_bk[,"SBS1"])
      fit_res <- fit_to_signatures(as.matrix(mut_mat[,i],1), signs)
      temp2Cos <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
      if(temp2Cos > c){
        signaturesNotPresent <- signaturesNotPresent[,colnames(signaturesNotPresent)!="SBS1",drop=FALSE]
        cancer_signaturesPresent <- signs
        incr <- temp2Cos - c
        c <- temp2Cos
        print(paste0("-- case ", colnames(mut_mat[,i,drop=FALSE]), "... added increase cosine SBS1 ",incr))
      }
    }
    if(! "SBS5" %in% colnames(cancer_signaturesPresent)){
      signs <- cbind(cancer_signaturesPresent, "SBS5"=cancer_signaturesPresent_bk[,"SBS5"])
      fit_res <- fit_to_signatures(matrix(mut_mat[,i], ncol=1), signs)
      temp2Cos <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
      if(temp2Cos > c){
        signaturesNotPresent <- signaturesNotPresent[,colnames(signaturesNotPresent)!="SBS5",drop=FALSE]
        cancer_signaturesPresent <- signs
        incr <- temp2Cos - c
        c <- temp2Cos
        print(paste0("-- case ", colnames(mut_mat[,i,drop=FALSE]), "... added increase cosine SBS5 ",incr))
      }
    }
    
    # Append Signature.9 if it has been removed, as we know all mutated cases have it
    if(!is.null(mutated_cases) & ! "SBS9" %in% colnames(cancer_signaturesPresent) & colnames(mut_mat)[i] %in% mutated_cases){
      signs <- cbind(cancer_signaturesPresent, "SBS9"=cancer_signaturesPresent_bk[,"SBS9"])
      fit_res <- fit_to_signatures(as.matrix(mut_mat[,i],1), signs)
      temp2Cos <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
      if(temp2Cos > c){
        signaturesNotPresent <- signaturesNotPresent[,colnames(signaturesNotPresent)!="SBS9",drop=FALSE]
        cancer_signaturesPresent <- signs
        incr <- temp2Cos - c
        c <- temp2Cos
        print(paste0("-- case ", colnames(mut_mat[,i,drop=FALSE]), "... added increase cosine SBS9 MCLL ",incr))
      }
    }
    }
    
    # Add signatures not present in this sample/clone but identified in other CLL/RS increases threshold_add
    if(!is.null(threshold_add)){
      add <- TRUE
      while(add){
        cosIncrease <- NULL
        for(j in colnames(signaturesNotPresent)){
          signs <- cbind(cancer_signaturesPresent, signaturesNotPresent[, colnames(signaturesNotPresent) == j,drop=FALSE])
          colnames(signs) <- c(colnames(cancer_signaturesPresent), j)
          fit_res <- fit_to_signatures(matrix(mut_mat[,i], ncol=1), as.matrix(signs))
          temp2Cos <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
          cosIncrease <- c(cosIncrease, temp2Cos - c)
        }
        names(cosIncrease) <- colnames(signaturesNotPresent)

        if(max(cosIncrease) > threshold_add){
          clnames <- colnames(cancer_signaturesPresent)
          cancer_signaturesPresent <- cbind(cancer_signaturesPresent, signaturesNotPresent[, colnames(signaturesNotPresent) == names(cosIncrease)[which.max(cosIncrease)]])
          colnames(cancer_signaturesPresent) <- c(clnames, names(cosIncrease)[which.max(cosIncrease)])
          signaturesNotPresent <- signaturesNotPresent[,colnames(signaturesNotPresent) != names(cosIncrease)[which.max(cosIncrease)]]
          fit_res <- fit_to_signatures(matrix(mut_mat[,i], ncol=1), as.matrix(cancer_signaturesPresent))
          c <- cosineSimilarity(mut_mat[,i], fit_res$reconstructed[,1])
          print(paste0("-- case ", colnames(mut_mat[,i,drop=FALSE]), "... added PCAWG signature (",
                       names(cosIncrease)[which.max(cosIncrease)],") increase cosine > threshold_add (",threshold_add,")"))
        } else{  
          add <- FALSE
        }
      }
    }
    
    # last fit
    fit_res <- fit_to_signatures(as.matrix(mut_mat[,i], ncol=1), as.matrix(cancer_signaturesPresent))
    finalCos <- cosineSimilarity(mut_mat[,i],fit_res$reconstructed[,1])
    contribution <- fit_res$contribution
    
    contributionMatrix[i, c(1:2, match(rownames(contribution), colnames(contributionMatrix)))] <- c(colnames(mut_mat)[i], round(finalCos,4), round(contribution))
  }
  
  contributionMatrix[,2:ncol(contributionMatrix)] <- sapply(contributionMatrix[,2:ncol(contributionMatrix)], as.numeric)
  contributionMatrix
  
  contributionMatrix[is.na(contributionMatrix)] <- 0
  
  smallContributionMatrix <- contributionMatrix[, c(T, T, colSums(contributionMatrix[, 3:ncol(contributionMatrix)]) > 0)]
  smallContributionMatrix <- as.data.frame(smallContributionMatrix)
  return(smallContributionMatrix)
}

# Function based on maftools function get_kataegis (https://github.com/PoisonAlien/maftools)
get_kat <- function(muts,nmuts,win){
  # print(nrow(muts))
  # Get clustered (cc) and non-clustered (nn) SNVs
  cc <- data.frame()
  kat_summary <- data.frame()
  # for each sample...
  for (sample in unique(muts$SAMPLE)){
    # print(sample)
    aux <- muts[muts$SAMPLE==sample,]
    # for each chromosome...
    for (chr in unique(aux$CHROM)){
      aux2 <- aux[aux$CHROM==chr,]
      # at least nmuts
      if (nrow(aux2)>=nmuts) {
        start <- 1
        end <- nmuts
        # print(paste0("start:end ",start,":",end))
        candidates <- aux2[start:end]
        # print(candidates$POSITION)
        while(end <= nrow(aux2)){
          if(mean(diff(candidates[, POSITION], na.rm = TRUE), na.rm = TRUE) > win){
            # Not complying kataegis definition, check next positions
            start <- start+1
            end <- end+1
            # print(paste0("NOPE start:end ",start,":",end))
            candidates <- aux2[start:end]
            # print(candidates$POSITION)
          }else{
            # While finding more kataegis mutations...
            while( mean(diff(candidates[, POSITION], na.rm = TRUE), na.rm = TRUE) <= win & end <= nrow(aux2) ){
              end <- end+1
              # print(paste0("YES start:end ",start,":",end))
              candidates <- aux2[start:end]
              # print(candidates$POSITION)
            }
            # Kataegis
            cc <- rbind(cc,aux2[start:(end-1)])
            
            x = aux2[(start):c(end-1)]  # start_idx not incremented after kat detected
            ycp = data.table::data.table(Chromosome = unique(x[,CHROM]),
                                         Start_Position = x[,min(POSITION)],
                                         End_Position = x[,max(POSITION)],
                                         nMuts = nrow(x),
                                         Avg_intermutation_dist = mean(x[,diff(POSITION)]))
            ycp[,Size := End_Position - Start_Position]
            
            x$con.class <- paste(x$context,x$mutType,sep="_")
            ycp = cbind(ycp,
                        data.table::dcast(data = x[,.N,.(con.class, SAMPLE)],
                                          SAMPLE ~ con.class, value.var = 'N'))
            
            kat_summary = data.table::rbindlist(l = list(kat_summary, ycp), fill = TRUE, use.names = TRUE)
            
            start <- end
            end <- start + nmuts -1
            candidates <- aux2[start:end]
          }
        }
      }
    }
  }
  nrow(cc)
  nrow(unique(cc))
 
  kat_summary <- as.data.frame(kat_summary)
  caid_changes <- c("ACT_C>T","ACC_C>T","GCT_C>T","GCC_C>T","ACT_C>G","ACC_C>G","GCT_C>G","GCC_C>G")
  
  kat_summary$cAID <- rowSums(kat_summary[,colnames(kat_summary) %in% caid_changes],na.rm = TRUE)
  kat_summary$other <- rowSums(kat_summary[,grepl(">",colnames(kat_summary)) & !colnames(kat_summary) %in% caid_changes],na.rm = TRUE)
  
  return(list(cc,kat_summary))
}

cosineSimilarity <- function(x,y){
  cos <- x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
  return(cos)
}

plot_extracted_sigs <- function(sigs,mut_types,title){
  signatures2 <- sigs
  signatures2$type <- mut_types
  signatures2 <- melt(signatures2)
  head(signatures2)
  signatures2$type2 <- unlist(lapply(signatures2$type, function(x) substr(x,3,5)))
    
  p<-ggplot(data=signatures2, aes(x=type, y=value, fill=type2)) +
        geom_bar(stat="identity", color="black",width=0.8) + facet_grid(variable~type2,scales="free")+ theme_classic() +
        theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 0, size = 4)) +
        scale_fill_manual(values=c("#14BAEC", "#010202", "#E42B28", "#989999", "#9FC866", "#EAC6C5")) +
        ggtitle(title)+ylab("Probability")+xlab("Context")+
        theme(strip.background = element_blank())

  return(p)
}

plotSignature2 <- function(sig,title,mut_types){
  sig <- as.data.frame(sig)
  colnames(sig) <- "count"
  sig_toplot <- data.frame(count=sig,context=mut_types,type=substr(mut_types,3,5))
  p<-ggplot(data=sig_toplot, aes(x=context, y=count, fill=type)) +
      geom_bar(stat="identity", color="black") + facet_grid(.~type,scales="free")+ theme_classic() +
      theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 0, size = 4)) +
      scale_fill_manual(values=c("#14BAEC", "#010202", "#E42B28", "#989999", "#9FC866", "#EAC6C5")) +
      ggtitle(title)+ylab("Probability")+xlab("Context")+
      theme(strip.background = element_blank()) + 
      geom_hline(aes(yintercept = Inf), color = "white", size=2) + # white space
      geom_hline(aes(yintercept = Inf, color = type),size = 8)+
      scale_colour_manual(values=c("#14BAEC", "#010202", "#E42B28", "#989999", "#9FC866", "#EAC6C5"))
  return(p)
}

# Function adapted from https://rdrr.io/github/sjdlabgroup/MutSigTools/src/R/mergeSignature.R
mergeSignature2 <-function(sigmatrix,sig,weights,mut_types,title="mergesig") {
  rownames(sigmatrix) <- mut_types
  mergesig <- rowSums(t(t(sigmatrix[,sig])*weights))
  mergesig_toplot <- data.frame(count=mergesig,context=rownames(sigmatrix),type=substr(rownames(sigmatrix),3,5))
  p<-ggplot(data=mergesig_toplot, aes(x=context, y=count, fill=type)) +
        geom_bar(stat="identity", color="black") + facet_grid(.~type,scales="free")+ theme_bw() +
        theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 0, size = 4)) +
        scale_fill_manual(values=c("#14BAEC", "#010202", "#E42B28", "#989999", "#9FC866", "#EAC6C5")) +
        ggtitle(title)
  return(list(mergesig,p))
}
