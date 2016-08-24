printStructureBarChartNoMetaData <- function(gene.name,...) {

   wh = genesymbol[gene.name]
   gr.txdb <- crunch(txdb,which=wh)
   colnames(values(gr.txdb))[4] <- "model"
   grl <- split(gr.txdb, gr.txdb$tx_id)
   p0 <- autoplot(grl, aes(type=model))
   grepname = paste(paste("^",gene.name,sep=""),"$",sep="")
   #-----------------------------------
   # lin-
   #-----------------------------------
   pb.53A.neg.ga <- readGAlignmentsFromBam(pb.53A.neg.bamfile,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   neg.len = length(pb.53A.neg.ga)
   if (neg.len != 0) {
     pb.ICE.53A.neg.data <- pb.ICE.53A.lin.neg[grep(grepname,pb.ICE.53A.lin.neg[,"gene"]),]

     pbgenes = pb.ICE.53A.neg.data[,"pbgene.x"]
     u.pbgenes = unique(pbgenes)
     majority = 0
     major.pbgene = ""
     for (i in 1:length(u.pbgenes)) {
       m = match(pb.ICE.53A.neg.data[,"pbgene.x"],u.pbgenes[i])
       if (sum(!is.na(m)) > majority) {
          majority = sum(!is.na(m))
          major.pbgene = u.pbgenes[i]
       }
     }
     m = match(pb.ICE.53A.neg.data[,"pbgene.x"],major.pbgene)
     pb.ICE.53A.neg.data = pb.ICE.53A.neg.data[!is.na(m),]

     name = paste(pb.ICE.53A.neg.data[,"id"],paste(pb.ICE.53A.neg.data[,"genomeLoc"],
          paste(pb.ICE.53A.neg.data[,"cluster"],paste(pb.ICE.53A.neg.data[,"fp"],paste(pb.ICE.53A.neg.data[,"pblength"],
          sep="/"),sep="/"),sep="/"),sep="|"),sep="|")

     m = match(names(pb.53A.neg.ga),name)
     # there is sometimes a mixup regarding matches -- clean up 
     if (sum(!is.na(m)) != length(pb.53A.neg.ga)) {
        pb.53A.neg.gar = pb.53A.neg.ga[!is.na(m),]
        pb.53A.neg.ga = pb.53A.neg.gar
     }
     neg.id = as.vector(as.character(unlist(pb.ICE.53A.neg.data[,c("id")])))
     t <- data.frame(id    = neg.id,
                  mate_id= neg.id,
                  TPM    = rep(0,length(name)),
                  logTPM = rep(0,length(name)),
                  length = rep(0,length(name)),
                  color  = as.character("blue"),
                  median = 0.0,
                  shared = rep(FALSE,length(name)))
     mcols(pb.53A.neg.ga) <- t
     ga.len = length(name)
     for (j in 1:ga.len) {
        m = match(name,names(pb.53A.neg.ga[j]))
        r = sum(!is.na(m))
        if (r > 0) {
          tpm = pb.ICE.53A.neg.data[!is.na(m),"TPM"]
          id = as.character(unlist(pb.ICE.53A.neg.data[!is.na(m),"id"]))
          length = pb.ICE.53A.neg.data[!is.na(m),"pblength"]
          mcols(pb.53A.neg.ga)$TPM[j]    = tpm
          mcols(pb.53A.neg.ga)$logTPM[j] = log10(tpm)
          mcols(pb.53A.neg.ga)$id[j]     = as.character(unlist(id))
          mcols(pb.53A.neg.ga)$length[j] = length
          n = match(mcols(pb.53A.neg.ga)$id[j],shared.df[,"mate_id"])
          s = sum(!is.na(n))
          if (s >0) {
             mcols(pb.53A.neg.ga)$shared[j] = TRUE
             mcols(pb.53A.neg.ga)$color[j]  = as.character("yellow")
          }
        }
     }
   }
   #-----------------------------------
   # lin+
   #-----------------------------------
   pb.53A.pos.ga <- readGAlignmentsFromBam(pb.53A.pos.bamfile,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   pos.len = length(pb.53A.pos.ga)
   if (pos.len != 0) {
     pb.ICE.53A.pos.data <- pb.ICE.53A.lin.pos[grep(grepname,pb.ICE.53A.lin.pos[,"gene"]),]
     pbgenes = pb.ICE.53A.pos.data[,"pbgene.x"]
     u.pbgenes = unique(pbgenes)
     majority = 0
     major.pbgene = ""
     for (i in 1:length(u.pbgenes)) {
       m = match(pb.ICE.53A.pos.data[,"pbgene.x"],u.pbgenes[i])
       if (sum(!is.na(m)) > majority) {
          majority = sum(!is.na(m))
          major.pbgene = u.pbgenes[i]
       }
     }
     m = match(pb.ICE.53A.pos.data[,"pbgene.x"],major.pbgene)
     pb.ICE.53A.pos.data = pb.ICE.53A.pos.data[!is.na(m),]


     name = paste(pb.ICE.53A.pos.data[,"id"],paste(pb.ICE.53A.pos.data[,"genomeLoc"],
            paste(pb.ICE.53A.pos.data[,"cluster"],paste(pb.ICE.53A.pos.data[,"fp"],paste(pb.ICE.53A.pos.data[,"pblength"],
            sep="/"),sep="/"),sep="/"),sep="|"),sep="|")

     m = match(names(pb.53A.pos.ga),name)
     # there is sometimes a mixup regarding matches -- clean up 
     if (sum(!is.na(m)) != length(pb.53A.pos.ga)) {
        pb.53A.pos.gar = pb.53A.pos.ga[!is.na(m),]
        pb.53A.pos.ga = pb.53A.pos.gar
     }
     pos.id = as.vector(as.character(unlist(pb.ICE.53A.pos.data[,c("id")])))
     t <- data.frame(id    = pos.id,
                  mate_id= pos.id,
                  TPM    = rep(0,length(name)),
                  logTPM = rep(0,length(name)),
                  length = rep(0,length(name)),
                  color  = as.character("brown"),
                  median = 0.0,
                  shared = rep(FALSE,length(name)))
     mcols(pb.53A.pos.ga) <- t
     ga.len = length(name)
     for (j in 1:ga.len) {
        m = match(name,names(pb.53A.pos.ga[j]))
        r = sum(!is.na(m))
        if (r > 0) {
          tpm = pb.ICE.53A.pos.data[!is.na(m),"TPM"]
          id = as.character(unlist(pb.ICE.53A.pos.data[!is.na(m),"id"]))
          length = pb.ICE.53A.pos.data[!is.na(m),"pblength"]
          mcols(pb.53A.pos.ga)$TPM[j]    = tpm
          mcols(pb.53A.pos.ga)$logTPM[j] = log10(tpm)
          mcols(pb.53A.pos.ga)$id[j]     = as.character(unlist(id))
          mcols(pb.53A.pos.ga)$length[j] = length
          n = match(mcols(pb.53A.pos.ga)$id[j],shared.df[,"mate_id"])
          s = sum(!is.na(n))
          if (s >0) {
             mcols(pb.53A.pos.ga)$shared[j] = TRUE
             mcols(pb.53A.pos.ga)$color[j]  = as.character("yellow")
          }
        }
     }
   }
   #-------------------------------------
   # find shared ids
   #-------------------------------------
   m.shared.ids = match(shared.lin.neg.lin.pos.chained.pb.ids[,4],pos.id)   
   shared.id = shared.lin.neg.lin.pos.chained.pb.ids[!is.na(m.shared.ids),]

   len = dim(shared.id)[1]
   # update the pos information
   if (pos.len != 0) {
     for (i in 1:len) {
     # get the mate info
         t = pb.ICE.53A.neg.data[,"id"] == shared.id[i,"negative"]
       tpm = pb.ICE.53A.neg.data[t,"TPM"]
       logtpm = log10(tpm)
       # add it to the pos info
       t = pb.ICE.53A.pos.data[,"id"] == shared.id[i,"positive"]
       name = paste(pb.ICE.53A.pos.data[t,"id"],paste(pb.ICE.53A.pos.data[t,"genomeLoc"],
          paste(pb.ICE.53A.pos.data[t,"cluster"],paste(pb.ICE.53A.pos.data[t,"fp"],
          paste(pb.ICE.53A.pos.data[t,"pblength"],
          sep="/"),sep="/"),sep="/"),sep="|"),sep="|")
       m = match(names(pb.53A.pos.ga),name)
       mcols(pb.53A.pos.ga)$mate_TPM[!is.na(m)]   = tpm
       mcols(pb.53A.pos.ga)$mate_logTPM[!is.na(m)] = logtpm
       mcols(pb.53A.pos.ga)$shared[!is.na(m)] =  TRUE
     }
   }

   # update the negative information
   if (length(pb.53A.pos.ga) != 0) {
     for (i in 1:len) {
       # get the mate info  
       t = pb.ICE.53A.pos.data[,"id"] == shared.id[i,"positive"]
       tpm = pb.ICE.53A.pos.data[t,"TPM"]
       logtpm = log10(tpm)
       # add it to the neg info
       t = pb.ICE.53A.neg.data[,"id"] == shared.id[i,"negative"]
       name = paste(pb.ICE.53A.neg.data[t,"id"],paste(pb.ICE.53A.neg.data[t,"genomeLoc"],
          paste(pb.ICE.53A.neg.data[t,"cluster"],paste(pb.ICE.53A.neg.data[t,"fp"],
          paste(pb.ICE.53A.neg.data[t,"pblength"],
          sep="/"),sep="/"),sep="/"),sep="|"),sep="|")
       m = match(names(pb.53A.neg.ga),name)
       mcols(pb.53A.neg.ga)$mate_TPM[!is.na(m)]   = tpm
       mcols(pb.53A.neg.ga)$mate_logTPM[!is.na(m)] = logtpm
       mcols(pb.53A.neg.ga)$shared[!is.na(m)] =  TRUE
     }
   }

   #---------------------------------------------------------
   # create a bar chart pos and neg
   #   around the median of the values for this gene
   #   we have here the TPMs -- and associate them with the 
   #   relative to total?  or just to their associated median
   #----------------------------------------------------------
   if (length(pb.53A.pos.ga) != 0) {
      gene.name = pb.ICE.53A.pos.data[1,"gene"]
   }
   else {
      gene.name = pb.ICE.53A.neg.data[1,"gene"]
   }

   if (length(pb.53A.neg.ga) != 0) {
      neg <- data.frame ( 
                      id     = mcols(pb.53A.neg.ga)$id,
                      pop    = c(rep("Lin-",length(pb.53A.neg.ga))),
                      gene   = c(rep(gene.name,length(pb.53A.neg.ga))),
                      TPM    = mcols(pb.53A.neg.ga)$logTPM,
                      shared = mcols(pb.53A.neg.ga)$shared,
                      col    = "darkblue",
                      median = 0.0)
     n.i = order(neg$TPM,decreasing=T)
     #----------------------------------------------------------------------
     # this held me up for a long time!  you have to transform the variable 
     # back to a character vector in order to use it
     #----------------------------------------------------------------------
     neg[,"col"]  = sapply(neg[,"col"], as.character)
   }
   
   if (length(pb.53A.pos.ga) != 0) {
      pos <- data.frame ( 
                      id     = mcols(pb.53A.pos.ga)$id,
                      pop    = c(rep("Lin+",length(pb.53A.pos.ga))),
                      gene   = c(rep(gene.name,length(pb.53A.pos.ga))),
                      TPM    = mcols(pb.53A.pos.ga)$logTPM,
                      shared = mcols(pb.53A.pos.ga)$shared,
                      col    = "brown",
                      median = 0.0)
      p.i = order(pos$TPM,decreasing=T)
      pos[,"col"]  = sapply(pos[,"col"], as.character)
   }
   #----------------------------------------------------------------------
   #
   #  only color green that set of isoforms which is greater --
   #  that is in which population more is occuring
   #   be that lin- or lin+
   #
   #----------------------------------------------------------------------
   pos.len = ifelse(length(pb.53A.pos.ga) != 0,dim(pos)[1],0)
   neg.len = ifelse(length(pb.53A.neg.ga) != 0,dim(neg)[1],0)
   if (pos.len > neg.len) {
      for (i in 1:pos.len) {
          if (pos$shared[i] == TRUE) {
             pos$col[i] = "darkgreen"
          }
          else {
             pos$col[i] = "brown"
          }
      }
   }
   if (neg.len > pos.len) {
      for (i in 1:neg.len) {
          if (neg$shared[i] == TRUE) {
             neg$col[i] = "darkgreen"
          }
      }
   }
   
   #----------------------------------------------------------------------
   # if else requires the syntax below -- else on line with return "}"
   #----------------------------------------------------------------------
   if (pos.len !=0 & neg.len != 0) {
     both <- rbind(neg[n.i,],pos[p.i,])
   } else if (pos.len !=0) {
     both <- pos[p.i,]
   } else  {
     both <- neg[n.i,]
   }
   both[,"col"] = sapply(both[,"col"],as.character)

   #----------------------------------------------------------
   #
   # in adobe we will put the barcharts together with the structure diagrams
   #
   #----------------------------------------------------------
   filename = paste(paste("../jpg/",gene.name,sep=""),".both.barplot.pdf",sep="")
   pdf(filename, width = 5,height = 10)
      barplot((both$TPM),col=as.character(both$col),horiz=T,xlim=c(0,5))
   dev.off()

   #----------------------------------------------------------
   # reread alignment to have no meta data - screws up join
   #----------------------------------------------------------   
   pb.53A.neg.ga <- readGAlignmentsFromBam(pb.53A.neg.bamfile,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   pb.53A.pos.ga <- readGAlignmentsFromBam(pb.53A.pos.bamfile,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   pos.len = length(pb.53A.pos.ga)
   neg.len = length(pb.53A.neg.ga)
   if (neg.len != 0) {
     name = paste(pb.ICE.53A.neg.data[,"id"],paste(pb.ICE.53A.neg.data[,"genomeLoc"],
          paste(pb.ICE.53A.neg.data[,"cluster"],paste(pb.ICE.53A.neg.data[,"fp"],paste(pb.ICE.53A.neg.data[,"pblength"],
          sep="/"),sep="/"),sep="/"),sep="|"),sep="|")

     m = match(names(pb.53A.neg.ga),name)
     # there is sometimes a mixup regarding matches -- clean up
     if (sum(!is.na(m)) != length(pb.53A.neg.ga)) {
        pb.53A.neg.gar = pb.53A.neg.ga[!is.na(m),]
        pb.53A.neg.ga = pb.53A.neg.gar
     }
     n.i = order(pb.ICE.53A.neg.data[,"TPM"],decreasing=T)
   }
   
   if (pos.len != 0) {
     name = paste(pb.ICE.53A.pos.data[,"id"],paste(pb.ICE.53A.pos.data[,"genomeLoc"],
          paste(pb.ICE.53A.pos.data[,"cluster"],paste(pb.ICE.53A.pos.data[,"fp"],paste(pb.ICE.53A.pos.data[,"pblength"],
          sep="/"),sep="/"),sep="/"),sep="|"),sep="|")

     m = match(names(pb.53A.pos.ga),name)
     # there is sometimes a mixup regarding matches -- clean up
     if (sum(!is.na(m)) != length(pb.53A.pos.ga)) {
        pb.53A.pos.gar = pb.53A.pos.ga[!is.na(m),]
        pb.53A.pos.ga = pb.53A.pos.gar
     }
     p.i = order(pb.ICE.53A.pos.data[,"TPM"],decreasing=T)
   }
   #----------------------------------------------------------
   # join and reorder based upon level of expression
   #----------------------------------------------------------
   if ((pos.len != 0) & (neg.len != 0)) {
      pb.neg.pos = c(pb.53A.pos.ga[p.i],pb.53A.neg.ga[n.i])
   } else if (pos.len != 0) {
      pb.neg.pos = pb.53A.pos.ga[p.i]
   } else {
      pb.neg.pos = pb.53A.neg.ga[n.i]
   }

   p3 <- autoplot(pb.neg.pos,fill="brown",color="brown")

   directory = paste("../jpg/",gene.name,sep="")
   filename = paste(directory,".pos.neg.gene.model.isoform.pdf",sep="")
   pdf(file = filename,width= 5, height = 10,title=gene.name,fonts="Times",paper="letter",colormodel="cmyk",compress=FALSE)
       tracks(p0,p3,heights=c(.15,.85))
   dev.off()
}
