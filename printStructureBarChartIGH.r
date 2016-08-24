printStructureBarChartIGH <- function(gene.name,...) {

# IGH - 
# lin.neg - isoforms on Chr14 r = IRanges(106053231,107170413)
#           pbgene.x = "PB.2080"
# total   - isoforms on chr14 r = IRanges(106053224,107199471)
#           pbgene.x = "PB.2666"
# lin.pos - isoforms on chr14 r = IRanges(106053230,107179279)
#           pbgene.x = "PB.1184"
#
   range= (106053224,107199471)
   r = IRanges(range)
   wh = GRanges(seqnames=Rle("chr14"),r)
   gr.txdb <- crunch(txdb,which=wh,type="all")
   colnames(values(gr.txdb))[4] <- "model"
   grl <- split(gr.txdb, gr.txdb$tx_id)
   
   lin.neg = paste(paste("^",PB.2080,sep=""),"$",sep="")
   #-----------------------------------
   # lin-
   #-----------------------------------
   pb.53A.neg.ga <- readGAlignmentsFromBam(pb.53A.neg.bamfile,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   neg.len = length(pb.53A.neg.ga)
   if (neg.len != 0) {

      pb.ICE.53A.neg.data <- pb.ICE.53A.lin.neg[grep("PB.2080",pb.ICE.53A.lin.neg[,c("pbgene.x")]),]
      pb.53A.neg.ga <- updateAlignmentData(pb.ICE.53A.neg.data, grepname, pb.53A.neg.ga)
      neg.len = length(pb.53A.neg.ga)

      name = paste(pb.ICE.53A.neg.data[,"id"],paste(pb.ICE.53A.neg.data[,"genomeLoc"],
        paste(pb.ICE.53A.neg.data[,"cluster"],paste(pb.ICE.53A.neg.data[,"fp"],paste(pb.ICE.53A.neg.data[,"pblength"],
        sep="/"),sep="/"),sep="/"),sep="|"),sep="|")

#
# weird bug -- there is a repeat of one of the alignments.
#              remove the last one....? -- neg is 31 alignments long
      pb.53A.neg.ga = pb.53A.neg.ga[1:31]
      neg.id = as.vector(as.character(unlist(pb.ICE.53A.neg.data[,c("id")])))

      t <- data.frame(id    = neg.id,
              TPM    = rep(0,length(neg.id)),
              logTPM = rep(0,length(neg.id)),
              length = rep(0,length(neg.id)),
              color  = as.character("blue"))
      mcols(pb.53A.neg.ga) <- t
      ga.len = length(name)
      for (k in 1:ga.len) {
         m = match(name,names(pb.53A.neg.ga[k]))
         r = sum(!is.na(m))
         if (r > 0) {
            tpm = pb.ICE.53A.neg.data[!is.na(m),"TPM"]
            id = as.character(unlist(pb.ICE.53A.neg.data[!is.na(m),"id"]))
            length = pb.ICE.53A.neg.data[!is.na(m),"pblength"]
            mcols(pb.53A.neg.ga)$TPM[k]    = tpm
            mcols(pb.53A.neg.ga)$logTPM[k] = log10(tpm)
            mcols(pb.53A.neg.ga)$id[k]     = as.character(unlist(id))
            mcols(pb.53A.neg.ga)$length[k] = length
         }
      }
   } else {
      neg.len = 0
   }
   
   #-----------------------------------
   # lin+
   #-----------------------------------
   pb.53A.pos.ga <- readGAlignmentsFromBam(pb.53A.pos.bamfile,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   pos.len = length(pb.53A.pos.ga)
   if (pos.len != 0) {

      pb.ICE.53A.pos.data <- getIsoformData(pb.ICE.53A.lin.pos, grepname, pb.53A.pos.ga)
      pb.53A.pos.ga <- updateAlignmentData(pb.ICE.53A.pos.data, grepname, pb.53A.pos.ga)
      pos.len = length(pb.53A.pos.ga)
      name = paste(pb.ICE.53A.pos.data[,"id"],paste(pb.ICE.53A.pos.data[,"genomeLoc"],
        paste(pb.ICE.53A.pos.data[,"cluster"],paste(pb.ICE.53A.pos.data[,"fp"],paste(pb.ICE.53A.pos.data[,"pblength"],
        sep="/"),sep="/"),sep="/"),sep="|"),sep="|")

#
# weird bug -- there is a repeat of one of the alignments.
#              remove the last one....? -- pos is 61 alignments long
      pb.53A.pos.ga = pb.53A.pos.ga[1:61]
      pos.id = as.vector(as.character(unlist(pb.ICE.53A.pos.data[,c("id")])))

      t <- data.frame(id    = pos.id,
              TPM    = rep(0,length(pos.id)),
              logTPM = rep(0,length(pos.id)),
              length = rep(0,length(pos.id)),
              color  = as.character("blue"))
      mcols(pb.53A.pos.ga) <- t
      ga.len = length(name)
      for (k in 1:ga.len) {
         m = match(name,names(pb.53A.pos.ga[k]))
         r = sum(!is.na(m))
         if (r > 0) {
            tpm = pb.ICE.53A.pos.data[!is.na(m),"TPM"]
            id = as.character(unlist(pb.ICE.53A.pos.data[!is.na(m),"id"]))
            length = pb.ICE.53A.pos.data[!is.na(m),"pblength"]
            mcols(pb.53A.pos.ga)$TPM[k]    = tpm
            mcols(pb.53A.pos.ga)$logTPM[k] = log10(tpm)
            mcols(pb.53A.pos.ga)$id[k]     = as.character(unlist(id))
            mcols(pb.53A.pos.ga)$length[k] = length
         }
      }
   } else {
      pos.len = 0
   }
   
   #---------------------------------------------------------
   # create a bar chart pos and neg
   #  because of certain events -- it can be that both
   # pos.len and neg.len are zero -- in this case we exit without
   #   printing anything.
   #----------------------------------------------------------
   if (neg.len != 0) {
      neg <- data.frame ( 
                   id     = mcols(pb.53A.neg.ga)$id,
                   pop    = c(rep("Lin-",length(pb.53A.neg.ga))),
                   gene   = c(rep(gene.name,length(pb.53A.neg.ga))),
                   logTPM = mcols(pb.53A.neg.ga)$logTPM,
                   TPM    = mcols(pb.53A.neg.ga)$TPM,
                   col    = "darkblue")
      n.i = order(neg$TPM,decreasing=T)
      nr.i = order(neg$TPM,decreasing=F)
      #----------------------------------------------------------------------
      # this held me up for a long time!  you have to transform the variable 
      # back to a character vector in order to use it
      #----------------------------------------------------------------------
      neg[,"col"]  = sapply(neg[,"col"], as.character)
   }
   if (pos.len != 0) {
     pos <- data.frame ( 
                   id     = mcols(pb.53A.pos.ga)$id,
                   pop    = c(rep("Lin+",length(pb.53A.pos.ga))),
                   gene   = c(rep(gene.name,length(pb.53A.pos.ga))),
                   logTPM = mcols(pb.53A.pos.ga)$logTPM,
                   TPM    = mcols(pb.53A.pos.ga)$TPM,
                   col    = "brown")
     p.i = order(pos$TPM,decreasing=T)
     pr.i = order(pos$TPM,decreasing=F)
     pos[,"col"]  = sapply(pos[,"col"], as.character)
   }
   #----------------------------------------------------------------------
   if (pos.len !=0 & neg.len != 0) {
      both <- rbind(neg[nr.i,],pos[pr.i,])
   } else if (pos.len !=0) {
     both <- pos[pr.i,]
   } else if (neg.len !=0) {
     both <- neg[nr.i,]
   }
   both[,"col"] = sapply(both[,"col"],as.character)
   #----------------------------------------------------------
   # in adobe we will put the barcharts together with the structure diagrams
   # first the abundance information.
   # 2014-12-3 - there are approximately 80 isoforms in the case of MPO - this is the highest
   #             number for this experiment. 
   #             1 structural canonical isoform for MPO 
   #             6 for lin+
   #             71 for lin-
   # 
   #----------------------------------------------------------
   filename = paste(paste("../jpg/",gene.name,sep=""),".both.barplot.pdf",sep="")
   max_isoforms = 80
   bar.width = 1/80
   pdf(filename, width = 5,height = 10)
     par(las=2,cex=0.7,mar=c(5,10,5,5))
     barplot((both$logTPM),col=as.character(both$col),horiz=T,xlim=c(0,5),
             names.arg=c(as.character(both$id)),xlab="log10(TPM)",main=as.character(gene.name))
   dev.off()

   #----------------------------------------------------------
   # print out the structure information 
   # height is normalized by the protein with the most isoforms
   # canonical = 3
   # pos.len = 6
   # neg.len = 71
   # total = 80
   #----------------------------------------------------------
   if ((pos.len != 0) & (neg.len != 0)) {
      pr.i = reverse(p.i)
      nr.i = reverse(n.i)
#----------------------------------------------------------
# to print out labels as well -- for debugging. -- unfortunately -- order is not preserved unless group.selfish is kept
#  and this causes labels to be printed -- requiring editing later
      p1.r <-autoplot((pb.53A.pos.ga[pr.i]),geom="alignment",group.selfish=T,fill="brown",color="brown")
      n1.r <-autoplot((pb.53A.neg.ga[nr.i]),geom="alignment",group.selfish=T,fill="darkblue",color="darkblue")
#----------------------------------------------------------
      directory = paste("../jpg/",gene.name,sep="")
      filename = paste(directory,".pos.neg.gene.model.isoform.pdf",sep="")
      neg.span = max(end(pb.53A.neg.ga)) - min(start(pb.53A.neg.ga))
      pos.span = max(end(pb.53A.pos.ga)) - min(start(pb.53A.pos.ga))
      if (neg.span > pos.span) {
          span=neg.span
      } else {
          span=pos.span
      }
      can.length = length(grl)
      pos.length = length(pb.53A.pos.ga)
      neg.length = length(pb.53A.neg.ga)
      tot.length = pos.length + neg.length + can.length
      ph = pos.length/tot.length
      nh = neg.length/tot.length
      ch = can.length/tot.length
      p0 <- autoplot(grl, gap.geom="chevron",aes(type=model))
      main.title = paste(gene.name,as.character(span),sep=".")
      
      tracks(p1.r,n1.r,
             heights=c(ph,nh),
             main=main.title)
      ggsave(filename,width=8,height=10)
   } else if (pos.len != 0) {
      directory = paste("../jpg/",gene.name,sep="")
      filename = paste(directory,".pos.neg.gene.model.isoform.pdf",sep="")
      pos.span = max(end(pb.53A.pos.ga)) - min(start(pb.53A.pos.ga))
      pr.i = reverse(p.i)
      p1.r <-autoplot((pb.53A.pos.ga[pr.i]),geom="alignment",group.selfish=T,fill="brown",color="brown")
      can.length = length(grl)
      pos.length = length(pb.53A.pos.ga)
      tot.length = pos.length + can.length
      ph = pos.length/tot.length
      ch = can.length/tot.length
      p0 <- autoplot(grl, gap.geom="chevron",aes(type=model))

      main.title = paste(gene.name,as.character(pos.span),sep=".")         
      tracks(p0,p1.r,
             heights=c(ch,ph),
             main=main.title)
      ggsave(filename,width=8,height=10)
   } else if (neg.len !=0) {
      directory = paste("../jpg/",gene.name,sep="")
      filename = paste(directory,".pos.neg.gene.model.isoform.pdf",sep="")
      neg.span = max(end(pb.53A.neg.ga)) - min(start(pb.53A.neg.ga))
      nr.i = reverse(n.i)
      n1.r <-autoplot((pb.53A.neg.ga[nr.i]),geom="alignment",group.selfish=T,fill="darkblue",color="darkblue")
      can.length = length(grl)
      neg.length = length(pb.53A.neg.ga)
      tot.length = neg.length + can.length
      nh = neg.length/tot.length
      ch = can.length/tot.length
      p0 <- autoplot(grl, gap.geom="chevron",aes(type=model))
      
      main.title = paste(gene.name,as.character(neg.span),sep=".")
      tracks(p0,n1.r,
             heights=c(ch,nh),
             main=main.title)

      ggsave(filename,width=8,height=10)
   }
}
