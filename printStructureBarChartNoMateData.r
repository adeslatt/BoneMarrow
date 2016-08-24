printStructureBarChartNoMateData <- function(start=start,end=end,chr=chr,grepgene.1=grepgene.1,gene.name=gene.name,bam.file.1=bam.file.1,isoform.data.1=isoform.data.1,gene.name.2=gene.name.2,grepgene.2=grepgene.2,bam.file.2=bam.file.2,isoform.data.2=isoform.data.2,col.1=col.1,col.2=col.2,filename.preface=filename.preface,ht=ht, wt=wt) {
   if (start > 0) {
      r = IRanges(start,end)
      wh = GRanges(seqnames=Rle(chr),r)
   } else {
      wh = genesymbol[genesymbol$symbol == gene.name]
   }
   wh=wh[1]
   #-----------------------------------
   # lin-
   #-----------------------------------
   pb.1.ga <- readGAlignments(bam.file.1,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   neg.len = length(pb.1.ga)
   if (neg.len != 0) {

      ICE.data.1.data <- getIsoformData(isoform.data.1, gene.name, gene.name.2,grepname, grepgene.1, pb.1.ga)
      pb.n = dim(ICE.data.1.data)[1]
      if (pb.n >0) {
         pb.1.ga <- updateAlignmentData(isoform.data=ICE.data.1.data, pb.ga=pb.1.ga, grepgene=grepgene.1)
         ICE.data.1.data <- updateICEData(ICE.data.1.data, pb.1.ga)
         ga.len = neg.len = length(pb.1.ga)
         n = strsplit(names(pb.1.ga),"\\|")
         m = unlist(n)
         ss = 1
         e = ga.len*3
         l = seq(from=ss,to=e,by=3)

         t <- data.frame(id    = c(m[l]),
              TPM    = rep(0,ga.len),
              logTPM = rep(0,ga.len),
              length = rep(0,ga.len),
              color  = as.character(col.1))
         mcols(pb.1.ga) <- t
         for (p in 1:ga.len) {
              tpm = ICE.data.1.data[ICE.data.1.data[,"id"] == m[l[p]],"TPM"]
              id = m[l[p]]
              length = ICE.data.1.data[ICE.data.1.data[,"id"] == m[l[p]],"pblength"]
              mcols(pb.1.ga)$TPM[p]    = tpm
              mcols(pb.1.ga)$logTPM[p] = log10(tpm)
              mcols(pb.1.ga)$id[p]     = as.character(unlist(id))
              mcols(pb.1.ga)$length[p] = length
         }
      } else {
        neg.len = 0
      }
   } else {
      neg.len = 0
   }
   
   #-----------------------------------
   # lin+
   #-----------------------------------
   pb.2.ga <- readGAlignments(bam.file.2,
                                 param = ScanBamParam(which = wh),
                                 use.names = TRUE)
   pos.len = length(pb.2.ga)
   if (pos.len != 0) {
      
      ICE.data.2.data <- getIsoformData(isoform.data.2, gene.name, gene.name.2,grepname, grepgene.2, pb.2.ga)
      pb.n = dim(ICE.data.2.data)[1]
      if (pb.n > 0) {
         pb.2.ga <- updateAlignmentData(isoform.data=ICE.data.2.data, pb.ga=pb.2.ga, grepgene=grepgene.2)
         ICE.data.2.data <- updateICEData(ICE.data.2.data, pb.2.ga)
         ga.len = pos.len = length(pb.2.ga)

         n = strsplit(names(pb.2.ga),"\\|")
         m = unlist(n)

         ss = 1
         e = ga.len*3
         l = seq(from=ss,to=e,by=3)

         t <- data.frame(id    = c(m[l]),
              TPM    = rep(0,ga.len),
              logTPM = rep(0,ga.len),
              length = rep(0,ga.len),
              color  = as.character(col.2))
         mcols(pb.2.ga) <- t
         for (p in 1:ga.len) {
              tpm = ICE.data.2.data[ICE.data.2.data[,"id"] == m[l[p]],"TPM"]
              id = m[l[p]]
              length = ICE.data.2.data[ICE.data.2.data[,"id"] == m[l[p]],"pblength"]
              mcols(pb.2.ga)$TPM[p]    = tpm
              mcols(pb.2.ga)$logTPM[p] = log10(tpm)
              mcols(pb.2.ga)$id[p]     = as.character(unlist(id))
              mcols(pb.2.ga)$length[p] = length
         }
      } else {
         pos.len = 0
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
                   id     = mcols(pb.1.ga)$id,
                   logTPM = mcols(pb.1.ga)$logTPM,
                   TPM    = mcols(pb.1.ga)$TPM,
                   col    = as.character(col.1))
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
                   id     = mcols(pb.2.ga)$id,
                   logTPM = mcols(pb.2.ga)$logTPM,
                   TPM    = mcols(pb.2.ga)$TPM,
                   col    = col.2)
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
   if (pos.len !=0 | neg.len !=0 ){
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
      setEPS()
      endfile  = paste(as.character(filename.preface),".both.barplot.eps",sep="")
      filename = paste(paste("../jpg/",gene.name,sep=""),endfile,sep="")
      max_isoforms = 80
      bar.width = 1/80
      bar.page.width = (wt/2)
      postscript(filename, height = ht+.5, width = bar.page.width)
        par(las=2,cex=0.7,mar=c(5,10,5,5))
        barplot((both$logTPM),col=as.character(both$col),horiz=T,xlim=c(0,5),
             names.arg=c(as.character(both$id)),xlab="log10(TPM)",main=as.character(gene.name))
      dev.off()
      rm(both)
   }
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
      #  print a shorter name
      #----------------------------------------------------------
      name = names(pb.2.ga[pr.i])
      nn = strsplit(name,"\\|")
      nnn = unlist(nn)
      e = length(nnn)
      ll = seq(from=1,to=e,by=3)
      pn = nnn[ll]
      pl = pb.2.ga[pr.i]
      names(pl) = pn
      #----------------------------------------------------------
      #  print a shorter name
      #----------------------------------------------------------
      name = names(pb.1.ga[nr.i])
      nn = strsplit(name,"\\|")
      nnn = unlist(nn)
      e = length(nnn)
      ll = seq(from=1,to=e,by=3)
      pn = nnn[ll]
      nl = pb.1.ga[nr.i]
      names(nl)<- pn

#----------------------------------------------------------
# to print out labels as well -- for debugging. -- unfortunately -- order is not preserved unless group.selfish is kept
#  and this causes labels to be printed -- requiring editing later
      p1.r <-autoplot(pl,geom="alignment",group.selfish=T,fill=col.2,color=col.2)
      n1.r <-autoplot(nl,geom="alignment",group.selfish=T,fill=col.1,color=col.1)
#----------------------------------------------------------
      directory = paste("/Users/adeslat/Desktop/newcode/code/Georgetown/boneMarrowRNASeq/jpg/",gene.name,sep="")
      endfile  = paste(filename.preface,".gene.model.isoform.eps",sep="")
      filename = paste(directory,endfile,sep="")
      neg.span = max(end(pb.1.ga)) - min(start(pb.1.ga))
      pos.span = max(end(pb.2.ga)) - min(start(pb.2.ga))
      if (neg.span > pos.span) {
          span=neg.span
      } else {
          span=pos.span
      }
      pos.length = length(pb.2.ga)
      neg.length = length(pb.1.ga)
      tot.length = pos.length + neg.length
      ph = pos.length/tot.length
      nh = neg.length/tot.length
      main.title = paste(gene.name,as.character(span),sep=".")
      setEPS()
      par(las=2,cex=0.7,mar=c(5,10,5,5))
      tracks(p1.r,
             n1.r,
             heights=c(ph,nh),
             main=main.title
          )
      ggsave(filename,width=wt,height=ht,units="in")
      rm (neg)
      rm (pos)
   } else if (pos.len != 0) {
      pr.i = reverse(p.i)
      directory = paste("../jpg/",gene.name,sep="")
      setEPS()
      endfile  = paste(filename.preface,".gene.model.isoform.eps",sep="")
      filename = paste(directory,endfile,sep="")
      pos.span = max(end(pb.2.ga)) - min(start(pb.2.ga))

      #----------------------------------------------------------
      #  print a shorter name
      #----------------------------------------------------------
      name = names(pb.2.ga[pr.i])
      nn = strsplit(name,"\\|")
      nnn = unlist(nn)
      e = length(nnn)
      ll = seq(from=1,to=e,by=3)
      pn = nnn[ll]
      pl = pb.2.ga[pr.i]
      names(pl) = pn
#----------------------------------------------------------
# to print out labels as well -- for debugging. -- unfortunately -- order is not preserved unless group.selfish is kept
#  and this causes labels to be printed -- requiring editing later
      p1.r <-autoplot(pl,geom="alignment",group.selfish=T,fill=col.2,color=col.2)
#----------------------------------------------------------
      pos.length = length(pb.2.ga)
      tot.length = pos.length
      ph = pos.length/tot.length
      main.title = paste(gene.name,as.character(pos.span),sep=".")
      setEPS()
      par(las=2,cex=0.7,mar=c(5,10,5,5))

      tracks(p1.r,
             heights=c(ph),
             main=main.title)
      ggsave(filename,width=wt,height=ht,units="in")
      rm(pos)
   } else if (neg.len !=0) {
      directory = paste("../jpg/",gene.name,sep="")
      endfile  = paste(filename.preface,".gene.model.isoform.eps",sep="")
      filename = paste(directory,endfile,sep="")
      neg.span = max(end(pb.1.ga)) - min(start(pb.1.ga))
      nr.i = reverse(n.i)
      #----------------------------------------------------------
      #  print a shorter name
      #----------------------------------------------------------
      name = names(pb.1.ga[nr.i])
      nn = strsplit(name,"\\|")
      nnn = unlist(nn)
      e = length(nnn)
      ll = seq(from=1,to=e,by=3)
      pn = nnn[ll]
      nl = pb.1.ga[nr.i]
      names(nl)<- pn
#----------------------------------------------------------
# to print out labels as well -- for debugging. -- unfortunately -- order is not preserved unless group.selfish is kept
#  and this causes labels to be printed -- requiring editing later
      n1.r <-autoplot(nl,geom="alignment",group.selfish=T,fill=col.1,color=col.1)
#----------------------------------------------------------
      neg.length = length(pb.1.ga)
      tot.length = neg.length
      nh = neg.length/tot.length
      main.title = paste(gene.name,as.character(neg.span),sep=".")
      setEPS()
      par(las=2,cex=0.7,mar=c(5,10,5,5))
      tracks(n1.r,
             heights=c(nh),
             main=main.title)
      ggsave(filename,width=wt,height=ht,units="in")
      rm(neg)
   }
}
