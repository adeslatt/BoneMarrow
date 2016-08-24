printGvizSetupAbundance <- function (gene.name = gene.name, chr=chr,from =from, to=to, font=font, width=width, height=height) {

   SR.total.gtf_track<-GeneRegionTrack(
      ill.bm.total.gtf_txdb,
      name = "SR Total gtf",
      genome="hg19",
      chromosome=chr, 
      start=from, 
      end=to)

   SR.lin.neg.gtf_track<-GeneRegionTrack(
      ill.lin.neg.gtf_txdb,
      name = "SR Lin- gtf",
      genome="hg19",
      chromosome=chr, 
      start=from, 
      end=to)

   LR.total_track<-GeneRegionTrack(
      pb.53AB.total.gtf_txdb,
      name = "LR Total",
      genome="hg19",
      chromosome=chr, 
      start=from, 
      end=to)

   LR.lin.neg_track<-GeneRegionTrack(
      pb.53A.lin.neg.gtf_txdb,
      name = "LR Lin-",
      genome="hg19",
      chromosome=chr,
      start=from,
      end=to)

   aSR.total.gtf.track<-AlignmentsTrack("../briefCommunications/ill.bm.total.gtf.accepted_hits.bam",
      isPaired=TRUE,
      genome="hg19",
      chromosome=chr, 
      coverageHeight=coverageHeight,
      sashimiHeight=sashimiHeight,
      type=type,
      name="SR total gtf",
      start=from, 
      end=to)

   aSR.lin.neg.gtf.track<- AlignmentsTrack("../briefCommunications/ill.lin.neg.no.mm.no.u.gtf.fa.accepted_hits.bam",
      isPaired=TRUE,
      genome="hg19",
      chromosome=chr, 
      coverageHeight=coverageHeight,
      sashimiHeight=sashimiHeight,
      type=type,
      name="SR Lin- gtf",
      start=from, 
      end=to)


   aLR.flnc.lin.neg_track<- AlignmentsTrack("../briefCommunications/53A_lineage_negative_allsize.flnc.fasta.sorted.bam",
      genome="hg19",
      chromosome=chr, 
      coverageHeight=coverageHeight,
      sashimiHeight=sashimiHeight,
      type=type,
      name="LR Lin- flnc",
      start=from, 
      end=to)

   aLR.flnc.total_track  <- AlignmentsTrack("../briefCommunications/53AB_total_allsize.flnc.fasta.sorted.bam",
      genome="hg19",
      chromosome=chr, 
      coverageHeight=coverageHeight,
      sashimiHeight=sashimiHeight,
      type=type,
      name="LR total flnc",
      start=from, 
      end=to)


}