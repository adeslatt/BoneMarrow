
printGviz <- function (gene.name = gene.name, chr=chr,from =from, to=to, font=fontsans, width=wt, height=ht) {

   # GenomeAxis Track, ideoTrack

   axisTrack <- GenomeAxisTrack()
   ideoTrack <- IdeogramTrack(genome="hg19", chromosome=chr,showID = TRUE, showBandId = TRUE)
   sTrack <- SequenceTrack(Hsapiens)

   # GeneRegionTrack(ill.bm.total.gtf_txdb)
   SR.total.gtf_track         <-GeneRegionTrack(ill.bm.total.gtf_txdb,         name = "SR Total gtf",        chromosome=chr, start=from, end=to)
   SR.total.no.gtf_track      <-GeneRegionTrack(ill.bm.total.no.gtf_txdb,      name = "SR Total no gtf",     chromosome=chr, start=from, end=to)
   SR.total.masked.gtf_track  <-GeneRegionTrack(ill.bm.total.no.AZU1.gtf_txdb, name = "SR Total masked gtf", chromosome=chr, start=from, end=to)
   SR.lin.neg.gtf_track       <-GeneRegionTrack(ill.lin.neg.gtf_txdb,          name = "SR Lin- gtf",         chromosome=chr, start=from, end=to)
   SR.lin.neg.no.gtf_track    <-GeneRegionTrack(ill.lin.neg.no.gtf_txdb,       name = "SR Lin- no gtf",      chromosome=chr, start=from, end=to)
   SR.lin.neg.masked.gtf_track<-GeneRegionTrack(ill.lin.neg.no.AZU1.gtf_txdb,  name = "SR Lin- masked gtf",  chromosome=chr, start=from, end=to)
   LR.total_track             <-GeneRegionTrack(pb.53AB.total.gtf_txdb,        name = "LR Total",            chromosome=chr, start=from, end=to)
   LR.lin.neg_track           <-GeneRegionTrack(pb.53A.lin.neg.gtf_txdb,       name = "LR Lin-",             chromosome=chr, start=from, end=to)

   aSR.total.no.gtf.track     <- AlignmentsTrack("../briefCommunications/ill.bm.total.tophat.2.0.14.fr-unstranded.no.gtf.accepted_hits.bam")
   aSR.total.gtf.track        <- AlignmentsTrack("../briefCommunications/ill.bm.total.gtf.accepted_hits.bam",isPaired=TRUE)
   aSR.total.masked.gtf.track <- AlignmentsTrack("../briefCommunications/ill.bm.total.tophat.no.AZU1.ELANE.CFD.accepted_hits.bam",isPaired=TRUE)

   aSR.lin.neg.no.gtf.track     <- AlignmentsTrack("../briefCommunications/ill.lin.neg.no.mm.no.u.no.gtf.accepted_hits.bam",isPaired=TRUE)
   aSR.lin.neg.gtf.track        <- AlignmentsTrack("../briefCommunications/ill.lin.neg.no.mm.no.u.gtf.fa.accepted_hits.bam",isPaired=TRUE)
   aSR.lin.neg.masked.gtf.track <- AlignmentsTrack("../briefCommunications/ill.lin.neg.no.AZU1.ELANE.CFD.accepted_hits.bam",isPaired=TRUE)

   bmt <- BiomartGeneRegionTrack(genome="hg19",chromosome=chr,start = from, end=to, 
       filter=list(with_ox_refseq_mrna=TRUE),stacking="dense")

   filename = paste(paste("../jpg",gene.name,sep="/"),"aSR.SR.LR.total.no.gtf.gviz.eps",sep=".")

# 
#               aSR.total.no.gtf.track,

   postscript(filename, fonts=c("sans"), width = width, height= height)
   plotTracks(c( ideoTrack,
               bmt,
               axisTrack,
               SR.total.no.gtf_track,
               LR.total_track,
               sTrack),
               from=from, to = to,chromosome=chr,type = c("coverage"),cex=0.5,min.height=8)
   dev.off()

   filename = paste(paste("../jpg",gene.name,sep="/"), "aSR.SR.total.gtf.gviz.eps",sep=".")

   postscript(filename, fonts=c(font), width = width, height = height)
   plotTracks(c( ideoTrack,
               axisTrack,
               SR.total.gtf_track,
               LR.total.gtf_track),
               from=from, to = to,chromosome=chr,type = c("coverage"),cex=0.5,min.height=8)
   dev.off()

   filename = paste(paste("../jpg",gene.name,sep="/"), "aSR.SR.LR.lin.neg.no.gtf.gviz.pdf")


   plotTracks(c( ideoTrack,
               bmt,
               axisTrack,
               SR.lin.neg.no.gtf_track,
               LR.lin.neg_track,
               sTrack),
               from=from, to = to,chromosome=chr,type = c("coverage"),cex=0.5,min.height=8)
    dev.off()

   filename = paste(paste("../jpg",gene.name,sep="/"), "aSR.SR.LR.lin.neg.gtf.gviz.pdf")

   plotTracks(c( ideoTrack,
               bmt,
               axisTrack,
               SR.lin.neg.gtf_track,
               LR.lin.neg_track,
               sTrack),
               from=from, to = to,chromosome=chr,type = c("coverage"),cex=0.5,min.height=8)
    dev.off()

    return()
}
