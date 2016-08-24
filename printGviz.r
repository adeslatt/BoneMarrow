printGviz <- function (gene.name = gene.name, chr=chr,from =from, to=to, font=font, width=width, height=height) {

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"SR.no.gtf.LR.lin.neg.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               SR.lin.neg.no.gtf_track,
               LR.lin.neg_track),
               from=from, to = to,chromosome=chr)
   dev.off()

#   setEPS()
#   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aSR.SR.no.gtf.LR.lin.neg.gviz.eps",sep=".")
#   postscript(filename, fonts=c("sans"), width = wt, height= ht)
#   plotTracks(list(ideoTrack,
#               bmt,
#               axisTrack,
#               aSR.lin.neg.no.gtf.track,
#               SR.lin.neg.no.gtf_track,
#               LR.lin.neg_track),
#               coverageHeight=coverageHeight,
#               sashimiHeight=sashimiHeight,
#               from=from, to = to,chromosome=chr)
#   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"SR.gtf.LR.lin.neg.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               SR.lin.neg.gtf_track,
               LR.lin.neg_track),
               from=from, to = to,chromosome=chr)
   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aSR.SR.gtf.lin.neg.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               aSR.lin.neg.gtf.track,
               SR.lin.neg.gtf_track),
               coverageHeight=coverageHeight,
               from=from, to = to,chromosome=chr)
   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aLR.flnc.lin.neg.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               aLR.flnc.lin.neg_track),
               type = "coverage",
               from=from, to = to,chromosome=chr)
   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aLR.flnc.lin.neg.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               aLR.flnc.lin.neg_track),
               type = c("coverage","sashimi"),
               coverageHeight=coverageHeight,
               sashimiHeight=sashimiHeight,
               from=from, to = to,chromosome=chr)
   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aLR.flnc.total.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               aLR.flnc.total_track),
               type = c("coverage","sashimi"),
               coverageHeight=coverageHeight,
               sashimiHeight=sashimiHeight,
               from=from, to = to,chromosome=chr)
   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aSR.gtf.lin.neg.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               aSR.lin.neg.gtf.track),
               coverageHeight=coverageHeight,
               sashimiHeight=sashimiHeight,
               from=from, to = to,chromosome=chr)
   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"SR.no.gtf.LR.total.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               SR.total.no.gtf_track,
               LR.total_track),
               from=from, to = to,chromosome=chr)
   dev.off()

#   setEPS()
#   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aSR.SR.no.gtf.total.gviz.eps",sep=".")
#   postscript(filename, fonts=c("sans"), width = wt, height= ht)
#   plotTracks(list(
#               ideoTrack,
#               bmt,
#               axisTrack,
#               aSR.total.no.gtf.track,
#               SR.total.no.gtf_track),
#               coverageHeight=coverageHeight,
#               sashimiHeight=sashimiHeight,
#               from=from, to = to,chromosome=chr)
#   dev.off()

   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"aLR.LR.flnc.total.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               aLR.flnc.total_track,
               LR.total_track),
               coverageHeight=coverageHeight,
               sashimiHeight=sashimiHeight,
               from=from, to = to,chromosome=chr)
   dev.off()


   setEPS()
   filename = paste(paste("../briefCommunications",gene.name,sep="/"),"SR.gtf.LR.total.gviz.eps",sep=".")
   postscript(filename, fonts=c("sans"), width = wt, height= ht)
   plotTracks(list(
               ideoTrack,
               bmt,
               axisTrack,
               SR.total.gtf_track,
               LR.total_track),
               from=from, to = to,chromosome=chr)
   dev.off()

    

   return()
}
