printGvizSetup <- function (gene.name = gene.name, chr=chr,from =from, to=to, font=font, width=width, height=height) {

   # gene model track
   type = c("coverage","sashimi")
   axisTrack <- GenomeAxisTrack()
   stackHeight = 0.2
   # IdeogramTrack
   ideoTrack <- IdeogramTrack(genome="hg19", chromosome=chr,showID = TRUE, showBandId = TRUE)

   data(geneModels)
   uniprot.gtf_track <- GeneRegionTrack(
      uniprot.gtf_txdb,
      name="uniprot",
      genome="hg19",
      chromosome=chr, 
      start=from, 
      end=to)

   ucsc.genes <- GeneRegionTrack(
      genes.gtf.txdb,
      name = gene.name,
      genome="hg19",
      chromosome=chr, 
      start=from, 
      end=to)


   return()
}