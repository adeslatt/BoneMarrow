printGvizSetupFirstCall <- function () {

   library(BSgenome.Hsapiens.UCSC.hg19)
   library(GenomicFeatures)
   library(Gviz)
   library(BSgenome.Hsapiens.UCSC.hg19)

seqdb    <-BSgenome.Hsapiens.UCSC.hg19

chrominfo<-data.frame(chrom=seqnames(seqdb), 
            length=seqlengths(seqdb), 
            is_circular=FALSE)

genes.gtf.txdb <- makeTxDbFromGFF(file="../briefCommunications/genes.gtf",
   format="gtf", 
   dataSource = "genes.gtf",
   organism = "Homo sapiens",
   chrominfo = chrominfo)

genes.hg38.gtf.txdb <- makeTxDbFromGFF(file="../briefCommunications/genes.hg38.gtf",
   format="gtf", 
   dataSource = "genes.hg38.gtf",
   organism = "Homo sapiens",
   chrominfo = chrominfo)

uniprot.gtf_txdb <- makeTxDbFromGFF(file = "../briefCommunications/uniprot-taxonomy.gff.gz",
   format="gff",
   dataSource = "Uniprot",
   organism = "Homo sapiens",
   chrominfo = chrominfo)

ill.bm.total.gtf_txdb <-makeTxDbFromGFF(file="../briefCommunications/ill.bm.total.gtf.m.transcripts.5.19.gtf",
   format="gtf", 
   dataSource = "Illumina Total gtf",
   organism = "Homo sapiens",
   chrominfo = chrominfo)

ill.lin.neg.gtf_txdb <-makeTxDbFromGFF(file="../briefCommunications/ill.lin.neg.gtf.transcripts.5.19.gtf",
   format="gtf",
   dataSource = "Illumina Lin- gtf",
   organism = "Homo sapiens",
   chrominfo = chrominfo)

pb.53AB.total.gtf_txdb<-makeTxDbFromGFF(file="../briefCommunications/53AB_total.good.5merge.collapsed.longest_rep.fa.gtf",
   format="gtf",
   dataSource = "PB Total",
   organism = "Homo sapiens",
   chrominfo = chrominfo)

pb.53A.lin.neg.gtf_txdb<-makeTxDbFromGFF(file="../briefCommunications/53A_negative.good.5merge.collapsed.longest_rep.fa.gtf",
   format="gtf",
   dataSource = "PB Lin-",
   organism = "Homo sapiens",
   chrominfo = chrominfo)


   return()
}