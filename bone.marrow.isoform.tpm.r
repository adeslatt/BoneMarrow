# 2016 August 25
# File:  bone.marrow.isoform.tpm.r
# Purpose:  to read in the abundance.csv file
#           the csv file in this case had been externally modified from the
#           count file provided as output from ToFU.
#           The modification was just the splitting out of the isoform number 
#           from the pacbio id.
#           This was so isoform informaiton could be captured.
#           I merged this with blast output information. -- I will change this to 
#           work with matchAnnot output for the future.
# Author:   Anne Deslattes Mays, PhD
#
#------------------------------------------------------------------------------------------
# PacBio ICE Sample Specific Determined Abundance 43A Lineage Negative
#------------------------------------------------------------------------------------------
pb.43A.lin.neg.abundance = read.csv("../sailfishv5/43A_negative.all.merge.collapsed.abundance.csv",
                         colClasses=c(transcript        = "character",
                                      pbgene            = "character",
                                      isoform           = "numeric",
                                      count_fl          = "numeric",
                                      count_nfl         = "numeric",
                                      count_nfl_amb     = "numeric",
                                      abundance_fl      = "numeric",
                                      abundance_nfl     = "numeric",
                                      abundance_nfl_amb = "numeric"),
                                       as.is=T,header=T)
#> dim(pb.43A.lin.neg.abundance)
#[1] 8167    9
#------------------------------------------------------------------------------------------
# PacBio ICE Sample Specific Determined Abundance 53AB_total
#------------------------------------------------------------------------------------------
pb.53AB.total.abundance = read.csv("../sailfishv5/53AB_total.good.5merge.collapsed.abundance.csv",
                         colClasses=c(transcript        = "character",
                                      pbgene            = "character",
                                      isoform           = "numeric",
                                      count_fl          = "numeric",
                                      count_nfl         = "numeric",
                                      count_nfl_amb     = "numeric",
                                      abundance_fl      = "numeric",
                                      abundance_nfl     = "numeric",
                                      abundance_nfl_amb = "numeric"),
                                       as.is=T,header=T)
#> dim(pb.53AB.total.abundance)
#[1] 16856     9
#> dim(pb.53AB_total.blast)
#[1] 16178    12
#------------------------------------------------------------------------------------------
# PacBio ICE Sample Specific Determined Abundance 53A Lineage Negative
#------------------------------------------------------------------------------------------
pb.53A.lin.neg.abundance = read.csv("../sailfishv5/53A_negative.good.5merge.collapsed.abundance.csv",
                         colClasses=c(transcript        = "character",
                                      pbgene            = "character",
                                      isoform           = "numeric",
                                      count_fl          = "numeric",
                                      count_nfl         = "numeric",
                                      count_nfl_amb     = "numeric",
                                      abundance_fl      = "numeric",
                                      abundance_nfl     = "numeric",
                                      abundance_nfl_amb = "numeric"),
                                       as.is=T,header=T)
#> dim(pb.53A.lin.neg.abundance)
#[1] 12225     9
# dim(pb.53A_negative.blast)
#[1] 11953    12
# ... missing several hundred hits....
#------------------------------------------------------------------------------------------
# PacBio ICE Sample Specific Determined Abundance 53A Lineage Positive
#------------------------------------------------------------------------------------------
pb.53A.lin.pos.abundance = read.csv("../sailfishv5/53A_positive.good.5merge.collapsed.abundance.csv",
                         colClasses=c(transcript        = "character",
                                      pbgene            = "character",
                                      isoform           = "numeric",
                                      count_fl          = "numeric",
                                      count_nfl         = "numeric",
                                      count_nfl_amb     = "numeric",
                                      abundance_fl      = "numeric",
                                      abundance_nfl     = "numeric",
                                      abundance_nfl_amb = "numeric"),
                                       as.is=T,header=T)
#> dim(pb.53A.lin.pos.abundance)
#[1] 5865    9
#> dim(pb.53A_positive.blast)
#[1] 5700   12
#------------------------------------------------------------------------------------------
# PacBio Only blast Annotation Dataframes
#------------------------------------------------------------------------------------------
pb.43A.lin.neg.blast.df     <- data.frame(id=pb.43A_negative.blast[,1],pb.43A_negative.blast)
pb.53AB.total.blast.df      <- data.frame(id=pb.53AB_total.blast[,1],  pb.53AB_total.blast)
pb.53A.lin.neg.blast.df     <- data.frame(id=pb.53A_negative.blast[,1],pb.53A_negative.blast)
pb.53A.lin.pos.blast.df     <- data.frame(id=pb.53A_positive.blast[,1],pb.53A_positive.blast)

#------------------------------------------------------------------------------------------
# PacBio ICE quantification dataframes
#------------------------------------------------------------------------------------------
pb.43A.lin.neg.abundance.df     <- data.frame(id=pb.43A.lin.neg.abundance[,1],pb.43A.lin.neg.abundance)
pb.53AB.total.abundance.df      <- data.frame(id=pb.53AB.total.abundance[,1], pb.53AB.total.abundance)
pb.53A.lin.neg.abundance.df     <- data.frame(id=pb.53A.lin.neg.abundance[,1],pb.53A.lin.neg.abundance)
pb.53A.lin.pos.abundance.df     <- data.frame(id=pb.53A.lin.pos.abundance[,1],pb.53A.lin.pos.abundance)

#------------------------------------------------------------------------------------------
# PacBio ICE quantification and Blast merge
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg <- merge(pb.43A.lin.neg.blast.df,pb.43A.lin.neg.abundance.df,by="id")
pb.ICE.53AB.total  <- merge(pb.53AB.total.blast.df, pb.53AB.total.abundance.df,by="id")
pb.ICE.53A.lin.neg <- merge(pb.53A.lin.neg.blast.df,pb.53A.lin.neg.abundance.df,by="id")
pb.ICE.53A.lin.pos <- merge(pb.53A.lin.pos.blast.df,pb.53A.lin.pos.abundance.df,by="id")

#------------------------------------------------------------------------------------------
# Creating the TPM Metric for ease of comparisons -- adding it to the above datastructure
#------------------------------------------------------------------------------------------
#
# Creating the TPM metric.    Here is the translation of the Bo Li paper
# for ease to apply it to our data structure
#
# Bo Li Bioinformatics 2010 doi:10.1093/bioinformatics/btp692 
# "RNA-Seq gene expression estimation with read mapping uncertainty"
#
# fraction of nucleotides vi of the transcriptome
#
#      ti*li       the length of the transcript* that transcript total
# vi =_______  ==  ___________________________________________________
#
#     sum(tj*lj)  sum(length of all transcripts * quantity of all transcripts)
#
# fraction of transcript ti of the transcriptome
#
#       vi/li          fraction of nucleotides vi of the transcriptome / transcript ti length
# ti= __________==___________________________________________________________________________
#     sum(vj/lj)   sum( fraction of all nucleotides vj of the transcriptome / transcript tj length)
#
#
# Translating to our data structure for pac bio we have for each isoform
#
#     
#    pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]      pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]
# vi=_____________________________________________  =___________________________________________
#    sum(pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"])      transcriptome_total_nucleotides
#
#
# and
#
#        vi/li       (pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]/sum(pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"]))/pb.ln[i,"length"]
#  ti = ________ =  ________________________________________________________________________________________________________________   
#      sum(vj/lj)  sum((pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"]/sum(pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"]))/pb.ln[i,"length"])
#
#
#                     ((pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]/ttn) /pb.ln[i,"length"]
#                  =_________________________________________________________________________________________________________________
#                    sum(((pb.ln[,"count_nfl_amb"]*pb.ln[,"length"]/ttn) /pb.ln[,"length"]))
#
#
#                     ((pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]/ttn) /pb.ln[i,"length"]
#                  =_________________________________________________________________________________________________________________
#                                     ttt
#                  
#  ttn = Transcriptome_total_nucleotides
#  ttt = transcriptome_total_transcripts
#
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# now make TPM for pac bio provided abundances  (post merge)
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg.ttn = sum(pb.ICE.43A.lin.neg[,"pblength"]*pb.ICE.43A.lin.neg[,"count_nfl_amb"])
pb.ICE.43A.lin.neg.ttt = sum(((pb.ICE.43A.lin.neg[,"count_nfl_amb"]*pb.ICE.43A.lin.neg[,"pblength"]/pb.ICE.43A.lin.neg.ttn) /pb.ICE.43A.lin.neg[,"pblength"]))

pb.ICE.53A.lin.neg.ttn = sum(pb.ICE.53A.lin.neg[,"pblength"]*pb.ICE.53A.lin.neg[,"count_nfl_amb"])
pb.ICE.53A.lin.neg.ttt = sum(((pb.ICE.53A.lin.neg[,"count_nfl_amb"]*pb.ICE.53A.lin.neg[,"pblength"]/pb.ICE.53A.lin.neg.ttn) /pb.ICE.53A.lin.neg[,"pblength"]))

pb.ICE.53A.lin.pos.ttn = sum(pb.ICE.53A.lin.pos[,"pblength"]*pb.ICE.53A.lin.pos[,"count_nfl_amb"])
pb.ICE.53A.lin.pos.ttt = sum(((pb.ICE.53A.lin.pos[,"count_nfl_amb"]*pb.ICE.53A.lin.pos[,"pblength"]/pb.ICE.53A.lin.pos.ttn) /pb.ICE.53A.lin.pos[,"pblength"]))

pb.ICE.53AB.total.ttn = sum(pb.ICE.53AB.total[,"pblength"]*pb.ICE.53AB.total[,"count_nfl_amb"])
pb.ICE.53AB.total.ttt = sum(((pb.ICE.53AB.total[,"count_nfl_amb"]*pb.ICE.53AB.total[,"pblength"]/pb.ICE.53AB.total.ttn) /pb.ICE.53AB.total[,"pblength"]))

#------------------------------------------------------------------------------------------
#  Add the TPM to the PacBio data structure
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg <- data.frame(pb.ICE.43A.lin.neg, 
   TPM = (((pb.ICE.43A.lin.neg[,"count_nfl_amb"]*pb.ICE.43A.lin.neg[,"pblength"]/pb.ICE.43A.lin.neg.ttn) /pb.ICE.43A.lin.neg[,"pblength"])/pb.ICE.43A.lin.neg.ttt)*10^6)

pb.ICE.53A.lin.neg <- data.frame(pb.ICE.53A.lin.neg, 
   TPM = (((pb.ICE.53A.lin.neg[,"count_nfl_amb"]*pb.ICE.53A.lin.neg[,"pblength"]/pb.ICE.53A.lin.neg.ttn) /pb.ICE.53A.lin.neg[,"pblength"])/pb.ICE.53A.lin.neg.ttt)*10^6)

pb.ICE.53A.lin.pos <- data.frame(pb.ICE.53A.lin.pos, 
   TPM = (((pb.ICE.53A.lin.pos[,"count_nfl_amb"]*pb.ICE.53A.lin.pos[,"pblength"]/pb.ICE.53A.lin.pos.ttn) /pb.ICE.53A.lin.pos[,"pblength"])/pb.ICE.53A.lin.pos.ttt)*10^6)

pb.ICE.53AB.total <- data.frame(pb.ICE.53AB.total, 
   TPM = (((pb.ICE.53AB.total[,"count_nfl_amb"]*pb.ICE.53AB.total[,"pblength"]/pb.ICE.53AB.total.ttn) /pb.ICE.53AB.total[,"pblength"])/pb.ICE.53AB.total.ttt)*10^6)


#------------------------------------------------------------------------------------------
# PacBio TPM derived from the abundancie and count_nfl_amb data (as opposed to Sailfish)
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg.i <- order(pb.ICE.43A.lin.neg[,"TPM"],decreasing=T)
pb.ICE.53A.lin.neg.i <- order(pb.ICE.53A.lin.neg[,"TPM"],decreasing=T)
pb.ICE.53A.lin.pos.i <- order(pb.ICE.53A.lin.pos[,"TPM"],decreasing=T)
pb.ICE.53AB.total.i  <- order(pb.ICE.53AB.total [,"TPM"],decreasing=T)

#------------------------------------------------------------------------------------------
# Creating the TPM Metric for ease of comparisons -- adding it to the above datastructure
#------------------------------------------------------------------------------------------
#
# Creating the TPM metric.    Here is the translation of the Bo Li paper
# for ease to apply it to our data structure
#
# Bo Li Bioinformatics 2010 doi:10.1093/bioinformatics/btp692 
# "RNA-Seq gene expression estimation with read mapping uncertainty"
#
# fraction of nucleotides vi of the transcriptome
#
#      ti*li       the length of the transcript* that transcript total
# vi =_______  ==  ___________________________________________________
#
#     sum(tj*lj)  sum(length of all transcripts * quantity of all transcripts)
#
# fraction of transcript ti of the transcriptome
#
#       vi/li          fraction of nucleotides vi of the transcriptome / transcript ti length
# ti= __________==___________________________________________________________________________
#     sum(vj/lj)   sum( fraction of all nucleotides vj of the transcriptome / transcript tj length)
#
#
# Translating to our data structure for pac bio we have for each isoform
#
#     
#    pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]      pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]
# vi=_____________________________________________  =___________________________________________
#    sum(pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"])      transcriptome_total_nucleotides
#
#
# and
#
#        vi/li       (pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]/sum(pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"]))/pb.ln[i,"length"]
#  ti = ________ =  ________________________________________________________________________________________________________________   
#      sum(vj/lj)  sum((pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"]/sum(pb.ln[j,"count_nfl_amb"]*pb.ln[j,"length"]))/pb.ln[i,"length"])
#
#
#                     ((pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]/ttn) /pb.ln[i,"length"]
#                  =_________________________________________________________________________________________________________________
#                    sum(((pb.ln[,"count_nfl_amb"]*pb.ln[,"length"]/ttn) /pb.ln[,"length"]))
#
#
#                     ((pb.ln[i,"count_nfl_amb"]*pb.ln[i,"length"]/ttn) /pb.ln[i,"length"]
#                  =_________________________________________________________________________________________________________________
#                                     ttt
#                  
#  ttn = Transcriptome_total_nucleotides
#  ttt = transcriptome_total_transcripts
#
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# now make TPM for pac bio provided abundances  (post merge)
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg.ttn = sum(pb.ICE.43A.lin.neg[,"pblength"]*pb.ICE.43A.lin.neg[,"count_nfl_amb"])
pb.ICE.43A.lin.neg.ttt = sum(((pb.ICE.43A.lin.neg[,"count_nfl_amb"]*pb.ICE.43A.lin.neg[,"pblength"]/pb.ICE.43A.lin.neg.ttn) /pb.ICE.43A.lin.neg[,"pblength"]))

pb.ICE.53A.lin.neg.ttn = sum(pb.ICE.53A.lin.neg[,"pblength"]*pb.ICE.53A.lin.neg[,"count_nfl_amb"])
pb.ICE.53A.lin.neg.ttt = sum(((pb.ICE.53A.lin.neg[,"count_nfl_amb"]*pb.ICE.53A.lin.neg[,"pblength"]/pb.ICE.53A.lin.neg.ttn) /pb.ICE.53A.lin.neg[,"pblength"]))

pb.ICE.53A.lin.pos.ttn = sum(pb.ICE.53A.lin.pos[,"pblength"]*pb.ICE.53A.lin.pos[,"count_nfl_amb"])
pb.ICE.53A.lin.pos.ttt = sum(((pb.ICE.53A.lin.pos[,"count_nfl_amb"]*pb.ICE.53A.lin.pos[,"pblength"]/pb.ICE.53A.lin.pos.ttn) /pb.ICE.53A.lin.pos[,"pblength"]))

pb.ICE.53AB.total.ttn = sum(pb.ICE.53AB.total[,"pblength"]*pb.ICE.53AB.total[,"count_nfl_amb"])
pb.ICE.53AB.total.ttt = sum(((pb.ICE.53AB.total[,"count_nfl_amb"]*pb.ICE.53AB.total[,"pblength"]/pb.ICE.53AB.total.ttn) /pb.ICE.53AB.total[,"pblength"]))

#------------------------------------------------------------------------------------------
#  Add the TPM to the PacBio data structure
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg <- data.frame(pb.ICE.43A.lin.neg, 
   TPM = (((pb.ICE.43A.lin.neg[,"count_nfl_amb"]*pb.ICE.43A.lin.neg[,"pblength"]/pb.ICE.43A.lin.neg.ttn) /pb.ICE.43A.lin.neg[,"pblength"])/pb.ICE.43A.lin.neg.ttt)*10^6)

pb.ICE.53A.lin.neg <- data.frame(pb.ICE.53A.lin.neg, 
   TPM = (((pb.ICE.53A.lin.neg[,"count_nfl_amb"]*pb.ICE.53A.lin.neg[,"pblength"]/pb.ICE.53A.lin.neg.ttn) /pb.ICE.53A.lin.neg[,"pblength"])/pb.ICE.53A.lin.neg.ttt)*10^6)

pb.ICE.53A.lin.pos <- data.frame(pb.ICE.53A.lin.pos, 
   TPM = (((pb.ICE.53A.lin.pos[,"count_nfl_amb"]*pb.ICE.53A.lin.pos[,"pblength"]/pb.ICE.53A.lin.pos.ttn) /pb.ICE.53A.lin.pos[,"pblength"])/pb.ICE.53A.lin.pos.ttt)*10^6)

pb.ICE.53AB.total <- data.frame(pb.ICE.53AB.total, 
   TPM = (((pb.ICE.53AB.total[,"count_nfl_amb"]*pb.ICE.53AB.total[,"pblength"]/pb.ICE.53AB.total.ttn) /pb.ICE.53AB.total[,"pblength"])/pb.ICE.53AB.total.ttt)*10^6)


#------------------------------------------------------------------------------------------
# PacBio TPM derived from the abundancie and count_nfl_amb data (as opposed to Sailfish)
#------------------------------------------------------------------------------------------
pb.ICE.43A.lin.neg.i <- order(pb.ICE.43A.lin.neg[,"TPM"],decreasing=T)
pb.ICE.53A.lin.neg.i <- order(pb.ICE.53A.lin.neg[,"TPM"],decreasing=T)
pb.ICE.53A.lin.pos.i <- order(pb.ICE.53A.lin.pos[,"TPM"],decreasing=T)
pb.ICE.53AB.total.i  <- order(pb.ICE.53AB.total [,"TPM"],decreasing=T)

