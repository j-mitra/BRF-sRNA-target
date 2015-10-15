#Author: Joydeep Mitra
#Email: jmitra@uga.edu
#Date: 21st September, 2015

library(seqinr)
library(randomForest)

#Read in FASTA files for sRNAs and corresponding mRNAs
srna.seqs<-read.fasta("./srna_fasta.txt", as.string = TRUE, set.att = FALSE,seqonly=TRUE)
mrna.seqs<-read.fasta("./mrna_fasta.txt", as.string = TRUE, set.att = FALSE,seqonly=TRUE)

#Stop if number of sequences in sRNA and target mRNA files do not match
if(length(srna.seqs)!=length(mrna.seqs))
{
  stop("Number of sequences in sRNA and mRNA files do not match!")
}


#Function that translates sequence into alternate alphabets
twolettercode<-function(s,code)
{
  if(code=="RY") #R/Y
  {
    tls<-gsub("a","r",s)
    tls<-gsub("g","r",tls)
    tls<-gsub("c","y",tls)
    tls<-gsub("u","y",tls)
    return(tls)
  }
  else if(code=="MK") #K/M
  {
    tls<-gsub("a","m",s)
    tls<-gsub("c","m",tls)
    tls<-gsub("g","k",tls)
    tls<-gsub("u","k",tls)
    return(tls)
  }
  else if(code=="WS") #S/W
  {
    tls<-gsub("g","s",s)
    tls<-gsub("c","s",tls)
    tls<-gsub("a","w",tls)
    tls<-gsub("u","w",tls)
    return(tls)
  }
  else
  {
    return(s)
  } 
}

#Store sRNA Patterns and names in a list
srna.patterns<-read.table("./srnaPatterns.txt",header=T,sep="\t")

srna.pattern.namelist<-list("AUGC"=c(),"RY"=c(),"MK"=c())
srna.pattern.patlist<-list("AUGC"=c(),"RY"=c(),"MK"=c())
for(i in 1:nrow(srna.patterns))
{
  if(srna.patterns[i,4]==0)
  {
    srna.pattern.namelist$AUGC<-c(srna.pattern.namelist$AUGC,as.character(srna.patterns[i,1]))
    srna.pattern.patlist$AUGC<-c(srna.pattern.patlist$AUGC,as.character(srna.patterns[i,2]))
  }
  else if(srna.patterns[i,4]==1)
  {
    srna.pattern.namelist$RY<-c(srna.pattern.namelist$RY,as.character(srna.patterns[i,1]))
    srna.pattern.patlist$RY<-c(srna.pattern.patlist$RY,as.character(srna.patterns[i,2]))
  }
  else if(srna.patterns[i,4]==3)
  {
    srna.pattern.namelist$MK<-c(srna.pattern.namelist$MK,as.character(srna.patterns[i,1]))
    srna.pattern.patlist$MK<-c(srna.pattern.patlist$MK,as.character(srna.patterns[i,2]))
  }
}


#Store mRNA Patterns and names in a data structure
mrna.patterns<-read.table("./mrnaPatterns.txt",header=T,sep="\t")

mrna.pattern.namelist<-list("AUGC"=c(),"RY"=c(),"MK"=c(),"WS"=c())
mrna.pattern.patlist<-list("AUGC"=c(),"RY"=c(),"MK"=c(),"WS"=c())

for(i in 1:nrow(mrna.patterns))
{
  if(mrna.patterns[i,4]==0)
  {
    mrna.pattern.namelist$AUGC<-c(mrna.pattern.namelist$AUGC,as.character(mrna.patterns[i,1]))
    mrna.pattern.patlist$AUGC<-c(mrna.pattern.patlist$AUGC,as.character(mrna.patterns[i,2]))
  }
  else if(mrna.patterns[i,4]==1)
  {
    mrna.pattern.namelist$RY<-c(mrna.pattern.namelist$RY,as.character(mrna.patterns[i,1]))
    mrna.pattern.patlist$RY<-c(mrna.pattern.patlist$RY,as.character(mrna.patterns[i,2]))
  }
  else if(mrna.patterns[i,4]==2)
  {
    mrna.pattern.namelist$WS<-c(mrna.pattern.namelist$WS,as.character(mrna.patterns[i,1]))
    mrna.pattern.patlist$WS<-c(mrna.pattern.patlist$WS,as.character(mrna.patterns[i,2]))
  }
  else if(mrna.patterns[i,4]==3)
  {
    mrna.pattern.namelist$MK<-c(mrna.pattern.namelist$MK,as.character(mrna.patterns[i,1]))
    mrna.pattern.patlist$MK<-c(mrna.pattern.patlist$MK,as.character(mrna.patterns[i,2]))
  }
}



#Function to calculate frequencies for each kmer pattern
compute_patternFrequencies<-function(inseq,seq_type)
{
  inseq<-toupper(inseq)
  if(seq_type=="srna")
  {
    pattern.namelist<-srna.pattern.namelist
    pattern.patlist<-srna.pattern.patlist
  }
  else if(seq_type=="mrna")
  {
    pattern.namelist<-mrna.pattern.namelist
    pattern.patlist<-mrna.pattern.patlist
  }
  frequencies<-c()
  freq.names<-c()
  no.alphabets<-length(pattern.namelist)
  for(i in 1:no.alphabets)
  {
    current.alphabet<-names(pattern.patlist)[i]
    no.patterns<-length(pattern.patlist[[i]])
    
    #Translate the sequence to current alphabet
    current.sequence<-twolettercode(inseq,current.alphabet)
    
    freq.names<-c(freq.names,pattern.namelist[[i]])
    for(j in 1:no.patterns)
    {
      current.pattern<-pattern.patlist[[i]][j]
      
      current.freq<-length(gregexpr(current.pattern,current.sequence)[[1]])/(nchar(current.sequence)-(nchar(current.pattern)-1))
      frequencies<-c(frequencies,current.freq)
      
      
      
    }
    
  }
  
  names(frequencies)<-freq.names
  return(frequencies)
  
}


#Calculate the frequencies of all sequences in the file and prepare input to sRNA-BRF
for(k in 1:length(srna.seqs))
{
    srna.frequencies<-compute_patternFrequencies(srna.seqs[[k]],"srna")
    mrna.frequencies<-compute_patternFrequencies(mrna.seqs[[k]],"mrna")
    
    if(k==1)
    {
      all.srna.freqs<-srna.frequencies
      all.mrna.freqs<-mrna.frequencies
    }
    else if(k>1)
    {
      all.srna.freqs<-rbind(all.srna.freqs,srna.frequencies)
      all.mrna.freqs<-rbind(all.mrna.freqs,mrna.frequencies)
    }
}

all.frequencies<-cbind(all.srna.freqs,all.mrna.freqs)
row.names(all.frequencies)<-NULL

#Load saved BRF model and make predictions
load("./BRF_model.Rdata")
target.predictions<-predict(best.model.mcc,all.frequencies,type="prob")
write.table(target.predictions[,2],"./prediction_probabilities.txt",row.names=F,quote=F,col.names=F)
#write.table(all.frequencies,"./predict_data.txt",sep="\t",quote=F,row.names=F)