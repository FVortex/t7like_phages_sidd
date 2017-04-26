#all in one directory

wd <- "/home/hp/t7like_phages_sidd/"

library(RCurl)
library(seqinr)

myfile <- getURL('ftp://ftp.genome.jp/pub/db/virushostdb/taxid2lineage_full_VH.tsv', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
viralhostdb <- read.table(textConnection(myfile), header=T, sep = '\t', fill = T)

viralhostdb_Autographivirinae <- viralhostdb[grepl('Autographivirinae', viralhostdb$cellular.organisms),]
taxids <- as.vector(viralhostdb_Autographivirinae$X2)
library(rentrez)
#find sequenced linked to a taxid


#loop for all phages
all_Autographivirinae_seqs <- list()
all_Autographivirinae_names <- c()
for (i in seq_along(taxids)) {
  tax_seqs <- entrez_link(db = "nuccore", dbfrom = "taxonomy", id=taxids[i])
#elink result with ids from 2 databases:
#[1] taxonomy_nuccore        taxonomy_nucleotide_exp

#grab them

tmp <- tempfile()
recs <- entrez_fetch(db="nuccore", id=tax_seqs$links$taxonomy_nuccore, rettype="fasta")
cat(recs, file=tmp)
seqs <- as.list(ape::read.dna(tmp, format="fasta")) #since some are store as list, some as matrix coercion is needed
#print(str(seqs)) 
###seqs <- as.character(seqs[which.max(as.vector(sapply(seqs, length))),])

seqs <- as.character(seqs[which.max(as.vector(sapply(seqs, length)))]) # for 1st etc
all_Autographivirinae_names <- c(all_Autographivirinae_names, names(seqs))
## if (TRUE%in%grepl('complete', labels(seqs), ignore.case=TRUE)&length(seqs)!=1) {
###seqs <- seqs[which(grepl('complete', labels(seqs), ignore.case=TRUE))[1],]
#print(str(seqs))
###  ###} else {
##  ######seqs <- seqs[1,]    
#print(str(seqs))  
##  }
print(str(seqs))  
##  ###if (length(which.max(lapply(as.character(seqs), length)))!=1) {
#assign(unlist(as.character(seqs)[which.max(lapply(as.character(seqs), length))]), 
##      #  paste0(viralhostdb_Autographivirinae[i,2], 'seq_from_nuccore'))
seqs <- as.character(unlist(seqs))
all_Autographivirinae_names <- c(all_Autographivirinae_names,  names(seqs)[1])
all_Autographivirinae_seqs[[i]]<- seqs
###names(all_Autographivirinae_seqs[i]) <- tmp2
#####names(all_Autographivirinae_seqs[[i]]) <- names(all_Autographivirinae_seqs[[i]])[1]
#  #unlist(as.character(seqs)[which.max(lapply(as.character(seqs), length))])
}
all_Autographivirinae_names <- gsub(' ', '_', all_Autographivirinae_names)
names(all_Autographivirinae_seqs) <- gsub(' ', '_', all_Autographivirinae_names)
all_Autographivirinae_seqs <- lapply(all_Autographivirinae_seqs, toupper)
#This is really simple with Entrez Direct
#lapply(all_Autographivirinae_seqs, length)
##system('epost -db taxonomy -id 10759 | elink -target nuccore | efetch -format uid')


####Cutting genomes into pieces

#fucntion to split
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}


for (i in seq_along(all_Autographivirinae_seqs)) {
  assign(paste0(all_Autographivirinae_names[[i]], '_by_chunks'), splitWithOverlap(unlist(all_Autographivirinae_seqs[[i]]), 10000, 2000))
}

phages_chunks<-sort(grep(pattern = '*_by_chunks$', ls(), value = T))
#t7chunks<-splitWithOverlap(t7gb, 10000, 2000)
lapply(phages_chunks, function(x) {str(get(x))})


#t7chunks<-splitWithOverlap(seq_along(t7gb), 10000, 2000)

dirs <- c('phages_by_10kbp_chunks', 'phages_by_10kbp_sidd')

for (i in dirs) {
  if (file.exists(paste0(wd, '/', i))) {
    # #unlink(paste0(wd, '/', i), recursive = T)
  } else {
    dir.create(paste0(wd, '/', i))  
  }
}


for (i in phages_chunks) {
  for (j in seq_along(get(i))) {
    write.fasta(get(i)[[j]], names = NULL, file.out = paste0(wd, 'phages_by_10kbp_chunks/', i, '_', j))
  }
}

phages_chunks

for (i in phages_chunks) {
  for (j in seq_along(get(i))) {
    # #aa<-as.character(read.fasta(paste0(wd, '/t7_genome_parts_string/', i), as.string = T, set.attributes = F))
    # #print(substr(aa, nchar(aa)-1000-10, nchar(aa)-1000+1))
    system(paste0('cd ', wd, 'sist/
                  ',
                  'perl -X master.pl -a M -f ', wd, 'phages_by_10kbp_chunks/', i, '_', j, ' -o' , wd, 'phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'))
  }
}
#phages_chunks <- phages_chunks[-c(2,5)]

for (i in phages_chunks) {
  interm<-c()
  print(i)
  for (j in seq_along(get(i))) {
    print(j)
    if (j==1) {
      touch<-read.csv(paste0(wd, '/phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'), sep = '\t', skip = 1)
      touch<-touch[-(9001:10000),]
    }
    else if (j==length(get(i))){
      
      if(i=="k11gbchunks") { touch<-read.csv(paste0(wd, '/phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'), sep = '\t', skip = 1)  #skipping needed since warnong message 'sequence if too short' is present
      } else {
        touch<-read.csv(paste0(wd, '/phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'), sep = '\t')
        touch<-touch[-(1:1000),]
      }
    }
    else{
      touch<-read.csv(paste0(wd, '/phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'),sep = '\t', skip=1)
      touch<-touch[-c(1:1000, 9001:10000),]
    } #warnings omitting
    interm<-rbind(interm, as.matrix(touch))
    
  } 
  assign(paste0(i, '_sidd'), interm)
}

#letters check-up
letters_tables <- lapply(all_Autographivirinae_seqs, function (x) {table(unlist(x))})
phages_with_spec_letters <- all_Autographivirinae_seqs[which(sapply(letters_tables, length)>4)]
letters_tables[which(sapply(letters_tables, length)>4)]
'%!in%' <- function(x,y)!('%in%'(x,y))
lapply(phages_with_spec_letters, function(x) {which(x%!in%c('A', 'C', 'G', 'T'))})
