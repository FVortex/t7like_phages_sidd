#all in one directory
rm(list = ls())
wd <- "/home/hp/t7like_phages_sidd/"
setwd(wd)

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
phages_with_phiols <- c('ba14', 'k11', '13a', 'phiYeO3',  't3', 't7', 'yepe2')
all_Autographivirinae_names[grepl(paste(phages_with_phiols, collapse="|"), all_Autographivirinae_names, ignore.case = T)]
#need to substitute T7, T3 and phiYeO3-12 genomes with used previously
##phages that DO contain phiols
library(ape)

t7gb<-read.GenBank(access.nb = 'NC_001604', as.character = T) #Group I: dsDNA viruses Order: 	Caudovirales Family: 	Podoviridae Subfamily: 	Autographivirinae Genus: T7likevirus
t7gb<-t7gb$NC_001604
#substitution
no <- which(grepl(paste('t7,', collapse="|"), all_Autographivirinae_names, ignore.case = T))
all_Autographivirinae_seqs[[no]] <- t7gb
names(all_Autographivirinae_seqs[no]) <- "NC_001604 Enterobacteria phage T7, complete genome"


t3gb<-read.GenBank(access.nb = 'NC_003298', as.character = T) #wikipedia Group: 	Group I (dsDNA) Order: 	CaudoviralesFamily: 	PodoviridaGenus: 	T7-like virusesSpecies: 	T3 phage
t3gb<-t3gb$NC_003298

no <- which(grepl(paste('t3,', collapse="|"), all_Autographivirinae_names, ignore.case = T))
all_Autographivirinae_seqs[[no]] <- t3gb
names(all_Autographivirinae_seqs[no]) <- "NC_003298 Enterobacteria phage T3, complete genome" 

phiYeO3_12gb<-read.GenBank(access.nb = 'NC_001271', as.character = T) # Based on its morphology, φYeO3-12 belongs to the family Podoviridae (25) and to type C in Bradley's classification (12); furthermore, it resembles a typical member of the T7 group (H.-W. Ackermann, personal communication) (1). Other Y. enterocolitica phages characterized to date by electron microscopy have been of type A in Bradley's classification (22) or have been classified into the families Myoviridae or Podoviridae (2). ------------ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC94659/
phiYeO3_12gb<-phiYeO3_12gb$NC_001271 


no <- which(grepl(paste('phiYeO3-12,', collapse="|"), all_Autographivirinae_names, ignore.case = T))
all_Autographivirinae_seqs[[no]] <- phiYeO3_12gb
names(all_Autographivirinae_seqs[no]) <- "NC_003298_Enterobacteria_phage_T3,_complete_genome" 

#and adding 4 more phages
phage13a_gb<-read.GenBank(access.nb = 'NC_011045.1', as.character = T) # Enterobacteria phage 13a [TAX:532076]Lineage	Viruses; dsDNA viruses,noRNAstage;Caudovirales;Podoviridae;Autographivirinae; T7virus; unclassified T7-like viruses
phage13a_gb<-phage13a_gb$NC_011045.1
all_Autographivirinae_seqs[49] <- list(phage13a_gb)
names(all_Autographivirinae_seqs[49]) <- "NC_011045.1_Enterobacteria_phage_13a,_complete_genome"

ba14gb<-read.GenBank(access.nb = 'NC_011040.1', as.character = T) # Enterobacteria phage BA14 [TAX:532074] 	Viruses; dsDNA viruses, no RNA stage; Caudovirales;Podoviridae; Autographivirinae; T7virus;unclassified T7-like viruses
ba14gb<-ba14gb$NC_011040.1

all_Autographivirinae_seqs[[50]] <- ba14gb
names(all_Autographivirinae_seqs[50]) <- "NC_011040.1_Enterobacteria_phage_BA14,_complete_genome"


yepe2gb<-read.GenBank(access.nb = 'NC_011038.1', as.character = T) # Yersinia phage Yepe2 [TAX:532078] 	Viruses; dsDNA viruses, no RNA stage; Caudovirales;Podoviridae; Autographivirinae; T7virus;unclassified T7-like viruses
yepe2gb<-yepe2gb$NC_011038.1

all_Autographivirinae_seqs[[51]] <- yepe2gb
names(all_Autographivirinae_seqs[51]) <- "NC_011038.1_Yersinia_phage_Yepe2,_complete_genome"

k11gb<-read.GenBank(access.nb = 'EU734173.1', as.character = T) #Klebsiella phage K11 [TAX:532077] 	Viruses; dsDNA viruses, no RNA stage; Caudovirales;Podoviridae; Autographivirinae; T7virus;unclassified T7-like viruses
k11gb<-k11gb$EU734173.1

all_Autographivirinae_seqs[[52]] <- k11gb
names(all_Autographivirinae_seqs[52]) <- "EU734173.1_Klebsiella_phage_K11,_complete_genome"


all_Autographivirinae_names <- c(all_Autographivirinae_names, "NC_011045.1_Enterobacteria_phage_13a,_complete_genome", "NC_011040.1_Enterobacteria_phage_BA14,_complete_genome","NC_011038.1_Yersinia_phage_Yepe2,_complete_genome", "EU734173.1_Klebsiella_phage_K11,_complete_genome")
names(all_Autographivirinae_seqs) <- all_Autographivirinae_names

####Cutting genomes into pieces

#fucntion to split
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

# replacement + adding of old genomes [c(1,2,5 ,49:52)]
for (i in seq_along(all_Autographivirinae_seqs)) {
  assign(paste0(all_Autographivirinae_names[[i]], '_by_chunks'), splitWithOverlap(unlist(all_Autographivirinae_seqs[[i]]), 10000, 2000))
}

phages_chunks<-sort(grep(pattern = '*_by_chunks$', ls(), value = T))[c(1,2,5 ,49:52)]
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
    print(paste('Processing', i, j))
    # #aa<-as.character(read.fasta(paste0(wd, '/t7_genome_parts_string/', i), as.string = T, set.attributes = F))
    # #print(substr(aa, nchar(aa)-1000-10, nchar(aa)-1000+1))
      system(paste0('cd ', wd, 'sist/
                    ',
             'perl -X master.pl -a M -f ', wd, 'phages_by_10kbp_chunks/', i, '_', j, ' -o ' , wd, 'phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'))
  }
}
#phages_chunks <- phages_chunks[-c(2,5)]
######file.rename(list.files(pattern="*1no_", paste0("water_", 1:700))
for (i in phages_chunks) {
  interm<-c()
  print(i)
  for (j in seq_along(get(i))) {
    print(j)
    if (j==1) {
      touch<-read.csv(paste0(wd, 'phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'), sep = '\t', skip = 1)
      touch<-touch[-(9001:10000),]
    }
    else if (j==length(get(i))){
      
      if(i=="k11gbchunks") { touch<-read.csv(paste0(wd, 'phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'), sep = '\t', skip = 1)  #skipping needed since warnong message 'sequence if too short' is present
      } else {
        touch<-read.csv(paste0(wd, 'phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'), sep = '\t', skip = 1)
        touch<-touch[-(1:1000),]
      }
    }
    else{
      touch<-read.csv(paste0(wd, 'phages_by_10kbp_sidd/Perl_sist_output_', i, 'no_', j, '.tsv'),sep = '\t', skip=1)
      touch<-touch[-c(1:1000, 9001:10000),]
    } #warnings omitting
    interm<-rbind(interm, as.matrix(touch))
    
  } 
  assign(paste0(i, '_sidd'), interm)
}

#plots

phages_sidds<-(grep(pattern = '*_by_chunks_sidd', ls(), value = T))


phiOLs_coords<-c(353+17, ##353..375#ba14
                 479+17, #479..501 #k11 
                 413+17, # 423..445 #13a 
                 
                 375+17, #nuccore says 375..397 fo phiYe03-12, TSS at 392 (?) 
                 NA, #SP6 lacks
                 366+17, #nuccore 366..388, TSS at 383 for t3
                # 405, #t7
                 644+17)   ##644..666) #TSS is explicit for t7#yepe2 
mean_phiOLs_coords <- mean(phiOLs_coords, na.rm = T)
range_phiOLs_coords <- range(phiOLs_coords, na.rm = T)

mean_phiORs_coords <- mean(sapply(grep('*gb$', ls(), value = T), function(x) {length(get(x))})-phiORs_coords)
range_phiORs_coords <- range(sapply(grep('*gb$', ls(), value = T), function(x) {length(get(x))})-phiORs_coords)



#ba14 353+17 ##353..375
#k11 479+17 #479..501
#13a 413+17 # 423..445

#yepe2 644+17   ##644..666

phiORs_coords<-c(39073+17,   ##39073..39095 #ba14 
                 
                 40503+17, #40503..40525 #k11  
                 38159+17, # 38159..38181 #13a
                 38797+17, #38797..38819 phi...
  #               NA, #sp6
                 37432+17, # 37432..37454 t3
                 39229, #t7
                 37717+17) #37717..37739)#yepe2

mean_phiOLs_coords <- mean(phiOLs_coords, na.rm = T)
range_phiOLs_coords <- range(phiOLs_coords, na.rm = T)

phages_with_phiols <- c('ba14', 'k11', '13a', 'phiYeO3',  't3', 't7', 'yepe')
names(all_Autographivirinae_seqs)[grepl(paste(phages_with_phiols, collapse="|"), names(all_Autographivirinae_seqs), ignore.case = T)]
#need to substitute T7, T3 and phiYeO3-12 genomes with used previously
#ba14 39073+17   ##39073..39095
#k11 40503+17 #40503..40525
#13a  38159+17 # 38159..38181

#yepe2 37717+17 #37717..37739

names(phiOLs_coords) <- paste0(phages, '_phiols_coord')
names(phiORs_coords) <- paste0(phages, '_phiors_coord')

svg('SIDD_for_48_phages.svg', height = 12, width = 12)
par(mar=c(1,1,1,1))
par(oma=c(1,1,1,1))
par(mfrow=c(10,5))
for (i in seq_along(phages_sidds)){
  #plot(get(phages_sidds[i])[,2], type='l', ylim=c(0,1), main=paste0('SIDD profile for complete ', ' ',  toupper(all_Autographivirinae_names[i]), ' ', ' DNA'), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  plot(get(phages_sidds[i])[,2], type='l', ylim=c(0,1), main = toupper(all_Autographivirinae_names[i]), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  abline(h=0.5, col='grey', lty=3)
  # #abline(v=c(phiOLs_coords[i], phiORs_coords[i]), col='red')
}
dev.off()
svg('SIDD_for_48_phages_left_flank.svg', height = 12, width = 12)
par(mar=c(1,1,1,1))
par(oma=c(1,1,1,1))
par(mfrow=c(10,5))
for (i in seq_along(phages_sidds)){
  #plot(get(phages_sidds[i])[,2], type='l', ylim=c(0,1), main=paste0('SIDD profile for complete ', ' ',  toupper(all_Autographivirinae_names[i]), ' ', ' DNA'), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  plot(get(phages_sidds[i])[1:1000,2], type='l', ylim=c(0,1), main = toupper(all_Autographivirinae_names[i]), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  abline(h=0.5, col='grey', lty=3)
  abline(v = mean_phiOLs_coords, lty = 3, col = 'red', lwd= 3)
  abline(v = range_phiOLs_coords, lty = 3, col = 'red', lwd= 1.5)
  # #abline(v=c(phiOLs_coords[i], phiORs_coords[i]), col='red')
}
dev.off()

svg('SIDD_for_48_phages_right_flank.svg', height = 12, width = 12)
par(mar=c(1,1,1,1))
par(oma=c(1,1,1,1))
par(mfrow=c(10,5))
for (i in seq_along(phages_sidds)){
  interv <- (nrow(get(phages_sidds[i]))-1000):nrow(get(phages_sidds[i]))
  #plot(get(phages_sidds[i])[,2], type='l', ylim=c(0,1), main=paste0('SIDD profile for complete ', ' ',  toupper(all_Autographivirinae_names[i]), ' ', ' DNA'), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  plot(get(phages_sidds[i])[interv,2], type='l', ylim=c(0,1), main = toupper(all_Autographivirinae_names[i]), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  abline(h=0.5, col='grey', lty=3)
  #abline(v = mean_phiOLs_coords, lty = 3, col = 'red', lwd= 3)
  #abline(v = range_phiOLs_coords, lty = 3, col = 'red', lwd= 1.5)
  # #abline(v=c(phiOLs_coords[i], phiORs_coords[i]), col='red')
}
dev.off()

svg('SIDD_for_48_phages_both_flanks.svg', height = 12, width = 12)
par(mar=c(1,1,1,1))
par(oma=c(1,1,1,1))
par(mfrow=c(10,10))
for (i in seq_along(phages_sidds)){
  #plot(get(phages_sidds[i])[,2], type='l', ylim=c(0,1), main=paste0('SIDD profile for complete ', ' ',  toupper(all_Autographivirinae_names[i]), ' ', ' DNA'), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  plot(get(phages_sidds[i])[1:1000,2], type='l', ylim=c(0,1), main = toupper(all_Autographivirinae_names[i]), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5, col = 'darkred')
  abline(h=0.5, col='grey', lty=3)
  abline(v = mean_phiOLs_coords, lty = 3, col = 'red', lwd= 3)
  abline(v = range_phiOLs_coords, lty = 3, col = 'red', lwd= 1.5)
  
  interv <- (nrow(get(phages_sidds[i]))-1000):nrow(get(phages_sidds[i]))
  #plot(get(phages_sidds[i])[,2], type='l', ylim=c(0,1), main=paste0('SIDD profile for complete ', ' ',  toupper(all_Autographivirinae_names[i]), ' ', ' DNA'), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5)
  plot(get(phages_sidds[i])[interv,2], type='l', ylim=c(0,1), main = toupper(all_Autographivirinae_names[i]), ylab='Opening probability', xlab='Sequence (nts)', lwd=1.5 , col = 'darkblue')
  abline(h=0.5, col='grey', lty=3)
  # #abline(v=c(phiOLs_coords[i], phiORs_coords[i]), col='red')
}
dev.off()




#letters check-up
letters_tables <- lapply(all_Autographivirinae_seqs, function (x) {table(unlist(x))})
phages_with_spec_letters <- all_Autographivirinae_seqs[which(sapply(letters_tables, length)>4)]
letters_tables[which(sapply(letters_tables, length)>4)]
'%!in%' <- function(x,y)!('%in%'(x,y))
lapply(phages_with_spec_letters, function(x) {which(x%!in%c('A', 'C', 'G', 'T'))})

all_Autographivirinae_hosts <- sub("(.*?) p.*", "\\1", all_Autographivirinae_names)
all_Autographivirinae_hosts <- sub(".* (.*?)", "\\1", all_Autographivirinae_hosts)
cols_hosts <- as.numeric(factor(all_Autographivirinae_hosts))
barplot((as.numeric(lapply(all_Autographivirinae_seqs, GC))), col = cols_hosts, ylim = c(-0.1, 0.65))
abline(h=0.5)
legend('bottom', ncol = 5, legend = unique(all_Autographivirinae_hosts), fill = cols_hosts, cex = 0.5)

#plots