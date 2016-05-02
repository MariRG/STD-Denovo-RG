# STD-Denovo-RG
#!/bin/bash

#clear console
clear
echo "Running de novo assembly pipeline on sequencer directory"

#set variables for program and reference locations
#TODO: define location of contaminant sequences to pass to bbmap. This is in-built in the commands, but can be cleaned up or passed

#QC stage
for i in {1..93}
do
  echo "Making subdirectories for assemblies"
  mkdir $i/assemblies $i/prodigal
  cd $i
  echo "Temporarily moving fastq gz's to directories for processing"
  cp *\_S$i\_L001_R2_001.fastq.gz R2.fastq.gz
  cp *\_S$i\_L001_R1_001.fastq.gz R1.fastq.gz
  echo "Unpacking fastq files for QC"
  gunzip *.gz
  echo "Quality trimming raw fastq reads"
  /software/bbmap/bbduk.sh in1=R1.fastq in2=R2.fastq out1=qtrim1.fq out2=qtrim2.fq qtrim=rl trimq=10 overwrite=true Xmx240g
  echo "Removing contaminant sequences"
  /software/bbmap/bbduk.sh in1=qtrim1.fq in2=qtrim2.fq out1=unmatched1.fq out2=unmatched2.fq outm1=matched1.fq outm2=matched2.fq ref=/databases/pweb.fasta k=31 hdist=1 stats=pweb_stats.txt overwrite=true Xmx240g
  /software/bbmap/bbduk.sh in1=unmatched1.fq in2=unmatched2.fq out1=cleaned1.fq out2=cleaned2.fq outm1=crypto_matched1.fq outm2=crypto_matched2.fq ref=/databases/reference_genomes/crypto.fasta k=31 hdist=1 stats=crypto_stats.txt overwrite=true Xmx240g
  /software/bbmap/bbduk.sh in1=cleaned1.fq in2=cleaned2.fq out1=cleaned3.fq out2=cleaned4.fq outm1=e_coli_matched1.fq outm2=e_coli_matched2.fq ref=/databases/reference_genomes/E_coli_K12.fasta k=31 hdist=1 stats=e_coli_stats.txt overwrite=true Xmx240g
  /software/bbmap/bbduk.sh in1=cleaned3.fq in2=cleaned4.fq out1=cleaned5.fq out2=cleaned6.fq outm1=adh_matched1.fq outm2=adh_matched2.fq ref=/databases/reference_genomes/adh.fasta k=31 hdist=1 stats=adh_stats.txt overwrite=true Xmx240g
  /software/bbmap/bbduk.sh in1=cleaned5.fq in2=cleaned6.fq out1=cleaned7.fq out2=cleaned8.fq ref=/databases/reference_genomes/phix.fasta k=31 hdist=1 stats=phix_stats.txt overwrite=true Xmx240g
  #echo "Removing duplicate sequences"
  #/software/bbmap/dedupe.sh in=
  echo "Merging overlapping sequences"
  /software/bbmap/bbmerge.sh in1=cleaned7.fq in2=cleaned8.fq out=merged.fq outu1=unmerged1.fq outu2=unmerged2.fq overwrite=true Xmx240g
  echo "Completed QC of files, removing temporary files"
  rm -rf cleaned*.fq qtrim*.fq adh*.fq crypto*.fq e_coli*.fq matched*.fq unmatched*.fq
  echo "Successfully removed temporary files"
  cd ..
  echo "Proceeding to sample $i"
done

#assembly stage
for i in {1..93}
do
  cd $i
  python /software/SPAdes-3.6.2-Linux/bin/spades.py -1 unmerged1.fq -2 unmerged2.fq -s merged.fq -t 40 --careful -o assemblies/
  echo "Completed assembly of sample $i"
  cd ..
done

#Removing bad assemblies from analysis
for i in {1..93}
do
  cd $i/assemblies
  python /github/snippets/fasta_trim.py -i contigs.fasta -o ../$i.contigs.long.fasta
  python /github/snippets/fasta_trim.py -i scaffolds.fasta -o ../$i.scaffolds.long.fasta
  cd ../../
  echo "Completed removing sequences less than 1kb from Sample $i"
done

#Renaming fasta headers to match experiment
for i in {106..109}
do
  cd $i/assemblies
  python /github/snippets/fasta_rename.py -i $i.contigs.long.fasta -o $i.contigs.renamed.fasta -n $i
  python /github/snippets/fasta_rename.py -i $i.scaffolds.long.fasta -o $i.scaffolds.renamed.fasta -n $i
  cd ../../
  echo "Completed renaming sequences according to mapping lookup table"
done



#ORF calling stage
for i in {1..4}
do
  cd combined$i/assemblies
  prodigal -a $i.prot.fasta -d $i.nuc.fasta -i $i.contigs.long.fasta -o $i.prodigal.out -p meta
  prodigal -a $i.scaff.prot.fasta -d $i.scaff.nuc.fasta -i $i.scaffolds.long.fasta -o $i.scaff.prodigal.out -p meta
  cd ../../
done
