wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/425/GCF_000005425.2_Build_4.0/GCF_000005425.2_Build_4.0_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/425/GCF_000005425.2_Build_4.0/GCF_000005425.2_Build_4.0_genomic.fna.gz
mkdir -p data/ref/ && mv *.gz data/ref/
echo -e "data.dir =/home/huangl/bin/bioinformatics/snpEFF/data\n#rice genome, version ref\nref.genome : ref" > snpEff.config
cd data/ref/
gunzip GCF_000005425.2_Build_4.0_genomic.gff.gz
gunzip GCF_000005425.2_Build_4.0_genomic.fna.gz
mv GCF_000005425.2_Build_4.0_genomic.fna sequences.fa
mv GCF_000005425.2_Build_4.0_genomic.gff genes.gff
cd ../../
java -jar $bin/snpEff.jar build -gff3 -v ref -c snpEff.config
java -jar $bin/snpEff.jar ann -c snpEff.config -v ref $vcf > anno.vcf
