# Structural Variant Calling Benchmark

This benchmark is based on the publicly available PacBio CCS 15kb dataset of
the Ashkenazim son HG002/NA24385. We provide each step how to reproduce the
final metrics with publicly available tools.

Please let us know, if you find mistakes or want your tool added.

# Comparison

|Run|F1^ %|Precision %|Recall %|FP|FN|FP+FN|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|[**pbsv**](https://github.com/PacificBiosciences/pbsv)|96.33|94.68|98.04|531|189|720|
|[**svim**](https://github.com/eldariont/svim)|94.12|93.76|94.48|606|532|1138|
|[**sniffles**](https://github.com/fritzsedlazeck/Sniffles)|93.56|93.30|93.83|650|595|1245|

# Get tools

Information how to install `conda` and add the `bioconda` channel is available
on https://bioconda.github.io/.

```sh
conda install pbmm2==0.12.0 pbsv==2.2.* svim==0.4.3 sniffles==1.0.10
```

# Get data
1) Create directory structure:
```sh
mkdir -p fastqs ref alns alns_merged tools/{pbsv,sniffles} giab
```

2) Download genome in a bottle annotations:
```sh
FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed > giab/HG002_SVs_Tier1_v0.6.bed
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz > giab/HG002_SVs_Tier1_v0.6.vcf.gz
```

3) Make sure that you have exact truvari version (check that all python dependencies are met, not described):
```sh
git clone https://github.com/spiralgenetics/truvari
(cd truvari; git reset --hard 600b4ed7)
```

4) Download hg19 reference with decoys and map non-ACGT characters to N:
```sh
curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > ref/human_hs37d5.fasta.gz
gunzip ref/human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' human_hs37d5.fasta
```

5) Download hg19 tandem repeat annotations:
```sh
curl -s https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed > ref/human_hs37d5.trf.bed
```

6) Download all `.fastq` files:
```sh
curl -s ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/HG002.15kb.Q20.tar.gz > fastqs/HG002.15kb.Q20.tar.gz
tar xvzfC fastqs/HG002.15kb.Q20.tar.gz fastqs/
```

# Alignment

6) Index reference:
```sh
pbmm2 index ref/human_hs37d5.fasta ref/human_hs37d5.mmi --preset CCS
```

7) Align each movie, can be done in parallel:
```sh
for i in fastqs/*.fastq; do
    FILENAME="${i#*fastqs/}"
    FILEPREFIX="${FILENAME%.*}"
    pbmm2 align ref/human_hs37d5.mmi $i "alns/hg19.${FILEPREFIX}.bam" --preset CCS \
                --sort --rg '@RG\tID:${FILEPREFIX}' --sample HG2
done
```

# Run pbsv

8a) Discover SV signatures for each alignment, can be done in parallel:
```sh
for i in alns/*.bam; do
    FILENAME="${i#*alns/}"
    FILEPREFIX="${FILENAME%.*}"
    pbsv discover $i "tools/pbsv/hg19.${FILEPREFIX}.svsig.gz" --tandem-repeats ref/human_hs37d5.trf.bed
done
```

8b) Call and polish SVs:
```sh
pbsv call ref/human_hs37d5.fasta tools/pbsv/*.svsig.gz tools/pbsv/hg2.pbsv.vcf --ccs -t INS,DEL
bgzip tools/pbsv/hg2.pbsv.vcf
tabix tools/pbsv/hg2.pbsv.vcf.gz
```

# Run svim

9a) Merge and re-sort alignments:
```sh
samtools merge -b alns_merged/bam_list.fn alns_merged/hs37d5.all.bam
samtools sort -n alns_merged/hs37d5.all.bam > alns_merged/hs37d5.all.querysorted.bam
```

9b) Run svim:
```sh
svim alignment --cluster_max_distance 1.4 tools/svim alns_merged/hs37d5.all.querysorted.bam
```

9c) Prepare for truvari:
```sh
cat tools/svim/final_results.vcf \
    | sed 's/INS:NOVEL/INS/g' \
    | sed 's/DUP:INT/INS/g' \
    | sed 's/DUP:TANDEM/INS/g' \
    | awk '{ if($1 ~ /^#/) { print $0 } else { if(($5=="<DEL>" || $5=="<INS>") && $6>40) { print $0 } } }' \
    > tools/svim/hg2.svim.vcf
bgzip tools/svim/hg2.svim.vcf
tabix tools/svim/hg2.svim.vcf.gz
```

# Run sniffles

10a) Add MD tag to BAM:
```sh
samtools calmd -bS alns_merged/hs37d5.all.bam ref/human_hs37d5.fasta > alns_merged/hs37d5.all.md.bam
```

10b) Run sniffles:
```sh
sniffles -s 3 --skip_parameter_estimation -m alns_merged/hs37d5.all.md.bam -v tools/sniffles/hg2.sniffles.vcf --report_seq
```
10c) Prepare for truvari:
```sh
cat <(cat tools/sniffles/hg2.sniffles.vcf | grep "^#") \
    <(cat tools/sniffles/hg2.sniffles.vcf | grep -vE "^#" | \
      grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) \
    | bgzip -c > tools/sniffles/hg2.sniffles.vcf.gz
tabix tools/sniffles/hg2.sniffles.vcf.gz
```

# Final comparison

11) Compare to ground truth:
```sh
truvari/truvari.py -f ref/human_hs37d5.fasta -b giab/HG002_SVs_Tier1_v0.6.vcf.gz\
                   --includebed giab/HG002_SVs_Tier1_v0.6.bed -o bench-pbsv --passonly\
                   --giabreport -r 1000 -p 0.00 -c tools/pbsv/hg2.pbsv.vcf.gz
truvari/truvari.py -f ref/human_hs37d5.fasta -b giab/HG002_SVs_Tier1_v0.6.vcf.gz\
                   --includebed giab/HG002_SVs_Tier1_v0.6.bed -o bench-svim --passonly\
                   --giabreport -r 1000 -p 0.00 -c tools/svim/final_results.filtered.vcf.gz
truvari/truvari.py -f ref/human_hs37d5.fasta -b giab/HG002_SVs_Tier1_v0.6.vcf.gz\
                   --includebed giab/HG002_SVs_Tier1_v0.6.bed -o bench-sniffles --passonly\
                   --giabreport -r 1000 -p 0.00 -c tools/sniffles/hg2.sniffles.vcf.gz
```

12) Parse results:
```sh
function sumsv() { cat $1 | grep ':' | tr -d ',' |sed "s/^[ \t]*//"| tr -d '"' |\
                   tr -d ' ' | tr ':' '\t' | awk '{ for (i=1; i<=NF; i++)  {
                     a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++)
                     { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' |\
                   tail -n 1 | awk '{ printf "%1.4f\t%1.4f\t%1.4f\t%10.0f\t%10.0f\t%10.0f\n", $2,$4,$11,$1,$8,$1+$8 }';}
cat <(echo -e "Run\tF1\tPrecision\tRecall\tFP\tFN\tFP+FN")\
    <(cat <(for i in bench-*; do printf $i"\t";sumsv $i/summary.txt;done) |\
    sed 's/bench-//g;' | sort -k 2 -n -r) | column -t
```
Or for markdown:
```sh
cat <(echo -e "Run\tF1\tPrecision\tRecall\tFP\tFN\tFP+FN")\
    <(echo -e ":-:\t:-:\t:-:\t:-:\t:-:\t:-:\t:-:")\
    <(cat <(for i in bench-*; do printf $i"\t";sumsv $i/summary.txt;done) |\
    sed 's/bench-//g;' | sort -k 2 -n -r) | tr '\t' ' ' | tr -s ' ' | tr ' ' '|' | awk '{print "|"$0"|"}'
```
