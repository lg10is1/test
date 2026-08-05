#ULTRA PER HAPLOTYPE
file_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ~/sczassembly)
base_name=$(basename "$file_path")
id="${base_name//.fasta/}"

mkdir -p repeatmasker/$id
cd repeatmasker/$id

~/ULTRA/ultra -o ultra.tsv $file_path --max_split -1 -t 5
tail -n +2 ultra.tsv | awk 'BEGIN{OFS="\t"}  $4 > 1 && $5 > 10 {print $1,$2,$3,$6,$4,$7}' > ultra.bed
rb liftover --bed ultra.bed ./$id.v0.9.align.paf -q > chm13.bed
#Rscript ~/repeatmasker/reformatbed.r
awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$5,$6,$8,$9,$13}' chm13.bed > int.bed
mv int.bed chm13.bed
sed -i 's/id:Z://g' chm13.bed
bedtools subtract -a chm13.bed -b ~/trgt38/reference/adotto_chm13_noflank.bed -A > int.bed
bedtools subtract -a int.bed -b ~/repeatmasker/censat.bed -A > int1.bed
bedtools subtract -a int1.bed -b ./chm13v2.0_SD.bed -A > int.bed
bedtools subtract -a int.bed -b ~/trgt38/reference/chm13_rm.bed -A > chm13_unannot.bed

#merge unannotated repeat regions across samples
while IFS= read -r dirname; do
    bedfile="${dirname}/chm13_unannot.bed"

    awk -v dir="$dirname" -v OFS="\t" '{
        seq = $8

        a = gsub(/A/, "&", seq);
        t = gsub(/T/, "&", seq);
        c = gsub(/C/, "&", seq);
        g = gsub(/G/, "&", seq);
        print $1,$2,$3,$8,$5,$6,$7,dir,a "_" t "_" c "_" g;
    }' "$bedfile" >> merged_chm13_unannot.bed
done < <(ls -d */)


awk '{print > "group_"$9".bed"}' merged_chm13_unannot.bed
echo "processing group file"
for group in group_*.bed; do
    bedtools merge -i <(
        awk '{
            len = length($4);
            val = ($7 - $6) / len;               
            rounded = sprintf("%.0f", val);
            print $0 "\t" rounded;
        }' "$group" | sort -k1,1 -k2,2n
    ) -c 4,8,10 -o distinct,distinct,distinct >> int.bed
done

rm group_*.bed

sort -k1,1 -k2,2n int.bed > sorted_collapsed_chm13_unannot.bed

echo "processing collpsed file"
bedtools cluster -i sorted_collapsed_chm13_unannot.bed > int.bed
awk '{
     first_char = substr($4, 1, 1);
     all_same = 1;
     for (i=2; i<=length($4); i++) {
         if (substr($4, i, 1) != first_char) {
             all_same = 0;
             break;
         }
     }
     if (!all_same) print;
}' int.bed > clustered_filtered_collapsed_TR.bed

#collapse each cluster into one trgt entry
REF_GENOME="./GCF_009914755.1_T2T-CHM13v2.0_genomic_chr.fasta"
INPUT="clustered_filtered_collapsed_TR.bed"
OUTPUT="chm13_CNTR.bed"
TRGT="chm13_CNTR_trgt.bed"

> "$OUTPUT"
> "$TRGT"
awk -F'\t' '{
    n = split($5, a, ",")
    print $0 "\t" n
}' "$INPUT" | sort -k8,8 -k9,9nr | awk -F'\t' '{
    if (prev_group != $8) {
        if (prev_group != "") {
            print prev_group, motif, occurance, chr, start, end, min, max
        }
        prev_group = $8
        motif = $4
        occurance = 0
        chr = $1
        start = $2 - 25
        min = 999999999
        max = 0
        count = 0
    }
    
    occurance += $9
    if ($9 > max_n) {
        max_n = $9
        motif = $4
    }
    end = $3 + 25
    
    split($6, reps, ",")
    for (i in reps) {
        val = reps[i] + 0
        if (val < min) min = val
        if (val > max) max = val
    }
}
END {
    if (prev_group != "") {
        print prev_group, motif, occurance, chr, start, end, min, max
    }
}' | while read -r group motif occurance chr start end min max; do
    leftseq=$(samtools faidx "$REF_GENOME" "${chr}:${start}-$((start+25))" | tail -n +2 | tr -d '\n' | tr '[:lower:]' '[:upper:]')
    rightseq=$(samtools faidx "$REF_GENOME" "${chr}:$((end-25))-${end}" | tail -n +2 | tr -d '\n' | tr '[:lower:]' '[:upper:]')
    
    echo -e "${chr}\t${start}\t${end}\t${motif}\t${occurance}\t${leftseq}\t${rightseq}\t${min}\t${max}" >> "$OUTPUT"
    echo -e "${chr}\t${start}\t${end}\tID=${chr}_${start}_${end};MOTIFS=${motif};STRUC=${leftseq}(${motif})n${rightseq}" >> "$TRGT"   
    echo "Processed group: $group"
done

#run trgt on candidate tr loci
~/trgt/trgt genotype --genome ./GCF_009914755.1_T2T-CHM13v2.0_genomic_chr.fasta \
      --repeats chm13_CNTR_trgt.bed \
      --reads $file \
      --output-prefix CNTR/$id -t 5

#mask and merge trgt result
out="MC_matrix_all_samples.tsv"

for v in *.vcf.gz; do
  s=${v%.vcf.gz}
  bcftools query -f '[%AP\t%MC\n]' "$v" |
  awk '{
    split($1,a,","); split($2,m,",");
    h1=(a[1]+0<0.5)?"NA":m[1];
    h2=(a[2]+0<0.5)?"NA":m[2];
    print h1"\t"h2
  }' > "${s}.tmp"
done
paste *.tmp > MC.txt

rm *.tmp

#print all sites that 1)non-NA>50%; 2) at least two different genotypes
awk -F'\t' '{
    na=0; delete val;
    for(i=1;i<=NF;i++){
        if($i=="NA"||$i=="") na++;
        else val[$i]=1;
    }
    if(na < (NF-1)/2 && length(val)>=2) print NR
}' MC.txt > CNTR_final_rowidx.txt

