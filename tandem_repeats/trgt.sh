
path=/path/to/$id.align.bam

#run trgt
~/trgt/trgt genotype --genome human_GRCh38_no_alt_analysis_set.fasta \
      --repeats ./adotto_hg38.bed\
      --reads $path \
      --output-prefix ~/trgt38/sczcall/$id -t 5

#reformat MC field
cd ~/trgt38/sczcall
for vcf in *.vcf.gz;do txt_file="${vcf%.vcf.gz}.txt"; bcftools query -f '%TRID\t[%MC]\n' "$vcf" > "$txt_file"; done
for file in *.txt;do  sed -i 's/,/\t/g' $file; awk '{

    delete rows
    max_rows=0
    
    for(col=2;col<=NF;col++){
        split($col, elements, "_")
        for(row=1;row<=length(elements);row++){
            rows[row][col]=elements[row]
            if(row>max_rows) max_rows=row
        }
    }
    
    for(row=1;row<=max_rows;row++){
        printf "%d %s ", row, $1  
        for(col=2;col<=NF;col++){
            printf "%s%s", rows[row][col], (col==NF?"\n":" ")
        }
    }
}' $file > tmp; mv tmp $file; echo $file; done #columns: motif id, TR id, hap1 TR count, hap2 TR count

