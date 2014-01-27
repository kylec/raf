# create allele counts file for each patient
for a in `ls *.vcf | cut -d- -f1-3 | sort -u`; do 
    echo $a; python processVcf.py -i "$a*full*.vcf" -o $a.txt; 
done