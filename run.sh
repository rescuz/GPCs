sh $dir_to_GPCs/GPCs.sh \
        -v $dir_to_vcflist/tmp.vcf.list \
        -b $dir_to_bed/chr.bed.test \
        -i 1 \
        -n test11w \
        -d  $dir_to_result \
        -V $dir_to_vcf \
        -S $dir_to_singularity/singularity \
        -t 1 \
        -s 4 \
        -r 0.2 \
        -f $dir_to_phenotype_file/phenotype_file \
        -T $dir_to_GPCs \
        -p $dir_to_pheno/pheno.txt \
        -M 1 \
        -m 0
echo "Still_waters_run_deep" > run.sh.sign
