#!/bin/bash
set -e
set -o pipefail
echo `hostname`
echo ==========start at : `date "+%Y-%m-%d %H:%M:%S"` ==========

func(){
	echo "Usage:"
	echo "GPCs can be used to perform correlation analysis between genotype and phenotype data."
	echo -e "Version:beta1.0"
	echo -e "Update:2022-05-12\n"
	echo -e "Design:Lin Li  Email:lilin2021@sjtu.edu.cn"
	echo -e "Coding:Bing Zeng   Email:zengbing1@genomics.com"
	echo -e "Testing:Chenyang Yu, Xiaohui Sun et al\n"
	echo "######------#######-------######"
	echo "$0 [-v vcf list] [-b bed] [-i thread] [-n name] [-d dir] [-t depth] [-s reads_support]   [-r reads_ratio] [-p pheno] [-f pheno_file] [-S singularity] [-a anno_base] [-m remove_tmp] [-T script_dir] [-V vcf_dir] [-M mode] [-P Pvalue] [-H  genome_version]"
	echo -e  "Description:\nthe parameters labeled star means mandatory parameters."
	echo -e  "According to inupt vcf file list, variant information, qc thresold and phenotype file to output annotation file about population frequency, group information in gene level and mutation level, and boxplots with two-side t.test.\nCautions: 
1. The phenotype inputed into our program must be quantitative phenotype.\n"
	echo     "* -v vcf, the input vcf file list, eg., SNP\tchr1\tvcf file path, INDEL\tall chromosomes\tvcf file path."
	echo -e "\n"
	echo    "* -b bed, the vcf-like bed file of your all variant, eg., chr1\t1023\t.\tA\tT.\n"
	echo -e "\n"
	echo -e "  -i thread, thread information used for annotating vcf file by means of annovar. Default is 2.\n"
	echo -e "* -n name, the prefix name of output files.\n"
	echo -e "  -d dir, the directory of output files. Default is current directory.\n"
	echo -e "  -t depth, the minimum threshold value to support calling variant. Default is 10.\n"
	echo -e "  -s reads_support, the minimum threshold value of reads to support effect variant. Default is 4.\n"
	echo -e "  -r reads_ratio, the minimum threshold value to support effect variant. Default is 0.2.\n"
	echo -e "* -T script_dir, the absolute path of scripts.\n"
	echo -e "*  -p pheno, the phenotype names that used to execute statistic analysis.The colum name must be "phenoname".\n"
	echo    "* -f pheno_file, the phenotype file, eg., Sample\tTC\tTG\n1023\t10\t50."
	echo -e "\n"
	echo -e "* -S singularity, the absolute path of singularity. We use the default path, and you can change as you wish.\n"
	echo -e "* -a anno_base, the annotation databases of Annovar. We use the default path, and you can change as you wish.\n"
	echo -e "* -V vcf_dir, the absolute parent path of vcf files. \n"
	echo -e "  -m remove_tmp, remove tmp directory generated in analysis, eg -m 1(remove tmp directory), -m 0 (retain tmp directory). Removing tmp directory as default.\n"
	echo -e "  -P Pvalue, the pvalue threshold that used to keep figures. Default: retaining figures with pvalue no greater than  0.05/Nmut and 0.05/Ngene for mut_level and gene_level, respectively.\n"
	echo -e "  -M mode, group mode used in analysis, eg -M 1(common vcf with depths of variants, perform quality control), -M 2(vcf after inputation, without depths information, don't perform quality control).\n"
	echo -e "  -H genome_version, the genome version used in this program. If your choose hg38, then you should download  database of hg38, and all your input files should be hg38. Default is hg38, and hg19 version also supported."
	exit -1
}
thread=2
dir=`pwd`
depth=10
reads_support=4
reads_ratio=0.2
Nonsig_retain=0
singularity=/share/app/singularity-3.2.0/bin/singularity
anno_base=/zfssz5/BC_PS/zengbing1/INSTALL/workdir/AIVAR/annovar/humandb/hg38
remove_tmp=1
Pvalue=strict
genome_version=hg38
while getopts 'v:b:i:n:d:t:s:r:p:f:h:a:m:S:T:V:M:P:H:' OPT;do
	case $OPT in
		v) vcf="$OPTARG";;
		b) bed="$OPTARG";;
		i) thread="$OPTARG";;
		n) name="$OPTARG";;
		d) dir="$OPTARG";;
		t) depth="$OPTARG";;
		s) reads_support="$OPTARG";;
		r) reads_ratio="$OPTARG";;
		p) pheno="$OPTARG";;
		f) pheno_file="$OPTARG";;
		a) anno_base="$OPTARG";;
		m) remove_tmp="$OPTARG";;
		S) singularity="$OPTARG";;
		T) script_dir="$OPTARG";;
		V) vcf_dir="$OPTARG";;
		M) mode="$OPTARG";;
		P) Pvalue="$OPTARG";;
		H) genome_version="$OPTARG";;
		h) func;;
		?) func;;
	esac

done
### annotation by means of transvar ##### 如果不需要规范变异数据一致性，此步骤可不执行
#/home/zengbing1/.local/bin/transvar canno  -l $inputfile   --refseq --longestcoding --gseq  >$inputfile.out
#### select variants according to input bed file   #####
rm -fr $dir/gene_level &&\
rm -fr $dir/mut_level &&\
if [ ! -d $dir/tmp  ] ; then mkdir -p $dir/tmp ; fi   && \
if [ ! -d $dir/gene_level  ] ; then mkdir -p $dir/gene_level ; fi && \
if [ ! -d $dir/mut_level  ] ; then mkdir -p $dir/mut_level ; fi && \
#if [ ! -d $dir/gene_level  ] ; then mkdir -p $dir/gene_level/Significant/Strict; mkdir -p $dir/gene_level/Significant/Loose; mkdir -p $dir/gene_level/Non-significant   ; fi   && \
#if [ ! -d $dir/mut_level  ] ; then mkdir -p $dir/mut_level/Significant/Strict ; mkdir -p $dir/mut_level/Significant/Loose; mkdir -p $dir/mut_level/Non-significant ; fi   
cat $bed |awk 'BEGIN{OFS="\t"}{if($1 !~/^chr/) $1="chr"$1;print $1,$2}' >$dir/tmp/$name.input.region   && \
bed_path=${bed%/*}
pheno_path=${pheno%/*}
pheno_file_path=${pheno_file%/*}
cat $vcf |while read type chr  file 
do
	$singularity exec --cleanenv \
	-B $script_dir:$script_dir \
	-B $dir:$dir \
	-B $dir/tmp:$dir/tmp \
	-B $vcf_dir:$vcf_dir \
	-B $bed_path:$bed_path \
	$script_dir/bin/GPCs_v1.sif \
	bcftools view -R $dir/tmp/$name.input.region $file -O z -o $dir/tmp/Select.$chr.$type.vcf.gz 
done
cat $vcf |while read type chr  file
do
	$singularity exec --cleanenv  \
	-B $script_dir:$script_dir \
        -B $dir:$dir \
        -B $dir/tmp:$dir/tmp \
	 $script_dir/bin/GPCs_v1.sif \
	bcftools index $dir/tmp/Select.$chr.$type.vcf.gz && \
	echo "$file variants extract succeed !!!"
done
echo " All variants extract succeed !!!"
vcflist=($dir/tmp/Select.*.vcf.gz) 
$singularity exec --cleanenv \
-B $script_dir:$script_dir \
-B $dir:$dir \
-B $dir/tmp:$dir/tmp \
$script_dir/bin/GPCs_v1.sif \
bcftools concat ${vcflist[@]}  -O z -o $dir/merge.$name.select.vcf.gz && \
echo " Merged vcf file was achieved !!!" 
###annoation by means of annovar######
$singularity exec --cleanenv \
        -B $script_dir:$script_dir \
        -B $dir:$dir \
        -B $dir/tmp:$dir/tmp \
        -B $vcf_dir:$vcf_dir \
        -B $anno_base:$anno_base \
        -B $bed_path:$bed_path \
        $script_dir/bin/GPCs_v1.sif \
perl /$script_dir/bin/table_annovar.pl -buildver $genome_version $dir/merge.$name.select.vcf.gz $anno_base --thread $thread --outfile $dir/tmp/$name -protocol refGene,ChinaMAP,gnomad211_exome,gnomad211_genome,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,1000g2015aug_sas,exac03,bravo-dbsnp-all,dbnsfp31a_interpro,dbscsnv11,avsnp150,dbnsfp42a,clinvar_20190305,gwasCatalog -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r  -nastring .  -remove --onetranscript --vcfinput   && \
echo "Vcf file annotation succeed !!! "
#####qc, generate data for plotting  and generate script about plot###
if [ $mode == 1 ] ; then
$singularity exec --cleanenv \
        -B $script_dir:$script_dir \
        -B $dir:$dir \
        -B $dir/tmp:$dir/tmp \
        -B $vcf_dir:$vcf_dir \
        -B $anno_base:$anno_base \
        -B $bed_path:$bed_path \
        -B $pheno_path:$pheno_path \
        -B $pheno_file_path:$pheno_file_path \
        $script_dir/bin/GPCs_v1.sif \
perl  $script_dir/s1.get_dp_ingroup.pl -vcf $dir/merge.$name.select.vcf.gz  -bed  $bed  -info $dir/tmp/${name}.${genome_version}_multianno.txt -name $name   -dir $dir -depth $depth -reads_ratio $reads_ratio -pheno  $pheno -pheno_file $pheno_file -reads_support $reads_support   -Pvalue $Pvalue
elif [ $mode == 2  ] ; then
$singularity exec --cleanenv \
        -B $script_dir:$script_dir \
        -B $dir:$dir \
        -B $dir/tmp:$dir/tmp \
        -B $vcf_dir:$vcf_dir \
        -B $anno_base:$anno_base \
        -B $bed_path:$bed_path \
        -B $pheno_path:$pheno_path \
        -B $pheno_file_path:$pheno_file_path \
        $script_dir/bin/GPCs_v1.sif \
perl  $script_dir/s1.get_dp_ingroup.inputation.pl -vcf $dir/merge.$name.select.vcf.gz  -bed  $bed  -info $dir/tmp/${name}.${genome_version}_multianno.txt -name $name   -dir $dir -depth $depth -reads_ratio $reads_ratio -pheno  $pheno -pheno_file $pheno_file -reads_support $reads_support -Pvalue $Pvalue
fi && \
echo "Data  pre-processing succeed  !!!"
#### execute R script#####
$singularity exec --cleanenv \
        -B $script_dir:$script_dir \
        -B $dir:$dir \
        -B $dir/tmp:$dir/tmp \
        -B $vcf_dir:$vcf_dir \
        -B $anno_base:$anno_base \
        -B $bed_path:$bed_path \
        -B $pheno_path:$pheno_path \
        -B $pheno_file_path:$pheno_file_path \
        $script_dir/bin/GPCs_v1.sif \
Rscript $dir/gene_level/${name}.allgene_phenotype.plot.R   &
$singularity exec --cleanenv \
        -B $script_dir:$script_dir \
        -B $dir:$dir \
        -B $dir/tmp:$dir/tmp \
        -B $vcf_dir:$vcf_dir \
        -B $anno_base:$anno_base \
        -B $bed_path:$bed_path \
        -B $pheno_path:$pheno_path \
        -B $pheno_file_path:$pheno_file_path \
        $script_dir/bin/GPCs_v1.sif \
Rscript $dir/mut_level/${name}.allmut_phenotype.plot.R &
wait
echo "Plot succeed !!!"
if [ $remove_tmp == "1"  ] ; then rm -fr $dir/tmp ; else echo "tmp directory was retained !!!"  ; fi  && \
echo "$dir/tmp was removed !!!" && \
if [ $JOB_ID ]; then echo "jobid=$JOB_ID" && qstat -j $JOB_ID | grep usage ; fi && \
echo ==========end at : `date "+%Y-%m-%d %H:%M:%S"` ==========   && \
echo "GPCs program is done!!!"
