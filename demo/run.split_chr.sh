dir=`pwd`
if [ -e $dir/qsu.sh  ] ; then rm -f $dir/qsu.sh ; fi
mkdir -p  {$dir/tmp,$dir/shell}
### the full file of variant #####
file=/path/to/bed_like_file
dir2=$dir/test
name=${file##*/}
for i in {1..22} X Y M 
do
	cat $file |sort -V |grep -w chr$i > $dir/tmp/chr$i.${name} && \
	if [ -s  $dir/tmp/chr$i.${name}   ] 
	then
		if [ ! -d $dir2/chr$i ] ; then  mkdir -p $dir2/chr$i ; fi  
		echo -n "sh /*/GPC/GPCs.sh \\
		-v /*/tmp.vcf.list \\
		-b $dir/tmp/chr$i.${name} \\
		-i 20 \\
		-n ChinaMAP.GPC.$name.chr$i \\
		-d  $dir2/chr$i \\
		-V /path/to/vcf \\
		-a /path/to/database \\
		-S /share/app/singularity-3.2.0/bin/singularity \\
		-t 10 \\
		-s 4 \\
		-r 0.2 \\
		-f /path/to/phenotype_file \\
		-T /path/to/result \\
		-p /path/to/pheno.txt \\
		-M 1 \\
		-m 1 
		echo "Still_waters_run_deep" > $dir/shell/ChinaMAP.GPC.$name.chr$i.sh.sign "  >$dir/shell/ChinaMAP.GPC.$name.chr$i.sh
		### change according to your cluster management system   ####
		echo "qsub -cwd -l vf=20g,p=20 -q mgi_ruijin.q -P MGI_RUIJIN -binding linear:20  $dir/shell/ChinaMAP.GPC.$name.chr$i.sh "  >>$dir/qsu.sh 
	fi
done
