dir=`pwd`
if [ -e $dir/qsu.sh  ] ; then rm -f $dir/qsu.sh ; fi
mkdir -p  {$dir/tmp,$dir/shell}
### the full file of variant #####
file=/path/to/bed_like_file
dir2=$dir/test
name=${file##*/}
### $n means the size to split phenotype ####
n=4
lines=$((`cat  /path/to/pheno.txt |wc -l` -1))
x=$((lines%n))
if [ $x == 0 ] ; then index=$((lines/n)) ;  elif [ $x -gt 0 ]  ; then  index=$((lines/n+1))  ; fi
index=$((index-1))
for ((i=0;i<=${index};i++))
do
        i1=$((n*i+2))
        i2=$((n*i+n-1+2))
	 if [ ! -d $dir2/tmp/  ] ; then mkdir -p $dir2/tmp/ ;fi
        sed -ne "1p" -e "${i1},${i2}p" /path/to/pheno.txt  >$dir2/tmp/$name.pheno.$i
	if [ ! -d $dir2/pheno.$i  ] ; then mkdir -p $dir2/pheno.$i ;fi
	echo "sh /*/GPC/GPCs.sh \\
                -v /path/to/tmp.vcf.list \\  
                -b $file \\
                -i 20 \\
                -n ChinaMAP.GPC.$name.pheno.$i \\
                -d $dir2/pheno.$i  \\
                -V /path/to/database/ \\
                -S /share/app/singularity-3.2.0/bin/singularity \\
                -t 10 \\
                -s 4 \\
                -r 0.2 \\
                -f /path/to/phenotype_file \\
		-T /path/to/result \\
                -p $dir2/test/tmp/$name.pheno.$i \\
                -M 1 \\
                -m 1 \\
                -N 1
                echo "Still_waters_run_deep" > $dir/shell/ChinaMAP.GPC.$name.pheno.$i.sh.sign "  >$dir/shell/ChinaMAP.GPC.$name.pheno.$i.sh
                ### change according to your cluster management system   ####
                echo "qsub -cwd -l vf=20g,p=20 -q mgi_ruijin.q -P MGI_RUIJIN -binding linear:20  $dir/shell/ChinaMAP.GPC.$name.pheno.$i.sh "  >>$dir/qsu.sh
done

