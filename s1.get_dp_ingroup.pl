#!/user/bin/perl -w
use strict;
use Getopt::Long;
my($vcf,$indel,$bed,$info,$dir,$name,$pheno,$ratio,$depth,$pheno_file,$reads_support,$Pvalue);
GetOptions(
	"vcf:s"=>\$vcf,
	"bed:s"=>\$bed,
	"info:s"=>\$info,
	"dir:s" =>\$dir,
	"name:s"=>\$name,
	"depth:s"=>\$depth,
	"reads_ratio:s"=>\$ratio,
	"pheno:s"=>\$pheno,
	"pheno_file:s"=>\$pheno_file,
	"Pvalue:s"=>\$Pvalue,
	"reads_support:s"=>\$reads_support
);
die $! if ( ! defined $vcf || ! defined $bed  || ! defined $info || ! defined $dir || ! defined $name || ! defined $depth || ! defined $ratio || ! defined $reads_support || ! defined $Pvalue);

open IN,"gzip -dc $vcf|"  or die $!;
open BED,"$bed" or die $!;
open HEAD,"gzip -dc $vcf|"  or die "$vcf file is non-existed!";
my (%hash,%g_hash);
my %ge_double=();
my %mut_double=();
open PHE,"$pheno_file" or die $!;
my (%hash_pheno,$head_pheno,$gene_mut_set,$pheno_set,%pheno_array,@num_arr,$Rcol_type,$Ngene,$Nmut);
if ( $pheno ){
	open F,$pheno or die $!;
	while(<F>){
		chomp;
		next if (/phenoname/);
		$pheno_array{$_}++;
	}
	close F;
}
while(<PHE>){
	chomp;
	my @tmp_pheno=(split /\t/,$_);
	if (/Sample/){
		if($pheno){
			my (@tmp2_array,@head_arr,@tmp3_Rcol);
			for my $n (1..((scalar @tmp_pheno)-1)){
				if (exists $pheno_array{$tmp_pheno[$n]}){
					push @num_arr,$n;
					push @head_arr,$tmp_pheno[$n];
					push @tmp2_array,"\"$tmp_pheno[$n]\"";
					push @tmp3_Rcol,"$tmp_pheno[$n]=col_number()";
				}
			}
			my $tmp1=join "\t",@head_arr;
			$hash_pheno{$tmp_pheno[0]}=$tmp1;
			$pheno_set=join ",",@tmp2_array;
			$Rcol_type=join ",",@tmp3_Rcol;
		}		
		else{
			my $tmp1=join "\t",@tmp_pheno[1..((scalar @tmp_pheno)-1)];
			$hash_pheno{$tmp_pheno[0]}=$tmp1;			
			my (@tmp2_array,@tmp3_Rcol);
			for my $n2 (@tmp_pheno[1..((scalar @tmp_pheno)-1)]){
				push @tmp2_array,"\"$n2\"";
				push @tmp3_Rcol,"$n2=col_number()";
			}
			$pheno_set=join ",",@tmp2_array;
			$Rcol_type=join ",",@tmp3_Rcol;
		}
	}
	else{  
		if($pheno){
			my @tmp1_array;
			for my $n1 (@num_arr){
				push @tmp1_array,$tmp_pheno[$n1]; 
			}
			my $tmp2=join "\t",@tmp1_array;
			$hash_pheno{$tmp_pheno[0]}=$tmp2;
		}
		else {
			my $tmp2=join "\t",@tmp_pheno[1..((scalar @tmp_pheno)-1)];
                        $hash_pheno{$tmp_pheno[0]}=$tmp2;
		}
	}	
}
close PHE;
while (<BED>){
	chomp;
	my @tmp=(split/\t/,$_);
	if($tmp[0] !~/^chr/ ){
		$tmp[0]="chr" . $tmp[0];
	}
	my $newk=join "\t",($tmp[0],$tmp[1],".",$tmp[3],$tmp[4]);
        $hash{$newk}=$tmp[0] . "-" . $tmp[1] ."-". $tmp[3] ."-" . $tmp[4];
	$tmp[3] =~tr/actguACTGU/tgacaTGACA/; ##reverse DNA/RNA sequence
	$tmp[4] =~tr/actguACTGU/tgacaTGACA/;
	my $newk1=join "\t",($tmp[0],$tmp[1],".",$tmp[3],$tmp[4]);
        $hash{$newk1}=$tmp[0] . "-" . $tmp[1] ."-". $tmp[3] ."-" . $tmp[4];
}
close BED;
open INFO,"$info" or die $!;
my $infoID;
while (<INFO>){
	chomp;
	my @tmp1=(split/\t/,$_);
	if (/input|^Chr/){
		for my $num3 (0..($#tmp1)){
			if($tmp1[$num3] =~/^Otherinfo4$/ ){
				$infoID = $num3;
			}
		}			
	}
	else{
		$g_hash{"$tmp1[$infoID]-$tmp1[$infoID+1]-$tmp1[$infoID+3]-$tmp1[$infoID+4]"}=$tmp1[6];  ##match gene to chr postion
	}
}
close INFO;
$Nmut=scalar(keys %g_hash);
$Ngene=scalar(values %g_hash);
my @sample;
while(<HEAD>){
	chomp;
	next if (/^##/);
        my @arr=(split/\t/,$_);
	if (/^#/){
		@sample=@arr[9..$#arr];
	}
	last;
}
close HEAD;
my (%u_na,%u_pos,%u_neg,%sam_ge,%sam_mut,%umut_na,%umut_hetpos,%umut_hompos,%umut_neg,@tmp_set);
%u_na=%u_pos=%u_neg=%sam_ge=%sam_mut=%umut_na=%umut_hetpos=%umut_hompos=%umut_neg=();
while(<IN>){
		chomp;
		next if (/^#/);
		my @arr=(split/\t/,$_);
        	my $key=join "\t",($arr[0],$arr[1],".",$arr[3],$arr[4]);
		if (exists $hash{$key}){
			for my $i (9..$#arr){
				my $mutM=$g_hash{$hash{$key}}.":".$arr[0]."-" .$arr[1]. "-" .$arr[3]."-".$arr[4]; 
				$mut_double{$mutM}{$sample[$i-9]}{$key}=$arr[$i];
				$ge_double{$g_hash{$hash{$key}}}{$sample[$i-9]}{$key}=$arr[$i];
				$sam_ge{$sample[$i-9]}{$g_hash{$hash{$key}}}++;
				$sam_mut{$sample[$i-9]}{$mutM}++;
			}
		}
}
close IN;
for my $i (sort (keys %mut_double)){
	push @tmp_set, "\"$i\"";
}
$gene_mut_set =join ",",@tmp_set;
for my $gmut (keys %mut_double){
	for my $sam (keys %{$mut_double{$gmut}}  ){
		for my $chrmut(keys %{$mut_double{$gmut}{$sam}}){
			my ($tmp0,$tmp1,$tmp2,$tmp3)=(split/:|,/,$mut_double{$gmut}{$sam}{$chrmut})[0,1,2,3];
			if (($tmp2+$tmp1) <$depth  || $tmp0 =~/\./  ){
				$umut_na{$gmut}{$sam}{$chrmut}="NA";
			}
			elsif(($tmp2+$tmp1)>= $depth && $tmp0 =~/0\/1|1\/0/ && ($tmp2/($tmp2+$tmp1) >=$ratio  &&  $tmp2 >=$reads_support ) ){
				 $umut_hetpos{$gmut}{$sam}{$chrmut}="Het_Mutant";	
			}
			elsif(($tmp2+$tmp1)>= $depth && $tmp0 =~/1\/1/ && ($tmp2/($tmp2+$tmp1) >=$ratio  &&  $tmp2 >=$reads_support ) ){
				 $umut_hompos{$gmut}{$sam}{$chrmut}="Hom_Mutant";	
			}
			elsif(($tmp2+$tmp1)>= $depth && ($tmp0 =~/0\/0/ || ($tmp2/($tmp2+$tmp1) <$ratio  ||  $tmp2 <$reads_support )) ){
                                 $umut_neg{$gmut}{$sam}{$chrmut}="Wild_type";
                        }
			else {  print "NA\t$sam\t$chrmut\t$mut_double{$gmut}{$sam}{$chrmut}\n";
				$umut_na{$gmut}{$sam}{$chrmut}="NA";
			}
		} 
	}
}
for my $gene (keys %ge_double){
	for my $sam (keys %{$ge_double{$gene}}){
		my $flag=0;
		my @mutset;
		for my $chrmut ( keys %{$ge_double{$gene}{$sam}}   ){
			my ($tmp0,$tmp1,$tmp2,$tmp3)=(split/:|,/,$ge_double{$gene}{$sam}{$chrmut})[0,1,2,3];		
			if(   $tmp2+$tmp1 <$depth || $tmp2<$reads_support  ||  $tmp2/($tmp2+$tmp1) <$ratio ||   $tmp0 !~/[01]\/[01]/  ) {
				$flag++;
				#print "$gene\t$sam\t$tmp2\n";
			}
			elsif ( $flag ==0 && $tmp0 =~/\d\/1|1\/0/ && ($tmp2/($tmp2+$tmp1) >=$ratio &&  $tmp2 >=$reads_support)){
				$chrmut =~ s/\.\t//g;
				$chrmut =~ s/\t/-/g;
				$u_pos{$gene}{$sam}=$u_pos{$gene}{$sam}.";"."$chrmut";
			}		
			elsif ( $flag ==0 &&    ( $tmp0 =~/0\/0/ || ($tmp2/($tmp2+$tmp1) <$ratio && ($tmp2+$tmp1) >=$depth )  || (($tmp2+$tmp1) >=$depth && $tmp2 <$reads_support ) )){
				$u_neg{$gene}{$sam}++;
			}

		}
	}
}
open OUT3, ">$dir/gene_level/$name.category";
print OUT3 "Sample";
open OUT4,">$dir/mut_level/$name.mut.category";
print OUT4 "Sample";
open MERGE2,">$dir/mut_level/merge_pheno.$name.mut.txt";
print MERGE2 "Sample";
open MERGE,">$dir/gene_level/merge_pheno.$name.txt" or die $!;
print MERGE "Sample";
for my $gene (keys %ge_double){
	print OUT3 "\t$gene";
	print MERGE "\t$gene";
}
for my $gmut (keys %mut_double){
	print OUT4 "\t$gmut";
	print MERGE2 "\t$gmut";
}
print OUT3 "\n";
print OUT4 "\n";
print MERGE qq(\t$hash_pheno{"Sample"}\n);
print MERGE2 qq(\t$hash_pheno{"Sample"}\n);
for my $sam (keys %sam_ge){
	print OUT3 "$sam";
	print MERGE "$sam";
	for my $gene (keys %ge_double){
		if (exists $u_pos{$gene}{$sam}){
			print OUT3 "\tMutant" .  "$u_pos{$gene}{$sam}";
			print MERGE "\tMutant";
		}
		elsif (  exists $u_neg{$gene}{$sam}){
			if ($u_neg{$gene}{$sam} == scalar(keys %{$ge_double{$gene}{$sam}} ) ){
				print OUT3 "\tWild_type";
				print MERGE "\tWild_type";
				}
		}
		else{
			print OUT3 "\tNA";
			print MERGE "\tNA";
		}
	}
	print OUT3 "\n";
	print MERGE "\t$hash_pheno{$sam}\n";
}
close OUT3;
close MERGE;
for my $sam (keys %sam_mut){
	print OUT4 "$sam";
	print MERGE2 "$sam";
	for my $gmut (keys %mut_double  ){
		if(exists $umut_hompos{$gmut}{$sam}){
			print OUT4 "\tHom_Mutant";
			print MERGE2 "\tHom_Mutant";
		}		
		elsif(exists $umut_hetpos{$gmut}{$sam}){
			print OUT4 "\tHet_Mutant";
			print MERGE2 "\tHet_Mutant";
		}		
		elsif(exists $umut_neg{$gmut}{$sam}){
			print OUT4 "\tWild_type";
			print MERGE2 "\tWild_type";
		}
		elsif(exists $umut_na{$gmut}{$sam}){
			print OUT4 "\tNA";
			print MERGE2 "\tNA";
		}
		
	}
	print OUT4 "\n";
	print MERGE2 "\t$hash_pheno{$sam}\n";	
}
close OUT4;
close MERGE2;
`cat $info|cut -f7 >$dir/tmp/${name}.multianno.gene`;
open R,">$dir/gene_level/$name.allgene_phenotype.plot.R" or die $!;
print R <<CODE;
library(ggpubr)
library(tidyverse)
options(warn=-1)
setwd("$dir/gene_level")
tmp_info_file<-read_table2("$dir/tmp/${name}.multianno.gene")
vec_out2<-vector()
for ( i in unique(tmp_info_file\$Gene.refGene)  ){
	for (j in c($pheno_set)){ ### Column name of pheno must be "phenoname"
		data1 <-read_table2("$dir/gene_level/merge_pheno.$name.txt",na=" ",col_types=cols($Rcol_type))
		data1<-data1[,c(i,j)]
		colnames(data1)=c("gene","phe")
		vec_out<-vector()
		for (k in c("Mutant","Wild_type","NA")){
			if (k %in% sort(unique(data1\$gene )) ){
				rawsample <- length(data1[which(data1[1] == k),]\$gene)
				da_pos <-data1[which(data1[1] == k),][2] %>% filter_all(all_vars( . != "NA"))
				colnames(da_pos)<-c("group")
				if(length(da_pos\$group)>0){
					five<-fivenum(as.numeric(da_pos\$group))
					vec_out<-c(vec_out,c(paste(i,k,sep="_"),rawsample,nrow(da_pos),j,mean(as.numeric(da_pos\$group)),five[1],five[2],five[3],five[4],five[5]))
				}
				else{ vec_out<-c(vec_out,c(paste(i,k,sep="_"),rawsample,"0",j,"NA","NA","NA","NA","NA","NA"))    }
			}
			else{vec_out<-c(vec_out,c(paste(i,k,sep="_"),"0","0",j,"NA","NA","NA","NA","NA","NA"))}
		}
		for (l in (1:10)){
			vec_out2<-c(vec_out2,paste(vec_out[l],vec_out[l+10],vec_out[l+2*10],sep=";"))
		}
		data1<-data1[which(data1[1] != "NA"),]  %>% filter_all(all_vars(. != "NA")) 
  		my_comparisons <- list(c("Mutant", "Wild_type"))###组别
  		colnames(data1)<-c("group","value")
  		if( length(data1[which(data1\$group =="Mutant"),]\$value) >=2 &&  length(data1[which(data1\$group =="Wild_type"),]\$value) >=2 ){
  			means<-aggregate(value ~ group,data1,mean)
  			means\$value<-round(means\$value,3)
  			outcome<-t.test(data1[which(data1\$group =="Mutant"),]\$value,data1[which(data1\$group =="Wild_type"),]\$value,alternative="two.sided")
  			vec_out2<-c(vec_out2,outcome\$p.value)
  			p <- ggboxplot(data1, x="group", y="value",color = "group",palette = c("#00AFBB", "#E7B800"),shape="group")+xlab(i)+ylab(j)
  			q<-p+stat_compare_means(comparisons = my_comparisons,method="t.test",alternative = "two.sided")+stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE)+geom_text(data = means, aes(label = paste("mean= ",value,sep=""), y = value + 0.08))
			if ("$Pvalue" =="strict"  ){
				if( outcome\$p.value < 0.05/$Ngene  && length(outcome\$p.value)>0  ){
					dir.create("$dir/gene_level/Significant/Strict",recursive=TRUE)
					ggsave(paste("$dir/gene_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
				}
			}
			else if( "$Pvalue" != "strict" && "$Pvalue" >=0.05 && length(outcome\$p.value)>0 ){
				if(outcome\$p.value >= 0.05  && "$Pvalue" >= outcome\$p.value  ){
					dir.create("$dir/gene_level/Non-significant",recursive=TRUE)
					ggsave(paste("$dir/gene_level/Non-significant/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
				}
				else if( outcome\$p.value >=0.05/$Ngene &&   outcome\$p.value< 0.05 ){
					dir.create("$dir/gene_level/Significant/Loose",recursive=TRUE)
					ggsave(paste("$dir/gene_level/Significant/Loose/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
				}		
				else if (outcome\$p.value < 0.05/$Ngene){
					dir.create("$dir/gene_level/Significant/Strict",recursive=TRUE)
					ggsave(paste("$dir/gene_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
				}
			}
			else if(  "$Pvalue" != "strict" && "$Pvalue" <0.05  && "$Pvalue" >=0.05/$Ngene && length(outcome\$p.value)>0  ){
				if (outcome\$p.value >= 0.05/$Ngene  && "$Pvalue" >= outcome\$p.value){
					dir.create("$dir/gene_level/Significant/Loose",recursive=TRUE)
					ggsave(paste("$dir/gene_level/Significant/Loose/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
				}
				else if  (outcome\$p.value < 0.05/$Ngene){
					dir.create("$dir/gene_level/Significant/Strict",recursive=TRUE)
                                        ggsave(paste("$dir/gene_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
			}
			else if ( "$Pvalue" != "strict" &&  "$Pvalue" <0.05/$Ngene && length(outcome\$p.value)>0 ){
				if( outcome\$p.value < "$Pvalue" ){
				dir.create("$dir/gene_level/Significant/Strict",recursive=TRUE)
				ggsave(paste("$dir/gene_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
				}

			}
  		}
  		else{
			vec_out2<-c(vec_out2,"NA")
        		print(paste("the t-test of ",i,"_",j," group couldn't perform because of the observation was not enough  !!!",sep=""))
  		}
	}
}
result_out<-matrix(vec_out2,ncol=11,byrow=T)
result_col<-c("gene_group","raw_sample","sample_num","pheno","mean","minimum", "lower-hinge", "median", "upper-hinge", "maximum","P value(two-sided t.test between Mutant and Wild type)")
colnames(result_out)<-result_col
write.table(result_out,paste("$name","all_gene_phenotype",sep="."),quote=F,sep="\\t",row.names=F,col.names=T)
CODE
close R;
open R1,">$dir/mut_level/$name.allmut_phenotype.plot.R" or die $!;
print R1 <<CODE;
library(ggpubr)
library(tidyverse)
options(warn=-1)
setwd("$dir/mut_level")
vec_out2<-vector()
for ( i in c($gene_mut_set)  ){
	for (j in c($pheno_set)){ ### Column name of pheno must be "phenoname"
		data1 <-read_table2("$dir/mut_level/merge_pheno.$name.mut.txt",na=" ",col_types=cols($Rcol_type))
		data1<-data1[,c(i,j)]
		colnames(data1)=c("gmut","phe")
		vec_out<-vector()
		for (k in c("Het_Mutant","Hom_Mutant","NA","Wild_type")){
			if (k %in% sort(unique(data1\$gmut )) ){
				rawsample <- length(data1[which(data1[1] == k),]\$gmut)
				da_pos <-data1[which(data1[1] == k),][2] %>% filter_all(all_vars(. != "NA"))
				colnames(da_pos)<-c("group")
				if(length(da_pos\$group)>0){
					five<-fivenum(as.numeric(da_pos\$group))
					vec_out<-c(vec_out,c(paste(i,k,sep="_"),rawsample,nrow(da_pos),j,mean(as.numeric(da_pos\$group)),five[1],five[2],five[3],five[4],five[5]))
}
				else{vec_out<-c(vec_out,c(paste(i,k,sep="_"),rawsample,"0",j,"NA","NA","NA","NA","NA","NA"))}
			}
			else{ 
				vec_out<-c(vec_out,c(paste(i,k,sep="_"),"0","0",j,"NA","NA","NA","NA","NA","NA"))
			}
		}
		for (l in (1:10)){
			vec_out2<-c(vec_out2,paste(vec_out[l],vec_out[l+10],vec_out[l+2*10],vec_out[l+3*10],sep=";"))
		}
		data1<-data1[which(data1[1] != "NA"),]  %>% filter_all(all_vars(. != "NA")) 
  		colnames(data1)<-c("group","value")
  		if (length(unique(data1\$group)) ==3 &&  length(data1[which(data1\$group =="Hom_Mutant"),]\$value)>=2 &&  length(data1[which(data1\$group =="Wild_type"),]\$value)>=2 &&  length(data1[which(data1\$group =="Het_Mutant"),]\$value)>=2  ){
			outcome1<-t.test(data1[which(data1\$group =="Hom_Mutant"),]\$value,data1[which(data1\$group =="Wild_type"),]\$value,alternative="two.sided")
			outcome2<-t.test(data1[which(data1\$group =="Het_Mutant"),]\$value,data1[which(data1\$group =="Wild_type"),]\$value,alternative="two.sided")
			outcome3<-t.test(data1[which(data1\$group =="Hom_Mutant"),]\$value,data1[which(data1\$group =="Het_Mutant"),]\$value,alternative="two.sided")
			vec_out2<-c(vec_out2,paste(outcome1\$p.value,outcome2\$p.value,outcome3\$p.value,sep="/"))
			color=c("#00AFBB", "#E7B800","red")
			my_comparisons <- list(c("Hom_Mutant","Wild_type"),c("Het_Mutant","Wild_type"),c("Hom_Mutant","Het_Mutant"))
	 		means<-aggregate(value ~ group,data1,mean)
        		means\$value<-round(means\$value,3)
        		p <- ggboxplot(data1, x="group", y="value",color = "group",palette = color,shape="group")+xlab(i)+ylab(j)
        		q<-p+stat_compare_means(comparisons = my_comparisons,method="t.test",alternative = "two.sided")+stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE)+geom_text(data = means, aes(label = paste("mean= ",value,sep=""), y = value + 0.08))
        		pvector<-c(outcome1\$p.value,outcome2\$p.value,outcome3\$p.value)
			if("$Pvalue" =="strict"){
				if(any(pvector<0.05/$Nmut)){
					dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
					ggsave(paste("$dir/mut_level/Significant/Strict/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
				}

			}
			else if( "$Pvalue" !="strict" && "$Pvalue" >=0.05  ){
				if(any(pvector <"$Pvalue")  && all(pvector >=0.05)){
					dir.create("$dir/mut_level/Non-significant",recursive=TRUE)
					ggsave(paste("$dir/mut_level/Non-significant/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
				}
				else if ( any(pvector<0.05) &&  all(pvector>=0.05/$Nmut)  ){
					dir.create("$dir/mut_level/Significant/Loose",recursive=TRUE)
					ggsave(paste("$dir/mut_level/Significant/Loose/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
				}
				else if ( any(pvector<0.05/$Nmut) ){
					dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
					ggsave(paste("$dir/mut_level/Significant/Strict/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
				}
			}
			else if ("$Pvalue" !="strict" && "$Pvalue" <0.05  && "$Pvalue">=0.05/$Nmut){
				if( any(pvector <"$Pvalue") && all(pvector>=0.05/$Nmut )){
					dir.create("$dir/mut_level/Significant/Loose",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Loose/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
				}
				else if ( any(pvector<0.05/$Nmut) ){
					dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Strict/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
				}
			}
			else if ("$Pvalue" !="strict" && "$Pvalue" <0.05/$Nmut  ){
				if ( any(pvector<"$Pvalue") ){
                                        dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Strict/${name}",gsub(":","-",i),j,"pdf",sep="."),q,width=5,height=5)
                                }
			}
  		}
  		else if(length(unique(data1\$group)) ==2){
			color=c("#00AFBB", "#E7B800")
			if(all(sort(unique(data1\$group)) == c("Hom_Mutant","Wild_type")) &&  length(data1[which(data1\$group =="Hom_Mutant"),]\$value)>=2 &&  length(data1[which(data1\$group =="Wild_type"),]\$value)>=2 ){
				outcome11<-t.test(data1[which(data1\$group =="Hom_Mutant"),]\$value,data1[which(data1\$group =="Wild_type"),]\$value,alternative="two.sided")
				vec_out2<-c(vec_out2,paste(outcome11\$p.value,"NA","NA",sep="/"))
				my_comparisons <- list(c("Hom_Mutant","Wild_type"))
			}
			else if(all(sort(unique(data1\$group)) == c("Het_Mutant","Wild_type")) && length(data1[which(data1\$group =="Het_Mutant"),]\$value)>=2 &&  length(data1[which(data1\$group =="Wild_type"),]\$value)>=2   ){
				outcome11<-t.test(data1[which(data1\$group =="Het_Mutant"),]\$value,data1[which(data1\$group =="Wild_type"),]\$value,alternative="two.sided")
				vec_out2<-c(vec_out2,paste("NA",outcome11\$p.value,"NA",sep="/"))
				my_comparisons <- list(c("Het_Mutant","Wild_type"))
			}
			else if(all(sort(unique(data1\$group)) == c("Hom_Mutant","Het_Mutant")) && length(data1[which(data1\$group =="Hom_Mutant"),]\$value)>=2 &&  length(data1[which(data1\$group =="Het_Mutant"),]\$value)>=2   ){
				outcome11<-t.test(data1[which(data1\$group =="Hom_Mutant"),]\$value,data1[which(data1\$group =="Het_Mutant"),]\$value,alternative="two.sided")
				vec_out2<-c(vec_out2,paste("NA","NA",outcome11\$p.value,sep="/"))
				my_comparisons <- list(c("Hom_Mutant","Het_Mutant"))
			}
			else{
				vec_out2<-c(vec_out2,paste("NA","NA","NA",sep="/"))
				print(paste("the t-test of ",i,"_",j," group couldn't perform because of the observation was not enough  !!!",sep=""))
				next
			}
			means<-aggregate(value ~ group,data1,mean)
       			means\$value<-round(means\$value,3)
        		p <- ggboxplot(data1, x="group", y="value",color = "group",palette = color,shape="group")+xlab(i)+ylab(j)
        		q<-p+stat_compare_means(comparisons = my_comparisons,method="t.test",alternative = "two.sided")+stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE)+geom_text(data = means, aes(label = paste("mean= ",value,sep=""), y = value + 0.08))
			if ("$Pvalue" =="strict"  ){
                                if( outcome1\$p.value < 0.05/$Nmut  && length(outcome1\$p.value)>0  ){
                                        dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
                        }
                        else if( "$Pvalue" != "strict" && "$Pvalue" >=0.05 && length(outcome1\$p.value)>0 ){
                                if(outcome1\$p.value >= 0.05  && "$Pvalue" >= outcome1\$p.value  ){
                                        dir.create("$dir/mut_level/Non-significant",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Non-significant/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
                                else if( outcome1\$p.value >=0.05/$Nmut &&   outcome1\$p.value< 0.05 ){
                                        dir.create("$dir/mut_level/Significant/Loose",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Loose/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
                                else if (outcome1\$p.value < 0.05/$Nmut){
                                        dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
                        }
                        else if(  "$Pvalue" != "strict" && "$Pvalue" <0.05  && "$Pvalue" >=0.05/$Nmut && length(outcome1\$p.value)>0  ){
                                if (outcome1\$p.value >= 0.05/$Nmut  && "$Pvalue" >= outcome1\$p.value){
                                        dir.create("$dir/mut_level/Significant/Loose",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Loose/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
                                else if  (outcome1\$p.value < 0.05/$Nmut){
                                        dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
                                        ggsave(paste("$dir/mut_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }
                        }
                        else if ( "$Pvalue" != "strict" &&  "$Pvalue" <0.05/$Nmut && length(outcome1\$p.value)>0 ){
                                if( outcome1\$p.value < "$Pvalue" ){
                                dir.create("$dir/mut_level/Significant/Strict",recursive=TRUE)
                                ggsave(paste("$dir/mut_level/Significant/Strict/${name}",i,j,"pdf",sep="."),q,width=4,height=5)
                                }

                        }
  		}
  		else {  vec_out2<-c(vec_out2,paste("NA","NA","NA",sep="/"))
 			print(paste("the t-test of ",i,"_",j," group couldn't perform because of the observation was not enough  !!!",sep=""))
  		}
	}
}
result_out<-matrix(vec_out2,ncol=11,byrow=T)
result_col<-c("mutation_group","raw_sample","sample_num","pheno","mean","minimum", "lower-hinge", "median", "upper-hinge", "maximum","P value of two-sided t.test(Hom_Mutant vs Wild_type/Het_Mutant vs Wild_type/Hom_Mutant vs Het_Mutant)")
colnames(result_out)<-result_col
write.table(result_out,paste("$name","all_mut_phenotype",sep="."),quote=F,sep="\\t",row.names=F,col.names=T)
CODE
close R1
