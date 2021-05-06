import vcf
from sys import argv 

# cf https://github.com/jdoughertyii/PyVCF for vcf line objetc attributs

tumour_vcf = vcf.VCFReader(open(argv[1])) # /home/jrudewicz/Documents/DREAM/data/Tumour2/Tumour2.mutect.vcf
tumour_battenberg = open(argv[2]).readlines() # /home/jrudewicz/Documents/DREAM/data/Tumour2/Tumour2.battenberg.txt
meta_table = open(argv[3],"w") # /home/jrudewicz/Documents/DREAM/data/Tumour2/meta_table_T2.tsv


head = "chr\tpos\tid\tref\talt\tqual\tdb\tsomatic\tvariant_type\tnorm_genotype\tnorm_ref_allele_depth\tnorm_alt_allele_depth\tnorm_base_quality\tnorm_depth\tnorm_allele_fraction\tnorm_variant_statut\ttum_genotype\ttum_ref_allele_depth\ttum_alt_allele_depth\ttum_base_quality\ttum_depth\ttum_allele_fraction\ttum_variant_statut\t"
head += "startpos\tendpos\tBAF\tpval\tLogR\tntot\tnMaj1_A\tnMin1_A\tfrac1_A\tnMaj2_A\tnMin2_A\tfrac2_A\tSDfrac_A\tSDfrac_A_BS\tfrac1_A_0.025\tfrac1_A_0.975\tnMaj1_B\tnMin1_B\tfrac1_B\tnMaj2_B\tnMin2_B\tfrac2_B\tSDfrac_B\tSDfrac_B_BS\tfrac1_B_0.025\tfrac1_B_0.975\tnMaj1_C\tnMin1_C\tfrac1_C\tnMaj2_C\tnMin2_C\tfrac2_C\tSDfrac_C\tSDfrac_C_BS\tfrac1_C_0.025\tfrac1_C_0.975\tnMaj1_D\tnMin1_D\tfrac1_D\tnMaj2_D\tnMin2_D\tfrac2_D\tSDfrac_D\tSDfrac_D_BS\tfrac1_D_0.025\tfrac1_D_0.975\tnMaj1_E\tnMin1_E\tfrac1_E\tnMaj2_E\tnMin2_E\tfrac2_E\tSDfrac_E\tSDfrac_E_BS\tfrac1_E_0.025\tfrac1_E_0.975\tnMaj1_F\tnMin1_F\tfrac1_F\tnMaj2_F\tnMin2_F\tfrac2_F\tSDfrac_F\tSDfrac_F_BS\tfrac1_F_0.025\tfrac1_F_0.975\n"

meta_table.write(head)

# construct a dict of battenberg file infos
dict_batt = list()
for i_batt_line in range(1,len(tumour_battenberg)):
	split_batt_line = tumour_battenberg[i_batt_line].strip().split()
	batt_line_split = "	".split(tumour_battenberg[i_batt_line])
	dict_batt.append([split_batt_line[0],int(split_batt_line[1]),int(split_batt_line[2]),"\t".join(split_batt_line[1:len(split_batt_line)])])

batt_abs = str("NA\t1\tNA\tNA\t1\t1\t1\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")

# parse vcf file
for vcf_line in tumour_vcf:
	
	# .FILTER is empty (in tumour2)

	statut = "absent"
	if len(vcf_line.ALT) > 1 :
		print("Warning : more than one ALT")
		continue
	vcf_line.ALT = str(vcf_line.ALT)
	
	info = ""
	if not 'DB' in vcf_line.INFO:
		vcf_line.INFO['DB'] = "False"
	info += str(vcf_line.INFO['DB'])+"\t"+str(vcf_line.INFO['SOMATIC'])+"\t"+str(vcf_line.INFO['VT'])+"\t"

	samples = "" 
	for sample in vcf_line.samples:
		samples += sample['GT']+"\t"+str(sample['AD'][0])+"\t"+str(sample['AD'][1])+"\t"+str(sample['BQ'])+"\t"+str(sample['DP'])+"\t"+str(sample['FA'])+"\t"+str(sample['SS'])+"\t"
	
	for batt_info in dict_batt:
		if vcf_line.CHROM == batt_info[0]:
			if vcf_line.POS >= batt_info[1] and vcf_line.POS <= batt_info[2]:
				meta_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s\n"%
					(vcf_line.CHROM,
					vcf_line.POS,
					vcf_line.ID,
					vcf_line.REF,
					vcf_line.ALT[1:len(vcf_line.ALT)-1],
					vcf_line.QUAL,
					info,
					samples,
					batt_info[3]
					))
				statut = "present"
				continue

	# if there is no CN data for this alteration
	if statut == "absent" :
	    	meta_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s\t%s\t%s\n"%
					(vcf_line.CHROM,
					vcf_line.POS,
					vcf_line.ID,
					vcf_line.REF,
					vcf_line.ALT[1:len(vcf_line.ALT)-1],
					vcf_line.QUAL,
					info,
					samples,
					vcf_line.POS,
					vcf_line.POS,
					str(batt_abs)
					))
		#print ("Warning : no CN state for the alteration %s %s %s %s"%
		#	(vcf_line.CHROM,
		#	vcf_line.POS,
		#	vcf_line.REF,
		#	vcf_line.ALT[1:len(vcf_line.ALT)-1]))