# bgzip and tabix index filtered Manta, Delly, CNVnator, and Lumpy outputs for inheritance checking
# Leave files as they are and save compressed versions to a folder here
# These files are small, so it should be alright to keep both compressed and uncompressed versions
# Manta
lines=`wc -l ../../2022_02_14/manta/manta_files.csv | cut -f 1 -d ' '`
for (( line=2; line<=$lines; line++ ));
do
	FILE=`head -n $line ../../2022_02_14/manta/manta_files.csv | tail -1 | cut -f 1 -d ,`
	SAMPLE=`head -n $line ../../2022_02_14/manta/manta_files.csv | tail -1 | cut -f 2 -d ,`

	bgzip -c $FILE > compressed_caller_outputs/manta/$SAMPLE'.manta.filtered.vcf.gz'
	tabix compressed_caller_outputs/manta/$SAMPLE'.manta.filtered.vcf.gz'
done

# CNVnator
lines=`wc -l ../../2022_02_14/cnvnator/cnvnator_files.csv | cut -f 1 -d ' '`
for (( line=2; line<=$lines; line++ ));
do
	FILE=`head -n $line ../../2022_02_14/cnvnator/cnvnator_files.csv | tail -1 | cut -f 1 -d ,`
	SAMPLE=`head -n $line ../../2022_02_14/cnvnator/cnvnator_files.csv | tail -1 | cut -f 2 -d ,`

	bgzip -c $FILE > compressed_caller_outputs/cnvnator/$SAMPLE'.cnvnator.filtered.vcf.gz'
	tabix compressed_caller_outputs/cnvnator/$SAMPLE'.cnvnator.filtered.vcf.gz'
done

# Delly
lines=`wc -l ../../2022_02_14/delly/delly_files.csv | cut -f 1 -d ' '`
for (( line=2; line<=$lines; line++ ));
do
	FILE=`head -n $line ../../2022_02_14/delly/delly_files.csv | tail -1 | cut -f 1 -d ,`
	SAMPLE=`head -n $line ../../2022_02_14/delly/delly_files.csv | tail -1 | cut -f 2 -d ,`

	bgzip -c $FILE > compressed_caller_outputs/delly/$SAMPLE'.delly.filtered.vcf.gz'
	tabix compressed_caller_outputs/delly/$SAMPLE'.delly.filtered.vcf.gz'
done

# Lumpy
lines=`wc -l ../lumpy/lumpy_files.csv | cut -f 1 -d ' '`
for (( line=2; line<=$lines; line++ ));
do
	FILE=`head -n $line ../lumpy/lumpy_files.csv | tail -1 | cut -f 1 -d ,`
	SAMPLE=`head -n $line ../lumpy/lumpy_files.csv | tail -1 | cut -f 2 -d ,`

	bgzip -c $FILE > compressed_caller_outputs/lumpy/$SAMPLE'.lumpy.filtered.vcf.gz'
	tabix compressed_caller_outputs/lumpy/$SAMPLE'.lumpy.filtered.vcf.gz'
done
