import gzip
import sys

vcf = sys.argv[1]

# Removes vcf records where both alleles are REF to significantly reduce file size

gts_to_skip = ['0/0', '0|0', './.', '.|.', '.'] 
def process_rec(rec):
	# if gt == '0|0' returns empty string
	# else return original string
	
	recs = rec.strip()
	recs = recs.split('\t')
	
	gt_info = recs[-1]
	gt_info = gt_info.split(':')[0]
	
	if gt_info in gts_to_skip:
		recn = ''
	else:
		recn = rec
	
	return recn


with gzip.open(vcf, 'r') as f:
	for rec in f:
		rec = rec.decode("utf-8")
		if rec.startswith('#'):
			sys.stdout.write(rec)
		else:
			recn = process_rec(rec)
			sys.stdout.write(recn)
