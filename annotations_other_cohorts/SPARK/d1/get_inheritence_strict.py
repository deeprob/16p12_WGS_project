#!/bin/python3

# this script adds inheritence patterns

import sys
import vcf
import pandas as pd
import os

infilename = 'vcfs/formatted/ssc_d1.select_annotations.txt'
sample_map_filename = '/data4/SPARK/Data1_Regeneron/SPARK.27K.mastertable.20190501.txt'


def get_record_at(vcfr, chrom, pos, alt_allele, ref_allele):
    for rec in vcfr.fetch(chrom, int(pos)-1, int(pos)+1):
        alt_alleles = rec.ALT
        # if start positions don't match go to next one
        if int(pos) != rec.POS:
            continue
        # if ref alleles don't match go to next one
        if ref_allele != rec.REF:
            continue
        # check that alt allele in anno record
        if alt_allele in alt_alleles:
            return rec

def quality_score(rec, data):
    # returns True if PASS and False if FAIL
    if rec.QUAL < 50:
        return False
    # check if GT is empty
    if data.GT == './.':
        return False
    # check if DP is none
    if data.DP is None:
        return False
    # check if FORMAT/DP>=8
    if data.DP < 8:
        return False
    # if Reference allele no need to check further
    if data.GT == '0/0':
        return True
    # check if FORMAT/AD[:1]>0
    if data.AD[1] == 0:
        return False
    # check if (FORMAT/AD[:1])/(FMT/DP)>=0.25
    if data.AD[1] / data.DP < 0.25:
        return False
    # check if (FMT/AD[:1])/(FMT/DP)<=0.75
    if data.AD[1] / data.DP > 0.75:
        return False
    # check if QUAL/(FMT/AD[:1])>=1.5
    if rec.QUAL / data.AD[1] < 1.5:
        return False
    return True



# load in samp mapp
sampmapp = pd.read_csv(sample_map_filename, sep='\t')
sampmapp = sampmapp.set_index('spid')

f = open(infilename, 'r')

for line in f:
    #
    # header
    if line.startswith('Family'):
        new_line = line.strip()
        new_line = new_line + '\tinheritence_strict\n'
        sys.stdout.write(new_line)
        continue
    #
    sline = line.strip().split('\t')
    #
    family      = sline[0]
    role        = sline[1]
    sample      = sline[2]
    chrom       = sline[3]
    pos         = sline[4]
    ref_allele  = sline[5]
    alt_allele  = sline[6]
    #
    # if sample is a parent then skip
    if role == 'Mother' or role == 'Father':
        new_line = line.strip()
        new_line = new_line + '\t.\n'
        sys.stdout.write(new_line)
        continue
    #
    # find out who is father and who is mother
    mother = sampmapp.loc[sample, 'mother']
    father = sampmapp.loc[sample, 'father']
    #
    # open mother vcf (if doesn't exist then continue)
    mother_vcf_filename = '/data4/SPARK/Data1_Regeneron/wes.gVCF/{}/{}.acmg56excluded.g.vcf.gz'.format(mother[-1], mother)
    if not os.path.exists(mother_vcf_filename):
        new_line = line.strip()
        new_line = new_line + '\t.\n'
        sys.stdout.write(new_line)
        continue
    mother_vcfr = vcf.Reader(filename=mother_vcf_filename)
    #
    # open father vcf (if doesn't exist then continue)
    father_vcf_filename = '/data4/SPARK/Data1_Regeneron/wes.gVCF/{}/{}.acmg56excluded.g.vcf.gz'.format(father[-1], father)
    if not os.path.exists(father_vcf_filename):
        new_line = line.strip()
        new_line = new_line + '\t.\n'
        sys.stdout.write(new_line)
        continue
    father_vcfr = vcf.Reader(filename=father_vcf_filename)
    #
    # get mother vcf record
    mother_record = get_record_at(mother_vcfr, chrom, pos, alt_allele, ref_allele)
    #
    # get father vcf record
    father_record = get_record_at(father_vcfr, chrom, pos, alt_allele, ref_allele)
    #
    # if either mother or father don't have a record then skip
    if mother_record is None or father_record is None:
        new_line = line.strip()
        new_line = new_line + '\t.\n'
        sys.stdout.write(new_line)
        continue
    #
    # get call for mother and father
    mother_call = mother_record.genotype(mother).data
    father_call = father_record.genotype(father).data
    #
    # get quality for mother and father
    mother_quality = quality_score(mother_record, mother_call)
    father_quality = quality_score(father_record, father_call)
    #
    # if either mother or the father do not pass quality then skip
    if (not mother_quality) or (not father_quality):
        new_line = line.strip()
        new_line = new_line + '\t.\n'
        sys.stdout.write(new_line)
        continue
    #
    # finally get inheritence
    if '1' in mother_call.GT and '1' in father_call.GT:
        inheritence = 'both'
    elif '1' in mother_call.GT:
        inheritence = 'mother'
    elif '1' in father_call.GT:
        inheritence = 'father'
    else:
        inheritence = 'de novo'
    #
    # construct new line
    new_line = line.strip()
    new_line = new_line + '\t' + inheritence + '\n'
    sys.stdout.write(new_line)


f.close()











