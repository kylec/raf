# extract chr, pos, allele counts from vcfs
# kyle chang

import glob
import vcf
import re
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('-o', dest='output')
args = parser.parse_args()

def main():
    vcfs = glob.glob(args.input)
    outputFile = open(args.output, 'w')
    outputFile.write('sample\tchr\tpos\tdbsnp\trefCount\tvarCount\ttotal\n')

    if vcfs:
        for vcf in vcfs:
            print 'processing %s ...' % vcf
            parseVcf(vcf, outputFile)
    else:
        raise Exception('Error: .vcf does not exist')
            
def parseVcf(inputFile, outputFile):
    i = 0
    try:
        # open vcf file
        vcfFile = open(inputFile, 'r')
        vcfReader = vcf.Reader(vcfFile)
        
        for row in vcfReader:
            i += 1
            chrom = re.sub('chr', '', row.CHROM)
            start = row.POS
        
            # is dbsnp = 1 else 0
            dbsnp = '0' if row.ID == '.' else '1'
        
            # TODO: may have multiple samples in these VCFs?
            # first sample is het
            if row.samples[0].is_het:
            
                # row.sample[0] (call object)
                sample = row.samples[0].sample
                shortSample = re.search(r'^\w+-(\w+-\w+)', sample).group(1)

                # return genotype info (call object)
                refCount = row.genotype(sample)['AD'][0]
                varCount = row.genotype(sample)['AD'][1]
                depth = row.genotype(sample)['DP']
            
                outputFile.write('\t'.join([shortSample, chrom, str(start), dbsnp, str(refCount), str(varCount), str(depth)]))
                outputFile.write('\n')
    except Exception:
        print 'Error: line %d in %s' % (i, inputFile)
        sys.exit(1)
            
if __name__ == '__main__':
    main()
    