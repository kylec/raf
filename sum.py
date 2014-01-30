import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input')
parser.add_argument('-o', dest='output')
args = parser.parse_args()

sites = dict()

def main():
    # hash all sites and counts
    files = glob.glob(args.input)
    for file in files:
        f = open(file,'r')
        # skip header
        f.readline() 
        for line in f:
            fields = line.rstrip('\n').split('\t')
            #dict[chrpos]
            key = fields[1] + '\t' + fields[2]
            ref = fields[3]
            var = fields[4]
            refCount = int(fields[6])
            depth = int(fields[8])
            
            if key in sites:
                sites[key][0] += 1
                sites[key][1] += int(refCount)
                sites[key][2] += int(depth)
            else:
                sites[key] = []
                sites[key].append(1)
                sites[key].append(refCount)
                sites[key].append(depth)
                sites[key].append(ref)
                sites[key].append(var)
            
    # print sites and counts
    w = open(args.output, 'w')
    w.write('chr\tpos\tsample\trefCount\ttotal\tref\tvar\n')
    for key, vals in sites.iteritems():
        w.write(key + '\t' + '\t'.join(str(val) for val in vals))
        w.write('\n')
    
if __name__ == '__main__':
    main()