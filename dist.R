setwd('~/Projects/raf')

files = system('ls TCGA*.txt | head -10', intern=T)
names = system('ls TCGA*.txt | cut -d- -f2,3 | sed \'s/.txt//\'| head -10', intern=T)
# matrix of samples af
d=c()
d_tr=c()
d_tv=c()
stats = c()

# read af 
for (i in 1:length(files)) {
  a = read.table(files[i], header=TRUE, sep='\t')
  tr = a[which( ((a$ref=='G' & a$var=='A') | (a$ref=='A' & a$var=='G') | (a$ref=='C' & a$var=='T') | (a$ref=='T' & a$var=='C')) ), ]
  tv = a[which( (!(a$ref=='G' & a$var=='A') & !(a$ref=='A' & a$var=='G') & !(a$ref=='C' & a$var=='T') & !(a$ref=='T' & a$var=='C')) ), ]
  head(tv)
  
  af = sample(a$refCount/a$total, 10000)
  af_tr = sample(tr$refCount/tr$total, 1000)
  af_tv = sample(tv$refCount/tv$total, 1000)
  d = cbind(d,af)
  d_tr = cbind(d_tr, af_tr)
  d_tv = cbind(d_tv, af_tv)
  stats =rbind(stats, summary(af))
}

#plot
colnames(d) = names
boxplot(d, names=names, las=2)
boxplot(d_tr, names=names, las=2)
boxplot(d_tv, names=names, las=2)
stats

