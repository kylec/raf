# sum af across samples
setwd('~/Projects/raf')
d = read.table('out.txt', header=T)
head(d)
dim(d)
# summary of af dist
hist(d$refCount/d$total, breaks=seq(0,1,.005))
summary(d$refCount/d$total)

# #samples vs #sites
#hist(log2(d$sample), breaks=seq(0,42,1))
hist(d$sample, breaks=seq(0,42,1))
summary(d$sample)

# # samples vs #af
plot(d$sample, d$refCount/d$total)

# # cov vs #af
plot(d$total, d$refCount/d$total, xlim=c(0,10000))

##------------------------------------------------##
# af distribution of each samples
files = system('ls TCGA*.txt | head -10', intern=T)
names = system('ls TCGA*.txt | cut -d- -f2,3 | sed \'s/.txt//\'| head -10', intern=T)
# matrix of samples af
d=c()
d_tr=c()
d_tv=c()
stats = c()

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

