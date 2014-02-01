# sum af across samples
setwd('~/Projects/raf')
d = read.table('out.txt', header=T)
head(d)
dim(d)

# summary of af dist
png('af-hetnum.png')
#par(mfrow=c(2,1))
af=d$refCount/d$total
h=hist(af, breaks=seq(0,1,.01))
plot(h, freq=FALSE, axes=F, main='Het sites ref allele distribution', xlab='Ref allele fraction', 
     ylab="Fraction of total het sites")
axis(2, pretty(h$density), pretty(h$density)/100)
axis(1, at=seq(0,1,.1))
s=summary(af)
legend('topright', legend=c(paste('med=',s[3]), paste('avg=',s[4])))
dev.off()
max(h$counts)/sum(h$counts)

dd = d[d$total>100 & d$total< 1000,]
d=dd
# #samples vs #sites
png('samples-hetnum.png')
h=hist(d$sample, breaks=seq(0,42,1))
plot(h, freq=FALSE, yaxt='n', main='How common is a het site?', xlab='# Samples with common het sites', 
     ylab="Fraction of total het sites")
axis(2, pretty(h$density), pretty(h$density))
s=summary(d$sample)
legend('topright', legend=c(paste('med=',s[3]), paste('avg=',s[4])))
dev.off()
max(h$counts)/sum(h$counts)

# # samples vs #af
png('samples-af.png')
plot(d$sample, d$refCount/d$total, main='Do common het sites have less variance in af? (100x-1000x)', 
     xlab="# samples with common het sites", ylab="ref allele fraction")
dev.off()

# # cov vs #af
png('cov-af.png')
plot(d$total, d$refCount/d$total, main = 'depth vs af', xlab='depth', ylab='ref allele fraction',
     xlim=c(0,10000))
dev.off()
summary(d$total)


##------------------------------------------------##
# af distribution for ti and tv
names = c('all', 'ti', 'tv')
d=c()
stats = c()

a=d
tr = a[which( ((a$ref=='G' & a$var=='A') | (a$ref=='A' & a$var=='G') | (a$ref=='C' & a$var=='T') | (a$ref=='T' & a$var=='C')) ), ]
tv = a[which( (!(a$ref=='G' & a$var=='A') & !(a$ref=='A' & a$var=='G') & !(a$ref=='C' & a$var=='T') & !(a$ref=='T' & a$var=='C')) ), ]
head(tv)
  
af = sample(a$refCount/a$total, 10000)
af_tr = sample(tr$refCount/tr$total, 10000)
af_tv = sample(tv$refCount/tv$total, 10000)
d = cbind(d,af, af_tr, af_tv)
stats =rbind(stats, summary(af), summary(af_tr), summary(af_tv))

#plot
colnames(d) = names
boxplot(d, names=names, las=2)
legend('topleft', paste('ti/tv=', signif(dim(tr)[1]/dim(tv)[1],digits=3)))
stats

