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
#d=dd
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
     xlim=c(0,5000))
dev.off()
summary(d$total)

#### example desc stats
library(pastecs)
options(digits=3)
stat.desc(d$refCount/d$total, basic=F) 
summary(d$refCount/d$total)

# filter by coverage
d = read.table('out.txt', header=T)
e=d[which(d$total>20 & d$total< 10000),]
d=e

# sample rows in data
prob = signif(median(d$refCount/d$total), 3)
probs =c(prob-.2,prob-.15,prob-.1, prob-.05, prob, prob+.05, prob+.1, prob+.15, prob+.2)

# plot expected vs emmperical
#for (i in 1:1) {
for (i in 1:length(probs)) {
  d_sample = d[sample(nrow(d), size=50000),]
  ref_sim = rbinom(length(d_sample$total), d_sample$total, prob=probs[i])
  png(paste("exp",i,"png",sep="."))
  plot(d_sample$refCount/d_sample$total, d_sample$total, type='p',cex=.5, col=rgb(100,0,0,50,maxColorValue=200), pch=16,xlim=c(0,1), xlab="", ylab="")
  par(new=T)
  plot(ref_sim/d_sample$total, d_sample$total, type='p',cex=.5, xlim=c(0,1),  col=rgb(0,0,0,50,maxColorValue=200), pch=16, xlab = "Reference Allele Fraction", ylab = "Depth")
  par(new=F)
  legend("topright", c(paste("prob=",probs[i])))
  dev.off()
}

# overall median .54 , take rbinom of a depth given by a red point
library(pastecs)
stat.desc(d1$refCount/d1$total, basic=F)
head(ref_sim)

# rbinom( length(depth, depth, prob = .54)

#### qqplot
qqplot(ref_sim, ref)

#### af distribution for ti and tv
names = c('all', 'ti', 'tv')
d=c()
stats = c()

a=dd
tr = a[which( ((a$ref=='G' & a$var=='A') | (a$ref=='A' & a$var=='G') | (a$ref=='C' & a$var=='T') | (a$ref=='T' & a$var=='C')) ), ]
tv = a[which( (!(a$ref=='G' & a$var=='A') & !(a$ref=='A' & a$var=='G') & !(a$ref=='C' & a$var=='T') & !(a$ref=='T' & a$var=='C')) ), ]
head(tv)
  
af = sample(a$refCount/a$total, 10000)
af_tr = sample(tr$refCount/tr$total, 10000)
af_tv = sample(tv$refCount/tv$total, 10000)
d = cbind(d,af, af_tr, af_tv)
stats =rbind(stats, summary(af), summary(af_tr), summary(af_tv))

#### plot
colnames(d) = names
boxplot(d, names=names, las=2)
legend('topleft', paste('ti/tv=', signif(dim(tr)[1]/dim(tv)[1],digits=3)))
stats

#### individual sample plots
files = system("ls TCGA*.txt", intern=T)
for (i in 1:length(files)){
  # plot individual af
  a = read.table(files[i], header=T)
  png(paste(files[i],"png", sep="."))
  plot(a$total, a$refCount/a$total, xlab ="depth", ylab = "ref allele fraction")
  dev.off()
}

#### sample vs sample af
k=length(files)/2
k=1
for (i in 1:k) {
  png(paste("af", i, "png", sep="."))
  a = read.table(files[i], header=T)
  b= read.table(files[i+1], header=T)
  i = i + 1
  c = merge(a, b, by=c('chr','pos'))
  
  #filtered
  cc= c[which( (c$refCount.x+c$varCount.x)>20 & (c$refCount.x+c$varCount.x) < 1000 & 
                          (c$refCount.y+c$varCount.y)>20 & (c$refCount.y+c$varCount.y) < 1000 ),]
  #head(c_filtered)
  cat(dim(cc)[1]/dim(c)[1])
  depth_x = cc$refCount.x+cc$varCount.x
  depth_y = cc$refCount.y+cc$varCount.y
  af_x = cc$refCount.x/depth_x
  af_y = cc$refCount.y/depth_y
  plot(af_x,af_y, xlab=c$sample.x[1], ylab=c$sample.y[1], xlim=c(0,1) , ylim=c(0,1), col=rgb(100,0,0,50,maxColorValue=200), pch=16)
  dev.off()
}
dim(c)
dim(cc)
c$sample.x[1]
x = c$refCount.x/(c$refCount.x+c$varCount.x)
y = c$refCount.y/(c$refCount.y+c$varCount.y)


# look at the allelels deviated. and see how it looks another person