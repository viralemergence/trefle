## trefle phylofactor of missing viruses
## danbeck@ou.edu
## last update 4/28/2021

## clean
rm(list=ls()) 
graphics.off()

## libraries
library(magrittr)
library(tidyverse)
library(ape)
library(caper)
library(phylofactor)
library(readxl)
library(treeio)
library(MASS)
library(ggtree)
library(plotrix)

## load trefle data
#setwd("~/Github/trefle")
setwd("~/Desktop/trefle")
tre <- read_csv("artifacts/trefle.csv")

## reorganize
tre %<>% rename(Virus = 'virus', Host = 'host')
tre %<>% mutate(Trefle = 1)

## read in clover and organize
clo <- read_csv("data/clover.csv")
clo %<>% dplyr::select(Virus, Host)
clo %<>% mutate(Clover = 1)

## merge
full_join(clo, tre) %>% 
  mutate(Clover = replace_na(Clover, 0)) %>%
  mutate(Trefle = replace_na(Trefle, 0)) -> df

## reorganize
df %<>% dplyr::mutate(Missing = as.numeric((Trefle-Clover)==1))

## aggregate missing
df1=df
df1 %<>% group_by(Host) %>% summarize(Missing = sum(Missing))

## aggregate clover
df2=df
df2 %<>% group_by(Host) %>% summarize(Clover = sum(Clover))

## combine
df=merge(df1,df2,by="Host")
rm(df1,df2)

## clean
rm(clo,tre)

## visualize
ggplot(df,aes(log1p(Clover),log1p(Missing)))+
  geom_point(alpha=0.5)+
  geom_smooth(method="gam",formula=y~s(x,bs="ts"))

## load Upham phylogeny
setwd("~/Desktop/clover/data/phylogenies")
tree=read.nexus("upham_tree_666.nex")

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))

## load translation table
trans=read.csv("mammal_phylo_translations.csv")

## translate
df=merge(df,trans,by="Host",all.x=T)
rm(trans)

## use upham names
df$tips=df$Host_Upham

## load taxonomy
setwd("~/Desktop/clover/data/phylogenies/Upham_S1_Data")
taxa=read_excel('masterTaxonomy.xlsx',3)

## get tips
taxa$tips=taxa$SciName
taxa$tips=sapply(strsplit(taxa$tips,'_'),function(x) paste(x[1],x[2],sep=' '))

## uppercase first
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## fix
taxa$Order=firstup(tolower(taxa$Order))
taxa$Family=firstup(tolower(taxa$Family))

## trim
taxa=data.frame(taxa[c("tips","Order","Family","Genus")])

## merge
df=merge(df,taxa,by="tips",all.x=T)
rm(taxa)

## duplicate tips
dups=table(df$tips)>1
dups=names(which(dups))
dups=df[df$tips%in%dups,]
length(unique(dups$tips))

## remove rows from df
df=df[!df$tips%in%dups$tips,]

## save metadata
meta=dups[!duplicated(dups$tips),]
meta$Missing=NULL
meta$Clover=NULL

## aggregate
dups1=aggregate(Missing~tips,data=dups,mean)
dups2=aggregate(Clover~tips,data=dups,mean)
dups=merge(dups1,dups2,by="tips")
rm(dups1,dups2)
dups$Missing=round(dups$Missing,0)
dups$Clover=round(dups$Clover,0)

## merge with metadata
dups=merge(dups,meta,by="tips",all.x=T)
rm(meta)

## rebind
df$type="single"
dups$type="averaged"
dups=dups[names(df)]
df=rbind.data.frame(df,dups)
rm(dups)

## double check
table(table(df$tips)>1)

## clean
df$Host_Upham=NULL

## trim the tree
ptree=keep.tip(tree,df$tips)
rm(tree)

## merge phylo and data
df$label=df$tips
df=df[match(ptree$tip.label,df$tips),]
cdata=comparative.data(phy=ptree,data=df,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## taxonomy
cdata$data$Species=cdata$data$tips
cdata$data$taxonomy=with(cdata$data,paste(Order,Family,Genus,Species,sep="; "))

## pgls lambda
# mod=pgls(Missing~1,data=cdata,lambda="ML")
# summary(mod)

## save tree
cdata$data$label=rownames(cdata$data)
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## negbin model fcn
model.fcn2 <- function(formula,data,...){
  fit <- tryCatch(MASS::glm.nb(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}

## negbin objective function
obj.fcn2 <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
  #if (!'negbin' %in% class(fit) & !'glm' %in% class(fit) & !'lm' %in% class(fit))
  if (!'negbin' %in% class(fit))
  {
    return(0)
  }
  else 
  {
    #fit2 <- MASS::glm.nb(Z.poisson~1,data = fit$model)
    fit$null.deviance-fit$deviance %>% return()
    #fit$twologlik %>% return()
  }
}

## phylofactor negative binomial
set.seed(1)
pf=gpf(Data=cdata$data,tree=cdata$phy,
            frmla.phylo=Missing~phylo,
            model.fcn = model.fcn2,
            objective.fcn = obj.fcn2,
            cluster.depends='library(MASS)',
            algorithm='phylo',nfactors=5,
            min.group.size=10)

## Holm procedure
hp=HolmProcedure(pf)
  
## set key
setkey(pf$Data,'Species')
  
## make data
dat=data.frame(pf$Data)
  
## make clade columns in data
for(i in 1:hp){
  
  dat[,paste0("Missing",'_pf',i)]=ifelse(dat$tip%in%cladeget(pf,i),'factor','other')
  
}

## make data frame to store taxa name, response, mean, and other
results=data.frame(matrix(ncol=6, nrow = hp))
colnames(results)=c('factor','taxa','tips','node',"clade",'other')

## set taxonomy
taxonomy=dat[c('Species','taxonomy')]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## loop
for(i in 1:hp){
  
  ## get taxa
  tx=pf.taxa(pf,taxonomy,factor=i)$group1
  
  ## get tail
  tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
  
  ## combine
  tx=paste(tx,collapse=', ')
  
  # save
  results[i,'factor']=i
  results[i,'taxa']=tx
  
  ## get node
  tips=cladeget(pf,i)
  node=ggtree::MRCA(pf$tree,tips)
  results[i,'tips']=length(tips)
  results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                           ifelse(is.null(node) & length(tips)!=1,NA,node))
  
  ## get means
  ms=(tapply(dat[,"Missing"],dat[,paste0("Missing",'_pf',i)],mean,na.rm=T))
  
  ## round
  ms=round(ms,0)
  
  ## add in
  results[i,'clade']=ms['factor']
  results[i,'other']=ms['other']
  
}

## plot base tree
pbase=ggtree(dtree,layout="fan",branch.length="none",size=0.25)

## get tree data
tdata=pbase$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=rescale(tdata$Missing,c(max(tdata$x),xmax)),
                tips=tdata$label)

## set x max
plus=4
pplus=plus+0.5

## plot clades
gg=pbase
for(i in 1:nrow(results)){
  
  ## add clade and label
  gg=gg+
    geom_hilight(node=results$node[i],
                 fill=ifelse(results$clade[i]<results$other[i],"blue","red"),
                 alpha=0.2)+
    geom_cladelabel(node=results$node[i],
                    label=results$factor[i],
                    offset=plus,
                    offset.text=pplus)
}

## add missing
gg=gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.5,alpha=0.5)

## export
setwd("~/Desktop/trefle/figures")
png("phylofactor-missing.png",width=5,height=5,units="in",res=600)
gg
dev.off()

## clean table
results2=results
results2$node=NULL
write.csv(results2,"ED table phylofactor.csv")

## organize clade
dat$pf=ifelse(dat$Missing_pf1=="factor","1",
              ifelse(dat$Missing_pf2=="factor","2",
                     ifelse(dat$Missing_pf3=="factor","3",
                            ifelse(dat$Missing_pf4=="factor","4","other"))))
dat$pf=factor(dat$pf)

## separate
dat2=dat[dat$pf=="other",]
dat3=dat[!dat$pf=="other",]

## plot scaling
ggplot(dat,aes(log1p(Clover),Missing))+
  
  ## plot nonpf points first
  geom_point(data=dat2,alpha=0.5,colour="grey80")+

  ## gam
  geom_smooth(method="gam",formula=y~s(x,bs="ts"),colour="black",alpha=0.5)+
  
  ## plot pf
  geom_point(data=dat3,alpha=0.75,size=2,aes(colour=pf))+
  scale_colour_viridis_d()+
  scale_fill_viridis_d()+
  
  ## pf smooths
  geom_smooth(data=dat3,method="gam",formula=y~s(x,bs="ts",k=5),
              aes(group=pf,colour=pf,fill=pf),alpha=0.25)+
  
  ## theme
  theme_bw()+
  theme(legend.position="top")+
  guides(colour=guide_legend(title="phylofactor clade"),
         fill=guide_legend(title="phylofactor clade"))
