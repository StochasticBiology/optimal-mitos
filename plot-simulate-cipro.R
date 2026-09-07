library(ggplot2)
library(igraph)
library(ggraph)
library(ggbeeswarm)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggforce)
library(sp)

# ensure pipeline.sh has been run before this code!

# Layout:
#### TASK 1: plot visualisations for a collection of example parameterisations
#### TASK 2: parameter scan through simulated dynamics
#### TASK 3: load experimental data
#### TASK 4: illustration of simulated and experimental data
#### TASK 5: compare experiments to similar-geometry simulations
#### TASK 6: visualise the morphospace of possible behaviours
#### TASK 7: quantify proximity to Pareto front
#### TASK 8: consider set-valued optimisation picture
#### TASK 9: "inference without data"
#### TASK 10: parameter determinants of Pareto front
#### TASK 11: density and other considerations across experiments


#### TASK 1: plot visualisations for a collection of example parameterisations

# before this, ensure that simulate.c has been run (see pipeline.sh)

# simulation code outputs trajectories to test[X].csv and adj mats to test[X].csv-am.csv, 
# where X labels a particular example parameterisation

# read these in and visualise
exptset = ""
picset = 10
g1 = g1a = g3 = list()
for(expt in 1:12) {
  # read trajectories
  fname = paste(c("test", expt-1, exptset, ".csv"), collapse="")
  traj.df = read.csv(fname)
  # subset last few steps of trajectories and plot
  sub = traj.df[traj.df$frame > max(traj.df$frame)-picset,]
  g1[[expt]] = ggplot(sub, aes(x=y,y=x,color=factor(mito)))  + 
    geom_path() +theme_void() +theme(legend.position="none") 
  
  # subset just final point and plot these in green with trajectories
  sub1 = traj.df[traj.df$frame == max(traj.df$frame),]
  #g1a[[expt]] = ggplot()  + geom_rect(aes(xmin = -2, xmax = 32, ymin = -2, ymax = 102), fill="black") +
    
  g1a[[expt]] = ggplot()  + geom_rect(aes(xmin = -17, xmax = 47, ymin = -2, ymax = 102), fill="black") +
    geom_path(data=sub, aes(x=y,y=x,color=factor(mito)), alpha=0.5) + 
    geom_point(data=sub1, aes(x=y,y=x), size=0.5, color="#66FF66") +theme_void() +theme(legend.position="none")  + 
    #theme(plot.background = element_rect(fill = "black")) + 
    theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))
  
  # pull adj mat and plot
  am.df = read.csv(paste(c("test", expt-1, exptset, ".csv-am.csv"), collapse=""))
  edgelist.frame = as.matrix(data.frame(t1 = as.character(am.df$mito1), t2 = as.character(am.df$mito2))) 
  amg = graph_from_edgelist(edgelist.frame, directed=F)
  g3[[expt]] = ggraph(amg, layout="nicely") + geom_edge_link(alpha=0.4, color="#AAAAAA") + geom_node_point(size=0.1) + theme_void() +
    theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))
}

# various layouts for talks
sf = 2
png("example-traces-nets.png", width=400*sf, height=700*sf, res=72*sf)
ggarrange(g1a[[4]], g3[[4]],
             g1a[[10]], g3[[10]],
             g1a[[12]], g3[[12]], nrow=3, ncol=2,
             widths=c(1,3))
dev.off()

sf = 2
png("example-traces-nets-2.png", width=700*sf, height=400*sf, res=72*sf)
ggarrange(g1a[[4]], g3[[4]],
             g1a[[10]], g3[[10]],
             g1a[[12]], g3[[12]], nrow=2, ncol=4,
             widths=c(1,2,1,2))
dev.off()

ga = ggarrange(plotlist = g1a, nrow=2, ncol=length(g1a)/2)
gc = ggarrange(plotlist = g3, nrow=2, ncol=length(g1a)/2)

gac = ggarrange(ga, gc, nrow=2)

ggarrange(g1a[[4]], g3[[4]],
             g1a[[10]], g3[[10]],
             g1a[[8]], g3[[8]],
             g1a[[5]], g3[[5]])

# simulation code also outputs summary statistics (mean min distance, mean edges) to outstats.csv
stats.df = read.csv(paste(c("outstats", exptset, ".csv"), collapse=""))

# plot summary stats
gb = ggplot(stats.df, aes(x=log(1/meanmin), y=log(1/meanedges),color=factor(params))) + geom_point(alpha=0.5) +
  labs(x="Physical clumping (log 1/mean min dist)", y = "Exchange isolation (log 1/mean degree)") + theme_light() +
  theme(legend.position="none")

png("examples-scatter.png", width=400*sf, height=250*sf)
ggarrange(gac, gb, widths=c(1.4,1))
dev.off()

myres = 2
png(paste(c("examples", exptset, ".png"), collapse=""), width=1000*myres, height=400*myres, res=72*myres)
ggarrange(ga, gb, nrow=1, widths=c(2,1))  
dev.off()

#### TASK 2: parameter scan through simulated dynamics

# pull summary statistics for large parameter scan from the full simulation code
stats.df = data.frame()
for(i in 1:5) {
  tmp.df = read.csv(paste0(c("outstatsscan-", i, ".csv"), collapse=""))
  stats.df = rbind(stats.df, tmp.df)
}

# plot of theoretical behaviours
if(FALSE) {
gx = ggplot(stats.df, aes(x=log(1/meanmin), y=log(1/meanedges),color=expt)) + 
  geom_point(alpha=0.5, size=0.25) + theme_classic() + theme(legend.position="none") +
  xlab("Physical clumping (log 1/meanmin)") + ylab("Exchange isolation (log 1/meanedges)")

myres = 2
png("fullset.png", width=400*myres, height=300*myres, res=72*myres)
gx 
dev.off()
}

# downsample for plotting convenience
stats.df$x = log(1/stats.df$meanmin)
stats.df$y = log(1/stats.df$meanedges)

set.seed(1)
samples = stats.df[runif(10000, min=1, max=nrow(stats.df)),]
samples = samples[!is.na(samples$x) & !is.na(samples$y) & samples$y < 1000,]

# assign categories to samples based on trading performance
samples$class = ""
samples$class[samples$y > 0] = "0"
samples$class[samples$y > 2] = "2"
samples$class[samples$y > 4] = "4"

# assemble plots of parameter distributions corresponding to the different classes
p.D = ggplot(samples[samples$class != "",], aes(x=D, fill=class)) + 
  geom_histogram(data=samples[samples$class == "0",], position="dodge", aes(x=D-0.005, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "2",], position="dodge", aes(x=D, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "4",], position="dodge", aes(x=D+0.005, y = after_stat(count / sum(count)))) 
p.kon = ggplot(samples[samples$class != "",], aes(x=kon, fill=class)) + 
  geom_histogram(data=samples[samples$class == "0",], position="dodge", aes(x=kon-0.005, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "2",], position="dodge", aes(x=kon, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "4",], position="dodge", aes(x=kon+0.005, y = after_stat(count / sum(count)))) 
p.koff = ggplot(samples[samples$class != "",], aes(x=koff, fill=class)) + 
  geom_histogram(data=samples[samples$class == "0",], position="dodge", aes(x=koff-0.005, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "2",], position="dodge", aes(x=koff, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "4",], position="dodge", aes(x=koff+0.005, y = after_stat(count / sum(count)))) 
p.V = ggplot(samples[samples$class != "",], aes(x=V, fill=class)) + 
  geom_histogram(data=samples[samples$class == "0",], position="dodge", aes(x=V-0.005, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "2",], position="dodge", aes(x=V, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "4",], position="dodge", aes(x=V+0.005, y = after_stat(count / sum(count)))) 
p.dmito = ggplot(samples[samples$class != "",], aes(x=dmito, fill=class)) + 
  geom_histogram(data=samples[samples$class == "0",], position="dodge", aes(x=dmito-0.005, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "2",], position="dodge", aes(x=dmito, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "4",], position="dodge", aes(x=dmito+0.005, y = after_stat(count / sum(count)))) 
p.kmito = ggplot(samples[samples$class != "",], aes(x=kmito, fill=class)) + 
  geom_histogram(data=samples[samples$class == "0",], position="dodge", aes(x=kmito-0.005, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "2",], position="dodge", aes(x=kmito, y = after_stat(count / sum(count)))) + 
  geom_histogram(data=samples[samples$class == "4",], position="dodge", aes(x=kmito+0.005, y = after_stat(count / sum(count)))) 

ggarrange(p.D, p.kon, p.koff, p.V, p.dmito, p.kmito, nrow=2)

# characterise Pareto front
mins = data.frame()
for(i in (-45:50)/10) {
  sub = samples[samples$y > i & samples$y < i+0.1,]
  if(nrow(sub) > 0) {
    # find best x for this y
    best = which(sub$x == min(sub$x))
    # pull the parameters for this best x
    mins = rbind(mins, sub[best,])
  }
}

chull(samples$x, samples$y)

# plot with Pareto front highlighted
ghull = ggplot(samples, aes(x=x, y=y)) + geom_point() + geom_point(data=mins, aes(x=x, y=y), color="red")

# plot distributions of parameters that lie on the Pareto front
gposts = ggarrange(ggplot(mins, aes(x=D)) + geom_histogram(),
                      ggplot(mins, aes(x=kon)) + geom_histogram(),
                      ggplot(mins, aes(x=koff)) + geom_histogram(),
                      ggplot(mins, aes(x=V)) + geom_histogram(),
                      ggplot(mins, aes(x=dmito)) + geom_histogram(),
                      ggplot(mins, aes(x=kmito)) + geom_histogram(),
                      nrow=2, ncol=3)

png("inference-expt.png", width=800*myres, height=300*myres, res=72*myres)
ggarrange(ghull, gposts, nrow=1)
dev.off()

#### TASK 3: load experimental data

n.wt = 18
n.msh1 = 28
n.friendly = 19
n.cipro = 33
n.msh1.cipro = 40

# read in experimental samples
# before doing this, ensure that the analysis script has been run across all experimental observations (see pipeline.sh)
flabels = c("plant-mito-dynamics-main/mtgfp-rawtrajectories/mtGFP-",
            "plant-mito-dynamics-main/msh1-rawtrajectories/MSH-",
            "plant-mito-dynamics-main/friendly-rawtrajectories/Friendly-",
            "cipro-rawtrajectories/cipro-",
            "cipro-rawtrajectories/msh1Cip-")
expt.ns = c(n.wt,
            n.msh1,
            n.friendly,
            n.cipro,
            n.msh1.cipro)
exptlabels = c("WT",
               "msh1",
               "friendly",
               "cipro",
               "msh1-cipro")

mypng = function(protocol, fstr, width, height, res) {
  fname = paste0(c("pr-", protocol, "-", fstr), collapse="")
  png(fname, width=width, height=height, res=res)
}

protocol = 5
expt.set = 1:5

expts.df = expts.traj.df = ams.df = data.frame()
for(expt in expt.set) {
  for(i in 1:expt.ns[expt]) {
    fname = paste0(flabels[expt], i, ".xml-stats.csv", collapse="")
    tmp = read.csv(fname)
    # because we are putting some files into the same glabel set,
    # figure out if we already have records for this glabel
    # if so, use the biggest current elabel+1, otherwise use elabel=1
    r = length(which(expts.df$glabel == exptlabels[expt]))
    if(r > 0) {
      this.elab = max(expts.df$elabel[expts.df$glabel==exptlabels[expt]])+1
    } else {
      this.elab = 1
    }
    tmp$elabel = this.elab
    tmp$glabel = exptlabels[expt]
    tmp$mean.halo.n = NULL
    expts.df = rbind(expts.df, tmp)
    fname = paste0(flabels[expt], i, ".xml-rawtrajs.csv", collapse="")
    tmp = read.csv(fname)
    tmp$elabel = this.elab
    tmp$glabel = exptlabels[expt]
    expts.traj.df = rbind(expts.traj.df, tmp)
    fname = paste0(flabels[expt], i, ".xml-amlist.csv", collapse="")
    tmp = read.csv(fname)
    tmp$elabel = this.elab
    tmp$glabel = exptlabels[expt]
    ams.df = rbind(ams.df, tmp)
  }
}

#### TASK 4: illustration of simulated and experimental data

# visualise (previously read) simulations and experiments together
# a nice example set of experiments
nice.g = c("WT", "WT", "msh1", "friendly")
nice.e = c(3, 5, 2, 1)
nice.g = c("WT", "msh1", "friendly", "msh1-cipro")
nice.e = c(5, 2, 4, 4)
#nice.g = unique(expts.df$glabel)
#nice.e = rep(1, length(nice.g))

# construct cell and adjacency matrix visualisations for each
cell.vis = am.vis = list()
md = 0.475
for(i in 1:length(nice.g)) {
  this.g = nice.g[i]
  this.e = nice.e[i]
  
  tsamp = 10
  sub = expts.traj.df[expts.traj.df$glabel==this.g & expts.traj.df$elabel == this.e & expts.traj.df$t > tsamp-10 & expts.traj.df$t <= tsamp,]
  sub10 = sub[sub$t==tsamp,]
  
  cell.vis[[i]] = ggplot()  + #geom_rect(aes(xmin = -2, xmax = 32, ymin = -2, ymax = 102), fill="black") +
    geom_path(data=sub, aes(x=y,y=x,color=factor(traj)), alpha=0.5) + 
    geom_point(data=sub10, aes(x=y,y=x), size=0.5, color="#66FF66") +theme_void() +theme(legend.position="none")  + 
    #theme(plot.background = element_rect(fill = "black")) + 
    theme(plot.margin = unit(c(md,md,md,md), "cm"),
          panel.background = element_rect(fill = "black", color = "black")) # + facet_wrap(~ elabel, scale="free")
  
  amg = graph_from_edgelist(as.matrix(
    ams.df[ams.df$glabel==this.g & ams.df$elabel==this.e &
             ams.df$t1 < 250 & ams.df$t2 < 250,2:3], directed=F))
  
  am.vis[[i]] = ggraph(amg, layout="nicely") + geom_edge_link(alpha=0.4, color="#AAAAAA") + geom_node_point(size=0.1) + theme_void() +
    theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))
}

# pull a collection together for visualisation
g.all.1 = ggarrange(g1a[[4]], g1a[[10]], g1a[[8]], g1a[[5]], nrow=1, labels=c("A", "B", "C", "D"))
g.all.2 = ggarrange(plotlist=cell.vis, nrow=1, labels=c("E", "F", "G", "H"))
g.all.3 = ggarrange(g3[[4]], g3[[10]], g3[[8]], g3[[5]], nrow=1)
g.all.4 = ggarrange(plotlist = am.vis, nrow=1)


ggarrange(g1a[[5]], cell.vis[[1]], g3[[5]], am.vis[[1]], nrow=2, ncol=2 )

sf = 2
mypng(protocol, "all-demo-image.png", width=1200*sf, height=450*sf, res=72*sf)
ggarrange(g.all.1, g.all.2, g.all.3, g.all.4, nrow=2, ncol=2, widths=c(1,1.1), heights=c(1.5,1))
dev.off()

#### TASK 5: compare experiments to similar-geometry simulations

# get mito counts from different experiments
nset.df = data.frame()
glab.set = unique(expts.traj.df$glabel)
# glab.set = c("WT", "msh1", "friendly", "cipro")
for(glab in glab.set) {
  for(elab in 1:max(expts.traj.df$elabel[expts.traj.df$glabel==glab])) {
    sub = expts.traj.df[expts.traj.df$glabel==glab & expts.traj.df$elabel == elab,]
    nset.df = rbind(nset.df, data.frame(glabel=glab, elabel=elab, n=nrow(sub[sub$t==1,])))
  }
}

# timeframe for comparing simulation and experiment
key.frame = 100

# loop through experiments and compare to subset of simulations with comparable geometry
complist = list()
for(glab in glab.set) {
  for(elab in 1:max(expts.traj.df$elabel[expts.traj.df$glabel==glab])) {
    df.sub = expts.traj.df[expts.traj.df$elabel==elab & 
                             expts.traj.df$glabel==glab,]
    data.mat = df.sub[,3:4]
    
    pca = princomp(data.mat)
    x.dim = max(pca$scores[,1])-min(pca$scores[,1])
    y.dim = max(pca$scores[,2])-min(pca$scores[,2])
    n.dim = nset.df$n[nset.df$glabel==glab & nset.df$elabel==elab]
    
    # get subset of simulations with comparable geometry
    subbing = stats.df[((stats.df$Cx > x.dim-10 & stats.df$Cx < x.dim+10) & 
                          (stats.df$Cy > y.dim-10 & stats.df$Cy < y.dim+10) & 
                          (stats.df$Cn > n.dim-10 & stats.df$Cn < n.dim+10) &
                          (stats.df$V <= 2 & stats.df$kmito*stats.df$D < 0.2)),]
    
    titlestr = paste0(glab, "-", elab, ": ", round(x.dim), "x", round(y.dim), "x", n.dim, collapse=" ")
    complist[[length(complist)+1]] = ggplot() + 
      geom_point(data=subbing, aes(x=x, y=y), alpha = 0.2) + 
      geom_point(data = expts.df[expts.df$glabel == glab & expts.df$elabel == elab &
                                   expts.df$frame==key.frame,], 
                 aes(x = log(1/mean.min.dist), y=log(1/mean.degree)), 
                 color="red") +
      theme_light() + theme(legend.position = "none") + 
      xlab("Physical clumping (log 1/meanmin)") + ylab("Exchange isolation (log 1/meanedges)") +
      ggtitle(titlestr)
    
    
  }
}

# construct megaplot of each of these
mypng(protocol, paste0("ind-trajs.png", collapse=""), width=1000*sf, height=2000*sf, res=72*sf)
ggarrange(plotlist = complist)
dev.off()

#### TASK 6: visualise the morphospace of possible behaviours

nlevel= length(unique(expts.df$glabel))
col.set = c("#FFFFFF", viridis::viridis(nlevel, option="inferno")[2:nlevel])
col.set.plus = c(col.set, "#AAAAAA")

# decide on a Pareto summary stat
samples$stat = samples$meanedges / samples$meanmin 
ggplot() +
  geom_point(data=samples, aes(x=x,y=y,color=log(stat))) 

# plot the pooled space with experimental data
expts.df$glabel = factor(expts.df$glabel, levels=glab.set)
g.all.1 = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=1) +
  scale_color_viridis() +
  geom_point(data=expts.df[expts.df$frame==key.frame,], 
             aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), size=3, shape=21) + 
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal()

g.all.1.alt = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=1) +
  scale_color_viridis() +
  geom_point(data=expts.df[expts.df$frame==key.frame,], 
             aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), stroke = 0.1, size=2, shape=21) + 
  scale_fill_manual(values=col.set) +
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal()
g.all.1.alt

# plot zoomed-in version
g.all.2 = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=0.25, size=10) +
  scale_color_viridis() +
  geom_point(data=expts.df[expts.df$frame==key.frame,], 
             aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), size=2, shape=21) + 
  labs(x="Clumping [log(1/min dist)]", y="Loneliness\n[log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal() + xlim(-1.5,0) + ylim(-2.5, -0.75) + theme(legend.position = "none")

g.all.2.alt = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=0.25, size=10) +
  scale_color_viridis() +
  geom_point(data=expts.df[expts.df$frame==key.frame,], 
             aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), stroke = 0.2, size=2, shape=21) + 
  scale_fill_manual(values=col.set) +
  labs(x="Clumping [log(1/min dist)]", y="Loneliness\n[log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal() + xlim(-1.5,0) + ylim(-2.5, -0.75) + theme(legend.position = "none")
g.all.2.alt

# plot showing differences in Pareto stat (plus followup stats)
sub.expts = expts.df[expts.df$frame==key.frame,]
sub.expts$stat = (sub.expts$mean.degree/sub.expts$mean.min.dist) #log(sub.expts$mean.min.dist / sub.expts$mean.degree)
g.all.3 = ggplot(sub.expts, aes(x=glabel, y=stat, fill=glabel)) + 
  geom_boxplot(color="black", alpha=0.7, outliers=FALSE, width=0.5) +
  geom_beeswarm(size=2,stroke=0.2, shape=21) +
  scale_color_manual(values=col.set) +
  scale_fill_manual(values=col.set) +
  labs(x="Experiment", y="Pareto statistic") +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1))

pairs = combn(length(unique(sub.expts$glabel)), 2)            # Generate all combinations of 2 integers from 1 to n
comps = lapply(1:ncol(pairs), function(i) pairs[, i])  # Convert columns to a list of vectors

if(protocol==5) {
  comps = list(c(1,2), c(1,3), c(1,4), c(1,5), c(2,5))
}

g.all.3.stats = g.all.3 + stat_compare_means(comparisons=comps, size=2)

mypng(protocol, "boxes-stats.png", width=400*sf, height=400*sf, res=72*sf)
print(g.all.3.stats)
dev.off()

kruskal.test(sub.expts$stat, sub.expts$glabel)
summary(aov(stat ~ glabel, data=sub.expts))
tukey.out = TukeyHSD(aov(stat ~ glabel, data=sub.expts))
tukey.out$glabel[tukey.out$glabel[,4] < 0.05,]
#wilcox.test(sub.expts$stat[sub.expts$glabel=="WT"], sub.expts$stat[sub.expts$glabel=="msh1"])
#wilcox.test(sub.expts$stat[sub.expts$glabel=="WT"], sub.expts$stat[sub.expts$glabel=="friendly"])
#wilcox.test(sub.expts$stat[sub.expts$glabel=="WT"], sub.expts$stat[sub.expts$glabel=="cipro"])
#wilcox.test(sub.expts$stat[sub.expts$glabel=="friendly"], sub.expts$stat[sub.expts$glabel=="msh1"])

# megaplot of the morphospace picture
ggarrange( g.all.1.alt, ggarrange(g.all.2.alt, g.all.3, nrow=2, labels=c("B", "C")), labels=c("A", ""), nrow=1, widths=c(1.5,1))

# for talks
g.all.1.a = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=1) +
  scale_color_viridis() +
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal()

g.all.1.a1 = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=1) +
  scale_color_viridis() +
  geom_point(data=mins, aes(x=x,y=y), color="red") +
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal()


g.all.1.b = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=1) +
  scale_color_viridis() +
geom_point(data=expts.df[expts.df$glabel=="WT" & expts.df$frame==key.frame,], 
           aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), size=3, shape=21) + 
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal()
  
mypng(protocol, "pareto-talks-1.png", width=640*sf, height=360*sf, res=72*sf)
print(g.all.1.a)
dev.off()
mypng(protocol, "pareto-talks-1a.png", width=640*sf, height=360*sf, res=72*sf)
print(g.all.1.a1)
dev.off()
mypng(protocol, "pareto-talks-2.png", width=640*sf, height=360*sf, res=72*sf)
print(g.all.1.b)
dev.off()
mypng(protocol, "pareto-talks-3.png", width=640*sf, height=360*sf, res=72*sf)
print(g.all.1)
dev.off()

#### TASK 7: quantify proximity to Pareto front

null.dists.df = expt.dists.df = data.frame()
sub = samples[1:100,]
for(i in 1:nrow(sub)) {
  this.x = log(1/sub$meanmin[i])
  this.y = log(1/sub$meanedges[i])
  dists = (mins$x-this.x)**2 + (mins$y-this.y)**2
  null.dists.df = rbind(null.dists.df, data.frame(label="All", d=min(dists)))
}
for(glab in glab.set) {
  sub = sub.expts[sub.expts$glabel==glab,]
  for(i in 1:nrow(sub)) {
    this.x = log(1/sub$mean.min.dist[i])
    this.y = log(1/sub$mean.degree[i])
    dists = (mins$x-this.x)**2 + (mins$y-this.y)**2
    expt.dists.df = rbind(expt.dists.df, data.frame(label=glab, d=min(dists)))
  }
}

plot.dists = rbind(null.dists.df, expt.dists.df)

plot.dists$label = factor(plot.dists$label, levels=c(glab.set, "All"))
plot.dists = plot.dists[plot.dists$d > 0,]

dist.colors = c(viridis::viridis(5, option="inferno"), "#AAAAAA")
g.dist.hyp = ggplot(plot.dists, aes(x=label, y=log(d), fill=label)) + 
  geom_boxplot(color="black", alpha=0.7, outliers=FALSE, width=0.5) +
  geom_beeswarm(size=2,stroke=0.2, shape=21) +
  scale_color_manual(values=col.set.plus) +
  scale_fill_manual(values=col.set.plus) + 
  theme_minimal() + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Subset", y="log (distance\nto Pareto front)") + ylim(-2, NA)

sf = 2
mypng(protocol, "multi-opt-cipro.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange( g.all.1.alt, 
           ggarrange(g.all.2.alt, 
                     ggarrange(g.dist.hyp, g.all.3.stats,
                               nrow=1, widths=c(1,1), labels=c("C", "D")), 
                     nrow=2, labels=c("B", "")), 
           nrow = 1, labels=c("A", ""), widths=c(1,1))
dev.off()

median(plot.dists$d[plot.dists$label!="All"])
max(plot.dists$d[plot.dists$label=="All"])

median(plot.dists$d[plot.dists$label=="All"])/median(plot.dists$d[plot.dists$label!="All"])
wilcox.test(plot.dists$d[plot.dists$label=="All"], plot.dists$d[plot.dists$label!="All"],)

wilcox.test(plot.dists$d[plot.dists$label=="All"], plot.dists$d[plot.dists$label=="WT"],)
wilcox.test(plot.dists$d[plot.dists$label=="All"], plot.dists$d[plot.dists$label=="msh1"],)
wilcox.test(plot.dists$d[plot.dists$label=="All"], plot.dists$d[plot.dists$label=="friendly"],)
wilcox.test(plot.dists$d[plot.dists$label=="msh1"], plot.dists$d[plot.dists$label=="WT"],)
wilcox.test(plot.dists$d[plot.dists$label=="friendly"], plot.dists$d[plot.dists$label=="WT"],)
wilcox.test(plot.dists$d[plot.dists$label=="msh1"], plot.dists$d[plot.dists$label=="friendly"],)

#### TASK 8: consider set-valued optimisation picture

#### still a WIP
nsamp = nrow(stats.df[stats.df$expt == 1,])

set.seed(1)
eref = sample(1:max(stats.df$expt), 30)
polys.df = ellipses.df = data.frame()
for(e in eref) {
  sub = stats.df[stats.df$expt == e,]
  sub$x = log(1/sub$meanmin)
  sub$y = log(1/sub$meanedges)
  sub = sub[is.finite(sub$x) & is.finite(sub$y),]
  hull = chull(sub$x, sub$y)
  tmp = data.frame(e = e, xs = sub$x[hull], ys = sub$y[hull])
  polys.df = rbind(polys.df, tmp)
  ellipses.df = rbind(ellipses.df, data.frame(e = e,
                                              x0=mean(sub$x), y0=mean(sub$y),
                                              a=sd(sub$x)/sqrt(nsamp), b=sd(sub$y)/sqrt(nsamp)))
}

g.set.1 = ggplot() + 
  geom_ellipse(data=ellipses.df, 
               aes(x0 = x0, y0 = y0, a = a, b = b, angle=0,
                   fill = factor(e)), alpha = 0.3, color="#FFFFFF") +
  scale_fill_viridis_d() + scale_color_viridis_d() + 
  theme_minimal() + theme(legend.position = "none") +
   labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]")

g.set.2 = ggplot() + 
  geom_polygon(data=polys.df, aes(x=xs, y=ys, fill=factor(e)), alpha = 0.4) + 
  geom_point(data=stats.df[stats.df$expt %in% eref & stats.df$meanedges > 0,], aes(x = log(1/meanmin), y=log(1/meanedges), color=factor(expt))) +
  theme_minimal() + theme(legend.position = "none") +
  scale_fill_viridis_d() + scale_color_viridis_d() +
     labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]")

sf = 2
mypng(protocol, "set-valued.png", width=600*sf, height=250*sf, res=72*sf)
ggarrange(g.set.1, g.set.2)
dev.off()

#### TASK 9: "inference without data"
# further inference, based on Pareto front and experimental observations

# get the mean wildtype behaviour
mean.x = mean(expts.df$mean.min.dist[expts.df$glabel=="WT"])
mean.y = mean(expts.df$mean.degree[expts.df$glabel=="WT"])

# subset simulations close to this
# previously we built "mins", a dataframe of simulations on the Pareto front
delta = 4
inv.set = mins[abs(mins$meanmin-mean.x) < delta & abs(mins$meanedges-mean.y) < delta,]
samples$hist.ref = "All"
mins$hist.ref = "Pareto"
inv.set$hist.ref = "Proximal"
mins$stat = 0
inv.set$stat = 0
### columns bug here XXX
hist.set = rbind(samples, mins, inv.set)
new.hist.set = hist.set %>% pivot_longer(cols=c("D", "kon", "koff", "V", "dmito", "kmito", "inter"))

# plot parameters in each class of closeness
sf = 2
mypng(protocol, "inference.png", width=600*sf, height=200*sf, res=72*sf)
ggplot(new.hist.set, aes(x=factor(value), y=..prop.., group =hist.ref, fill=hist.ref)) +
  geom_bar(position = "dodge", alpha=0.8) + 
  scale_fill_manual(values=c("#FFAAAA", "#AA5555", "#440000")) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  labs(x = "", y="", fill = "Subset of\nmorphospace") +
  theme_minimal()
dev.off()

#### TASK 10: parameter determinants of Pareto front

baseplot = ggplot(data=samples, aes(x=x,y=y)) + scale_color_gradient(low = "blue", high = "red") + theme_light() #+ theme(legend.position = "none") 

p.alpha = 0.2
p.size = 1
p.list = list(
  baseplot + geom_point(aes(color=log(Cn/(Cx*Cy))), alpha = 0.2), 
  baseplot + geom_point(aes(color=Cn), alpha = 0.2), 
  baseplot + geom_point(aes(color=Cx*Cy), alpha = 0.2), 
  baseplot + geom_point(aes(color=inter), alpha=p.alpha, size=p.size), 
  baseplot + geom_point(aes(color=D), alpha=p.alpha, size=p.size),
  baseplot + geom_point(aes(color=V), alpha=p.alpha, size=p.size), 
  baseplot + geom_point(aes(color=dmito), alpha=p.alpha, size=p.size), 
  baseplot + geom_point(aes(color=kmito), alpha=p.alpha, size=p.size),
  baseplot + geom_point(aes(color=kon), alpha=p.alpha, size=p.size), 
  baseplot + geom_point(aes(color=koff), alpha=p.alpha, size=p.size), 
  baseplot + geom_point(aes(color=kon/koff), alpha=p.alpha, size=p.size) )

mypng(protocol, "scatter-vars.png", width=900*sf, height=600*sf, res=72*sf)
ggarrange(plotlist=p.list)
dev.off()

# So the specific case koff = 0 splits the distribution. Removing this:
# - Higher kon/koff and V pushes us to the low-y Pareto bound
# - Higher kmito generally pushes us to the Parteo front
# - Density based sets up an x-gradient across the scatter

#### let's look at speeds
df_diff <- expts.traj.df %>%
  arrange(elabel, glabel, traj, t) %>%  # Ensure data is sorted by label and t
  group_by(elabel, glabel, traj) %>%    # Group by label
  mutate(
    x_diff = x - lag(x),  # Calculate difference in x
    y_diff = y - lag(y),  # Calculate difference in y
    t_diff = t - lag(t)   # Calculate difference in t
  ) %>%
  filter(t_diff == 1) %>%  # Keep only rows where t is contiguous
  select(elabel, glabel, traj, x_diff, y_diff, t_diff)  # Select relevant columns

df_diff$speed = sqrt(df_diff$x_diff**2 + df_diff$y_diff**2)
df_diff$glabel = factor(df_diff$glabel, levels = unique(expts.df$glabel))
mypng(protocol, "speed-distn.png", width=400*sf, height=300*sf, res=72*sf)
ggplot(df_diff[df_diff$speed > 0,], aes(x=glabel, y=log10(speed), fill=glabel)) + 
  geom_violin() +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values=col.set) + 
  labs(x = "Experiment", y = "log(Speed / um s-1)", fill="Experiment") +
         theme_light()
dev.off()

### trajectory counts
traj.df = expts.traj.df %>%
  group_by(glabel, elabel) %>%
  summarise(unique_traj_count = n_distinct(traj))
traj.df$glabel = factor(traj.df$glabel, levels = unique(expts.traj.df$glabel))
ggplot(traj.df, aes(x=glabel, y=unique_traj_count)) + geom_violin() + geom_beeswarm()

traj.count.aov = aov(unique_traj_count ~ glabel, data = traj.df)
summary(traj.count.aov)
TukeyHSD(traj.count.aov)

traj.df.len = expts.traj.df %>%
  group_by(glabel, elabel, traj) %>%
  summarise(traj_count = n()) %>%
  group_by(glabel, elabel) %>%
  summarise(mean_traj_count = mean(traj_count))

traj.df.len$glabel = factor(traj.df.len$glabel, levels = unique(expts.traj.df$glabel))
ggplot(traj.df.len, aes(x=glabel, y=mean_traj_count) ) + geom_violin() + geom_beeswarm()

#### DIFFERENT FRAMES
expts.1 = expts.df[expts.df$frame==1,]
expts.100 = expts.df[expts.df$frame==100,]

if(nrow(expts.1)==nrow(expts.100)) {
  expts.100$n1 = expts.1$num.vertices
} else {
  warning("Final frame - initial frame mismatch")
}

expts.100$glabel = factor(expts.100$glabel, levels = unique(expts.traj.df$glabel))

# Function to calculate the area of the convex hull
calculate_convex_hull_area <- function(data) {
  # Get the convex hull points (indexes)
  hull_indices <- chull(data$x, data$y)
  
  # Extract the coordinates of the convex hull points
  hull_points <- data[hull_indices, ]
  
  # Create a Polygon object from the hull points
  poly <- Polygon(as.matrix(hull_points[, c("x", "y")]))
  
  # Return the area of the polygon
  return(poly@area)
}

# Group the dataframe by the combination of label1 and label2
cell.areas <- expts.traj.df %>%
  group_by(elabel, glabel) %>%
  summarise(area = calculate_convex_hull_area(cur_data())) %>%
  ungroup()
cell.areas$glabel = factor(cell.areas$glabel, levels = unique(expts.traj.df$glabel))
ggplot(cell.areas, aes(x=glabel, y=area)) + geom_violin() + geom_boxplot()

if(FALSE) {
  
expts.traj.df.sub = expts.traj.df[expts.traj.df$traj < 10,]
# intractable on mac
traj.areas <- expts.traj.df.sub %>%
  group_by(elabel, glabel, traj) %>%
  summarise(area = calculate_convex_hull_area(cur_data())) %>%
  ungroup()

traj.areas$glabel = factor(traj.areas$glabel, levels = unique(expts.traj.df$glabel))
ggplot(traj.areas, aes(x=glabel, y=log(area+1))) + geom_violin() + geom_beeswarm()
}



# plot a couple of small ones
if(FALSE) {
sub.small = expts.traj.df[expts.traj.df$glabel=="WT-2" & expts.traj.df$elabel <= 2,]
ggplot(sub.small, aes(x=x,y=y)) + geom_point() + facet_wrap(~elabel, scales = "free")
ggarrange(
  ggplot(expts.100, aes(x=glabel, y=num.vertices)) + geom_violin() + geom_beeswarm() +ggtitle("N vertex"),
  ggplot(expts.100, aes(x=glabel, y=singletons/num.vertices)) + geom_violin() + geom_beeswarm() +ggtitle("Singleton proportion"),
  ggplot(expts.100, aes(x=glabel, y=num.edges)) + geom_violin() + geom_beeswarm() +ggtitle("N edges"),
  ggplot(expts.100, aes(x=glabel, y=betweenness)) + geom_violin() + geom_beeswarm() + ggtitle("Betweenness"),
  ggplot(expts.100, aes(x=glabel, y=mean.min.dist)) + geom_violin() + geom_beeswarm() + ggtitle("Mean min dist"),
   ggplot(expts.100, aes(x=glabel, y=mean.degree)) + geom_violin() + geom_beeswarm() + ggtitle("Mean degree"),
   ggplot(expts.100, aes(x=glabel, y=sd.degree)) + geom_violin() + geom_beeswarm() + ggtitle("SD degree"),
  ggplot(expts.100, aes(x=glabel, y=cc.num)) + geom_violin() + geom_beeswarm() + ggtitle("CC num"),
  ggplot(traj.df.len, aes(x=glabel, y=mean_traj_count) ) + geom_violin() + geom_beeswarm() + ggtitle("Traj length"),
  ggplot(cell.areas, aes(x=glabel, y=area)) + geom_violin() + geom_beeswarm() + ggtitle("Hull areas"),
  nrow=2, ncol=5
)

}

#### TASK 11: density and other considerations across experiments

# attempt some normalisations by cell size
expts.100.size = expts.100
expts.100.size$area = 0
for(i in 1:nrow(expts.100)) {
  ref = which(cell.areas$elabel == expts.100$elabel[i] & cell.areas$glabel == expts.100$glabel[i])
  expts.100.size$area[i] = cell.areas$area[ref]
}

#### CHOICE OF DENSITIES
expts.100.size$density = expts.100.size$num.vertices/expts.100.size$area
#expts.100.size$density = expts.100.size$n1/expts.100.size$area

ggplot(expts.100.size, aes(x=glabel, y=mean.min.dist/sqrt(area))) + geom_boxplot()
ggplot(expts.100.size, aes(x=glabel, y=mean.degree)) + geom_boxplot()

ggplot(expts.100.size, aes(x=glabel, y=mean.min.dist/density)) + geom_boxplot()
ggplot(expts.100.size, aes(x=glabel, y=mean.degree/density)) + geom_boxplot()

g.density = ggplot(expts.100.size, aes(x=glabel, y=density, fill=glabel)) + 
  geom_boxplot(color="black", alpha=0.7, outliers=FALSE, width=0.5) +
  geom_beeswarm(size=2,stroke=0.2, shape=21) +
  scale_color_manual(values=col.set.plus) +
  scale_fill_manual(values=col.set.plus) + 
  theme_minimal() + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Subset", y="Density /\nmitos per cell area") 

  
g.density

sf = 2
mypng(protocol, "multi-opt-cipro-density.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange( g.all.1.alt, 
           ggarrange(g.all.2.alt, g.density,
                    g.dist.hyp, g.all.3.stats,
                               nrow=2, ncol= 2, 
                    labels=c("B", "C", "D", "E")), 
           nrow = 1, labels=c("A", ""), widths=c(1,1))
dev.off()


ggplot(expts.100.size, aes(x=glabel, y=mean.degree)) + geom_violin() + geom_boxplot() 

expts.100.size$traj.len = traj.df.len$mean_traj_count

# these at least induce a shift in the right direction (although the results aren't in the right direction yet)
ggplot(expts.100.size, aes(x=glabel, y=density)) + geom_violin() + geom_boxplot()
ggplot(expts.100.size, aes(x=glabel, y=mean.degree/density)) + geom_violin() + geom_boxplot()
ggplot(expts.100.size, aes(x=glabel, y=mean.degree/traj.len)) + geom_boxplot()

ggarrange(
  ggplot(expts.100.size, aes(x=num.vertices, y=mean.degree, color=glabel)) + geom_point(),
  ggplot(expts.100.size, aes(x=num.vertices, y=1/mean.min.dist, color=glabel)) + geom_point(),
ggplot(expts.100.size, aes(x=density, y=mean.degree, color=glabel)) + geom_point(),
ggplot(expts.100.size, aes(x=density, y=1/mean.min.dist, color=glabel)) + geom_point()
)

ggarrange(
  ggarrange(
    ggplot(expts.100.size, aes(x=density, y=mean.degree, color=glabel)) + geom_point(),
    ggplot(expts.100.size, aes(x=density, y=(1/mean.min.dist)**2, color=glabel)) + geom_point(),
    nrow=1
  ),
  ggarrange(
    ggplot(expts.100.size, aes(x=glabel, y=density)) + geom_violin() + geom_boxplot(),
    ggplot(expts.100.size, aes(x=glabel, y=mean.degree)) + geom_violin() + geom_boxplot(),
    ggplot(expts.100.size, aes(x=glabel, y=(1/mean.min.dist)**2)) + geom_violin() + geom_boxplot(),
    nrow=1
  ),
  nrow=2
)

lm(mean.degree ~ density, data=expts.100.size)

ggplot(expts.100.size, aes(x=glabel, y=area)) + geom_violin() + geom_boxplot()


ggplot(expts.100.size, aes(x=glabel, y=mean.degree/density)) + geom_violin() + geom_boxplot()
ggplot(expts.100.size, aes(x=glabel, y=(1/mean.min.dist**2)/density)) + geom_violin() + geom_boxplot()

expts.100.size$norm.degree = expts.100.size$mean.degree/expts.100.size$density
expts.100.size$norm.mmd = expts.100.size$mean.min.dist/expts.100.size$density
tester = expts.100.size
tester$glabel[tester$glabel=="WT-2"] = "WT"
tester$glabel[tester$glabel=="msh1-2"] = "msh1"
aov.mmd = aov(norm.mmd ~ glabel, data=tester)
summary(aov.mmd)
TukeyHSD(aov.mmd)

aov.deg = aov(norm.degree ~ glabel, data=tester)
summary(aov.deg)
TukeyHSD(aov.deg)

expts.100.size$stat = expts.100.size$mean.degree/expts.100.size$density
ggplot(expts.100.size, aes(x=density, y=stat, color=glabel)) + geom_point()
ggplot(expts.100.size, aes(x=glabel, y=stat, color=glabel)) + geom_beeswarm()
t.test(expts.100.size$stat[expts.100.size$glabel=="WT"], expts.100.size$stat[expts.100.size$glabel=="friendly"])

ggplot(tester, aes(x=glabel, y=mean.degree/density)) + geom_violin() + geom_boxplot()
ggplot(tester, aes(x=glabel, y=(1/mean.min.dist**2)/density)) + geom_violin() + geom_boxplot()


expts.100.bak = expts.100
expts.100 = expts.100.size
png("comparison.png", width=800*sf, height=600*sf, res=72*sf)
ggarrange(
  ggplot(expts.100, aes(x=glabel, y=num.vertices)) + geom_violin() + geom_boxplot() +ggtitle("N vertex"),
  ggplot(expts.100, aes(x=glabel, y=singletons/num.vertices)) + geom_violin() + geom_boxplot() +ggtitle("Singleton proportion"),
  ggplot(expts.100, aes(x=glabel, y=num.edges)) + geom_violin() + geom_boxplot() +ggtitle("N edges"),
  ggplot(expts.100, aes(x=glabel, y=betweenness)) + geom_violin() + geom_boxplot() + ggtitle("Betweenness"),
  ggplot(expts.100, aes(x=glabel, y=mean.min.dist)) + geom_violin() + geom_boxplot() + ggtitle("Mean min dist"),
  ggplot(expts.100, aes(x=glabel, y=mean.degree)) + geom_violin() + geom_boxplot() + ggtitle("Mean degree"),
  ggplot(expts.100, aes(x=glabel, y=max.degree)) + geom_violin() + geom_boxplot() + ggtitle("Max degree"),
  ggplot(expts.100, aes(x=glabel, y=cc.num)) + geom_violin() + geom_boxplot() + ggtitle("CC num"),
  ggplot(traj.df.len, aes(x=glabel, y=mean_traj_count) ) + geom_violin() + geom_boxplot() + ggtitle("Traj length"),
  ggplot(cell.areas, aes(x=glabel, y=area)) + geom_violin() + geom_boxplot() + ggtitle("Hull areas"),
  ggplot(expts.100, aes(x=glabel, y=density)) + geom_violin() + geom_boxplot() + ggtitle("Mito density"),
  ggplot(df_diff, aes(x=glabel, y=log10(speed))) + geom_violin() + geom_boxplot() + ggtitle("Speeds"),
  nrow=3, ncol=4
)
dev.off()

