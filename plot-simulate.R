library(ggplot2)
library(gridExtra)
library(igraph)
library(ggraph)
library(ggbeeswarm)
library(dplyr)
library(tidyr)

######## first a collection of example parameterisations

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
  g1a[[expt]] = ggplot()  + geom_rect(aes(xmin = -2, xmax = 32, ymin = -2, ymax = 102), fill="black") +
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
grid.arrange(g1a[[4]], g3[[4]],
             g1a[[10]], g3[[10]],
             g1a[[12]], g3[[12]], nrow=3, ncol=2,
             widths=c(1,3))
dev.off()

sf = 2
png("example-traces-nets-2.png", width=700*sf, height=400*sf, res=72*sf)
grid.arrange(g1a[[4]], g3[[4]],
             g1a[[10]], g3[[10]],
             g1a[[12]], g3[[12]], nrow=2, ncol=4,
             widths=c(1,2,1,2))
dev.off()

ga = grid.arrange(grobs = g1a, nrow=2)
gc = grid.arrange(grobs = g3, nrow=2)

gac = grid.arrange(ga, gc, nrow=2)

grid.arrange(g1a[[4]], g3[[4]],
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
grid.arrange(gac, gb, widths=c(1.4,1))
dev.off()

myres = 2
png(paste(c("examples", exptset, ".png"), collapse=""), width=1000*myres, height=400*myres, res=72*myres)
grid.arrange(ga, gb, nrow=1, widths=c(2,1))  
dev.off()

######## now onto parameter scan

# pull summary statistics for large parameter scan from the full simulation code
stats.df = read.csv(paste(c("outstatsscan", exptset, ".csv"), collapse=""))

# plot of theoretical behaviours
gx = ggplot(stats.df, aes(x=log(1/meanmin), y=log(1/meanedges),color=expt)) + 
  geom_point(alpha=0.5, size=0.25) + theme_classic() + theme(legend.position="none") +
  xlab("Physical clumping (log 1/meanmin)") + ylab("Exchange isolation (log 1/meanedges)")

myres = 2
png("fullset.png", width=400*myres, height=300*myres, res=72*myres)
gx 
dev.off()

# downsample for plotting convenience
stats.df$x = log(1/stats.df$meanmin)
stats.df$y = log(1/stats.df$meanedges)

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

grid.arrange(p.D, p.kon, p.koff, p.V, p.dmito, p.kmito, nrow=2)

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
gposts = grid.arrange(ggplot(mins, aes(x=D)) + geom_histogram(),
                      ggplot(mins, aes(x=kon)) + geom_histogram(),
                      ggplot(mins, aes(x=koff)) + geom_histogram(),
                      ggplot(mins, aes(x=V)) + geom_histogram(),
                      ggplot(mins, aes(x=dmito)) + geom_histogram(),
                      ggplot(mins, aes(x=kmito)) + geom_histogram(),
                      nrow=2)

png("inference-expt.png", width=800*myres, height=300*myres, res=72*myres)
grid.arrange(ghull, gposts, nrow=1)
dev.off()

######## now engage with experiments

n.wt = 18
n.msh1 = 28
n.friendly = 19

# read in WT samples
expts.df = expts.traj.df = ams.df = data.frame()
for(i in 1:n.wt) {
  fname = paste(c("plant-mito-dynamics-main/mtgfp-rawtrajectories/mtGFP-", i, ".xml-stats.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "WT"
  tmp$mean.halo.n = NULL
  expts.df = rbind(expts.df, tmp)
  fname = paste(c("plant-mito-dynamics-main/mtgfp-rawtrajectories/mtGFP-", i, ".xml-rawtrajs.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "WT"
  expts.traj.df = rbind(expts.traj.df, tmp)
  fname = paste(c("plant-mito-dynamics-main/mtgfp-rawtrajectories/mtGFP-", i, ".xml-amlist.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "WT"
  ams.df = rbind(ams.df, tmp)
}

# read in msh1 samples
for(i in 1:n.msh1) {
  fname = paste(c("plant-mito-dynamics-main/msh1-rawtrajectories/MSH-", i, ".xml-stats.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "msh1"
  tmp$mean.halo.n = NULL
  expts.df = rbind(expts.df, tmp)
  fname = paste(c("plant-mito-dynamics-main/msh1-rawtrajectories/MSH-", i, ".xml-rawtrajs.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "msh1"
  expts.traj.df = rbind(expts.traj.df, tmp)
  fname = paste(c("plant-mito-dynamics-main/msh1-rawtrajectories/MSH-", i, ".xml-amlist.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "msh1"
  ams.df = rbind(ams.df, tmp)
}

# read in friendly samples
for(i in 1:n.friendly) {
  fname = paste(c("plant-mito-dynamics-main/friendly-rawtrajectories/Friendly-", i, ".xml-stats.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "friendly"
  tmp$mean.halo.n = NULL
  expts.df = rbind(expts.df, tmp)
  fname = paste(c("plant-mito-dynamics-main/friendly-rawtrajectories/Friendly-", i, ".xml-rawtrajs.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "friendly"
  expts.traj.df = rbind(expts.traj.df, tmp)
  fname = paste(c("plant-mito-dynamics-main/friendly-rawtrajectories/Friendly-", i, ".xml-amlist.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "friendly"
  ams.df = rbind(ams.df, tmp)
}

# visualise (previously read) simulations and experiments together
# a nice example set of experiments
nice.g = c("WT", "WT", "msh1", "friendly")
nice.e = c(3, 5, 1, 1)

# construct cell and adjacency matrix visualisations for each
cell.vis = am.vis = list()
for(i in 1:length(nice.g)) {
  this.g = nice.g[i]
  this.e = nice.e[i]
  
  sub = expts.traj.df[expts.traj.df$glabel==this.g & expts.traj.df$elabel == this.e & expts.traj.df$t<=10,]
  sub10 = sub[sub$t==10,]
  
  cell.vis[[i]] = ggplot()  + #geom_rect(aes(xmin = -2, xmax = 32, ymin = -2, ymax = 102), fill="black") +
    geom_path(data=sub, aes(x=y,y=x,color=factor(traj)), alpha=0.5) + 
    geom_point(data=sub10, aes(x=y,y=x), size=0.5, color="#66FF66") +theme_void() +theme(legend.position="none")  + 
    #theme(plot.background = element_rect(fill = "black")) + 
    theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"),
          panel.background = element_rect(fill = "black", color = "black")) + facet_wrap(~ elabel, scale="free")
  
  amg = graph_from_edgelist(as.matrix(ams.df[ams.df$glabel=="WT" & ams.df$elabel==i,2:3], directed=F))
  
  am.vis[[i]] = ggraph(amg, layout="nicely") + geom_edge_link(alpha=0.4, color="#AAAAAA") + geom_node_point(size=0.1) + theme_void() +
    theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))
}

# pull a collection together for visualisation
g.all.1 = ggarrange(g1a[[4]], g1a[[10]], g1a[[8]], g1a[[5]], nrow=1)
g.all.2 = ggarrange(plotlist=cell.vis, nrow=1)
g.all.3 = ggarrange(g3[[4]], g3[[10]], g3[[8]], g3[[5]], nrow=1)
g.all.4 = ggarrange(plotlist = am.vis, nrow=1)

sf = 2
png("all-demo-image.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(g.all.1, g.all.2, g.all.3, g.all.4, nrow=2, ncol=2, heights=c(2,1))
dev.off()

########## actual optimality analysis

# get mito counts from different experiments
nset.df = data.frame()
for(glab in c("WT", "msh1", "friendly")) {
  for(elab in 1:max(expts.traj.df$elabel[expts.traj.df$glabel==glab])) {
    sub = expts.traj.df[expts.traj.df$glabel==glab & expts.traj.df$elabel == elab,]
    nset.df = rbind(nset.df, data.frame(glabel=glab, elabel=elab, n=nrow(sub[sub$t==1,])))
  }
}

# timeframe for comparing simulation and experiment
key.frame = 100

# loop through experiments and compare to subset of simulations with comparable geometry
complist = list()
for(glab in c("WT", "msh1", "friendly")) {
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
png(paste0("ind-trajs-", trajmult, ".png", collapse=""), width=1000*sf, height=1200*sf, res=72*sf)
ggarrange(plotlist = complist)
dev.off()

######## pooled morphospace

# decide on a Pareto summary stat
samples$stat = samples$meanedges / samples$meanmin 
ggplot() +
  geom_point(data=samples, aes(x=x,y=y,color=log(stat))) 

# plot the pooled space with experimental data
expts.df$glabel = factor(expts.df$glabel, levels=c("WT", "msh1", "friendly"))
g.all.1 = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=1) +
  scale_color_viridis() +
  geom_point(data=expts.df[expts.df$frame==key.frame,], 
             aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), size=3, shape=21) + 
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal()

# plot zoomed-in version
g.all.2 = ggplot() +
  geom_point(data=samples[samples$V<=2 & samples$kmito*samples$D < 0.1 & samples$Cn < 150,], aes(x=x,y=y,color=log(stat)), alpha=0.25, size=10) +
  scale_color_viridis() +
  geom_point(data=expts.df[expts.df$frame==key.frame,], 
             aes(x = log(1/mean.min.dist), y=log(1/mean.degree), fill=glabel), size=2, shape=21) + 
  labs(x="Clumping [log(1/min dist)]", y="Loneliness [log(1/mean degree)]", color="Pareto statistic", fill="Experiment") +
  theme_minimal() + xlim(-1.5,0) + ylim(-2.5, -0.75) + theme(legend.position = "none")

# plot showing differences in Pareto stat (plus followup stats)
sub.expts = expts.df[expts.df$frame==key.frame,]
sub.expts$stat = (sub.expts$mean.degree/sub.expts$mean.min.dist) #log(sub.expts$mean.min.dist / sub.expts$mean.degree)
g.all.3 = ggplot(sub.expts, aes(x=glabel, y=stat, fill=glabel, color=glabel)) + 
  geom_boxplot(alpha=0.2, outliers=FALSE, width=0.2) +
  geom_beeswarm() +
  labs(x="Experiment", y="Pareto statistic") +
  theme_minimal() + theme(legend.position = "none")

kruskal.test(sub.expts$stat, sub.expts$glabel)
summary(aov(stat ~ glabel, data=sub.expts))
wilcox.test(sub.expts$stat[sub.expts$glabel=="WT"], sub.expts$stat[sub.expts$glabel=="msh1"])
wilcox.test(sub.expts$stat[sub.expts$glabel=="WT"], sub.expts$stat[sub.expts$glabel=="friendly"])
wilcox.test(sub.expts$stat[sub.expts$glabel=="friendly"], sub.expts$stat[sub.expts$glabel=="msh1"])

# megaplot of the morphospace picture
ggarrange( g.all.1, ggarrange(g.all.2, g.all.3, nrow=2, labels=c("B", "C")), labels=c("A", ""), nrow=1, widths=c(1.5,1))

############# set-valued optimisation
#### still a WIP
summary_df <- samples %>%
  group_by(group) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

########### "inference without data"
# further inference, based on Pareto front and experimental observations

# get the mean wildtype behaviour
mean.x = mean(expts.df$mean.min.dist[expts.df$glabel=="WT"])
mean.y = mean(expts.df$mean.degree[expts.df$glabel=="WT"])

# subset simulations close to this
# previously we built "mins", a dataframe of simulations on the Pareto front
delta = 2
inv.set = mins[abs(mins$meanmin-mean.x) < delta & abs(mins$meanedges-mean.y) < delta,]
samples$hist.ref = "All"
mins$hist.ref = "Pareto"
inv.set$hist.ref = "Proximal"
hist.set = rbind(samples, mins, inv.set)
new.hist.set = hist.set %>% pivot_longer(cols=c("D", "kon", "koff", "V", "dmito", "kmito", "inter"))

# plot parameters in each class of closeness
ggplot(new.hist.set, aes(x=value, y=after_stat(density), fill=hist.ref)) +
  geom_histogram(position = "dodge", width=2) + 
  facet_wrap(~name, scales = "free") +
  theme_minimal()

####### exploring parameter determinants of Pareto front
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

png("scatter-vars.png", width=900*sf, height=600*sf, res=72*sf)
grid.arrange(grobs=p.list)
dev.off()

# So the specific case koff = 0 splits the distribution. Removing this:
# - Higher kon/koff and V pushes us to the low-y Pareto bound
# - Higher kmito generally pushes us to the Parteo front
# - Density based sets up an x-gradient across the scatter

