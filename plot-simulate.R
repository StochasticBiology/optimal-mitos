library(ggplot2)
library(gridExtra)
library(igraph)
library(ggraph)

######## first a collection of example parameterisations

# simulation code outputs trajectories to test[X].csv and adj mats to test[X].csv-am.csv, 
# where X labels a particular example parameterisation

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

# read in WT samples
expts.df = expts.traj.df = data.frame()
for(i in 1:12) {
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
}

# read in msh1 samples
for(i in 1:12) {
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
}

# so revisiting this, when we compare experimental data to simulation WITH SAME DIMENSIONS as the experiment
# we often fall on the front. but there's a subtlety about # of mitos vs # of trajectories to resolve
complist = list()
for(glab in c("WT", "msh1")) {
  for(elab in 1:12) {
    df.sub = expts.traj.df[expts.traj.df$elabel==elab & 
                             expts.traj.df$glabel==glab,]
    data.mat = df.sub[,3:4]
    
    pca = princomp(data.mat)
    x.dim = max(pca$scores[,1])-min(pca$scores[,1])
    y.dim = max(pca$scores[,2])-min(pca$scores[,2])
    n.dim = max(df.sub$traj)
    
  #  subbing = samples[((samples$Cx > x.dim-10 & samples$Cx < x.dim+10) & 
  #                       (samples$Cy > y.dim-10 & samples$Cy < y.dim+10) & 
  #                       (samples$Cn > n.dim-10 & samples$Cn < n.dim+10)),]
    subbing = stats.df[((stats.df$Cx > x.dim-10 & stats.df$Cx < x.dim+10) & 
                      (stats.df$Cy > y.dim-10 & stats.df$Cy < y.dim+10) & 
                      (stats.df$Cn > n.dim-10 & stats.df$Cn < n.dim+10)),]
    
    titlestr = paste0(glab, "-", elab, ": ", round(x.dim), "x", round(y.dim), "x", n.dim, collapse=" ")
    complist[[length(complist)+1]] = ggplot() + 
      geom_point(data=subbing, aes(x=x, y=y), alpha = 0.2) + 
      geom_point(data = expts.df[expts.df$glabel == glab & expts.df$elabel == elab &
                                   expts.df$frame==120,], 
                 aes(x = log(1/mean.min.dist), y=log(1/mean.degree)), 
                 color="red") +
      theme_light() + theme(legend.position = "none") + 
      xlab("Physical clumping (log 1/meanmin)") + ylab("Exchange isolation (log 1/meanedges)") +
      ggtitle(titlestr)
      
    
  }
}

png("ind-trajs.png", width=1000*sf, height=800*sf, res=72*sf)
ggarrange(plotlist = complist)
dev.off()

# here we don't look very Pareto -- but of course this isn't conditioned on density
ggplot() +
  geom_point(data=samples, aes(x=x,y=y,color="color")) +
  geom_point(data=expts.df[expts.df$frame==120,], 
       aes(x = log(1/mean.min.dist), y=log(1/mean.degree), color=glabel)) 
 
           
data.mat = expts.traj.df[expts.traj.df$elabel==12 & expts.traj.df$glabel=="WT",3:4]
ggplot(data.mat, aes(x=x, y=y)) + geom_point()
pca = princomp(data.mat)
ggarrange(ggplot(as.data.frame(data.mat), aes(x=x, y=y)) + geom_point(),

ggplot(as.data.frame(pca$scores), aes(x=Comp.1, y=Comp.2)) + geom_point())
ggplot(expts.traj.df, aes(x=x,y=y,color=factor(traj))) + geom_line() + 
  facet_wrap(glabel~elabel) + theme(legend.position="none")

# plot summary of the experimental samples
ggplot() + 
  geom_path(data = expts.df[expts.df$frame>100,], aes(x = log(1/mean.min.dist), y=log(1/mean.degree), color=factor(glabel)), size=1)

# contruct piece-by-piece summary figure of experimental data overlaid on theoretical behaviours
g.pset = ggplot() + geom_point(data=samples, aes(x=x, y=y, color=expt), alpha = 0.2) + 
  geom_point() + theme_light() + theme(legend.position = "none") + #ylim(-4.5,2.5) + xlim(-1.5,4) +
  xlab("Physical clumping (log 1/meanmin)") + ylab("Exchange isolation (log 1/meanedges)")

g.pset.hull = g.pset + geom_point(data=mins, aes(x=x, y=y), color="red", size = 1) 

g.pset.hull.exp1 = g.pset.hull + 
  geom_point(data = expts.df[expts.df$glabel == "WT" & expts.df$frame>100,], aes(x = log(1/mean.min.dist), y=log(1/mean.degree)), color="black", size=1) 

g.pset.hull.exp2 = g.pset.hull.exp1 + 
   geom_point(data = expts.df[expts.df$glabel == "msh1" & expts.df$frame>100,], aes(x = log(1/mean.min.dist), y=log(1/mean.degree)), color="orange", alpha = 0.2, size=1) 

myres = 2
png("pset-expt-1.png", width=500*myres, height=300*myres, res=72*myres)
print(g.pset)
dev.off()
png("pset-expt-2.png", width=500*myres, height=300*myres, res=72*myres)
print(g.pset.hull)
dev.off()
png("pset-expt-3.png", width=500*myres, height=300*myres, res=72*myres)
print(g.pset.hull.exp1)
dev.off()
png("pset-expt-4.png", width=500*myres, height=300*myres, res=72*myres)
print(g.pset.hull.exp2)
dev.off()

#### pick a subset of samples with cell dimensions that correspond to a given experiment

# this is a reasonable range for experiment 3, for example
xr = c(80, 120); yr = c(20, 40); nr = c(100, 160)
#xr = c(50, 70); yr = c(20, 40); nr = c(60, 100)

subbing = samples[((samples$Cx > xr[1] & samples$Cx < xr[2]) & 
                     (samples$Cy > yr[1] & samples$Cy < yr[2]) & 
                     (samples$Cn > nr[1] & samples$Cn < nr[2])),]
gs.pset = ggplot() + geom_point(data=subbing, aes(x=x, y=y, color=expt), alpha = 0.2) + 
  geom_point() + theme_light() + theme(legend.position = "none") + #ylim(-4.5,2.5) + xlim(-1.5,4) +
  xlab("Physical clumping (log 1/meanmin)") + ylab("Exchange isolation (log 1/meanedges)")

gs.pset.hull = gs.pset #+ geom_point(data=mins, aes(x=x, y=y), color="red", size = 1) 

gs.pset.hull.exp1 = gs.pset.hull + 
  geom_point(data = expts.df[expts.df$glabel == "WT" & expts.df$frame>100,], aes(x = log(1/mean.min.dist), y=log(1/mean.degree)), color="black", size=1) 

gs.pset.hull.exp2 = gs.pset.hull.exp1 + 
  geom_point(data = expts.df[expts.df$glabel == "msh1" & expts.df$frame>100,], aes(x = log(1/mean.min.dist), y=log(1/mean.degree)), color="orange", alpha = 0.2, size=1) 

gs.pset.hull.exp2

#### exploring determinants of Pareto front
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

###### not clear from here onwards

posts = read.csv("outstatsinfer.csv")
posts = posts[posts$score < 2.1*2.1,]

gpostsabc = grid.arrange(ggplot(posts, aes(x=D)) + geom_histogram(),
                         ggplot(posts, aes(x=kon)) + geom_histogram(),
                         ggplot(posts, aes(x=koff)) + geom_histogram(),
                         ggplot(posts, aes(x=V)) + geom_histogram(),
                         ggplot(posts, aes(x=dmito)) + geom_histogram(),
                         ggplot(posts, aes(x=kmito)) + geom_histogram(),
                         ggplot(posts[posts$score < 10,], aes(x=score)) + geom_histogram(),
                         ggplot(posts, aes(x=log(1./meanmin), y=log(1./meanedges), color=score)) + geom_point(),
                         nrow=2)

posts$x = 0; posts$y = 0; posts$type = "ABC"
mins$score = 0; mins$type = "IWD"
all.infer = rbind(posts, mins)

gpostsall = grid.arrange(ggplot(all.infer, aes(x=D, fill=type)) + geom_histogram(position="dodge") + theme(legend.position="none"),
                         ggplot(all.infer, aes(x=kon, fill=type)) + geom_histogram(position="dodge")+ theme(legend.position="none"),
                         ggplot(all.infer, aes(x=koff, fill=type)) + geom_histogram(position="dodge")+ theme(legend.position="none"),
                         ggplot(all.infer, aes(x=V, fill=type)) + geom_histogram(position="dodge")+ theme(legend.position="none"),
                         ggplot(all.infer, aes(x=dmito, fill=type)) + geom_histogram(position="dodge")+ theme(legend.position="none"),
                         ggplot(all.infer, aes(x=kmito, fill=type)) + geom_histogram(position="dodge")+ theme(legend.position="none"),
                         nrow=2)


png("inference2-expt.png", width=800*myres, height=300*myres, res=72*myres)
grid.arrange(ghullexpt2, gpostsall, nrow=1)
dev.off()

gbexpt = gb + geom_path(data=expt.df[expt.df$frame>100,], aes(x = 1/mean.min.dist, y=1/mean.degree, color="a")) + scale_x_continuous(trans="log") + scale_y_continuous(trans="log")
# XXX also use mean.halo.n
grid.arrange(gac, gbexpt, nrow=1)

