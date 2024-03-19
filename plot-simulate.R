library(ggplot2)
library(gridExtra)
library(av)
library(gganimate)
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
expts.df = data.frame()
for(i in 1:12) {
  fname = paste(c("../plant-mito-dynamics-main\ 3/mtgfp-rawtrajectories/mtGFP-", i, ".xml-stats.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "WT"
  tmp$mean.halo.n = NULL
  expts.df = rbind(expts.df, tmp)
}

# read in msh1 samples
for(i in 1:12) {
  fname = paste(c("../plant-mito-dynamics-main\ 3/msh1-rawtrajectories/MSH-", i, ".xml-stats.csv"), collapse="")
  tmp = read.csv(fname)
  tmp$elabel = i
  tmp$glabel = "msh1"
  tmp$mean.halo.n = NULL
  expts.df = rbind(expts.df, tmp)
}

# plot summary of the experimental samples
ggplot() + 
  geom_path(data = expts.df[expts.df$frame>100,], aes(x = log(1/mean.min.dist), y=log(1/mean.degree), color=factor(glabel)), size=1)

# contruct piece-by-piece summary figure of experimental data overlaid on theoretical behaviours
g.pset = ggplot() + geom_point(data=samples, aes(x=x, y=y, color=expt), alpha = 0.2) + 
  geom_point() + theme_light() + theme(legend.position = "none") + ylim(-4.5,2.5) + xlim(-1.5,4) +
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

##############
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

