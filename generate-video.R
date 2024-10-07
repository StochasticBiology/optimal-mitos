library(ggplot2)
library(gridExtra)
library(igraph)
library(ggraph)
library(ggbeeswarm)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggforce)

######## now engage with experiments

n.wt = 18
n.msh1 = 28
n.friendly = 19
n.cipro = 33
n.msh1.cipro = 40
n.no.cipro = 32
n.msh1.no.cipro = 35

# read in WT samples
flabels = c("plant-mito-dynamics-main/mtgfp-rawtrajectories/mtGFP-",
            "plant-mito-dynamics-main/msh1-rawtrajectories/MSH-",
            "plant-mito-dynamics-main/friendly-rawtrajectories/Friendly-",
            "plant-mito-dynamics-main/cipro-rawtrajectories/cipro-",
            "plant-mito-dynamics-main/cipro-rawtrajectories/msh1Cip-",
            "plant-mito-dynamics-main/cipro-rawtrajectories/mtgfp-",
            "plant-mito-dynamics-main/cipro-rawtrajectories/msh1MS-")
expt.ns = c(n.wt,
            n.msh1,
            n.friendly,
            n.cipro,
            n.msh1.cipro,
            n.no.cipro,
            n.msh1.no.cipro)
exptlabels = c("WT",
               "msh1",
               "friendly",
               "cipro",
               "msh1-cipro",
               "WT-2",
               "msh1-2")

mypng = function(protocol, fstr, width, height, res) {
  fname = paste0(c("pr-", protocol, "-", fstr), collapse="")
  png(fname, width=width, height=height, res=res)
}

expt.set = c(1)
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


# visualise (previously read) simulations and experiments together
# a nice example experiment
this.g = "WT"
this.e = 3

# initialise storage
cell.vis = am.vis = list()
md = 0.475

# extract data for this experiment
sub.traj = expts.traj.df[expts.traj.df$glabel==this.g & expts.traj.df$elabel == this.e,]
sub.am = ams.df[ams.df$glabel==this.g & ams.df$elabel== this.e,]
sub.amg = graph_from_edgelist(as.matrix(sub.am[,2:3], directed=F))
# pull the frame in which each edge first appears
E(sub.amg)$frame = sub.am[,1]

# create a graph layout for the final social network
layout_data <- create_layout(sub.amg, layout = 'nicely')

# initialise list of filenames
fset = c()
# loop through timesteps
for(i in 1:max(sub.traj$t)) {
  # pull trajectories up to this time point, and plot
  sub.traj.t = sub.traj[sub.traj$t <= i,]
  cell.vis[[i]] = ggplot()  + 
    geom_path(data=sub.traj.t, aes(x=y,y=x,color=factor(traj)), alpha=0.5) + 
    geom_point(data=sub.traj.t[sub.traj.t$t == i,], aes(x=y,y=x), size=0.5, color="#66FF66") +theme_void() +theme(legend.position="none")  + 
    theme(plot.margin = unit(c(md,md,md,md), "cm"),
          panel.background = element_rect(fill = "black", color = "black"))
  
  # plot social network using previous layout. alpha reports whether we've passed the timepoint where each edge appears or not
  am.vis[[i]] = ggraph(layout_data) + 
    geom_edge_link(aes(alpha=-sign(frame-i)), color="#AAAAAA") + 
    geom_node_point(size=0.1) + 
    scale_edge_alpha_continuous(range=c(0,1)) + theme_void() + theme(legend.position="none") +
    theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))
  
  # output pair of images to file and append filename to list
  sf = 2
  fname = paste0(c("vis-", i, ".png"), collapse="")
  fset = c(fset, fname)
  png(fname, width=400*sf, height=200*sf, res=72*sf)
  print(ggarrange(cell.vis[[i]], am.vis[[i]]))
  dev.off()  
}

# command line using ImageMagick to create an animated gif
compile.str = paste0(c("convert -delay 20 -loop 0", fset, "vis-animated-smaller.gif"), collapse=" ")
system(compile.str)

