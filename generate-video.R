library(ggplot2)
library(igraph)
library(ggraph)
library(ggpubr)

system("mkdir tmp-img/")

for(expt.type in c("Empirical", "Simulation")) {
  if(expt.type == "Empirical") {
    
    # empirical studies
    src.titles = c("mtGFP 1",
                   "mtGFP 2",
                   "msh1",
                   "msh1 + cipro")
    src.traj.files = c("plant-mito-dynamics/mtgfp-rawtrajectories/mtGFP-3.xml-rawtrajs.csv",
                       "plant-mito-dynamics/mtgfp-rawtrajectories/mtGFP-16.xml-rawtrajs.csv",
                       "plant-mito-dynamics/msh1-rawtrajectories/MSH-2.xml-rawtrajs.csv",
                       "cipro-rawtrajectories/msh1Cip-38.xml-rawtrajs.csv")
    src.am.files = c("plant-mito-dynamics/mtgfp-rawtrajectories/mtGFP-3.xml-amlist.csv",
                     "plant-mito-dynamics/mtgfp-rawtrajectories/mtGFP-16.xml-amlist.csv",
                     "plant-mito-dynamics/msh1-rawtrajectories/MSH-2.xml-amlist.csv",
                     "cipro-rawtrajectories/msh1Cip-38.xml-amlist.csv")
    src.vid.files = c("plant-mito-dynamics/mtgfp-videos/GFP3.avi",
                      "plant-mito-dynamics/mtgfp-videos/GFP16.avi",
                      "plant-mito-dynamics/msh1-videos/MSH2.avi",
                      "cipro-videos/msh1Cip-38-msh1onCip-s3-3ddriftcorrect-cell2.avi")
  } else {
    # simulation studies
    src.titles = paste("Simulation", 1:12)
    src.traj.files = paste("test", 1:12, ".csv", sep="")
    src.am.files = paste("test", 1:12, ".csv-am.csv", sep="")
    src.vid.files = rep(NA, 12)
  }
  
  for(expt in 1:length(src.traj.files)) {
    # initialise storage
    cell.vis = am.vis = list()
    md = 0.475
    this.label = gsub(" ", "-", src.titles[expt])
    # extract data for this experiment
    sub.traj = read.csv(src.traj.files[expt])
    sub.am = read.csv(src.am.files[expt])
    if(colnames(sub.traj)[1] == "frame") {
      colnames(sub.traj) = c("t", "traj", "x", "y")
    }
    sub.amg = graph_from_edgelist(as.matrix(sub.am[,2:3]+1, directed=FALSE))
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
        geom_path(data=sub.traj.t, aes(x=x,y=y,color=factor(traj)), alpha=0.5) + 
        geom_point(data=sub.traj.t[sub.traj.t$t == i,], aes(x=x,y=y), size=0.5, color="#66FF66") +theme_void() +theme(legend.position="none")  + 
        theme(plot.margin = unit(c(md,md,md,md), "cm"),
              panel.background = element_rect(fill = "black", color = "black")) +
        xlim(0, max(c(sub.traj.t$x, sub.traj.t$y))) + 
        ylim(0, max(c(sub.traj.t$x, sub.traj.t$y)))
      
      # plot social network using previous layout. alpha reports whether we've passed the timepoint where each edge appears or not
      am.vis[[i]] = ggraph(layout_data) + 
        geom_edge_link(aes(alpha=-sign(frame-i)), color="#AAAAAA") + 
        geom_node_point(size=0.1) + 
        ggtitle(src.titles[expt]) +
        scale_edge_alpha_continuous(range=c(0,1)) + 
        theme_void() + theme(legend.position="none") +
        theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"),
              plot.title = element_text(hjust = 1, family = "serif")) 
      
      
      # output pair of images to file and append filename to list
      sf = 2
      fname = paste0(c("tmp-img/expt-", expt, "-vis-", i, ".png"), collapse="")
      fset = c(fset, fname)
      png(fname, width=400*sf, height=200*sf, res=72*sf)
      print(ggarrange(cell.vis[[i]], am.vis[[i]]))
      dev.off()  
    }
    compile.str = paste0(c("convert -delay 20 -loop 0 ", 
                           paste0(fset, collapse=" "), " ",
                           this.label, "-vis-animated-smaller.gif"), collapse="")
    system(compile.str)
    
    if(!is.na(src.vid.files[expt])) {
      Sys.sleep(1)
      compile.str = paste0(c("./concatenate-video.sh ", this.label,  
                             "-vis-animated-smaller.gif ", 
                             gsub(" ", "\\\\ ", src.vid.files[expt]), 
                             " ", this.label, 
                             "-combined.mp4"), collapse="")
      system(compile.str)
      
    } else {
      compile.str = paste0(c("convert ", this.label, 
                             "-vis-animated-smaller.gif ", this.label,
                             "-simulation.mp4"), collapse="")
      system(compile.str)
    }
  }
}

# stack several of these outputs for illustrations
system("./stacker.sh")
