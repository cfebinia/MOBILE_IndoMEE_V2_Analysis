library(ggplot2)
library(ggpubr)
library(cowplot)

axis_labels <- function(cap, prefix="PCo"){
  eig <- cap$CA$eig
  var <- eig/sum(eig)
  axis_labs <- paste(prefix,1:length(names(var)), " (", round(var*100,1),"%)", sep="")
  axis_labs
}

pretty_PCA <- function(cap, 
                       axes, 
                       metadata, 
                       colour.by, 
                       colour.scheme,
                       colour.order,
                       shape.by,
                       shape.scheme,
                       shape.order,
                       display=TRUE,
                       output.path=NULL){
  
  PC <- axes
  
  # Get PCoA coordinates
  xycoord <- summary(cap)
  xycoord <- data.frame(xycoord$sites[,PC])
  xycoord$sampleid <- rownames(xycoord)
  rownames(xycoord) <- NULL
  colnames(xycoord) <- c("Xaxis","Yaxis","sampleid")
  
  # Axis limits
  xlims <- c(min(xycoord$Xaxis), max(xycoord$Xaxis)) * 1.1
  ylims <- c(min(xycoord$Yaxis), max(xycoord$Yaxis)) * 1.1
  
  # Merge with metadata
  ggdata <- merge(xycoord, metadata, by="sampleid")
  ggdata <- ggdata[, c(colnames(xycoord), colour.by, shape.by)]
  colnames(ggdata) <- c(colnames(xycoord), "my_colour", "my_shape")
  
  # Order
  ggdata$my_colour <- factor(ggdata$my_colour, levels = colour.order)
  ggdata$my_shape <- factor(ggdata$my_shape, levels = shape.order)
  ggdata$sampleid <- factor(ggdata$sampleid, unique(ggdata$sampleid))
  
  # Eigenvalues and labels
  eig <- cap$CA$eig
  var <- eig/sum(eig)
  axis_labs <- paste("PCo", 1:length(var), " (", round(var*100,1), "%)", sep="")
  xtitle <- axis_labs[PC[1]]
  ytitle <- axis_labs[PC[2]]
  
  mycols <- colour.scheme[colour.order]
  myshapes <- shape.scheme[shape.order]
  
  # make base plot
  a <- ggplot() +
    geom_vline(xintercept = 0, lty="dashed")+
    geom_hline(yintercept=0, lty="dashed")
  
  # Scatter Plot
  a <- a + 
    theme_bw(base_size = 12) +
    geom_point(data=ggdata, aes(x=Xaxis, y=Yaxis, fill=my_colour, shape=my_shape, size=my_shape), col="black") + 
    scale_fill_manual(values=mycols) +
    scale_shape_manual(values=shape.scheme)+
    scale_size_manual(values=c(3.5,3.2,2.7)) +
    scale_x_continuous(expand = c(0,0),
                       limits = xlims,
                       breaks = seq(-10,10,1),
                       position = "top")+
    scale_y_continuous(expand = c(0,0),
                       limits = ylims,
                       breaks = seq(-10,10,1),
                       position = "left") +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size=12),
          plot.margin=margin(0,0,0,0,"cm"),
          legend.position = "none")+
    xlab(xtitle) + ylab(ytitle)
  
  
  # Boxplot X axis
  pop_cols <- mycols[pop.ord] # set color by population
  ggdata$my_colour <- factor(ggdata$my_colour, rev(pop.ord))
  ggdata <- ggdata[order(ggdata$my_colour),]
  b <- ggplot(ggdata,
              aes(x=my_colour, y=Xaxis, fill=my_colour))+
    theme_bw(base_size = 12) +
    geom_hline(yintercept = 0, lty="dashed") +
    geom_boxplot(width=0.8, outlier.size = 2,
                 outlier.shape=21, outlier.fill = "white") +
    scale_fill_manual(values=pop_cols[pop.ord]) +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_blank(),
          axis.title.y = element_text(colour = "white",size=12),
          panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin=margin(0,0,0,0,"cm"),
          legend.position="none")+
    scale_y_continuous(expand = c(0,0),
                       limits = xlims,
                       breaks = seq(-10,10,1),
                       position = "left") +
    xlab("sample distance") + ylab(" ")+
    coord_flip()
  
  # Boxplot Y axis
  ggdata$my_colour <- factor(ggdata$my_colour, pop.ord)
  c <- ggplot(ggdata,
              aes(x=my_colour, y=Yaxis, fill=my_colour)) +
    theme_bw(base_size = 12) +
    geom_hline(yintercept = 0, lty="dashed")+
    geom_boxplot(width=0.8, outlier.size = 2,
                 outlier.shape=21, outlier.fill = "white") +
    scale_fill_manual(values=pop_cols) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x=element_text(colour = "white",size=12),
          plot.margin=margin(0,0,0,0,"cm"),
          legend.position = "none")+
    scale_y_continuous(expand = c(0,0),
                       limits = ylims,
                       breaks = seq(-10,10,1),
                       position = "right")+
    scale_x_discrete(position="top")+
    xlab("sample distance") + ylab(" ")
  
  # fix labels and extract legend
  pop_count <- data.frame(table(metadata[,colour.by]))
  colnames(pop_count)  <- c(colour.by,"sample")
  pop_count[,colour.by] <- factor(pop_count[,colour.by], pop.ord)
  pop_count
  
  poplab <-paste0(pop_count[,colour.by]," (n=",pop_count$sample,")")
  names(poplab) <- pop_count[,colour.by]
  
  ggdata$label <- c()
  for(pop in names(poplab)){
    ggdata[ggdata$my_colour == pop, "label"] <- poplab[pop]
  }
  ggdata <- ggdata[order(ggdata$my_colour),]
  ggdata$label <- factor(ggdata$label, unique(ggdata$label))
  
  legend.col <-  mycols[as.character(unique(ggdata$my_colour))]
  names(legend.col) <- poplab[names(legend.col)]
  
  p <- ggplot(ggdata, aes(x=label, y=Xaxis, fill=label))+ theme_bw(base_size = 10) +
    geom_boxplot(width=0.8, outlier.size = 2, outlier.shape=21, outlier.fill = "white") +
    scale_fill_manual(values=legend.col) +
    theme(legend.text = element_text(size=10),
          legend.title = element_blank(),
          legend.key.size = unit(0.3,"cm"),
          plot.background = element_rect(colour = "white"))
  
  p <- ggplot(ggdata)+ theme_void(base_size = 5) +
    geom_point(aes(x=label, y=Xaxis, shape=my_shape))+
    geom_boxplot(aes(x=label, y=Xaxis, fill=label))+    
    scale_fill_manual(values=legend.col) +
    scale_shape_manual(values=shape.scheme)+
    theme(legend.text = element_text(size=5),
          legend.title = element_blank(),
          legend.key.size = unit(0.3,"cm"),
          legend.box.margin = margin(0,0,0,0,"cm"))
  
  legend <- cowplot::get_legend(p)
  legend <- as_ggplot(legend) + theme(plot.margin = margin(-1.5,1,0,0,"cm"))
  
  ## Legend
  p_legend <- ggplot(ggdata, aes(x=sampleid, y=Xaxis, fill=my_colour, shape=my_shape)) +
    geom_point() +
    theme_void() +
    scale_fill_manual(values = mycols) +
    scale_shape_manual(values = myshapes)
  
  legend <- cowplot::get_legend(p_legend)
  
  ## Final arrangement using ggarrange
  final_plot <- ggarrange(
    ggarrange(a, c, ncol = 2, widths = c(4, 1)),
    ggarrange(b, as_ggplot(legend), ncol = 2, widths = c(4, 1)),
    nrow = 2,
    heights = c(4, 1)
  )
  
  # Save or display
  if (is.null(output.path)) {
    print(final_plot)
  } else {
    filename <- paste0(output.path, paste0("PCo", PC[1], "_", PC[2]))
    ggsave(paste0(filename, ".png"), plot=final_plot, width=10, height=10, dpi=300)
    ggsave(paste0(filename, ".svg"), plot=final_plot, width=10, height=10)
    message("Output saved to: ", output.path)
  }
  
  # display
  if(display==T){
    print(final_plot)
  }
}
