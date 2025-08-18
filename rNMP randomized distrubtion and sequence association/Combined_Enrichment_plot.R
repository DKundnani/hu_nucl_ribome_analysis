#Necessary libraries
library(grid)
library(gridExtra)
library(ggpubr)
library(egg)

#Auxiliary function
corresponding_deoxy=function(ribo){
  if(ribo!='U'){
    return(ribo)
  }else{
    return('T')
  }
}




#Getting the necessary data
window_size=1e6
list_chrom=const.CHROMS
for(chrom_list in list_chrom){
  sample_list=const.get_sample_list('normal')
  
  alpha=0.5
  ticks_size=25
  ticks_size_y=18
  text_size=30
  text_size_y=30
  linewidth=1
  celltypes_colors=c('CD4T'='#FF7F0E','hESC'='#2CA02C','HEK293T'='#D62728',
                     'RNH2A-KO T3-8'='#E377C2','RNH2A-KO T3-17'='#9467BD')
  
  
  #Upload Data
  Data <- tidyr::expand_grid(sample=sample_list,chrom=chrom_list) |> 
    purrr::pmap(function(sample,chrom){
      cell=const.get_sample_info(sample)$cell
      fn.ribos_binned(sample,chrom,window_size) |> read_rds() |> 
        as_tibble() |>  mutate(cell=cell) |> 
        gather(key='riboNuc','count',contains('count')) |> 
        separate('riboNuc',c('common','riboNuc')) |> 
        select(-common) |> dplyr::rename(chrom=seqnames) |> 
        select(cell,sample,chrom,strand,start,end,width,riboNuc,count)
    }) |> purrr::list_rbind()
  
  Data$cell <- factor(Data$cell,levels=c('CD4T','hESC-H9','HEK293T-WT','HEK293T-RNASEH2A-KO-T3-17','HEK293T-RNASEH2A-KO-T3-8'),
                      labels = c('CD4T','hESC','HEK293T','RNH2A-KO T3-17','RNH2A-KO T3-8'))
  
  #Factor columns
  Data$sample <- Data$sample |> factor(levels=sample_list,labels=purrr::map(sample_list,const.get_label))
  Data <- Data |> mutate_if(is.character, as.factor)
  
  #Get Background Data
  Data_bg <- purrr::map(chrom_list,function(chrom){
    fn.kmer_binned(chrom,window_size,'1') |> read_rds() |> as_tibble() |> 
      gather(key='deoxyNuc','count',contains('count')) |> 
      separate('deoxyNuc',c('common','deoxyNuc')) |> 
      dplyr::rename(chrom=seqnames) |> select(chrom,start,end,width,strand,deoxyNuc,count)
  })|> purrr::list_rbind()
  
  #Factor columns
  Data_bg <- Data_bg |> mutate_if(is.character, as.factor)
  
  Data_bg_ef=Data_bg |> group_by(chrom,strand,start,end,width) |> 
    mutate(EF=count/sum(count)) |> ungroup()
  
  
  #Relative Frequency: Each Ribo, with respect to the Deoxys in its given window
  grouping_vars=c("cell","sample","chrom","strand","start","end","width")
  
  Data_ef <- Data |> rowwise() |> mutate(deoxyNuc=corresponding_deoxy(riboNuc)) |> 
    left_join(Data_bg,by=c('chrom','strand','start','end','width','deoxyNuc')) |> 
    group_by(cell,sample,chrom,strand,start,end,width) |> 
    mutate(EF=(count.x/sum(count.x)/(count.y/sum(count.y))))
  
  Data_ef$EF[is.nan(Data_ef$EF)] <- 1
  
  Data_ef_cell=Data_ef |> group_by(cell,chrom,strand,start,end,width,riboNuc) |> summarise(EF=mean(EF))
  
  
  p_up <- ggplot()+geom_line(data=Data_ef_cell |> filter(strand=='+'),
                             aes(x=start,y=EF,color=cell,group=interaction(cell,riboNuc)),linewidth=linewidth,alpha=alpha)+
    facet_grid(rows=vars(riboNuc))+geom_hline(yintercept=1,color='#212F3C',linetype='dashed')+
    geom_hline(yintercept=0,color='#212F3C',linetype='solid')
  
  
  p_down<- ggplot()+geom_line(data=Data_ef_cell |> filter(strand=='-'),
                              aes(x=start,y=-EF,color=cell,group=interaction(cell,riboNuc)),linewidth=linewidth, alpha=alpha)+
    facet_grid(rows=vars(riboNuc))+geom_hline(yintercept=-1,color='#212F3C',linetype='dashed')+
    geom_hline(yintercept=0,color='#212F3C',linetype='solid')
  
  
  p_up <- p_up+scale_y_continuous(label=c('',1,''),breaks=c(0,1,2))+
    ylab('W(+)\nRelative Fraction')+xlab('')+
    scale_x_continuous(label=utils.get_pretty_base_pairs)+
    theme_minimal()+theme(legend.position = 'none',axis.line.y = element_line(),
                          axis.line.x=element_blank(),
                          axis.text.x = element_text(size=ticks_size),
                          axis.text.y = element_text(size=ticks_size_y),
                          strip.text.y = element_text(angle=0,size=text_size_y),
                          axis.title = element_text(size=text_size),
                          panel.spacing = unit(1, "lines"),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank())+
    scale_color_manual(values=celltypes_colors)
  
  p_down <- p_down+scale_y_continuous(label=c('',1,''),breaks=c(0,-1,-2))+
    ylab('C(-)\nRelative Fraction')+xlab(paste('Chromosome',substring(chrom_list,4,5)))+
    scale_x_continuous(label=utils.get_pretty_base_pairs)+
    theme_minimal()+theme(legend.position = 'none',axis.line = element_line(),
                          axis.line.x=element_blank(),
                          axis.text.x = element_text(size=ticks_size),
                          axis.text.y = element_text(size=ticks_size_y),
                          strip.text.y = element_text(angle=0,size=text_size_y),
                          axis.title = element_text(size=text_size),
                          panel.spacing = unit(1, "lines"),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_color_manual(values=celltypes_colors)
  
  
  #Adding Extra label
  labels_df<- data.frame(riboNuc=c('A','C','G','U'), lab_up=c('','','','0'),
                         lab_down=c('0','','',''))
  
  min_x_coord=0
  max_x_coord=max(Data_ef$start)
  x_position=-max_x_coord*0.06
  size_extra_label=6
  
  p_down<-p_down + geom_text(data=labels_df, x =x_position, y = 0, aes(label = lab_down), size = size_extra_label) +
    coord_cartesian(clip = "off",xlim=c(min_x_coord,max_x_coord)) +
    theme(plot.margin = unit(c(1,1,1,3), "lines"))
  
  p_up<-p_up + geom_text(data=labels_df, x =x_position, y = 0, aes(label = lab_up), size = size_extra_label) +
    coord_cartesian(clip = "off",xlim=c(min_x_coord,max_x_coord)) +
    theme(plot.margin = unit(c(1,1,1,3), "lines"))
  
  #Uploading other info
  file_ef='deepali_data/REZ_100000_1.5_0.8_ef.tsv'
  
  Data_ef0 <- file_ef |> read_tsv() |> filter(Chromosome %in% chrom_list) 
  
  Data_ef2 <- Data_ef0 |> mutate(start=((Start)%/%window_size)*window_size,
                                 end=(End%/%window_size)*window_size,
                                 width_i=End-Start) |> 
    group_by(Chromosome,Strand,Celltype,Sample,start) |> 
    summarise(EF=sum(EF*width_i)/sum(width_i), width=sum(width_i), end=max(end)) |> 
    group_by(Chromosome,Strand,Celltype,start,end) |> 
    summarise(EF=mean(EF))
  Data_ef2$Celltype <- factor(Data_ef2$Celltype,levels = c('CD4T','hESC','HEK293T',"RNH2A-KO T3-17","RNH2A-KO T3-8"))

  
  p0 <- ggplot()+geom_line(data=Data_ef2 |> filter(Strand=='+'),aes(x=start,y=EF,color=Celltype),linewidth=linewidth, alpha=alpha)+
    geom_line(data=Data_ef2 |> filter(Strand=='-'),aes(x=start,y=-EF,color=Celltype),linewidth=linewidth, alpha=alpha)+
    geom_hline(yintercept=0)+geom_hline(yintercept = 1,color='#212F3C',linetype='dashed')+
    geom_hline(yintercept = -1,color='#212F3C',linetype='dashed')
  
  p0 <- p0 +scale_y_continuous(label=abs,breaks=c(-3,-2,-1,0,1,2,3))+
    ylab('rNMP EF\nC(-)     W(+)')+xlab('')+
    scale_x_continuous(label=utils.get_pretty_base_pairs)+
    guides(color=guide_legend(title='Cell Type'))+theme_minimal()+
    scale_color_manual(values=celltypes_colors)+
    theme(axis.text.x = element_text(size=ticks_size),
          axis.text.y = element_text(size=ticks_size_y),
          axis.line = element_line(),
          axis.title = element_text(size=text_size))
  
  p_legend=get_legend(p0)
  p0 <- p0+theme(legend.position='none',
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p1 <- ggplot()+geom_line(data=Data2 |> filter(strand=='+'),aes(x=start,y=EF,color=cell),linewidth=linewidth, alpha=alpha)+
    geom_line(data=Data2 |> filter(strand=='-'),aes(x=start,y=-EF,color=cell),linewidth=linewidth, alpha=alpha)
  
  p1 <- p1 +scale_y_continuous(label=abs)+
    ylab('rNMP EF')+xlab('')+
    scale_x_continuous(label=utils.get_pretty_base_pairs)+
    guides(color=guide_legend(title='Cell Type'))+theme_minimal()+
    scale_color_manual(values=celltypes_colors)+
    theme(axis.text = element_text(size=ticks_size),axis.line = element_line(),
          axis.title = element_text(size=text_size),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position='none')
  
  
  png(filename = paste('channagiri/plot_png/Combined_Enrichment_',chrom_list,'.png',sep=''), 
      width=15, 
      height = 12, 
      units = 'in',res = 300)
  grid.newpage()
  grid.draw(egg::ggarrange(p_up,p0,p_down, ncol = 1))
  dev.off()
}