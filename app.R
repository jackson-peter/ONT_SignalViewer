###############################
####### VERSION 2.1 ###########
###############################
#test
library(shiny)
library(tidyverse)
library(data.table)
library(DT)
library(bedr)
library(shinythemes)
library(shinydashboard)
library(shinyFiles)
library(ggsci)
library(ggpubr)
library(ggtext)
library(Biostrings)
library(BSgenome)

##################################################################################################################################################################### 
####################################################### GLOBALS #####################################################################################################
##################################################################################################################################################################### 
HTSLIB_PATH = "/biotools/htslib/1.9/bin/"

# If new files, just add symlinks in corresponding folder
## TO DO: modify architecture if several projects for same species.
kmer_dirs=list.dirs(path = "/home/jpeter/ShinyApps/app_f5c_v2/data", full.names = T, recursive = F)
names(kmer_dirs)=lapply(kmer_dirs, basename)

##################################################################################################################################################################### 
####################################################### Functions ###################################################################################################
##################################################################################################################################################################### 

# FUNCTION: get ROI string from chr, start, end
get.ROI.string <- function(chromosome, start, end){
  region_tabix_format <- paste0(chromosome, ":", start, "-", end)
  print(region_tabix_format)
  return(region_tabix_format)
}
# / get ROI string

# FUNCTION: Get chr, start & end from ROI string
get.coords.from.ROI.string <- function(ROI_string){
  split <- str_split(ROI_string,":")
  reslist <- list(chr = split[[1]][1],
                  start = as.numeric(str_split(split[[1]][2], '-')[[1]][1]),
                  end= as.numeric(str_split(split[[1]][2], '-')[[1]][2]))
  print(reslist)
  return(reslist)
}
# / get coords from ROI string

# FUNCTION: get the reverse complement of a sequence
get.ReverseComp <- function(sequence) {
  comp=stringi::stri_reverse(chartr("NATGCatgcn-", "NTACGTACGN-", sequence ))
  return(comp)
}
# / get ReverseComp

# FUNCTION: read bed
read.bed <- function(bedfile, header=FALSE) {
  bed_df <- read_tsv(bedfile,col_names = header)
  cnames <- colnames(bed_df)
  cnames[1] <- "chromosome"
  cnames[2] <- "start"
  cnames[3] <- "stop"
  colnames(bed_df) <- cnames
  return(bed_df)
}
# / read.bed

# FUNCTION: read tabix
read.polished.tabix <- function(file, ROI, ref_fasta, header = FALSE) {
  print(ROI)
  print(file)
  ROI_vals <- get.coords.from.ROI.string(ROI)
  ROI_chrom=ROI_vals$chr
  ROI_start=ROI_vals$start
  ROI_end=ROI_vals$end
  
  ref_seq = readDNAStringSet(ref_fasta)
  names(ref_seq) = gsub(" .*", "", names(ref_seq))
  
  ROI_seq <- data.table(ref_nucl = as.vector(subseq(ref_seq[[ROI_chrom]], ROI_start, ROI_end)),
                          chromosome = ROI_chrom,
                          ref_position = ROI_start:ROI_end)
  
  tabix_colnames <- c("chromosome", "ref_position", "ref_kmer", "read_index","strand", "event_index", "event_level_mean", "event_stdv", "event_length", "model_kmer",
                      "model_mean", "model_stdv", "standardized_level", "raw_signal_start", "raw_signal_end")
  
  if(header) {
    dt <- fread(cmd = paste(file.path(HTSLIB_PATH, "tabix"), file, ROI, "-h"))
  } else {
    dt <- fread(cmd = paste(file.path(HTSLIB_PATH, "tabix"), file, ROI), col.names = tabix_colnames)
  }

  dt <- dt %>%
    mutate(kmer_expected = case_when(get.ReverseComp(model_kmer) == ref_kmer ~ model_kmer, model_kmer == ref_kmer ~ model_kmer)) %>%
    filter(model_kmer != "NNNNNN") %>%
    mutate(zscore = (model_mean - event_level_mean) / model_stdv) %>%
    left_join(ROI_seq, by=c("ref_position", "chromosome"))
  
  return(dt)
}
# / Read tabix

############################################################################################################################################################## 
####################################################### UI ###################################################################################################
############################################################################################################################################################## 

ui <- dashboardPage(
  
  header <- dashboardHeader(title="Event Viewer"),
  
  sidebar <-  dashboardSidebar("Parameters",
                               #minified = FALSE,
                               collapsed = FALSE,
                               width="400px",
                               sidebarMenu( id="tabs",
                                            selectInput("species", "Select a species:", kmer_dirs),
                                            menuItem("ROI_selector", tabName = "Select ROIs", icon=icon("gear"),
                                                     h5("ROI coordinates"),
                                                     textInput(inputId = 'ROI',
                                                               label = "Enter ROI (Chr:start-end)",
                                                               value = ""),
                                                     hr()
                                            ), #/menuitem
                                            actionButton(inputId = "SubmitROIs",
                                                         label = "Submit ROIs",
                                                         style="color: #fffc300; background-color: #e95420; border-color: #c34113;
                                                         >border-radius: 10px; >border-width: 2px")
                                            
                               ) #/sidebarmenu
                               
  ), # dashboardsidebar
  body <-  dashboardBody(
    div(class = "span", tabsetPanel(
      tabPanel("Signal Viewer",plotOutput("Event_lvl_mean")),
      tabPanel("Signal Distribution",plotOutput("Event_lvl_bp")),
      tabPanel("Dwell Time", plotOutput("Event_dwell_t")),
      tabPanel("Zscore", plotOutput("Zscore"))
      ), # / tabsetpanel
    ) # /div
  ) # /dashboardbody
  
) # /dashboardpage

##################################################################################################################################################################
####################################################### SERVER ###################################################################################################
##################################################################################################################################################################

server <- function(input, output, session) {
  
  # f5c_data is the data table of the ROI
  f5c_data <- eventReactive(input$SubmitROIs,{
    print("building DF")
    ROI_vals <- get.coords.from.ROI.string(input$ROI)
    ROI_chrom=ROI_vals$chr
    ROI_start=ROI_vals$start
    ROI_end=ROI_vals$end
    
    ref_fasta = list.files(path=input$species, pattern="fa|fas|fasta", full.names = T)
    kmer_files <- list.files(path=file.path(input$species, "kmer_files"), pattern="*kmers_sorted.tsv.gz$", full.names = T, recursive = F)
    names(kmer_files)=lapply(kmer_files,basename)
    
    ROI_DF <- rbindlist(lapply(kmer_files, read.polished.tabix, ROI=input$ROI, ref_fasta), idcol = "origin")
    
    ROI_DF$ROI_name=input$ROI
    print(head(ROI_DF))
    return(ROI_DF)
  })
  
  ########## OUTPUTS GENERATION ##########
  ### Tab1 ###
  # Graph showing signal levels
  output$Event_lvl_mean <- renderPlot({
    req(input$SubmitROIs)

    dt <- f5c_data()
    labels_vect <-  pull(dt%>% distinct(chromosome, ref_position, ref_nucl)%>% select(ref_nucl), ref_nucl)

    ggplot(dt, aes(x=ref_position, y=event_level_mean, group=read_index, color=as.factor(read_index), alpha=0.2)) +
      facet_wrap(~origin, nrow = length(unique(dt$origin))) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks = seq(min(dt$ref_position),max(dt$ref_position)), labels=labels_vect)+
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_markdown())+
      ggtitle(unique(dt$ROI_name))

  })
  
  ### Tab2 ###
  # Graph showing signal distribution (with violin plots)
  output$Event_lvl_bp <- renderPlot({
    req(input$SubmitROIs)

    dt <- f5c_data()
    labels_vect <-  pull(dt%>% distinct(chromosome, ref_position, ref_nucl)%>% select(ref_nucl), ref_nucl)

    ggplot(dt, aes(x=as.factor(ref_position), y=event_level_mean, group=interaction(as.factor(ref_position), origin), fill=origin, alpha=0.2)) +
      geom_violin() +
      stat_compare_means(label =  "p.signif") +
      scale_x_discrete(breaks = seq(min(dt$ref_position),max(dt$ref_position)), labels = labels_vect)+
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_markdown())+
      ggtitle(unique(dt$ROI_name))

  })
  
  ### Tab3 ###
  # Graph showing dwell time (with violin plots)
  output$Event_dwell_t <- renderPlot({
    req(input$SubmitROIs)
    
    dt <- f5c_data() %>%
      mutate(dwell_time=raw_signal_end-raw_signal_start)
    
    labels_vect <-  pull(dt%>% distinct(chromosome, ref_position, ref_nucl)%>% select(ref_nucl), ref_nucl)
    
    ggplot(dt, aes(x=as.factor(ref_position), y=dwell_time, group=interaction(as.factor(ref_position), origin), fill=origin, alpha=0.2)) +
      geom_violin() +
      stat_compare_means(label =  "p.signif") +
      scale_x_discrete(breaks = seq(min(dt$ref_position),max(dt$ref_position)), labels = labels_vect)+
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_markdown())+
      ggtitle(unique(dt$ROI_name))
  })
  
  ### Tab4 ###
  # Z-scores
  output$Zscore <- renderPlot({
    req(input$SubmitROIs)
    dt <- f5c_data()
    
    labels_vect <-  pull(dt%>% distinct(chromosome, ref_position, ref_nucl)%>% select(ref_nucl), ref_nucl)
    ggplot(dt, aes(x=as.factor(ref_position), y=zscore, group=interaction(as.factor(ref_position), origin), fill=origin, color=origin, alpha=0.2)) +
      geom_violin() +
      stat_compare_means(label =  "p.signif") +
      geom_jitter() +
      scale_x_discrete(breaks = seq(min(dt$ref_position),max(dt$ref_position)), labels = labels_vect)+
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_markdown())+
      ggtitle(unique(dt$ROI_name))
  })

  
}

### Run Application ###
shinyApp(ui, server)
