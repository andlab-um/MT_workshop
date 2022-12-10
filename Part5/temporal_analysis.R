  title: "MTprocessing"
  author: "Xinyi Xu"
  output:
    pdf_document: 
    fig_height: 4
    fig_width: 6
---
  # Preparations
  ## Load libraries

library(abind)
library(mousetrap)
library(ggplot2)
library(dplyr)
library(readbulk)
library(afex)

## Custom ggplot2 theme

theme_set(theme_classic()+ 
            theme(
              axis.line = element_line(colour = "black"),
              axis.ticks = element_line(colour = "black"),
              axis.text = element_text(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA)
            ))
rm(list=ls())
options(width=90)


# Import new dataset from block4 and 7


setwd("/Users/orlacamus/Desktop/cityU_code/IAT")

iat <- read_bulk("sampledata", fun=read_mt, extension=".mt")


iat <- mt_import_wide(iat)
mt_1 <- mt_remap_symmetric(iat)
mt_1 <- mt_align_start_end(mt_1)
mt_1 <- mt_time_normalize(mt_1)
mt_data <- mt_measures(mt_1)




##add label congruent/incongruent as condition1
condition1 <- rep(1,400)
condition1[1:200] <- 'congruent'
condition1[201:400] <- 'incongruent'
condition1 <- rep(condition1,2)


mt_data[['data']]<-cbind(condition1,mt_data[['data']])
mt_data[['measures']]<-cbind(condition1,mt_data[['measures']])




## Filtering
# Only keep trials with correct answers and RT<5000
mt_data <- mt_subset(mt_data,RT<=5000,check='data')

mt_data <- mt_subset(mt_data,error==0,check='data')


########################################################data preprocess

# Analysis

## Aggregate trajectories

# Fig. 4
mt_plot_aggregate(mt_data, use="tn_trajectories",
                  x="xpos", y="ypos", 
                  color="condition1", subject_id="subjID")+
  scale_color_manual(values=c("darkorange","steelblue"))


## Calculate measures

# Calculate velocity and acceleration
mt_data <- mt_derivatives(mt_data)

# Calculate sample entropy
mt_data <- mt_sample_entropy(mt_data, use="tn_trajectories")


## Curvature

### Aggregate analyses

# Aggregate MAD values per participant and condition
agg_mad <- mt_aggregate_per_subject(mt_data, 
                                    use_variables="MAD", use2_variables="condition1",
                                    subject_id="subjID")

# Compare aggregated MAD values
t.test(MAD~condition1,data=agg_mad,paired=TRUE)

# Calculate descriptives
agg_mad %>% 
  group_by(condition1) %>% 
  summarise_at("MAD",.funs=c("mean","sd")) %>%
  as.data.frame()


### Trial level analyses

# Create data.frame that contains the
# trial variables and mouse-tracking indices
results <- merge(mt_data$data, mt_data$measures, by="mt_id")



## Temporal analyses

### Average x positions

# Plot aggregate time-normalized x-positions (Fig. 6)
mt_plot_aggregate(mt_data, use="tn_trajectories",
                  x="steps", y="xpos", color="condition1",
                  subject_id="subjID", points=TRUE)+
  scale_color_manual(values=c("darkorange","steelblue"))
# Aggregate time-normalized trajectories per condition
# separately per participant
av_tn_trajectories <- mt_aggregate_per_subject(mt_data,
                                               use="tn_trajectories", use2_variables="condition1",
                                               subject_id="subjID")
# Paired t-tests on coordinates
xpos_t_tests <- 
  with(av_tn_trajectories,
       sapply(unique(steps),function(i){
         t.test(xpos[condition1=="congruent" & steps==i],
                xpos[condition1=="incongruent" & steps==i],
                paired = TRUE)$p.value})
  )
# Retrieve all significant t-tests
which(xpos_t_tests<.05)





