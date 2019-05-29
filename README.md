#TEMPLATE MANUSCRIPT#
===============================================

# Folders #
### Note: please feel to change the organization of your manuscript folder to your preference; so long it's organized and clearly documented. ## 
## data ##	
This is where you would store intermediate data produced that is specific to this manuscript. All other more general data should be stored in the Huang_lab_data folder and sourced from there. 
## doc ##		
This is where you would store and update the manuscript, figure, and supplementary documents.
>fishermen_manuscript.docx: read and make sure we are clear on the key components of the "fishing trip" before we commit to it. Often you may find it useful to continue to have a draft manuscript document to keep track of your Methods/Key Findings/Flow of your discovery. 
## reference ##
If you download reference as PDFs this is a good place to store it.
## analysis	## 
This directory should be git version-controlled (https://github.com/Huang-lab). For an analysis manuscript one way of organization would be putting all scripts related to one kind of analysis into one folder. Software tools should always have their independent tool folder. 
# Example scripts that you may use or delete #
>global_aes_out.R: used to define some color codes for other R plotting codes. ggplot (https://ggplot2.tidyverse.org/index.html) coupled with RStudio is one top killer app for computational biology plotting and we highly encourage using that for plotting once you have tabular data. 
>some_stat_functions.R: some demo of stat functions you may or may not use. Multiple higher level scripting languages (ex. python/R/matlab) offer strong solutions so that choice is more up to you. 