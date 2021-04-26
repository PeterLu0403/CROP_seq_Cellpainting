# CROP_seq_Cellpainting
This repository describes the procedures of processing the image datasets of CROPseq experiments using Cellpainting methods.
## Data preparation
### Features extraction
 The image datasets were taken on ZEISS Celldiscoverer 7 system. The methods of [CellPainting](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223290/) were used to extract the morphological features of individual cells. The [pipelines](Pipelines/) were created in CellProfiler (Version 3.1.9), which contains 3 parts, illumination correction, quality control, and features extraction. Also, the pipelines of generating individual images of each cells were created from the features extraction parts.

To extract the features, the CellProfiler were deployed in Linux environments in the Nectar Cloud (The National eResearch Collaboration Tools and Resources project). [Here](NectarSetup/), 8 virtual machines were set up to process the image datasets, each of which has 32 virtual CPUs and 64 GB of RAM. The process took a few days to finish and the analysis output a set of SQLite files. In total, the experiment contains 3 batches and each batch contains 5 plates of cells. 

Here is the [script](R/Read_sql_database_CP.R) of generating and merging *.csv files from the sqlite database. The [metadata](metadata/) of which gene was knockout is also added to the exported *.csv files.  

### Data preparation





### Transformation, normalization and batch-effect evaluation
## Data analysis
###