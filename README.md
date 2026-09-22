# Brown_Spot_Project-
## Overview

This repository contains the scripts and organized analysis materials used for the Brown Spot Needle Blight (BSNB) project. The workflow supports the processing and analysis of fungal sequence data, with particular emphasis on the culturable internal transcribed spacer (ITS) region of symptomatic p. taeda needle. It moves from chromatogram-derived sequence inputs through sequence-quality review, ITS preparation and alignment, and downstream population-genetic summaries.

### Workflow at a glance

The analysis follows the sequence below. Individual scripts should be reviewed before execution because input names, directory paths, and parameters may have been tailored to the original dataset.
1. [Process chromatogram or sequence-read files](ab1processngscript.R) supports preparation of AB1-derived sequence data for later review and analysis.

2. [Prepare ITS data](WorkingscriptITS.R) and the (ITS_script/) directory contain the ITS-focused processing workflow, including preparation of sequence data for comparison across samples.

3.Organize sequence inputs. [Base-called sequence files](basecalled/.)

4.Review sequence quality. [Quality-control materials and intermediate outputs are organized in](QC/.)

5.[Create and curate alignments.](minmaalign.sh) is a shell-based alignment helper used within the sequence-processing workflow. Alignment settings and any required external software should be checked directly in the script before use.

6.[Run genome assembly with](flygeno.sh) preserves an additional command-line step used in the project workflow. 

7.[Calculate genetic-diversity and population-genetic metrics.](SydowiaGeneticdiversitymetrics.R) and the scripts in popgenscripts/ support downstream analysis of genetic variation.

8.[Generate project-level analyses and summaries] (BSNBdraft.R) contains additional R-based analysis, exploratory work, and/or figures used to synthesize the BSNB dataset.

### Software and requirements

The workflow uses R scripts and Unix shell scripts. A working installation of R and a POSIX-compatible shell environment are therefore required. The exact R packages and external command-line tools should be installed according to the library calls and commands in each script. Before running the workflow, confirm that all input and output paths match the local copy of the repository.

### Contact

For questions about the Brown Spot Needle Blight project, its data, or the intended use of the scripts, please contact the repository owner or the project research team
### lab project 
https://willoughbylab.auburn.edu/ 


