MSR
===
The MultiScale Representation (MSR) of Genomic Signals is collection of segmentations of a genomic signal at different spatial scales. Each segment in a scale-specific segmentation is scored for signal enrichment or depletion by the Significant Fold Change (SFC) score.

There are two main components to create the MSR. The third component is optional.
1. Creating the multiscale segmentation 
2. Computating of the SFC
3. Pruning to detect relevant yet non-redundant segment in the MSR

Below we describe these components using the Matlab msr_example_script.m. Additionally, we provide a script to create MSRs from ENCODE data and we provide Matlab Compiler Runtimes for Unix and Windows that enable users to run the MSR application without installing Matlab. Note that depending on the size of the data and parameters settings computing the MSR can be demanding in terms of memory and CPU time.


## msr_example_script.m
[msr_example_script.m](../master/msr_example_script.m)

The code consists of four parts:

1. The multiscale segmentation algorithm
2. Computation of the SFC
3. Pruning to detect relevant yet non-redundant segment in the MSR
4. A plotting routine to visualize the MSR

These four parts are run with an example genomic signal in the main script: msr_example_script.m.
Below we describe these four parts in more detail.


### 1. The multiscale segmentation algorithm
```bash
#create the multiscale segmentation using signal V1 and L scales
[KM,SegmentEnd] = MSS(V1,L);
#KM contains the hierarchy between segments and SegmentEnd contains the end positions of the segments
```

### 2. Computation of the SFC
```bash
#Computation of the SFC taking into account mappability map
Enrichment1 = SignificantFoldChange(SegmentEnd,V1,UM); 
#Computation of the SFC without mappability map or other background signal
Enrichment2 = SignificantFoldChange(SegmentEnd,V1); 
#Computation of the SFC comparing V2 versus V1 without mappability map
Enrichment3 = SignificantFoldChange(SegmentEnd,V2,UM,V1); 
#Computation of the SFC comparing V2 versus V1 with mappability map
Enrichment4 = SignificantFoldChange(SegmentEnd,V2,[],V1); 
#Enrichmentx contains the SFC scores for all the segments
```

### 3. Pruning to detect relevant yet non-redundant segment in the MSR
```bash
#Pruning the multiscale representation with parameters R and T, and only keeping segment with SFC larger than Eth
KeepP = Prune(Enrichment1,KM,SegmentEnd,T,R,Eth);
#KeepP contains binary values for all segments indicated whether they were kept (1) or pruned (0).
```
Additional information on the pruning:
The segment is only kept (i.e. not pruned) when it has a score higher that all of its children and all of its parents. To be able to detect segments at different scales a size constraint is introduced: The segment under investigation is not compared to all its children and parents, but only to those children with segment lengths larger than S∙R and parents with length smaller than S/R. Here, S is the length of the segment under investigation and R is parameter between 0 and 1. The default setting for R is 0.2, which means that children 5 times smaller than the segment and parents 5 times larger than the segment are not considered. There is one exception: when all children at scales 1 to n-1 have a lower score than the segment under investigation, all children are pruned, also those smaller than S∙R. When R is 0, there is effectively no size constraint, i.e. the segment under investigation is compared to all its children and parents. In that case, a genomic position is part of at most one segment that is not pruned. 
Further, a slack parameter, T, is introduced, which prevents higher-scale segments from breaking up into smaller ones. Specifically, it prevents segments from being pruned, because one of the children has a slightly better score. Specifically, the segment under investigation is kept (i.e. not pruned) when its score (denoted by X) multiplied by T, (i.e. X∙T) is larger than the scores of all its children, and its (unadjusted) score X is larger than the scores of all its parents. The default setting for T is 1.05.

### 4. A plotting routine to visualize the MSR
```bash
#Plotting the MSR for the genomic signal V1 in the range x from scale Lmin to Lmax using SFC scores from Enrichment2. The last argument can be empty or contain the pruning results (KeepP)
plotMSS(x,Lmin,Lmax,V1,Enrichment2,KM,SegmentEnd,[]);
```

## msr_encodedatapipeline.m
[msr_encodedatapipeline.m](../master/msr_encodedatapipeline.m)

This script takes a signal downloaded from ENCODE (in WIG format) and creates the MSR. To use the ENCODE signal unzip the WIG file in the [ENCODE data folder](../master/Data/ENCODE). The unique mappability map used in this pipeline is downloaded from [the Uniqueome website](http://grimmond.imb.uq.edu.au/uniqueome/downloads/). The MSR representation is output as a BED file, where each line represent a segment on a particular scale with the SFC score in the last column:
```bash
track name="MSR_H3K04ME3"				
chr1	52287910	52288310	scale10	0.2736
chr1	121996410	121996710	scale10	0.2909
chr1	153193710	153194110	scale10	0.0079
chr1	52820810	52821210	scale11	0.4899
chr1	63176410	63176910	scale11	1.3695
chr1	79773210	79773310	scale11	1.4212
chr1	90172910	90173010	scale11	2.9788
```

## Matlab Compiler Runtimes
For [Windows 64 bit](../master/MCR/Windows64) and [Unix 64 bit](../master/MCR/Unix64)

Here, we provide standalone executables for Unix (64bit) and Windows (64bit) that enable users to run the MSR  application without installing Matlab. To use these runtime components, first download the Matlab Compiler Runtime Release 2013a (8.1) from [the Matlab website](http://www.mathworks.com/products/compiler/mcr/). Install the MCR and set paths as described [here](http://www.mathworks.com/help/compiler/working-with-the-mcr.html) and alternatively [here] (../master/MCR/AdditionalFunctions/Working with the MCR - MATLAB & Simulink.pdf) (click 'View Raw').

Copy the appropriate folder ([this one for Windows](../master/MCR/Windows64),[this one for Unix](../master/MCR/Unix64)) to a local directory or server. There are executables to create a MSR for an individual signal as well as for a BedGraph file. (See [here](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) for information about this format. Other genomic signal formats, such as WIG can easily be transformed into a BedGraph using the binary utilities available via the [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/admin/exe/).) The paths to input and output data as well as parameter settings should be stored in a parameter file, which forms the input argument to the standalone executable. Below we provide a detailed description and examples.

### Standalone for an individual signal
Computing the MSR for a signal requires certain input data and (optional) parameters:
```bash
signal='signal.txt' #path to tab delimited text with numbers presenting the signal for which the MSR has to be computed. (required)
out='msr.txt'       #path to MSR output file. Columns: scale# segmentbegin segmentend SFC. (required)
background='bg.txt' #path to file with background mappability map (optional, default none)
versus='vs.txt'     #path to tab delimited text with numbers presenting the signal to which the original signal is compared with computing the SFC (in beta, optional, default none)
L=30                #number of scales (optional, default 25)
Pth=1e-6            #P-value threshold (optional, default 1e-6)
PruneFlag=-1        #Output Pruned MSR (1), complete MSR (-1) (big file) or complete MSR exluding segments with SFC=0 (optional, default 0)
DepletionFlag=0     #When pruning, also output depleted segments (1) or not (0) (optional, default 0)
Eth=0               #Threshold for including enriched or depleted segments in the pruned MSR (optional, default 0)
T=1                 #Pruning parameter T (optional, default 1.05)
R=0                 #Pruning parameter R (optional, default 0.2)
```
These parameters are stored in a parameter file, such as [this](../master/MCR/Windows64/data/parameters/parameterfile_signal_example1_Windows.txt). It is very important to make sure that the paths to the input data, i.e. the signals and background mappability map, as well as the output folder, exist and are correct. Make sure you change the paths in the parameter file to reflect the actual locations of these data on your local or network drive. 
The executable is called with the parameter file as an argument as explained below.

#### Windows
The exectuable is called [msr_runtime_SIGNAL.exe](../master/MCR/Windows64/msr_runtime_SIGNAL.exe). 
Run the executable from the command prompt like this:
```bash
msr_runtime_SIGNAL U:\matlab\MSRGIT\MCR\Windows64\data\parameters\parameterfile_signal_example1_Windows.txt
```
Specifically, the exectuable has one argument, which is the path to the parameter file. Make sure that the path to the parameter file is correct. It will be different from the example above. Run the executable from the directory containing the executable. Alternatively, make sure the executable is part of your Windows path, or specify the path to the exectable when calling it. 
See [here](../master/MCR/Windows64/RUNEXAMPLE_SIGNAL.txt) for the example file. See [here](../master/MCR/Windows64/data/parameters/parameterfile_signal_example1_Windows.txt) for an example parameter file. Example signals can be found ([here](../master/MCR/Windows64/data/in)). Make sure you unzip in.zip to use these example signals.


#### Unix
The exectuable is called [msr_runtime_SIGNAL](../master/MCR/Unix64/msr_runtime_SIGNAL). 
Use the the [wrapper](../master/MCR/Unix64/run_msr_runtime_SIGNAL.sh) to run the executable like this:
```bash
./run_msr_runtime_SIGNAL.sh /titan/cancerregulome9/workspaces/mcr/mcr/v81/ ./data/parameters/parameterfile_signal_example1_Unix.txt
```
Specifically, the wrapper shell sript has two arguments. First, the path to where the Matlab Compiler Runtimes is installed. Second, the path to the parameter file. Make sure that these two paths are correct. They will be different from the example above. Run the wrapper shell sript from the directory containing the wrapper shell sript. Alternatively, make sure the wrapper shell sript is part of your Unix path, or specify the path to the wrapper shell sript when calling it. 
See [here](../master/MCR/Unix64/RUNEXAMPLE_SIGNAL.txt) for the example file. Note that in this case the parameter file is the second argument and the first argument is the path to the MCR v81. See [here](../master/MCR/Unix64/data/parameters/parameterfile_signal_example1_Unix.txt) for an example parameter file. Example signals can be found ([here](../master/MCR/Unix64/data/in)). Make sure you unzip in.zip to use these example signals.

### Standalone for a BEDgraph file
Computing the MSR for a Bedgraph file requires certain input data and (optional) parameters:
```bash
folder='/tmp/'            #path to folder for temporary data (created if not existing, required) 
chrsizes='mm9cs.txt'      #path to tab delimited text file with chromosome sizes (required)
signal='signal.bedgraph'  #path to BEDgraph file for which the MSR should be computed (required)
out='out.bed'             #MSR output file in BED format. Columns: chr segmentbegin segmentend scale# SFC (required)
background='um.bedgraph'  #path to BEDgraph file with background mappability map (optional, default none)
versus='vs.bedgraph'      #path to BEDgraph file with signal to which the original signal is compared compared with computing the SFC (in beta, optional, default none)
IS=10                     #resolution in basepairs (optional, default 10)
L=30                      #number of scales (optional, default 25)
Pth=1e-6                  #P-value threshold (optional, default 1e-6)
PruneFlag=-1              #Output Pruned MSR (1), complete MSR (-1) (big file) or complete MSR exluding segments with SFC=0 (optional, default 0)
DepletionFlag=0           #When pruning, also output depleted segments (1) or not (0) (optional, default 0)
Eth=0                     #Threshold for including enriched or depleted segments in the pruned MSR (optional, default 0)
T=1                       #Pruning parameter T (optional, default 1.05)
R=0                       #Pruning parameter R (optional, default 0.2)
```
Note on the chromosome sizes file: Such a file can be easily generated with [UCSC fetchChromSizes utility](http://hgdownload.cse.ucsc.edu/admin/exe/). MSRs will be generated for chromosomes that are mentioned in the file with chromosome sizes and for which data is recorded in the BedGraph file. [Here](../master/MCR/Windows64/data/chrsizes/mm9cs.txt) is an example chromosome size file.

The parameters are stored in a parameter file, such as [this](../master/MCR/Windows64/data/parameters/parameterfile_BED_example1_Windows.txt). It is very important to make sure that the paths to the input data, i.e. the signals, chromosome sizes and background mappability map, as well as the output folder and temporary folder, exist and are correct. Make sure you change the paths in the parameter file to reflect the actual locations of these data on your local or network drive. 
The executable is called with the parameter file as an argument as explained below. 

#### Windows
The exectuable is called [msr_runtime_BED.exe](../master/MCR/Windows64/msr_runtime_BED.exe). 
Run the executable from the command prompt like this:
```bash
msr_runtime_BED U:\matlab\MSRGIT\MCR\Windows64\data\parameters\parameterfile_BED_example1_Windows.txt
```
Specifically, the exectuable has one argument, which is the path to the parameter file. Make sure that the path to the parameter file is correct. It will be different from the example above. Run the executable from the directory containing the executable. Alternatively, make sure the executable is part of your Windows path, or specify the path to the exectable when calling it.
See [here](../master/MCR/Windows64/RUNEXAMPLE_BED.txt) for the example file. See [here](../master/MCR/Windows64/data/parameters/parameterfile_BED_example1_Windows.txt) for an example parameter file. Example BedGraphs can be found ([here](../master/MCR/Windows64/data/in)). Make sure you unzip in.zip to use these example signals.

#### Unix
The exectuable is called [msr_runtime_BED](../master/MCR/Unix64/msr_runtime_BED). 
Use the the [wrapper](../master/MCR/Unix64/run_msr_runtime_BED.sh) to run the executable like this:
```bash
./run_msr_runtime_BED.sh /titan/cancerregulome9/workspaces/mcr/mcr/v81/ ./data/parameters/parameterfile_BED_example1_Unix.txt
```
Specifically, the wrapper shell sript has two arguments. First, the path to where the Matlab Compiler Runtimes is installed. Second, the path to the parameter file. Make sure that these two paths are correct. They will be different from the example above. Run the wrapper shell sript from the directory containing the wrapper shell sript. Alternatively, make sure the wrapper shell sript is part of your Unix path, or specify the path to the wrapper shell sript when calling it. 
See [here](../master/MCR/Unix64/RUNEXAMPLE_BED.txt) for the example file. Note that in this case the parameter file is the second argument and the first argument is the path to the MCR v81. See [here](../master/MCR/Unix64/data/parameters/parameterfile_BED_example1_Unix.txt) for an example parameter file. Example signals can be found ([here](../master/MCR/Unix64/data/in)). Make sure you unzip in.zip to use these example signals.

#### Creation of the runtime
The executables were created by running make_mcc.m ([Windows](../master/MCR/Windows64/make_mcc.m)) ([Unix](../master/MCR/Unix64/make_mcc.m)).

## Warranty Disclaimer and Copyright Notice
Theo Knijnenburg
Institute for Systems Biology
December 2013
Copyright (C) 2003-2013 Institute for Systems Biology, Seattle, Washington, USA.
The Institute for Systems Biology and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA



