MSR
===
The MultiScale Representation (MSR) of Genomic Signals is collection of segmentations of a genomic signal at different spatial scales. Each segment in a scale-specific segmentation is scored for signal enrichment or depletion by the Significant Fold Change (SFC) score.

There are two main components to create the MSR. The third component is optional.
1. Creating the multiscale segmentation 
2. Computating of the SFC
3. Pruning to detect relevant yet non-redundant segment in the MSR

Below we describe these components using the MATLAB msr_example_script.m. Additionally, we provide a script to create MSRs from ENCODE data and we provide Matlab Compiler Runtimes for Unix and Windows that enable users to run the MSR application without installing MATLAB. 


## msr_example_script.m
[msr_example_script.m](../msr_example_script.m)

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

## Warranty Disclaimer and Copyright Notice
Theo Knijnenburg
Institute for Systems Biology
August 2013
Copyright (C) 2003-2013 Institute for Systems Biology, Seattle, Washington, USA.
The Institute for Systems Biology and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA



