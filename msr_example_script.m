% Theo Knijnenburg
% Institute for Systems Biology
%
% August 2013
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warranty Disclaimer and Copyright Notice
% 
% Copyright (C) 2003-2013 Institute for Systems Biology, Seattle, Washington, USA.
% 
% The Institute for Systems Biology and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
% 
% This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
% 
% You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;
addpath(genpath(pwd));

%% Path to fast convolution
if ispc
    addpath([pwd '\Win_CONVNFFT_Folder']);rmpath([pwd '\Unix_CONVNFFT_Folder']);
elseif isunix
    addpath([pwd '/Unix_CONVNFFT_Folder']);rmpath([pwd '/Win_CONVNFFT_Folder']);
end

%% Load genomic signals
load ExampleGenomicSignals;
%V1      RNA Polymerase II ChIP-seq signal of unstimulated murine bone marrow macrophage (BMM) on chromosome 1, 1Mb samples at 10bp resolution
%V2      Same but for BMM after 2 hours of LPS stimulation
%B       Control ChIP-Seq signal of BMMs with immunoglobulin G derived from rabbits that were not immunized with specific target antigens.
%UM      Mappability landscape for ChIP-seq fragments (see paper) 

%% Segmentation
L = 35; %number of scales
[KM,SegmentEnd] = MSS(V1,L);

%% Enrichment computation
clc
Enrichment1 = SignificantFoldChange(SegmentEnd,V1,UM); %Enrichment taking into account mappability map
Enrichment2 = SignificantFoldChange(SegmentEnd,V1); %No mappability map or other background signal
Enrichment3 = SignificantFoldChange(SegmentEnd,V2,UM,V1); %V2 versus V1 with mappability map
Enrichment4 = SignificantFoldChange(SegmentEnd,V2,[],V1); %V2 versus V1 without mappability map

%% Pruning
Eth = 0;
T = 1.5;
R = 0.2;
KeepP = Prune(Enrichment1,KM,SegmentEnd,T,R,Eth);
T = 1.05;
R = 0;
KeepA = Prune(Enrichment1,KM,SegmentEnd,T,R,Eth);

%% Plotting
x = [400000 600000];%range to plot
Lmin = 5;
Lmax = 30;
plotMSS(x,Lmin,Lmax,V1,Enrichment2,KM,SegmentEnd,[]);








