% Theo Knijnenburg
% Institute for Systems Biology
%
% November 2013
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
rmpath(genpath([pwd '/MCR/']));

%% Path to fast convolution
if ispc
    addpath([pwd '\Win_CONVNFFT_Folder']);rmpath([pwd '\Unix_CONVNFFT_Folder']);
elseif isunix
    addpath([pwd '/Unix_CONVNFFT_Folder']);rmpath([pwd '/Win_CONVNFFT_Folder']);
end

%% load Chromosome information
A = importdata([pwd '/Data/MouseGenome/GenomeChromSizesM37.txt']);
ChrName = A.textdata;
ChrSize = A.data;
% removing random chromosomes
I = [];
for n = 1:length(ChrName);
    if ~isempty(strfind(ChrName{n},'random'));
        I = [I n];
    end
end
% exclude mitochondrial chromosome
I = [I find(strcmp(ChrName,'chrM'))];
ChrName(I) = [];
ChrSize(I) = [];
clear A I n

%% Parameters
IS = 10;        %sample size in bp's
L = 50;         %number of scales
Pth = 1e-6;     %P-value threshold for SFC computation
Eth = 0;        %SFC threshold in pruning
T = 1.05;       %Slack parameter in pruning
R = 0.2;        %Size constraint in pruning

%% Read in WIG file (BigWig downloaded from ENCODE project at UCSC, UCSC Accession: wgEncodeEM002659, GEO sample accession:GSM1000065, and transformed into wig files using UCSC’s BigWigToWig binary utility)
wigfile = [pwd '/Data/ENCODE/H3K04ME3.WIG'];
fid = fopen(wigfile);
E = textscan(fid,'%s%n%n%n','Delimiter','\t','HeaderLines',0,'CommentStyle','#');
fclose(fid);
for c = 1:length(ChrName);
    S = zeros(ChrSize(c),1);
    pos = find(strcmp(E{1},ChrName{c}));
    for p = 1:length(pos);
        S((E{2}(pos(p))+1):E{3}(pos(p))) = E{4}(pos(p));
    end
    S = filter(ones(1,IS)/IS,1,S);
    X = IS:IS:ChrSize(c);
    S = S(X);
    save([pwd '/Data/ENCODE/Signal_' ChrName{c} '.mat'],'S');
end

%% Read in uniqueome mappability signal (Wig.gz downloaded from http://grimmond.imb.uq.edu.au/uniqueome/downloads/)
wigfile = [pwd '/Data/Uniqueome/mm9_uniqueome.unique_starts.base-space.35.1.Wig'];
tmpfile = [pwd '/Data/Uniqueome/tmp.wig'];

for c = 1:length(ChrName);
    display(['Creating temp file for ' ChrName{c} '...']);
    [status,result] = unix(['tail -n +1 ' wigfile ' | awk ''$1=="' ChrName{c}    '"'' > ' tmpfile]);
    display(['Reading temp file for ' ChrName{c} '...']);
    fid = fopen(tmpfile);
    E = textscan(fid,'%s%n%n%n','Delimiter','\t','HeaderLines',0,'CommentStyle','#');
    fclose(fid);
    display(['Running signal for ' ChrName{c} '...']);
    S = zeros(ChrSize(c),1);
    pos = find(strcmp(E{1},ChrName{c}));
    for p = 1:length(pos);
        S((E{2}(pos(p))+1):E{3}(pos(p))) = E{4}(pos(p));
    end
    S = filter(ones(1,IS)/IS,1,S);
    X = IS:IS:ChrSize(c);
    S = S(X);
    display(['Saving signal for ' ChrName{c} '...']);
    save([pwd '/Data/Uniqueome/UnMap_' ChrName{c} '.mat'],'S');
    [status,result] = unix(['rm ' tmpfile]);
end


%% Loop chromosomes to create MSR
for c = 1:length(ChrName);
    load([pwd '/Data/Uniqueome/UnMap_' ChrName{c} '.mat'],'S');
    UM = S;
    load([pwd '/Data/ENCODE/Signal_' ChrName{c} '.mat'],'S');
    [KM,SegmentEnd] = MSS(S,L);
%     Enrichment = SignificantFoldChange(SegmentEnd,S,[],[],Pth); 
    Enrichment = SignificantFoldChange(SegmentEnd,S,UM,[],Pth); 
    Keep = Prune(Enrichment,KM,SegmentEnd,T,R,Eth);
    save([pwd '/Data/ENCODE/MSR_' ChrName{c} '.mat'],'KM','SegmentEnd','Enrichment','Keep');
end

%% Create output BED file
bedfile = [pwd '/Data/ENCODE/MSR_H3K04ME3.bed'];
PruneFlag = 1;   %output only enriched segments that are not pruned
fid = fopen(bedfile,'w+');
fprintf(fid,['track name="MSR_H3K04ME3"' '\n']);
for c = 1:length(ChrName);
    display(['Running ' ChrName{c} '...']);
    load([pwd '/Data/ENCODE/MSR_' ChrName{c} '.mat'],'SegmentEnd','Enrichment','Keep');
    str = ChrName{c};
    fs = [];
    for s = 1:length(str);
        fs = [fs '%c'];
    end
    fs = [fs '\t%d\t%d\tscale%d\t%6.4f\n'];
    for nn = 1:L;
        fprintf(1,'.');
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 IS*SegmentEnd{nn}(1:SL-1)+1],IS*SegmentEnd{nn});
        Segments(1,:) = Segments(1,:)-1;
        if PruneFlag;
            pos = find(Keep{nn}==1);
        else
            pos = 1:SL;
        end
        LP = length(pos);
        if LP>0;
            fprintf(fid,fs,[ones(LP,1)*str Segments(:,pos)' nn*ones(LP,1) Enrichment{nn}(pos)]');
        end
    end
end
fclose(fid);



    

