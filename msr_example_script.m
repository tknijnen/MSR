clear all;close all;clc;

%% Windows
addpath([pwd '\Win_CONVNFFT_Folder']);rmpath([pwd '\Unix_CONVNFFT_Folder']);

%% Unix
% addpath([pwd '/Unix_CONVNFFT_Folder']);rmpath([pwd '/Win_CONVNFFT_Folder']);

%% Signal
load V

%% Segmentation
L = 35; %number of scales
[KM,SegmentEnd,Enrichment] = MSS(V,L);

%% Plotting
x = [10000 20000];%range to plot
Lmin = 5;
Lmax = 30;
plotMSS(x,Lmin,Lmax,V,Enrichment,KM,SegmentEnd);