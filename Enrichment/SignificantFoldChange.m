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

function Enrichment = SignificantFoldChange(SegmentEnd,V1,UM,V2)

%% Scenarios
BackgroundFlag = 0;
VersusFlag = 0;
if exist('UM','var');
    if ~isempty(UM);
        BackgroundFlag = 1;
    end
end
if exist('V2','var');
    if ~isempty(V2);
        VersusFlag = 1;
    end
end

%% Parameters
minLsize1 = 0;
minLsize2 = 10;
th = 1e-5;
pth = 1e-6;
maxEFC = 10;                        %only in versus mode with backgroundflag
foldchange_or_additive = 1;         %only in versus mode with backgroundflag

%% Loop scales depending on scenario
if VersusFlag==0&BackgroundFlag==0;
    
    CV = cumsum(V1);
    UM = ones(size(V1));
    CL = cumsum(UM);
    
    mu = mean(V1);
    sigma2 = var(V1);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
        fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        
        %% Computing enrichment
        Enrichment{nn} = EFCn(O,Le.*mu,sqrt(G2r.*sigma2),pth);
        
    end
    
elseif VersusFlag==0&BackgroundFlag==1;
    
    MLV = max(UM);
    V1 = V1./MLV;
    UM = UM./MLV;
    CL = cumsum(UM);
    CV = cumsum(V1);
    p = CV(end)/CL(end);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
        fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        G1 = CV(end)*ones(size(O));
        T = CL(end)*ones(size(O));
        
        %% Computing enrichment
        %     Enrichment{nn} = EFC(T,G1,G2,O,pth,th)';
        Enrichment{nn} = EFC_inflatedvariance(T,G1,G2,G2r,O,pth,th);
        
    end
    
elseif VersusFlag==1&BackgroundFlag==0;
    
    CV = cumsum(V1);
    UM = ones(size(V1));
    CL = cumsum(UM);
    CR = cumsum(V2);
    
    mu = mean(V1);
    sigma2 = var(V1);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
        fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        Rr    = (CR(Segments(2,:))-CR(Segments(1,:))+V2(Segments(1,:)));
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        
        %% Computing enrichment
        Enrichment{nn} = EFCn(Rr,O,sqrt(G2r.*sigma2),pth);
        
    end
    
elseif VersusFlag==1&BackgroundFlag==1;
    
    MLV = max(UM);
    V1 = V1./MLV;
    UM = UM./MLV;
    V2 = V2./MLV;
    CL = cumsum(UM);
    CV = cumsum(V1);
    p = CV(end)/CL(end);
    CR = cumsum(V2);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
        fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        Rr    = (CR(Segments(2,:))-CR(Segments(1,:))+V2(Segments(1,:)));
        
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        
        %% Computing enrichment
        Enrichment{nn} = EFC_versus_iv(G2,G2r,O,Rr,pth,th,maxEFC,foldchange_or_additive);
    end
end






