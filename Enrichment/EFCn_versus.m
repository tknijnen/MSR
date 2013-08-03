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

%not finished
function S = EFCn_versus(O,m,s,Rr,G2r,pth,maxEFC,foldchange_or_additive);

if nargin<8;
    foldchange_or_additive = 1;
end
if nargin<7;
    maxEFC = 10;
end
if nargin<6;
    pth = 1e-6;
end

%% Enrichment
Zs = -norminv(pth);
Os = Zs.*s + m;
if foldchange_or_additive==1;
    E = Rr./Os;
    E(E<1|isnan(E))=1;
else
    E = Rr-Os;
    E = E./G2r;
    E(E<0|isnan(E))=0;
end

%% Depletion
Zs = norminv(pth);
Os = Zs.*s + m;
if foldchange_or_additive==1;
    D = Rr./Os;
    D(D>1|isnan(D)|Os<0)=1;
else
    D = Rr-Os;
    D = D./G2r;
    D(D>0|isnan(D))=0;
end


%% Score
if foldchange_or_additive==1;
    S = ones(size(E));
    S(E>1) = E(E>1);
    S(D<1) = D(D<1);
    S = log2(S);
    S(S>maxEFC) = maxEFC;
    S(S<-maxEFC) = -maxEFC;
else
    S = zeros(size(E));
    S(E>0) = E(E>0);
    S(D<0) = D(D<0);
end


%% Score =
S = ones(size(E));
S(E>1) = E(E>1);
S(D<1) = D(D<1);
S = log2(S);
