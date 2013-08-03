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

function [KM,SegmentEnd] = MSS(V,L);

time1 = clock;
display(['Running ']);

%% Gaussian scale space and window parameters
N = length(V);
deltatau = .5*log(2);
Pmin = 1e-3;                              %Cutoff of Gaussian window   -   will substantially effect window size in Gaussian Scale Space
kmin = (1-exp(-2*deltatau)).^(-1/2);
kn = ceil(kmin);

%% Starting nodes
Kids = unique([1+find(diff(V)~=0)' N]);   %default        %resolution of Kids will substantially effect the Linkage Loop            
% Kids = 1:N;                             %better/nicer   %resolution of Kids will substantially effect the Linkage Loop
DImax = max(V)-min(V);
w1 = 1;
nTZ = min(find(V~=0))-1;                %Number of trailing zeros
nEZ = min(find(flipud(V)~=0))-1;        %Number of ending zeros

%% Matrices
M = zeros(N,2);
M(:,2) = V;
sigma = ones(1,L);
sigma(1) = 1; 
r = zeros(1,L);
KM = cell(L,1);
SegmentEnd = cell(L,1);
SegmentEnd{1} = Kids;

for n = 2:L

    %% Gaussian scale space
    M(:,1) = M(:,2);
    sigma(n) = exp((n-1)*deltatau);
    X = norminv(Pmin,0,sigma(n));
    window = normpdf(round(X):-round(X),0,sigma(n));
    windowSize = length(window);
    window = window/sum(window);
    tic;
    if windowSize<1250;
        M(:,2) = conv2(V,window','same');
    else
        M(:,2) = convnfft(V',window,'same');
    end
    display(['Level : ' num2str(n) '   Scale : ' num2str(round(sigma(n))) '   WindowSize : ' num2str(windowSize) '   CPUtime : ' num2str(toc) ' seconds']);

    %% Search Volume
    if n==2;
        r(n) = sigma(n);
        r(n) = ceil(kmin*r(n));
    else
        r(n) = sqrt(sigma(n)^2-sigma(n-1)^2);
        r(n) = ceil(kmin*r(n));
    end
    Dcp = -r(n):r(n);
    D = exp((Dcp.^2)./(-2*(sigma(n)^2-sigma(n-1)^2)))./exp(((.5*sigma(n)).^2)./(-2*(sigma(n)^2-sigma(n-1)^2)));
    J = abs(Dcp)>.5*sigma(n);
    D(~J) = 1+1e-6*D(~J);
    
    %% Linkage Loop
    UK = unique(Kids);   %unique kids to linked to parent
    K = length(UK);
    P = zeros(size(UK));            %parent matrix
    lDcp = length(Dcp);
    rn = r(n);

    tic;
    if lDcp<50;

        % over window
        TempScore = zeros(1,K);
        Score = zeros(1,K);
        for w = 1:lDcp;
            J = UK+Dcp(w)>=1&UK+Dcp(w)<=N;
            TempScore(J) = D(w)*(1-abs(M(UK(J),1)-M(UK(J)+Dcp(w),2))/DImax);
            I = TempScore>Score;  %better parent was found
            I = J&I;
            Score(I) = TempScore(I); %updating score
            P(I) = UK(I)+Dcp(w);  %updating parent
        end
        P = sort(P); %prevent overlapping segments
        SegmentEnd{n} = [SegmentEnd{n-1}(diff(P)~=0) N];
        UP = [P(diff(P)~=0) P(end)];
        Gp = zeros(N,1);
        Gp(UP) = [SegmentEnd{n}(1) diff(SegmentEnd{n})];
        Gp(UP(1)) = Gp(UP(1))-nTZ;
        Gp(UP(end)) = Gp(UP(end))-nEZ;
        Gpmax = max(Gp);
        for w = 1:lDcp;
            J = UK+Dcp(w)>=1&UK+Dcp(w)<=N;
            TempScore(J) = D(w)*((1-abs(M(UK(J),1)-M(UK(J)+Dcp(w),2))/DImax)+(w1*Gp(UK(J)+Dcp(w))/Gpmax));
            I = TempScore>Score;  %better parent was found
            I = J&I;
            Score(I) = TempScore(I); %updating score
            P(I) = UK(I)+Dcp(w);  %updating parent
        end
        P = sort(P); %prevent overlapping segments
        SegmentEnd{n} = [SegmentEnd{n-1}(diff(P)~=0) N];

    else

        % over kids
        for k = 1:length(UK);
            px1 = max(1,UK(k)+Dcp(1));
            px2 = min(N,UK(k)+Dcp(lDcp));
            x1 = 1    + (rn - (UK(k)-px1));
            x2 = lDcp - (rn - (px2  -UK(k)));
            [i,j] = max(D(x1:x2).*(1-abs(M(UK(k),1)-M(px1:px2,2))/DImax)');
            P(k) = px1+j-1;
        end
        P = sort(P); %prevent overlapping segments
        SegmentEnd{n} = [SegmentEnd{n-1}(diff(P)~=0) N];
        UP = [P(diff(P)~=0) P(end)];
        Gp = zeros(N,1);
        Gp(UP) = [SegmentEnd{n}(1) diff(SegmentEnd{n})];
        Gp(UP(1)) = Gp(UP(1))-nTZ;
        Gp(UP(end)) = Gp(UP(end))-nEZ;
        Gpmax = max(Gp);
        for k = 1:length(UK);
            px1 = max(1,UK(k)+Dcp(1));
            px2 = min(N,UK(k)+Dcp(lDcp));
            x1 = 1    + (rn - (UK(k)-px1));
            x2 = lDcp - (rn - (px2  -UK(k)));
            [i,j] = max(D(x1:x2).*((1-abs(M(UK(k),1)-M(px1:px2,2))/DImax)+(w1*Gp(px1:px2)/Gpmax))');
            P(k) = px1+j-1;
        end
        P = sort(P); %prevent overlapping segments
        SegmentEnd{n} = [SegmentEnd{n-1}(diff(P)~=0) N];

    end
    display(['Level : ' num2str(n) '   Scale : ' num2str(round(sigma(n))) '   Number of Kids : ' num2str(K) '   CPUtime : ' num2str(toc) ' seconds']);
    KM{n} = cat(1,UK,P);
    Kids = P;

    display(['Level : ' num2str(n) '   Total Time : '  num2str(round(etime(clock,time1))) ' seconds']);
    display('___');
end
display(' ');
display(['Total Time : '  num2str(round(etime(clock,time1))) ' seconds']);





