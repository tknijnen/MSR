function [KM,SegmentEnd,Enrichment] = MSS(V,L);

time1 = clock;
display(['Running ']);

%% parameters
N = length(V);
deltatau = .5*log(2);
Pmin = 1e-3;            %will substantially effect window size in Gaussian Scale Space
kmin = (1-exp(-2*deltatau)).^(-1/2);
kn = ceil(kmin);
%resolution of Kids will substantially effect the Linkage Loop
Kids = unique([1+find(diff(V)~=0)' N]);   %OK              %removing trailing zeros is probably a good idea
% Kids = 1:N;                               %better/nicer    %resolution of Kids will substantially effect the Linkage Loop
DImax = max(V)-min(V);
w1 = 1;
nTZ = min(find(V~=0))-1;                %Number of trailing zeros
nEZ = min(find(flipud(V)~=0))-1;        %Number of ending zeros

%% matrices
M = zeros(N,2);
M(:,2) = V;
sigma = ones(1,L);
sigma(1) = 1; %exp(0*(deltatau))
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
        %         display('First Pass...');
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
        %         cat(1,SegmentEnd{n-1},UK,P)
        %         cat(1,UP,Gp(UP)')
        %         display('Second Pass...');
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
        %         display('First Pass...');
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
        %         display('Second Pass...');
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

%% Compute enrichment (Z-score)
display('Enrichment computation');
CV = cumsum(V);
Enrichment = cell(L,1);
mu = mean(V);
sigma2 = var(V);
minLsize=10;

for n = 1:L;
    fprintf(1,'.');
    SL = length(SegmentEnd{n});
    Segments = cat(1,[1 SegmentEnd{n}(1:SL-1)+1],SegmentEnd{n});
    Le = Segments(2,:)-Segments(1,:)+1;
    S = CV(Segments(2,:))-CV(Segments(1,:))+V(Segments(1,:));
    Enrichment{n} = ((S'-Le.*mu)./sqrt(max(Le,minLsize).*sigma2));
end
display(' ');
display(['Total Time : '  num2str(round(etime(clock,time1))) ' seconds']);





