function plotMSS(x,Lmin,Lmax,V,Enrichment,KM,SegmentEnd,K);

%% parameters
lw = 1;
mlw = 1;
fac2 = 1;       %Linewidths for segmentation tree
fac3 = 1.5;     %Linewidths for segmentation tree
XTG = 11;
facs = 0.0008;  %to plot very small (smaller than pixelwidth) peaks
facs = 0;
wsl = round(min(2e3,diff(x)/25)); %smoothing window of original signal for plotting
plotflag=0; %if not 1, skip the Gaussian scale and the tree view

%% Plotting original signal
N = length(V);
L = length(KM);
R = 1:N;
x1 = x(1);
x2 = x(2);

close all;
h = figure; set(h,'Color','w','Position',[1 29 1280 929]);
h1 = subplot('Position',[.08 .6 .91 .35]);hold on

if wsl>1;
    window = ones(1,wsl);
    window = window/sum(window);
    if wsl<1250;
        V = conv2(V,window','same');
    else
        V = convnfft(V',window,'same');
    end
end

plot(R(x1:x2),V(x1:x2),'k','LineWidth',lw);set(gca,'XLim',[R(x1) R(x2)]);
title([num2str(x1) ' - ' num2str(x2)]);set(gca,'XTick',[]);
ylabel('Signal');
h2 = subplot('Position',[.08 .05 .91 .55]);

%%
if plotflag==1;
    %Gaussian scale space
    h2 = subplot('Position',[.08 .2 .91 .40]);
    
    deltatau = .5*log(2);
    Pmin = 1e-3;

    M = zeros(x2-x1+1,Lmax);
    M(:,1) = V(x1:x2);

    sigma = ones(1,L);
    sigma(1) = 1; %exp(0*(deltatau))
    
    for n = 2:Lmax
        
        %% Gaussian scale space
        sigma(n) = exp((n-1)*deltatau);
        X = norminv(Pmin,0,sigma(n));
        window = normpdf(round(X):-round(X),0,sigma(n));
        windowSize = length(window);
        window = window/sum(window);
        if windowSize<1250;
            M(:,n) = conv2(M(:,1),window','same');
        else
            M(:,n) = convnfft(M(:,1)',window,'same');
        end
        
    end
    
    M = M(:,Lmin:Lmax);
    axes(h2);
    imagesc(log(M'+1));set(gca,'XTick',[]);XL = get(gca,'Xlim');YL = get(gca,'Ylim');
    ylabel('Scale');
    YT = get(h2,'YTick');
    ytl = Lmin:Lmax;
    YTL = round(ytl(YT));
    set(h2,'YTickLabel',YTL);
    hold on
    
    drawnow;pause;
    
    %Segmentation tree
    
    ED = zeros(1,0);
    pos = cell(Lmax,1);
    for mm = Lmin:Lmax;
        pos{mm} = find(SegmentEnd{mm}>x1&SegmentEnd{mm}<x2);
        if isempty(pos{mm});
            pos{mm} = length(SegmentEnd{mm});
        elseif length(pos{mm})<length(SegmentEnd{mm})
            pos{mm} = [pos{mm} pos{mm}(end)+1];
        end
        if any(pos{mm})&mm~=Lmax
            if mm==Lmin;
                plot(KM{mm+1}(:,pos{mm})-x1+1,[mm-1.5 mm]+2-Lmin,'w','LineWidth',fac3*(mlw*mm/Lmax).^fac2);
            else
                plot(KM{mm+1}(:,pos{mm})-x1+1,[mm-1 mm]+2-Lmin,'w','LineWidth',fac3*(mlw*mm/Lmax).^fac2);
            end
        end
        ED = cat(2,ED,Enrichment{mm}(:)');
        ED(ED==0) = [];
    end
    
    drawnow;pause;
    
else
   
    pos = cell(Lmax,1);
    ED = zeros(1,0);
    for mm = Lmin:Lmax;
        pos{mm} = find(SegmentEnd{mm}>x1&SegmentEnd{mm}<x2);
        if isempty(pos{mm});
            pos{mm} = min(find(SegmentEnd{mm}>x2));
            if isempty(pos{mm});
                pos{mm} = length(SegmentEnd{mm});
            end
        elseif length(pos{mm})<length(SegmentEnd{mm})
            pos{mm} = [pos{mm} pos{mm}(end)+1];
        end
        ED = cat(2,ED,Enrichment{mm}(:)');
        ED(ED==0) = [];
    end
    
end
    
%% Plot enrichment
YL = [.5 Lmax-Lmin+1.5];
XL = [.5 diff(x)+1.5];

N = 30;
flag = 0;
while flag==0;
    x21 = prctile(ED(ED<0),linspace(0,100,N));
    x22 = prctile(ED(ED>0),linspace(0,100,N));
    if any(isnan(x21)); x21 = linspace(-1,-0.1,N); end
    if any(isnan(x22)); x22 = linspace(.1,1,N); end
    X2 = [x21 0 x22];
    Y = cat(1,interp1(linspace(0,1,3),cat(1,[0 0 .5],[0 0 1],[.8 .8 1]),linspace(0,1,N)),...
        [1 1 1],...
        interp1(linspace(0,1,3),cat(1,[1 .8 .8],[1 0 0],[.5 0 0]),linspace(0,1,N)));
    if length(X2)==length(unique(X2));
        flag = 1;
    end
    N = N - 1;
end
X2
Y
N

fac = 1;

axes(h2);
cla;

if ~exist('Eth','var');
    Eth = 0;
end

if ~isempty(K);
    M2 = ones(Lmax-Lmin+1,x2-x1+1,3)*.7;
    M3 = cell(0,4);q = 0;
    for mm = Lmin:Lmax;
        if any(pos{mm})
            SE = SegmentEnd{mm}(pos{mm});
            E = Enrichment{mm}(pos{mm});
            E = E(:)';
            SL = length(SE);
            Segments = cat(1,[1 SE(1:SL-1)-1],SE);
            KP = K{mm}(pos{mm});
            Segments(1,1) = x1;
            Segments(2,end) = x2;
            for s = find(abs(E)>Eth);
                rgb = interp1(X2,Y,E(s)).^fac;
                if KP(s)==0;
                     rgb = .7*[1 1 1];
                end
                vecpos = (Segments(1,s):Segments(2,s))-x1+1;
                if (length(vecpos)/(x2-x1))<facs;
                    vecpos = [round(mean(vecpos)-(facs*(x2-x1))) round(mean(vecpos)+(facs*(x2-x1)))];
                    vecpos = min(max(1,vecpos),x2-x1+1);
                    vecpos = vecpos(1):vecpos(2);
                end
                M2(mm-Lmin+1,vecpos,1) = rgb(1);
                M2(mm-Lmin+1,vecpos,2) = rgb(2);
                M2(mm-Lmin+1,vecpos,3) = rgb(3);
                if length(vecpos)>1e4;
                    q = q + 1;
                    M3{q,1} = mm-Lmin+1;
                    M3{q,2} = mean(vecpos);
                    M3{q,3} = E(s);
                    M3{q,4} = mean(rgb);
                end
            end
        end
    end
    image(M2);
else
    M2 = ones(Lmax-Lmin+1,x2-x1+1,3);
    M3 = cell(0,4);q = 0;
    for mm = Lmin:Lmax;
        if any(pos{mm})
            SE = SegmentEnd{mm}(pos{mm});
            E = Enrichment{mm}(pos{mm});
            E = E(:)';
            SL = length(SE);
            Segments = cat(1,[1 SE(1:SL-1)-1],SE);
            Segments(1,1) = x1;
            Segments(2,end) = x2;
            for s = find(abs(E)>Eth);
                rgb = interp1(X2,Y,E(s)).^fac;
                vecpos = (Segments(1,s):Segments(2,s))-x1+1;
                if (length(vecpos)/(x2-x1))<facs;
                    vecpos = [round(mean(vecpos)-(facs*(x2-x1))) round(mean(vecpos)+(facs*(x2-x1)))];
                    vecpos = min(max(1,vecpos),x2-x1+1);
                    vecpos = vecpos(1):vecpos(2);
                end
                M2(mm-Lmin+1,vecpos,1) = rgb(1);
                M2(mm-Lmin+1,vecpos,2) = rgb(2);
                M2(mm-Lmin+1,vecpos,3) = rgb(3);
                if length(vecpos)>1e4;
                    q = q + 1;
                    M3{q,1} = mm-Lmin+1;
                    M3{q,2} = mean(vecpos);
                    M3{q,3} = E(s);
                    M3{q,4} = mean(rgb);
                end
            end
            
        end
    end
    image(M2);
end
set(gca,'Xlim',XL,'Ylim',YL,'XTick',linspace(XL(1),XL(2),XTG),'XTickLabel',[]);
ylabel('Scale');grid;
YT = get(gca,'YTick');
set(gca,'YTickLabel',Lmin+YT-1);

% for n = 1:q;
%     ht = text(M3{n,2},M3{n,1},num2str(M3{n,3},3));
%     set(ht,'FontSize',7,'HorizontalAlignment','center','Color',ones(3,1)*double(M3{n,4}<.3));
% end
% 
% 



