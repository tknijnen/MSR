function plotMSS(x,Lmin,Lmax,V,Enrichment,KM,SegmentEnd,X2);

x1 = x(1);
x2 = x(2);
deltatau = .5*log(2);
Pmin = 1e-3;

N = length(V);
L = length(KM);
R = 1:N;

M = zeros(x2-x1+1,Lmax);
M(:,1) = V(x1:x2);

sigma = ones(1,L);
sigma(1) = 1; %exp(0*(deltatau))

lw = 0.5;
mlw = 1;
fac2 = 1;
fac3 = 1.5;

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

h = figure; set(h,'Color','w','Position',[1 29 1280 929]);
h1 = subplot('Position',[.08 .6 .91 .35]);plot(R(x1:x2),M(:,1),'k','LineWidth',lw);set(gca,'XLim',[R(x1) R(x2)]);
title([num2str(x1) ' - ' num2str(x2)]);set(gca,'XTick',[]);
M = M(:,Lmin:Lmax);
h2 = subplot('Position',[.08 .05 .91 .55]);imagesc(log(M'+1));set(gca,'XTick',[]);XL = get(gca,'Xlim');YL = get(gca,'Ylim');
axes(h2);
ylabel('Scale');
YT = get(h2,'YTick');
sigma = sigma(Lmin:Lmax);
YTL = round(sigma(YT));
set(h2,'YTickLabel',YTL);
axes(h1);
ylabel('Signal');
axes(h2); hold on
clear pos

drawnow;pause;
ED = zeros(Lmax-Lmin+1,2);
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
    ED(mm,1) = max(Enrichment{mm});
    ED(mm,2) = prctile(abs(Enrichment{mm}),85);
end

drawnow;pause;

ED(ED==0)=NaN;

if ~exist('X2','var');
    X2 = [-realmax -prctile(ED(:,1),25) -prctile(ED(:,2),75) 0 prctile(ED(:,2),75) prctile(ED(:,1),25) realmax];
end

Y = [[0 0 .5]' [0 0 .5]' [0 0 1]' [1 1 1]' [1 0 0]' [.5 0 0]' [.5 0 0]']';
fac = 1;

axes(h2); 
cla;
M2 = zeros(Lmax-Lmin+1,x2-x1+1,3);
for mm = Lmin:Lmax;
    if any(pos{mm})
        SE = SegmentEnd{mm}(pos{mm});
        E = Enrichment{mm}(pos{mm});
        SL = length(SE);
        Segments = cat(1,[1 SE(1:SL-1)-1],SE);
        for s = 2:SL
            rgb = interp1(X2,Y,E(s)).^fac;
            M2(mm-Lmin+1,(Segments(1,s):Segments(2,s))-x1+1,1) = rgb(1);
            M2(mm-Lmin+1,(Segments(1,s):Segments(2,s))-x1+1,2) = rgb(2);
            M2(mm-Lmin+1,(Segments(1,s):Segments(2,s))-x1+1,3) = rgb(3);
        end
        
    end
end
image(M2);

set(gca,'Xlim',XL,'Ylim',YL,'XTick',[]);



