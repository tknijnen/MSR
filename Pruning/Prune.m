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


function SegmentsFlag = Prune(Enrichment,KM,SegmentEnd,T,R,Eth)

display('Pruning');

L = length(Enrichment);
SegmentsSI = cell(1,L);
SegmentsFlag = cell(1,L);

for n = 1:L;
    fprintf(1,'.');
    SL = length(Enrichment{n});
    SegmentsFlag{n} = double(Enrichment{n}>Eth);
    SegmentsSI{n} = zeros(4,SL);
    SE = SegmentEnd{n};
    Segments = cat(1,[1 SE(1:SL-1)+1],SE);
    SegmentsSI{n}(3,:) = diff(Segments,[],1)+1;
    if n>1;
       [~,b,c] = unique(KM{n}(2,:),'first');
       SegmentsSI{n}(1,:) = b;
       SegmentsSI{n-1}(4,:) = c;
       [~,b] = unique(KM{n}(2,:),'last');
       SegmentsSI{n}(2,:) = b;
    end
end

% Round 1
for n = fliplr(1:L);
    AS = find(SegmentsFlag{n});
    fprintf(1,'.');
    for s = 1:length(AS);
        Z = Enrichment{n}(AS(s));
        SS = SegmentsSI{n}(3,AS(s));
        flag = 0;
        nc = n;
        pos = zeros(n,2);
        pos(n,:) = AS(s);
        while flag==0&nc>1;
            nc = nc - 1;
            pos(nc,1) = SegmentsSI{nc+1}(1,pos(nc+1,1));
            pos(nc,2) = SegmentsSI{nc+1}(2,pos(nc+1,2));
            if pos(nc,1)==pos(nc,2); %same segment
                SegmentsFlag{nc}(pos(nc,1):pos(nc,2)) = 0; %PRUNE
            else
                vec = pos(nc,1):pos(nc,2);
                sp = find(Enrichment{nc}(vec)>=(Z*T));
                if ~isempty(sp); %better segment found
                    if any(SegmentsSI{nc}(3,vec(sp))>(R*SS)) %better segment is same size
                        SegmentsFlag{n}(AS(s)) = 0;
                    end
                    flag=1;
                else %no better segment found
                    sp2 = find(Enrichment{nc}(vec)<=Z);
                    if ~isempty(sp2); %prune worse segments
                        SegmentsFlag{nc}(vec(sp2)) = 0; %PRUNE
                    end
                end
            end
        end
        if flag==0;
            for nc = 1:n-1;
                SegmentsFlag{nc}(pos(nc,1):pos(nc,2)) = 0;
            end
        end
        if SegmentsFlag{n}(AS(s))==1;
            flag = 0;
            nc = n;
            pos = SegmentsSI{nc}(4,AS(s)); 
            while flag==0&nc<L;
               nc = nc + 1;
               if SegmentsSI{nc}(3,pos)<(SS/R); %small enough to be discarded
                   if SegmentsFlag{nc}(pos)==1; %prune larger yet small enough selected segments
                        SegmentsFlag{nc}(pos)=0; %PRUNE
                   end
               else
                   flag=1;
               end
               pos = SegmentsSI{nc}(4,pos);  
            end
        end
    end
end




