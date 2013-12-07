function msr_runtime_SIGNAL(parameterfile)

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

%% Read in parameter file
L = 25;         %number of scales
Pth = 1e-6;     %P-value threshold for SFC computation
PruneFlag=0;    %Output Pruned MSR (1), complete MSR (-1) or complete MSR exluding segments with SFC=0 (optional, default 0)
DepletionFlag=0;%When pruning, also output depleted segments (1) or not (0)
Eth = 0;        %SFC threshold in pruning
T = 1.05;       %Slack parameter in pruning
R = 0.2;        %Size constraint in pruning
UM = [];        %unique mappability map
S = [];         %signal
VS = [];        %versus signal
OF = [];        %output file

fid=fopen(parameterfile);
if fid==-1;error('input file not found');end

PF = cell(0,2);
p = 0;
rmc = {'''','"',';',' ','\t'};
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    pos = find(tline=='=');
    if ~isempty(pos);
        p=p+1;
        PF{p,1} = tline(1:pos(1)-1);
        tline = tline(pos(1)+1:end);
        pos = find(tline=='#');
        if ~isempty(pos);
            tline(pos(1):end)=[];
        end
        while any(strcmp(tline(1),rmc)); tline(1)=[]; end
        while any(strcmp(tline(end),rmc)); tline(end)=[]; end
        PF{p,2} = tline;
    end
end
fclose(fid);
display(PF);

for n=1:p;
    switch lower(PF{n,1});
        case 'signal'
            S = dlmread(PF{n,2});
        case 'out'
            OF = PF{n,2};
        case 'background'
            try            
                UM = dlmread(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: background not used');
            end
        case 'versus'
            try
                VS = dlmread(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: versus signal not used');
            end
        case 'l'
            try            
                L = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: L is illspecified, using default');
            end
        case 'pth'
            try
                Pth = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: Pth is illspecified, using default');
            end
        case 'pruneflag'
            try
                PruneFlag = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: PruneFlag is illspecified, using default');
            end
        case 'depletionflag'
            try
                DepletionFlag = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: DepletionFlag is illspecified, using default');
            end
        case 'eth'
            try
                Eth = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: Eth is illspecified, using default');
            end
        case 't'
            try
                T = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: T is illspecified, using default');
            end
        case 'r'
            try
                R = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: R is illspecified, using default');
            end
    end
end

%some checks
if isempty(OF);error('output file not specified');end
if min(size(S))~=1;error('signal should be ONE column or vector');end
S = S(:);
if ~isempty(UM);
    if min(size(UM))~=1;error('background signal (unique mappability map) should be ONE column or vector');end
    UM = UM(:);
    if ~all(size(UM)==size(S));
        error('signal and background do not match');
    end
end
if ~isempty(VS);
    if min(size(VS))~=1;error('versus signal should be ONE column or vector');end
    VS = VS(:);
    if ~all(size(VS)==size(S));
        error('signal and versus signal do not match');
    end
end

%% Create MSR
[KM,SegmentEnd] = MSS(S,L);
Enrichment = SignificantFoldChange(SegmentEnd,S,UM,VS,Pth);
if PruneFlag==1;
    Keep = Prune(Enrichment,KM,SegmentEnd,T,R,Eth);
    if DepletionFlag==1;
        E2 = Enrichment;
        for nn = 1:L;
            E2{nn} = -E2{nn};
        end
        KeepD = Prune(E2,KM,SegmentEnd,T,R,Eth);
        for nn = 1:L;
            Keep{nn} = double(Keep{nn}|KeepD{nn});
        end
    end
end
fid = fopen(OF,'w+');
fclose(fid);
for nn = 1:L;
    fprintf(1,'.');
    SL = length(SegmentEnd{nn});
    Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
    if PruneFlag==1;
        pos = find(Keep{nn}==1);
    elseif PruneFlag==-1
        pos = 1:SL;
    else
        pos = find(Enrichment{nn}~=0);
    end
    dlmwrite(OF,[nn*ones(length(pos),1) Segments(:,pos)' Enrichment{nn}(pos)],'delimiter','\t','-append');
end







    

