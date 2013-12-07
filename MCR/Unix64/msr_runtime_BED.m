function msr_runtime_BED(parameterfile)

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
IS = 10;        %resolution in bps
OF = [];        %output file
UM = [];        %unique mappability map file
S = [];         %signal file
VS = [];        %versus signal file
chrsizes = [];  %file with chromosome sizes
folder = [];    %folder for temporary data (created if not existing)

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
        case 'chrsizes'
            if ~exist(PF{n,2},'file');
                error('chromosome file not found');
            else
                chrsizes = PF{n,2};
            end
        case 'folder'
            if ~exist(PF{n,2},'dir');
                mkdir(PF{n,2});
            end
            folder = PF{n,2};
        case 'signal'
            if ~exist(PF{n,2},'file');
                error('signal file not found');
            else
                S = PF{n,2};
            end
        case 'out'
            OF = PF{n,2};
        case 'background'
            if ~exist(PF{n,2},'file');
                display('WARNING: background signal not found and thus used');
            else
                UM = PF{n,2};
            end
        case 'versus'
            if ~exist(PF{n,2},'file');
                display('WARNING: versus signal not found and thus used');
            else
                VS = PF{n,2};
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
        case 'is'
            try
                IS = str2num(PF{n,2});
            catch
                le =lasterror;display(le.message);
                display('WARNING: IS is illspecified, using default');
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
if isempty(S);error('no signal specified');end
if isempty(chrsizes);error('missing chromsome sizes file');end
if isempty(folder);error('folder for temp data missing');end
fidOF = fopen(OF,'w+');
fprintf(fidOF,['track type=bed name="MSR based on ' parameterfile '"\n']);

%Start parsing files
ChrName = createSIGNALfromBED('S',S,chrsizes,folder,IS);
if ~isempty(UM);
    createSIGNALfromBED('UM',UM,chrsizes,folder,IS);
end
if ~isempty(VS);
    createSIGNALfromBED('VS',VS,chrsizes,folder,IS);
end

for c=1:length(ChrName);
    load([folder '/S_' ChrName{c} '.mat']);
    display(['Running ' ChrName{c} '...']);
    if ~isempty(S);
        
        X=S;
        if ~isempty(UM);
            if exist([folder '/UM_' ChrName{c} '.mat'],'file');
                load([folder '/UM_' ChrName{c} '.mat']);
                UM=S;
            end
        end
        if ~isempty(VS);
            if exist([folder '/VS_' ChrName{c} '.mat'],'file');
                load([folder '/VS_' ChrName{c} '.mat']);
                VS=S;
            end
        end
        S=X;clear X;
        
        %some checks
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
        
        %% Choose writing mode
        str = ChrName{c};
        fs = [];
        for s = 1:length(str);
            fs = [fs '%c'];
        end
        fs = [fs '\t%d\t%d\tscale%d\t%6.4f\n'];
        
        for nn = 1:L;
            display(['Writing output for scale ' num2str(nn)]);
            fprintf(1,'.');
            SL = length(SegmentEnd{nn});
            Segments = cat(1,[1 IS*SegmentEnd{nn}(1:SL-1)+1],IS*SegmentEnd{nn});
            Segments(1,:) = Segments(1,:)-1;
            if PruneFlag==1;
                pos = find(Keep{nn}==1);
            elseif PruneFlag==-1
                pos = 1:SL;
            else
                pos = find(Enrichment{nn}~=0);
            end
            LP = length(pos);
            if LP>0;
                fprintf(fidOF,fs,[ones(LP,1)*str Segments(:,pos)' nn*ones(LP,1) Enrichment{nn}(pos)]');
            end
        end
    end
end
fclose(fidOF)




