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

function ChrName = createSIGNALfromBED(pp,bedfile,chrsizes,folder,IS)

%% Set default of sample size to 10 bp
if exist('IS','var');
    if isempty(IS);
        IS = 10;
    end
else
    IS = 10;
end

%% remove / or \ from folder
if strcmp(folder(end),'/')|strcmp(folder(end),'/');
    folder(end)=[];
end

%% Parse chromosome size file
fid = fopen(chrsizes);
E = textscan(fid,'%s%n','Delimiter','\t','HeaderLines',0,'CommentStyle','#');
fclose(fid);
ChrName = cell(length(E{1}),1);
ChrSize = zeros(length(E{1}),1);
for c = 1:length(ChrName);
    ChrName{c} = E{1}{c};
    ChrSize(c) = E{2}(c);
end

%% Parse BED(graph) file
if isunix
    try %reading bed file per chromosome
        for c = 1:length(ChrName);
            if exist([folder '/' pp '_' ChrName{c} '.mat'],'file');
               display(['File ' folder '/' pp '_' ChrName{c} '.mat already exists, skipping parsing']);
            else
                tmpfile = 'tmp.bed';
                display(['Creating temp file for ' ChrName{c} '...']);
                [status,result] = unix(['tail -n +1 ' bedfile ' | awk ''$1=="' ChrName{c}    '"'' > ' folder '/' tmpfile]);
                readBEDperChr(pp,[folder '/' tmpfile],folder,ChrName{c},ChrSize(c),IS);
                display(['Removing temp file for ' ChrName{c} '...']);
                [status,result] = unix(['rm ' folder '/' tmpfile]);
            end
        end
    catch %read as a whole
        readBED(pp,bedfile,folder,ChrName,ChrSize,IS);
    end
else %read be file as a whole (possible memory errors with large files)
    readBED(pp,bedfile,folder,ChrName,ChrSize,IS);
end

