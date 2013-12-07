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

function readBED(pp,bedfile,folder,CN,CS,IS);

flag=0;
for c = 1:length(CN);
    if ~exist([folder '/' pp '_' CN{c} '.mat'],'file');
        flag=1;
    end
end

if flag==1;
    
    display(['Reading ' bedfile]);
    fid = fopen(bedfile);
    E = textscan(fid,'%s%n%n%n','Delimiter','\t','HeaderLines',1); %assuming one header line
    fclose(fid);
    
    for c = 1:length(CN);
        display(['Running ' CN{c} '...']);
        pos = find(strcmp(E{1},CN{c}));
        if length(pos)<10;
            display('Less than 10 lines in BED file for this chromosome, skipping');
            S = [];
        else
            S = zeros(CS(c),1);
            for p = 1:length(pos);
                S((E{2}(pos(p))+1):E{3}(pos(p))) = E{4}(pos(p));
            end
            S = filter(ones(1,IS)/IS,1,S);
            S = S(IS:IS:CS(c));
        end
        save([folder '/' pp '_' CN{c} '.mat'],'S');
    end
else
    display('All chromosome files already exist, skipping');
end



