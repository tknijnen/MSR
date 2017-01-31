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

function_name = 'msr_runtime_SIGNAL.m';
fid = fopen([pwd '/make_mcc.sh'],'w+');
fprintf(fid,['/tools/matlab_R2013a/bin/mcc -R ''-logfile,log_msr_runtime_SIGNAL.txt'' -m -I ' pwd '/ -d ' pwd ' ' function_name ' -a ./mfiles \n']);
fclose(fid);
[status,result] = unix(['chmod 777 -f ' pwd '/make_mcc.sh']);
[status,result] = unix([pwd '/make_mcc.sh']);

% function_name = 'msr_runtime_BED.m';
% fid = fopen([pwd '/make_mcc.sh'],'w+');
% fprintf(fid,['/tools/matlab_R2013a/bin/mcc -R ''-logfile,log_msr_runtime_BED.txt'' -m -I ' pwd '/ -d ' pwd ' ' function_name ' -a ./mfiles \n']);
% fclose(fid);
% [status,result] = unix(['chmod 777 -f ' pwd '/make_mcc.sh'])
% [status,result] = unix([pwd '/make_mcc.sh'])