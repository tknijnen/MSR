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

function S = EFC(T,G1,G2,O,pth,th);

if nargin<6;
    th = 1e-3;
end
if nargin<5;
    pth = 1e-6;
end

%% Enrichment
Zs = -norminv(pth);
a = (Zs.^2).*G2+(G2.^2);
b = -(2.*O.*G2 + (Zs.^2).*G2);
c = (O.^2);
psp = (-b + sqrt((b.^2)-(4.*a.*c)))./(2.*a);
psn = (-b - sqrt((b.^2)-(4.*a.*c)))./(2.*a);
Pse = NaN*ones(size(psp));
Cmp = ((O-(psp.*G2))./sqrt(G2.*psp.*(1-psp)));
Cmn = ((O-(psn.*G2))./sqrt(G2.*psn.*(1-psn)));
Ip = abs(Cmp-Zs)<th;
In = abs(Cmn-Zs)<th;
if any(any(In&Ip));error('hfs');end
Pse(Ip) = psp(Ip);
Pse(In) = psn(In);
E = (G2.*Pse)./((G2.*G1)./T);
E(E<1|isnan(E))=1;


%% Depletion
Zs = norminv(pth);
a = (Zs.^2).*G2+(G2.^2);
b = -(2.*O.*G2 + (Zs.^2).*G2);
c = (O.^2);
psp = (-b + sqrt((b.^2)-(4.*a.*c)))./(2.*a);
psn = (-b - sqrt((b.^2)-(4.*a.*c)))./(2.*a);
Psd = NaN*ones(size(psp));
Cmp = ((O-(psp.*G2))./sqrt(G2.*psp.*(1-psp)));
Cmn = ((O-(psn.*G2))./sqrt(G2.*psn.*(1-psn)));
Ip = abs(Cmp-Zs)<th;
In = abs(Cmn-Zs)<th;
if any(any(In&Ip));error('hfs');end
Psd(Ip) = psp(Ip);
Psd(In) = psn(In);
D = (G2.*Psd)./((G2.*G1)./T);
D(D>1|isnan(D))=1;

%% Score =
S = ones(size(E));
S(E>1) = E(E>1);
S(D<1) = D(D<1);
S = log2(S);
