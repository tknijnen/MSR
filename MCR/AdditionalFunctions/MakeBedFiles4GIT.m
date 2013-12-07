%% File definitions
% sfile = '/titan/cancerregulome9/workspaces/users/tknijnen/CPD/ChipSeqData/20090611_1919_A_BMM_NoStim_0000_PolII.bed';
% vfile = '/titan/cancerregulome9/workspaces/users/tknijnen/CPD/ChipSeqData/20090529_1922_A_BMM_LPS_0240_PolII.bed';
bdir =  '/titan/cancerregulome9/workspaces/users/tknijnen/MSS/MSS_bin/UnMap/';
sdir =  '/titan/cancerregulome9/workspaces/users/tknijnen/MSS/MSS_bin/20090611_1919_A_BMM_NoStim_0000_PolII/';
vdir =  '/titan/cancerregulome9/workspaces/users/tknijnen/MSS/MSS_bin/20090529_1922_A_BMM_LPS_0240_PolII/';
chrs = [9 16 19];

sbedg = '/users/tknijnen/matlab/MSRGIT/MCR/Unix64/data/in/PolII_chr9-16-19.bedgraph';
vbedg = '/users/tknijnen/matlab/MSRGIT/MCR/Unix64/data/in/PolII_after2hLPS_chr9-16-19.bedgraph';
umbedg = '/users/tknijnen/matlab/MSRGIT/MCR/Unix64/data/in/UM_chr9-16-19.bedgraph';

if exist(sbedg,'file');[status,result] = unix(['rm ' sbedg]);end
if exist(vbedg,'file');[status,result] = unix(['rm ' vbedg]);end
if exist(umbedg,'file');[status,result] = unix(['rm ' umbedg]);end

%% SBED
% [status,result] = unix(['head -1 ' sfile ' > ' sbed]);
% for n = 1:length(chrs);
%     n
%     [status,result] = unix(['tail -n +1 ' sfile ' | awk ''$1=="' 'chr' num2str(chrs(n)) '"'' >> ' sbed]);
% end
% n+1
% [status,result] = unix(['cut -f ' '1,2,3,5' ' ' sbed ' > ' sbedg]);
% 
% %% VBED
% [status,result] = unix(['head -1 ' vfile ' > ' vbed]);
% for n = 1:length(chrs);
%     n
%     [status,result] = unix(['tail -n +1 ' vfile ' | awk ''$1=="' 'chr' num2str(chrs(n)) '"'' >> ' vbed]);
% end
% n+1
% [status,result] = unix(['cut -f ' '1,2,3,5' ' ' vbed ' > ' vbedg]);

%% UMBED
IS = 10;
fid = fopen(umbedg,'w+');
fprintf(fid,'track type=bedGraph name="Unique mappability at resolution 100bp"\n');
for n = 1:length(chrs);
    n
    load([bdir '158_chr' num2str(chrs(n)) '.mat']);
    LV = filter(ones(1,IS)/IS,1,LV);
    R = IS:IS:length(LV);
    if R(end)~=length(LV);
        R(end+1) = length(LV);
    end
    LV = LV(R);
    R = R*10; %stored as resolution 10 bp
    SL = length(R);
    Segments = cat(1,[1 R(1:SL-1)+1],R);
    Segments(1,:) = Segments(1,:)-1;
    Segments(:,LV==0)=[];
    LV(LV==0)=[];
    str = ['chr' num2str(chrs(n))];
    fs = [];
    for s = 1:length(str);
        fs = [fs '%c'];
    end
    fs = [fs '\t%d\t%d\t%6.4f\n'];
    fprintf(fid,fs,[ones(length(LV),1)*str Segments' LV]');
end
fclose(fid);
n+1;

%% SBED
IS = 10;
fid = fopen(sbedg,'w+');
fprintf(fid,'track type=bedGraph name="Pol II signal at resolution 100bp"\n');
for n = 1:length(chrs);
    n
    load([sdir 'chr' num2str(chrs(n)) '.mat'],'V');
    LV = V; clear V;
    LV = filter(ones(1,IS)/IS,1,LV);
    R = IS:IS:length(LV);
    if R(end)~=length(LV);
        R(end+1) = length(LV);
    end
    LV = LV(R);
    R = R*10; %stored as resolution 10 bp
    SL = length(R);
    Segments = cat(1,[1 R(1:SL-1)+1],R);
    Segments(1,:) = Segments(1,:)-1;
    Segments(:,LV==0)=[];
    LV(LV==0)=[];
    str = ['chr' num2str(chrs(n))];
    fs = [];
    for s = 1:length(str);
        fs = [fs '%c'];
    end
    fs = [fs '\t%d\t%d\t%6.4f\n'];
    fprintf(fid,fs,[ones(length(LV),1)*str Segments' LV]');
end
fclose(fid);
n+1;

%% VBED
IS = 10;
fid = fopen(vbedg,'w+');
fprintf(fid,'track type=bedGraph name="Pol II signal after 2hr LPS at resolution 100bp"\n');
for n = 1:length(chrs);
    n
    load([vdir 'chr' num2str(chrs(n)) '.mat'],'V');
    LV = V; clear V;
    LV = filter(ones(1,IS)/IS,1,LV);
    R = IS:IS:length(LV);
    if R(end)~=length(LV);
        R(end+1) = length(LV);
    end
    LV = LV(R);
    R = R*10; %stored as resolution 10 bp
    SL = length(R);
    Segments = cat(1,[1 R(1:SL-1)+1],R);
    Segments(1,:) = Segments(1,:)-1;
    Segments(:,LV==0)=[];
    LV(LV==0)=[];
    str = ['chr' num2str(chrs(n))];
    fs = [];
    for s = 1:length(str);
        fs = [fs '%c'];
    end
    fs = [fs '\t%d\t%d\t%6.4f\n'];
    fprintf(fid,fs,[ones(length(LV),1)*str Segments' LV]');
end
fclose(fid);
n+1;

