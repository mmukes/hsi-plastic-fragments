% extract_fragments.m
% Cut each 4x5 grid HSI cube (20 fragments) into one small HSI image per
% fragment and save it in PET, HDPE, PP, PS folders.
% Segmentation is same as FunX.m (PSEP 2023 paper).
% Label of fragment k of scanNNNN = row (NNNN-1)*20+k of 1hotcode.csv
% (fragment numbering is the bwlabel order)
% Mukesh Kumar Singh

clear;close all; clc;

myFolder=fileparts(mfilename('fullpath'));          % folder of this code
startFolder=fullfile(myFolder,'HSI_images');
if ~isfolder(startFolder)
    startFolder=myFolder;
end
inFolder=uigetdir(startFolder,'Select the folder with scan0001.bil ... and 1hotcode.csv');
if isequal(inFolder,0)
    error('No folder selected');
end
outFolder=fullfile(myFolder,'fragments');           % output folder

%% labels (one hot code, columns PET HDPE PP PS)
Y=readmatrix(fullfile(inFolder,'1hotcode.csv'));
if size(Y,2)~=4
    error('1hotcode.csv must have 4 columns');
end
classNames={'PET','HDPE','PP','PS'};
allNames=[classNames,{'Unlabelled','Unmatched'}];
nCount=zeros(1,6);                                  % fragments in each folder
pad=2;                                              % extra pixels around the fragment
if ~isfolder(outFolder)
    mkdir(outFolder);
end

%% list of cubes
theFiles=dir(fullfile(inFolder,'*.bil'));
tF=length(theFiles);
fprintf(1,'%d cubes found in %s\n',tF,inFolder);
idx=cell(0,10);                                     % rows of fragment_index.csv

for f=1:tF
    baseFileName=theFiles(f).name;
    fullFileName=fullfile(theFiles(f).folder,baseFileName);
    [~,base]=fileparts(baseFileName);
    hdrText=fileread(fullfile(theFiles(f).folder,[base '.hdr']));

    %% segmentation (same as FunX)
    [L,nOb]=segmentHSI(fullFileName);

    %% which rows of 1hotcode.csv
    tok=regexp(base,'^scan(\d+)$','tokens','once');
    blk=0;                                          % 0 --> no label
    if isempty(tok)
        warning('%s: name is not scanNNNN, labels not found',baseFileName);
    elseif str2double(tok{1})*20>size(Y,1)
        warning('%s: no rows for this cube in 1hotcode.csv',baseFileName);
    elseif nOb~=20
        warning('%s: %d fragments found instead of 20',baseFileName,nOb);
    else
        blk=str2double(tok{1});
    end
    if blk==0
        fprintf(1,'%-14s %2d fragments  -> Unmatched\n',baseFileName,nOb);
    else
        fprintf(1,'%-14s %2d fragments\n',baseFileName,nOb);
    end

    %% raw cube, all 168 bands, uint16 (for cropping)
    raw=multibandread(fullFileName,[400 320 168],'uint16=>uint16',0,'bil','ieee-le');

    for k=1:nOb
        if blk==0
            cls='Unmatched';
        else
            y=Y((blk-1)*20+k,:);
            if sum(y)==1
                [~,ci]=max(y);
                cls=classNames{ci};
            else
                cls='Unlabelled';                   % not exactly one class
            end
        end
        if ~isfolder(fullfile(outFolder,cls))
            mkdir(fullfile(outFolder,cls));
        end

        % rectangle around fragment k
        % (L is from the cropped cube --> +6 rows, +9 cols in raw cube)
        [r,c]=find(L==k);
        r1=max(min(r)-pad,1); r2=min(max(r)+pad,size(L,1));
        c1=max(min(c)-pad,1); c2=min(max(c)+pad,size(L,2));
        sub=raw(r1+6:r2+6,c1+9:c2+9,:);
        mask=L(r1:r2,c1:c2)==k;                     % white = fragment, black = belt

        % save .bil + .hdr + mask
        name=sprintf('%s_frag%02d',strrep(base,' ','_'),k);
        outBase=fullfile(outFolder,cls,name);
        multibandwrite(sub,[outBase '.bil'],'bil','precision','uint16','machfmt','ieee-le');
        desc=sprintf('%s fragment %d (%s); crop rows %d-%d, cols %d-%d of the original %s', ...
                     base,k,cls,r1+6,r2+6,c1+9,c2+9,baseFileName);
        writeHdr([outBase '.hdr'],hdrText,size(sub),desc);
        imwrite(mask,[outBase '_mask.png']);

        i=find(strcmp(allNames,cls));
        nCount(i)=nCount(i)+1;
        idx(end+1,:)={[cls '/' name '.bil'],cls,baseFileName,k,blk,r1+6,c1+9, ...
                      size(sub,1),size(sub,2),nnz(mask)}; %#ok<AGROW>
    end
end

%% fragment_index.csv
T=cell2table(idx,'VariableNames',{'file','class','source_cube','fragment', ...
             'scan_number','first_row','first_col','rows','cols','fragment_pixels'});
writetable(T,fullfile(outFolder,'fragment_index.csv'));

fprintf(1,'\nFragments written:\n');
for i=1:6
    fprintf(1,'  %-10s %4d\n',allNames{i},nCount(i));
end
fprintf(1,'  %-10s %4d\n','Total',height(T));


%% ----------------------------------------------------------------------
function [L,nOb]=segmentHSI(fullFileName)
% same steps as FunX.m, up to bwlabel
X= multibandread(fullFileName, [400 320 168], 'uint16', 0, 'bil', 'ieee-le');
X=X(7:393,10:313,15:160);
X(X==0)=1;                                  % replacing 0 with 1 (if any 0 is present)

%% probability of each pixel wrt to its spectral bands
P=sum(X,3);                                 % adding the hyper cube in 3rd dimension-->gives summation matrix
p=X./P;                                     % EW division --> probability cube

%% reference pixel window (belt)
a=[5,15];                                   % center pixel location of the window
b=4;                                        % half size of square window
c=a(1)-b:a(1)+b;
d=a(2)-b:a(2)+b;
qw=p(c(1):c(2*b+1),d(1):d(2*b+1),:);        % 9*9 sub cube

%% relative entropy p|q
qmean=mean(qw,[1,2]);                       % mean along each band
Et=sum(p.*(log(p./qmean)),3);               % relative entropy
nEt=(Et-min(Et(:)))/(max(Et(:))-min(Et(:)));% normalized relative entropy

%% power law & thresholding
g1=imadjust(nEt,[],[],0.05);                % power law transformation with gamma=0.05
binim=imbinarize(g1,graythresh(g1));        % Otsu thresholding

%% removing noise
BW1=bwareaopen(binim,150);                  % removes objects <150 pixels
BW3=imopen(BW1,strel('diamond',1));         % opening operation
BW3=bwareaopen(BW3,300);

%% labelling objects
[L,nOb]=bwlabel(BW3);
end

function writeHdr(fname,srcHdr,sz,desc)
% ENVI header, wavelengths and reflectance scale factor copied from the cube
wl=regexp(srcHdr,'wavelength\s*=\s*\{[^}]*\}','match','once');
rsf=regexp(srcHdr,'reflectance scale factor\s*=\s*[^\r\n]*','match','once');
fid=fopen(fname,'w');
fprintf(fid,'ENVI\n');
fprintf(fid,'description = {%s}\n',desc);
fprintf(fid,'samples = %d\nlines = %d\nbands = %d\n',sz(2),sz(1),sz(3));
fprintf(fid,'header offset = 0\nfile type = ENVI Standard\ndata type = 12\n');
fprintf(fid,'interleave = bil\nbyte order = 0\n');
if ~isempty(rsf), fprintf(fid,'%s\n',rsf); end
if ~isempty(wl),  fprintf(fid,'%s\n',wl);  end
fclose(fid);
end
