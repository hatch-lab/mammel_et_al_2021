function [IS, omeMeta, ImSize]=imreadLif2(liffile,idx)

reader=bfGetReader(liffile); %create a reader using the .lif file name as input
glob=reader.getGlobalMetadata();
ser=reader.getSeriesMetadata();
javaMethod('merge', 'loci.formats.MetadataTools', ...
    glob, ser, 'Global ');


omeMeta=reader.getMetadataStore();
nSeries=reader.getSeriesCount();
x = reader.getBitsPerPixel();
disp(nSeries)

%iSeries=input(['Which series out of ( ' int2str(nSeries) ' )? ']);
iSeries=idx;

% Get dimensions of the dataset
sx=reader.getSizeX();
sy=reader.getSizeY();
sz=reader.getSizeZ();
sC=reader.getSizeC();

ImSize.X=sx;
ImSize.Y=sy;
Imsize.Z=sz;
ImSize.C=sC;

if x==8
IS=zeros(sx,sy,sz,sC,'uint8');
elseif x==16
    IS=zeros(sx,sy,sz,sC,'uint16');
else
    disp('Unknown pixel depth...');
    return
end

reader.setSeries(iSeries-1);
iT=1;

for j=1:sC
    for i=1:sz
        iPlane=reader.getIndex(i-1,j-1,iT-1)+1; %nSeries is the 0 in getIndex
        IS(:,:,i,j)=bfGetPlane(reader,iPlane);
    end
end