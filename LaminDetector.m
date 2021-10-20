function data=LaminDetector(liffile,idx)

% This function segments the lamin network of micronuclei to 
% identify holes in the mesh from 3D image datasets acquired with
% super res STED technique using the AROS algorithm developped by Mark 
% Kittisopikul and colleagues.
%
% INPUT:        -liffile: string name of the image file (.lif or .tif, 8bits)
%               -idx: image index or series in multiseries datasets
%
%OUTPUT:        -data: structure containing global morphometrics of the
%whole micronucleus and morphometrics of single identified holes.
%
%               data.name: name of the file
%               data.nms:m-by-n-by-z array of the weigthed segmented lamin
%               mesh from SAROS detector
%               data.nmsMedianIntMaxInt: 1-by-2 vector containing median
%               and max intensities of nms
%               data.InnerRadii:1-by3 vector containing radii of the fitted
%               ellipsoid
%               data.InnerArea: Area of the fitted ellipsoid
%               data.InnerVolume: Volume of the fitted ellipsoid
%               data.OuterPrincipalAxisLength: length of the principal axis
%               of the convex hull containing the MN
%               data.OuterArea: Area of the convex hull
%               data.OuterVolume: Volume of the convex hull
%               data.StatsTop (data.StatsBot): table containing the
%               properties of identified holes in the top and bottom
%               hemispheres of the MN. One hole per row, with 7 columns:
%                       -Top (Botttom): unique hole ID
%                       -Area: Hole area
%                       -Eccentricity
%                       -Solidity
%                       -Perimeter
%                       -MeanIntensity
%                       -IntRatio: Intensity ratio between hole and non
%                       hole region
%               data.scale:1-by-2 vector of physical voxel size in xy and z
%               axis
%               data.TotalzPlanes: number of z planes in the image stack
%               data.BottomzPlanes:number of z planes of the bottom
%               hemisphere

% Examples:
%  myData = LaminDetector('mySTEDImage.lif',3);
%  myData = LaminDetector('STED_Example.tif',1);

% Julien Dubrulle, Ph.D.
% Fred Hutchinson Cancer Research Center
% 2021

%% Load imagedata and initialize
fname=[liffile(1:end-4) '_' int2str(idx) 'f.tif']; %filename of processed image for export

[IS, omeMeta]=imreadLif2(liffile,idx); %reads image file using bioformat converters

% the imageset typically contains 5 channels:

[sx, sy, sz, sc]=size(IS);
CENTRO=IS(:,:,:,1); %channel1, Centromere marker
CENTROP=max(CENTRO,[],3);
CHROM=IS(:,:,:,2); %channel2, chromosome paint marker
CHROMP=max(CHROM,[],3);
HISTONE=IS(:,:,:,3); %channel3, histone marker
HISTONEP=max(HISTONE,[],3);
DAPI=IS(:,:,:,4); %channel4, DAPI
DP=max(DAPI,[],3);
LAMIN=IS(:,:,:,5); %channel5, Lamin antibody labelling

LAMIN=medfilt3(LAMIN);
Lp=max(LAMIN,[],3);

data=struct;

Xvoxel=omeMeta.getPixelsPhysicalSizeX(0).value();
VoxX=Xvoxel.doubleValue();
Yvoxel=omeMeta.getPixelsPhysicalSizeY(0).value();
VoxY=Yvoxel.doubleValue();
Zvoxel=omeMeta.getPixelsPhysicalSizeZ(0).value();
VoxZ=Zvoxel.doubleValue();
scale=[VoxX VoxZ];

%% Select micronucleus of interest

% Display montage of markers to identify MN of interest
montage(cat(3,CENTROP,CHROMP,HISTONEP,Lp),'DisplayRange',[0 100]);pause;
close;

% Select ROI containing MN of interest
figure('Name','Select ROI'), imshow(Lp,[]);
roi=drawrectangle;pause;
RW=createMask(roi);close;

%% Crop, adjust and segment micronucleus and find best z-focus
tic


SRW=regionprops(RW,'BoundingBox','Area');
rect=SRW.BoundingBox;
RWdC=imcrop(RW,rect);
RWdCcat=repmat(RWdC,[1 1 sz]);
[ca, cb, cc]=size(RWdCcat);

cLp=imcrop(Lp,rect);
lowhigh=stretchlim(cLp);
T=minminT(imadjust(Lp,lowhigh));
LAMINC=zeros(ca, cb, cc,'Like',IS);
ms=zeros(sz,1);
for i=1:sz
    LAMINC(:,:,i)=imadjust(imcrop(LAMIN(:,:,i),rect),lowhigh,[]);
    ms(i)=fmeasure(LAMINC(:,:,i),'BREN'); %measure focus
end

[~, ims]=max(ms);
disp(ims)
if ims==1
    data=[];
    return;
end

LAMINTC=LAMINC-1*T;
CENTROPC=imcrop(CENTROP,rect);
CHROMPC=imcrop(CHROMP,rect);
HISTONEPC=imcrop(HISTONEP,rect);
DAPIC=imcrop(DP,rect);

% Save cropped channels for QC
figure, imshow(imadjust(medfilt2(CENTROPC)),'InitialMagnification',150);CEdataI=print('-RGBImage','-r0');close;
figure, imshow(imadjust(medfilt2(CHROMPC)),'InitialMagnification',150);CHdataI=print('-RGBImage','-r0');close;
figure, imshow(imadjust(medfilt2(HISTONEPC)),'InitialMagnification',150);HIdataI=print('-RGBImage','-r0');close;
figure, imshow(imadjust(medfilt2(cLp)),'InitialMagnification',150);LAdataI=print('-RGBImage','-r0');close;

% Segment MN boundaries based on DAPI signal
TD=minminT(DAPIC);
Db=DAPIC-TD;
DIM=medfilt2(imdiffusefilt(Db),[5 5]);

DW=DIM>0;
h=fspecial('average',10);
DWf=imfilter(DW,h);
DWf=imfill(DWf,'holes');
DWfss=imclearborder(DWf); %if MN touches main nucleus, expect imclearborder to clear whole image

hh=fspecial('average',20);
DWfss=imfilter(DWfss,hh);

% if micronucleus touching main, further segment by hand
if ~any(DWfss(:))
    figure('Name','Circle OOI'), imshow(cLp,[]);
    roi=drawfreehand(gca,'Closed',1,'Multiclick',1);pause; %contour the MN
    RWW=createMask(roi);close;
    DWfss=DWf & RWW;
end

%% SAROS algo to segment lamin network

NMS=zeros(size(LAMINTC));
for i=1:sz
    I=LAMINTC(:,:,i);
    if any(I(:))
        [~,~,nms]= steerableAdaptiveResolutionOrientationSpaceDetector...
            ( double(I));
        NMS(:,:,i)=nms;
    end
end


NMS(repmat(~DWfss,1,1,sz))=0; 

%save the segmentation results for QC
figure, imshow(max(NMS,[],3),[0 50],'InitialMagnification',200);
NMSDataI=print('-RGBImage','-r0');close


%% Binarize Lamin network
NMS(isnan(NMS))=0;
NMSt=imtophat(NMS,strel('disk',2));
NMStg=imgaussfilt3(NMSt,2);
NMStgval=NMStg(NMStg>0);
mN=max(NMStgval);
aN=mean(NMStgval);

%BBW=NMStg>mN/10; %shallow segmentation
%BBW=NMStg>10; %very shallow segmentation
BBW=NMStg>aN; %stringent segmentation

Bout=max(NMStg,[],3);
BoutW=Bout>mN/5; %to detect perimeter --Stringent--
%BoutW=Bout>mN/8; %to detect perimeter --shallow--

BLayer=max(max(BBW,[],2),[],1);
%firstPlane=find(BLayer,1,'First');
lastPlane=find(BLayer,1,'Last');
nPlanes=numel(find(BLayer));
BBWd=false(size(BBW));
for i=1:sz
    BBWd(:,:,i)=imdilate(BBW(:,:,i),strel('disk',2));
end

%% Fit ellipsoid to segmented lamin mesh

S=regionprops3(BBW,'VoxelList');
SF=regionprops3(single(BBWd),'ConvexVolume','ConvexImage');
ICV=SF.ConvexImage;
SSFF=regionprops3(ICV{1},'PrincipalAxisLength');
VL=vertcat(S.VoxelList{:});
VL(:,[1 2])=VL(:,[1 2])*scale(1);
VL(:,3)=VL(:,3)*scale(2);
[~,Radius_LSE,~,v] = ellipsoid_fit(VL);



 


%% Compute properties for top hemisphere

TopW=BBW(:,:,1:ims-1);

MTopW=logical(max(TopW,[],3)) | bwperim(bwconvhull(BoutW),8);

TUCK=imclearborder(~MTopW);
TUCKa=bwareaopen(TUCK,288); %0.12um^2 (minimum hole area)

LTCT=max(LAMINTC(:,:,1:ims-1),[],3);

StatsTop=regionprops('table',TUCKa,LTCT,'Area','Perimeter','Eccentricity','Solidity','MeanIntensity'); 
TSK=mean(LTCT(MTopW));
StatsTop.IntRatio=StatsTop.MeanIntensity/TSK;
StatsTop.Area=StatsTop.Area*scale(1)^2; %convert to um^2
StatsTop.Perimeter=StatsTop.Perimeter*scale(1); %convert to um

% display segmented holes on top of lamin signal
figure, imshow(imadjust(LTCT),[],'InitialMagnification',200);hold on;
B=bwboundaries(TUCKa);
for i=1:length(B)
    bb=B{i};
    plot(gca,bb(:,2),bb(:,1),'-r','LineWidth',1);
end

drawnow;
LTDataI=print('-RGBImage','-r0');close;

%% Compute properties for bottom hemisphere


if (lastPlane-ims)>3 %need at least 3 bottom zstacks


BotW=BBW(:,:,ims+1:end);

MBotW=logical(max(BotW,[],3)) | bwperim(bwconvhull(BoutW),8);

BUCK=imclearborder(~MBotW);
BUCKa=bwareaopen(BUCK,288);

LTCB=max(LAMINTC(:,:,ims+1:end),[],3);

StatsBot=regionprops('table',BUCKa,LTCB,'Area','Perimeter','Eccentricity','Solidity','MeanIntensity'); 
BSK=mean(LTCB(MBotW));
StatsBot.IntRatio=StatsBot.MeanIntensity/BSK;
StatsBot.Area=StatsBot.Area*scale(1)^2; %convert to um^2
StatsBot.Perimeter=StatsBot.Perimeter*scale(1); %convert to um


figure, imshow(imadjust(LTCB),[],'InitialMagnification',200);hold on;
B=bwboundaries(BUCKa);
for i=1:length(B)
    bb=B{i};
    plot(gca,bb(:,2),bb(:,1),'-r','LineWidth',1);
end

drawnow;
LTDataBI=print('-RGBImage','-r0');close;
else
    StatsBot=[];
    LTDataBI=uint8(zeros(size(cLp)));
    figure, imshow(LTDataBI,[0 255],'InitialMagnification',200);
    LTDataBI=print('-RGBImage','-r0');close;
end

%% Gather images and metrics


TopPan=cat(2,CEdataI,CHdataI,HIdataI,LAdataI);
[sa, sb, sc]=size(TopPan);
BotPan=cat(2,NMSDataI,LTDataI,LTDataBI);
[ba, bb, bc]=size(BotPan);

if sb~=bb
    BotPan=imresize3(BotPan,[ba,sb,sc]);
end
    
    
    FinalFigure=uint8(zeros(sa+ba,sb,3));
    FinalFigure(1:sa,:,:)=TopPan;
    FinalFigure(sa+1:end,:,:)=BotPan;
    
    figure, imshow(FinalFigure);
    


export_fig(fname,'-tif','-r150')
nelTop=size(StatsTop,1);
nelBot=size(StatsBot,1);
rname=[fname(1:end-4) '_'];
vtop=1:nelTop;vtop=vtop';
toplab=cellstr([repmat(rname,nelTop,1) int2str(vtop)]);


vbot=1:nelBot;vbot=vbot';
botlab=cellstr([repmat(rname,nelBot,1) int2str(vbot)]);

TTop=table(toplab,'VariableNames',{'Top'});
TBot=table(botlab,'VariableNames',{'Bottom'});


PAL=SSFF.PrincipalAxisLength/2;
PAL(1:2)=PAL(1:2)*scale(1);
PAL(3)=PAL(3)*scale(2);


data.name=fname;
data.nms=NMS;
data.nmsMedianIntMaxInt=[median(NMStgval) max(NMStgval)];
data.InnerRadii=Radius_LSE; %already in microns
data.InnerArea= 4*pi*( ((Radius_LSE(1)*Radius_LSE(2))^1.6 + (Radius_LSE(1)*Radius_LSE(3))^1.6 + (Radius_LSE(2)*Radius_LSE(3))^1.6)/3)^(1/1.6);
data.InnerVolume=4/3*pi*Radius_LSE(1)*Radius_LSE(2)*Radius_LSE(3); %already in microns^3
data.OuterPrincipalAxisLength=PAL;
data.OuterArea= 4*pi*( ((PAL(:,1).*PAL(:,2)).^1.6 + (PAL(:,1).*PAL(:,3)).^1.6 + (PAL(:,2).*PAL(:,3)).^1.6)/3).^(1/1.6);
data.OuterVolume=SF.ConvexVolume*VoxX^2*VoxZ;
if ~isempty(StatsTop)
data.StatsTop=[TTop StatsTop];
else
    data.StatsTop=StatsTop;
end
if ~isempty(StatsBot)
data.StatsBot=[TBot StatsBot];
else
    data.StatsBot=StatsBot;
end
data.scale=scale;
data.TotalzPlanes=nPlanes;
data.BottomzPlanes=(lastPlane-ims);
formatspec='The number of slices is %d and the number of bottom slices is %d';
sprintf(formatspec,nPlanes,(lastPlane-ims))

toc

end

function thresvalue=minminT(I)

thresvalue = max([min(max(I,[],1)) min(max(I,[],2))])	 ;

end

