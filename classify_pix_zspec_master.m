
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Takes CEST data and performs pixelwise classification into 1 of 3 categories (hydrogel, muscle, noise)
%   Pixelwise classification based on correlation coefficient and NMSE calculations
%   Subsequent Lorentzian fitting and figure generation performed based on classification result
%
%   INPUT(s):
%       None
%   OUTPUT(s):
%       None - automatically analyzes data, generates figures, and saves
%       results
%
%   REQUIRED SCRIPT(s):
%       -zspec_lorentzian_fit_hydrogel.m
%       -zspec_lorentzian_fit_muscle.m
%       -b0_correction.m
%
%   AUTHOR(S):
%       - Bonnie Lam (bonnie.lam@berkeley.edu)
%
%   DATE:
%       - 2023/02/26
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% make save directory
saveDir=[pwd,'\results'];
mkdir(saveDir);

%% image information
frame_size=128;
num_offsets=60;
num_frames=73;

%% load data
% select paths
selPath_full=uigetdir(pwd,'Please select folder with CEST data');
PathName_full=[selPath_full,'\pdata\1\2dseq'];
selPath_rare=uigetdir(pwd,'Please select folder with RARE data');
PathName_rare=[selPath_rare,'\pdata\1\2dseq'];

% read in data
fileID_full=fopen(PathName_full,'r');
fileID_rare=fopen(PathName_rare,'r');
precision='int16';
data_full=fread(fileID_full,Inf,precision);
data_rare=fread(fileID_rare,Inf,precision);
data_full=reshape(data_full,[frame_size frame_size num_frames]);
data_rare=reshape(data_rare,[frame_size frame_size]);

% extract cest offset images and reference images
load('idx_ref.mat');
load('idx_cest.mat');
data_ref=data_full(:,:,idx_ref);
data_cest=data_full(:,:,idx_cest);

%% make/load mask, convert indices
figure; imshow(data_rare,[0 15000]);
mask_leg=roipoly;
indices_leg=find(mask_leg);

save([saveDir,'\mask_leg'],'mask_leg');
save([saveDir,'\indices_leg'],'indices_leg');

% % load indices + masks (if already drawn)
% load([saveDir,'\indices_leg.mat']);
% load([saveDir,'\mask_leg.mat']);

%% calculate + display uncorrected z-spectra
% interpolate references
load('freq_aav.mat');
freq_interp=linspace(-10,10,1000);
freq_aav_ref=[freq_aav(1) freq_aav(5:5:end)];
[x,y,f]=ndgrid(-63.5:63.5,-63.5:63.5,freq_aav_ref);
[xi,yi,fi]=ndgrid(-63.5:63.5,-63.5:63.5,freq_aav);
data_ref_interp=interpn(x,y,f,data_ref,xi,yi,fi,'linear');

% calculate + save pixelwise zspec
cest_spec=data_cest./data_ref_interp;
cest_spec_reshape=reshape(cest_spec,[frame_size*frame_size num_offsets]);
save([saveDir,'\cest_spec.mat'],'cest_spec');
save([saveDir,'\data_rare.mat'],'data_rare');

%% B0 correction on pixelwise zspec
num_categories=2;
num_pix=length(indices_leg);
x_corr_pix=zeros(num_pix,num_offsets);
zspec_pix=cest_spec_reshape(indices_leg,:);

for iter=1:num_pix
    zspec_curr=cest_spec_reshape(indices_leg(iter),:)';
    [x_corr_pix(iter,:),b0_shift_curr]=b0_correction(freq_aav,zspec_curr);
end

save([saveDir,'\x_corr_pix.mat'],'x_corr_pix');
save([saveDir,'\zspec_pix.mat'],'zspec_pix');

%% classify pixels
num_pix=length(indices_leg);
noise_thresh=0.97;
class_results=strings(num_pix,1);    
corrVals=zeros(num_pix,num_categories);   % in order: corr coefficient of hydrogel and muscle

% load model zspec
load('zspec_ref_hydrogel.mat');
load('zspec_ref_muscle.mat');

for iter=1:num_pix
    % interpolate zspec
    zspec_test_curr=zspec_pix(iter,:)';
    x_corr_test_curr=x_corr_pix(iter,:);
    zspec_test_interp=interp1(x_corr_test_curr,zspec_test_curr,freq_interp,'linear','extrap')';
    
    % calculate correlation coefficient
    corrVal_hydrogel=corr2(zspec_test_interp,zspec_high_3uT_hydrogel_interp);
    corrVal_muscle=corr2(zspec_test_interp,zspec_high_3uT_muscle_interp);
    corrVals(iter,1)=corrVal_hydrogel;
    corrVals(iter,2)=corrVal_muscle;
    
    % determine pixel classification
    corrVal_max=max([corrVal_hydrogel corrVal_muscle]);
    if (corrVal_max<noise_thresh)
        class_results(iter)='noise';
    elseif (isequal(corrVal_max,corrVal_hydrogel))
        class_results(iter)='hydrogel';
    else
        class_results(iter)='muscle';
    end

end

corrVals_max=max(corrVals,[],2);

% save results
save([saveDir,'\class_results.mat'],'class_results');
save([saveDir,'\corrVals.mat'],'corrVals');
save([saveDir,'\corrVals_max.mat'],'corrVals_max');

%% perform Lorentzian fitting based on classification category
% initialize variables containing x_corr and zspec data for different categories
num_hydrogel=sum(class_results=='hydrogel');
num_muscle=sum(class_results=='muscle');

% extract corresponding frequency and zspec for different pixel categories
[idx_hydrogel,idy_hydrogel]=find(class_results=='hydrogel');
[idx_muscle,idy_muscle]=find(class_results=='muscle');

x_corr_hydrogel=x_corr_pix(idx_hydrogel,:)';
x_corr_muscle=x_corr_pix(idx_muscle,:)';

zspec_hydrogel=zspec_pix(idx_hydrogel,:)';
zspec_muscle=zspec_pix(idx_muscle,:)';

% perform Lorentzian fitting of z-spectra based on respective classification
F1_hydrogel=struct;
F1_muscle=struct;

if (num_hydrogel~=0)
    F1_hydrogel=zspec_lorentzian_fit_hydrogel(x_corr_hydrogel,zspec_hydrogel,'2stepLD');
end
if (num_muscle~=0)
    F1_muscle=zspec_lorentzian_fit_muscle(x_corr_muscle,zspec_muscle,'2stepLD');
end

%% NMSE exclusion of noisy pixels
NMSE_thresh=0.04;

indices_hydrogel=indices_leg(idx_hydrogel);
indices_muscle=indices_leg(idx_muscle);

[F1_hydrogel_final,indices_hydrogel_final]=NMSE_exclusion(num_hydrogel,indices_hydrogel,F1_hydrogel,NMSE_thresh);
[F1_muscle_final,indices_muscle_final]=NMSE_exclusion(num_muscle,indices_muscle,F1_muscle,NMSE_thresh);
num_hydrogel_final=size(indices_hydrogel_final,1);
num_muscle_final=size(indices_muscle_final,1);

% save results
save([saveDir,'\indices_hydrogel.mat'],'indices_hydrogel');
save([saveDir,'\indices_hydrogel_final.mat'],'indices_hydrogel_final');
save([saveDir,'\indices_muscle.mat'],'indices_muscle');
save([saveDir,'\indices_muscle_final.mat'],'indices_muscle_final');

save([saveDir,'\F1_hydrogel.mat'],'F1_hydrogel');
save([saveDir,'\F1_muscle.mat'],'F1_muscle');
save([saveDir,'\F1_hydrogel_final.mat'],'F1_hydrogel_final');
save([saveDir,'\F1_muscle_final.mat'],'F1_muscle_final');

%% display binary maps of different ROIs
% show original image for comparison
figure; imshow(data_rare,[0 15000]);

% display hydrogel map
mask_hydrogel_final=zeros(frame_size,frame_size);
if (num_hydrogel_final~=0)
    mask_hydrogel_final(indices_hydrogel_final)=1;
    figure; imshow(mask_hydrogel_final)
else
    disp('no pixels in this category')
end

% display muscle map
mask_muscle_final=zeros(frame_size,frame_size);
if (num_muscle_final~=0)
    mask_muscle_final(indices_muscle_final)=1;
    figure; imshow(mask_muscle_final)
else
    disp('no pixels in this category')
end

%% calculate ROI contrast values (averaged from pixelwise Lorentzian fitting)
if (num_hydrogel_final~=0)
    aav_contrast_hydrogel_final=mean(F1_hydrogel_final.contrasts.aav);
else
    aav_contrast_hydrogel_final=NaN;
end

if (num_muscle_final~=0)
    aav_contrast_muscle_final=mean(F1_muscle_final.contrasts.aav);
else
    aav_contrast_muscle_final=NaN;
end

aav_contrast=zeros(1,2);
aav_contrast(1)=aav_contrast_hydrogel_final;
aav_contrast(2)=aav_contrast_muscle_final;
save([saveDir,'\aav_contrast.mat'],'aav_contrast');

%% display histograms for CEST contrast values per category
if (num_hydrogel_final~=0)
    figure; histogram(F1_hydrogel_final.contrasts.aav,'BinWidth',1.0);
    xlabel('AAV CEST contrast (%)');
    ylabel('# of pixels');
    xlim([0 10]);
    set(gca,'fontsize',14);
    title('hydrogel ROI (corr + NMSE filtering)');
    savefig(gcf,[saveDir,'\hist_aav_hydrogel_final.fig']);
end

if (num_muscle_final~=0)    
    figure; histogram(F1_muscle_final.contrasts.aav,'BinWidth',1.0);
    xlabel('AAV CEST contrast (%)');
    ylabel('# of pixels');
    xlim([0 5]);
    set(gca,'fontsize',14);
    title('muscle ROI (corr + NMSE filtering)');
    savefig(gcf,[saveDir,'\hist_aav_muscle_final.fig']);
end

%% generate AAV contrast maps
if (num_hydrogel_final~=0)
    contrast_pix_hydrogel_final_map=zeros(frame_size,frame_size);
    contrast_pix_hydrogel_final_map(indices_hydrogel_final)=F1_hydrogel_final.contrasts.aav;
    figure;imshow(contrast_pix_hydrogel_final_map);
    caxis([0 8]);
    colormap(hot); colorbar;
    set(gca,'fontsize',20);
    savefig(gcf,[saveDir,'\contrast_pix_hydrogel_final_map.fig']);
    save([saveDir,'\contrast_pix_hydrogel_final_map.mat'],'contrast_pix_hydrogel_final_map');
end

if (num_muscle_final~=0)
    contrast_pix_muscle_final_map=zeros(frame_size,frame_size);
    contrast_pix_muscle_final_map(indices_muscle_final)=F1_muscle_final.contrasts.aav;
    figure;imshow(contrast_pix_muscle_final_map);
    caxis([0 8]);
    colormap(hot); colorbar;
    set(gca,'fontsize',20);
    savefig(gcf,[saveDir,'\contrast_pix_muscle_final_map.fig']);
    save([saveDir,'\contrast_pix_muscle_final_map.mat'],'contrast_pix_muscle_final_map');
end

%% median filter contrast maps
if (num_hydrogel_final~=0)
    contrast_pix_hydrogel_final_filt_map=medfilt2(contrast_pix_hydrogel_final_map);
    figure;imshow(contrast_pix_hydrogel_final_filt_map);
    caxis([0 8]);
    colormap(hot); colorbar;
    set(gca,'fontsize',20);
    savefig(gcf,[saveDir,'\contrast_pix_hydrogel_final_filt_map.fig']);
    save([saveDir,'\contrast_pix_hydrogel_final_filt_map.mat'],'contrast_pix_hydrogel_final_filt_map');
end

if (num_muscle_final~=0)
    contrast_pix_muscle_final_filt_map=medfilt2(contrast_pix_muscle_final_map);
    figure;imshow(contrast_pix_muscle_final_filt_map);
    caxis([0 8]);
    colormap(hot); colorbar;
    set(gca,'fontsize',20);
    savefig(gcf,[saveDir,'\contrast_pix_muscle_final_filt_map.fig']);
    save([saveDir,'\contrast_pix_muscle_final_filt_map.mat'],'contrast_pix_muscle_final_filt_map');
end
