% QA evaluation for coils included in the 2020 manuscript:
% "Multi-Site Benchmarking of Clinical 13C RF Coils at 3 T", 
% by S�nchez-Heredia JD*, Olin RB*, McLean MA, Laustsen C, Hansen AE,
% Hanson LG and Ardenkj�r-Larsen JH.
% *shared first authorship

% last edited: December 2019, Rie Beck Olin

clear; close all; clc

QAeval = false;     % set to true for QA evaluation output and figures  
est_T2s = false;    % set to true for estimation of T2* based on FWHM
smap = true;        % set to false for sum-of-squares coil combination

%% coil #1, cylindrical phantom
csi_path = 'data/coil1_cyl.mat'; noise_path = []; noiseBW = 0.793; 

[SNRim_coil1_cyl,~,~,T2s_coil1_cyl] = SNReval_carbon13coils(csi_path,...
    noise_path,noiseBW,QAeval,est_T2s);

%% coil #1, SAM phantom
csi_path = 'data/coil1_sam.mat'; noise_path = []; noiseBW = 0.793;

[SNRim_coil1_sam,~,~,T2s_coil1_sam] = SNReval_carbon13coils(csi_path,...
    noise_path,noiseBW,QAeval,est_T2s);

%% coil #2, SAM phantom
csi_path = 'data/coil2_sam.mat'; noise_path = []; noiseBW = 0.845;

[SNRim_coil2_sam,~,~,T2s_coil2_sam] = SNReval_carbon13coils(csi_path,...
    noise_path,noiseBW,QAeval,est_T2s);

%% coil #3, cylindrical phantom
csi_path = 'data/coil3_cyl.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil3_cyl,~,~,T2s_coil3_cyl] = SNReval_carbon13coils(csi_path,...
    noise_path,noiseBW,QAeval,est_T2s);

%% coil #3, SAM phantom
csi_path = 'data/coil3_sam.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil3_sam,~,~,T2s_coil3_sam] = SNReval_carbon13coils(csi_path,...
    noise_path,noiseBW,QAeval,est_T2s);

%% coil #4, cylindrical phantom
csi_path = 'data/coil4_cyl_8ch.mat'; noise_path = []; noiseBW = 0.845;

[SNRim_coil4_cyl,ncov_coil4_cyl,ncorr_coil4_cyl,T2s_coil4_cyl] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);

%% coil #4, SAM phantom
csi_path = 'data/coil4_sam_8ch.mat'; noise_path = []; noiseBW = 0.845;

[SNRim_coil4_sam,ncov_coil4_sam,ncorr_coil4_sam,T2s_coil4_sam] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);

%% coil #5, cylindrical phantom
csi_path = 'data/coil5_cyl.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil5_cyl,ncov_coil5_cyl,ncorr_coil5_cyl,T2s_coil5_cyl] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);

%% coil #6, cylindrical phantom
csi_path = 'data/coil6_cyl.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil6_cyl,ncov_coil6_cyl,ncorr_coil6_cyl,T2s_coil6_cyl] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);

%% coil #6, SAM phantom
csi_path = 'data/coil6_sam.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil6_sam,ncov_coil6_sam,ncorr_coil6_sam,T2s_coil6_sam] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);

%% coil #7, cylindrical phantom
csi_path = 'data/coil7_cyl.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil7_cyl,~,~,T2s_coil7_cyl] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);

%% coil #7, cylindrical phantom
csi_path = 'data/coil7_sam.mat'; noise_path = []; noiseBW = 0.845; 

[SNRim_coil7_sam,~,~,T2s_coil7_sam] = ...
    SNReval_carbon13coils(csi_path,noise_path,noiseBW,QAeval,est_T2s,smap);


%% display all SNR images
% shift SNRim_coil1_cyl and SNRim_coil6_cyl to align with other cylinder 
% phantom images
[nx,ny] = size(SNRim_coil1_cyl);
SNRim_coil1_cyl_shft = cat(1,SNRim_coil1_cyl(2:end,:),zeros(1,ny));
SNRim_coil1_cyl_shft = cat(2,SNRim_coil1_cyl_shft(:,2:end),zeros(nx,1));

imfontsize = 16;

figure; 
imagesc_row(cat(4,cat(3,SNRim_coil1_cyl_shft,SNRim_coil3_cyl,SNRim_coil4_cyl)...
    ,cat(3,SNRim_coil5_cyl,SNRim_coil6_cyl,SNRim_coil7_cyl)));
h = colorbar; caxis([0 250]); ylabel(h,'SNR')
hold on; text(1+(0:nx:(nx*3-1)),ones(1,3)*(ny-1),...
    {'Coil #1: Birdcage','Coil #3: Helmholtz','Coil #4: 8-CH array'},...
    'fontsize',imfontsize,'color',[1 1 1]);
text(1+(0:nx:(nx*3-1)),ones(1,3)*(2*ny-1),...
    {'Coil #5: 16-CH array','Coil #6: 14-CH array',...
    'Coil #7: Surface coil'},'fontsize',imfontsize,'color',[1 1 1]);
hold off; set(gca,'fontsize',imfontsize)

% shift to align with other SAM phantom images
SNRim_coil2_sam_shft = cat(2,zeros(nx,1),SNRim_coil2_sam(:,1:end-1));
SNRim_coil2_sam_shft = cat(1,zeros(2,ny),SNRim_coil2_sam_shft(1:end-2,:));
SNRim_coil1_sam_shft = cat(1,SNRim_coil1_sam(2:end,:),zeros(1,ny));
SNRim_coil6_sam_shft = cat(1,SNRim_coil6_sam(2:end,:),zeros(1,ny));
SNRim_coil7_sam_shft = cat(1,zeros(1,ny),SNRim_coil7_sam(1:end-1,:));

figure; 
imagesc_row(cat(4,cat(3,SNRim_coil1_sam_shft,SNRim_coil2_sam_shft,SNRim_coil3_sam)...
    ,cat(3,SNRim_coil4_sam,SNRim_coil6_sam_shft,SNRim_coil7_sam_shft)));
h = colorbar; caxis([0 250]); ylabel(h,'SNR')
hold on; text(1+(0:nx:(nx*3-1)),ones(1,3)*(ny-1),...
    {'Coil #1: Birdcage','Coil #2: Birdcage','Coil #3: Helmholtz'},...
    'fontsize',imfontsize,'color',[1 1 1]);
text(1+(0:nx:(nx*3-1)),ones(1,3)*(2*ny-1),...
    {'Coil #4: 8-CH array','Coil #6: 14-CH array',...
    'Coil #7: Surface coil'},'fontsize',imfontsize,'color',[1 1 1]);
hold off; set(gca,'fontsize',imfontsize)

% SNR profiles cylinder phantom
AP = -360/2+15:15:360/2; RL = -360/2+15:15:360/2;
APidx = 13-1:13+1; RLidx = 12-1:12+1;
profilelegend = {'Coil #1','Coil #2','Coil #3','Coil #4','Coil #5',...
    'Coil #6','Coil #7'};

SNRprofile_cyl_AP = squeeze(mean(cat(3,SNRim_coil1_cyl_shft(:,APidx),...
    SNRim_coil3_cyl(:,APidx),SNRim_coil4_cyl(:,APidx),...
    SNRim_coil5_cyl(:,APidx),SNRim_coil6_cyl(:,APidx),...
    SNRim_coil7_cyl(:,APidx)),2));
SNRprofile_cyl_RL = squeeze(mean(cat(3,SNRim_coil1_cyl_shft(RLidx,:),...
    SNRim_coil3_cyl(RLidx,:),SNRim_coil4_cyl(RLidx,:),...
    SNRim_coil5_cyl(RLidx,:),SNRim_coil6_cyl(RLidx,:),...
    SNRim_coil7_cyl(RLidx,:)),1)).';

figure; colorder = get(gca,'colororder');
subplot(121); 
set(gca,'colororder',colorder([1,3:7],:),'NextPlot','ReplaceChildren'); 
plot(AP,SNRprofile_cyl_AP,'linewidth',1.5)
axis([-inf inf 0 250]); ylabel 'SNR'; xlabel 'A-P [mm]'
set(gca,'fontsize',imfontsize)
subplot(122); 
set(gca,'colororder',colorder([1,3:7],:),'NextPlot','ReplaceChildren'); 
plot(RL,SNRprofile_cyl_RL,'linewidth',1.5)
axis([-inf inf 0 250]); ylabel 'SNR'; xlabel 'R-L [mm]'
legend(profilelegend([1,3:7]),'location','northwest')
set(gca,'fontsize',imfontsize)

% SNR profiles SAM phantom
APidx = 13-1:13+1; RLidx = 13-1:13+1;

SNRprofile_sam_AP = squeeze(mean(cat(3,SNRim_coil1_sam_shft(:,APidx),...
    SNRim_coil2_sam_shft(:,APidx),SNRim_coil3_sam(:,APidx),...
    SNRim_coil4_sam(:,APidx),SNRim_coil6_sam_shft(:,APidx),...
    SNRim_coil7_sam_shft(:,APidx)),2));
SNRprofile_sam_RL = squeeze(mean(cat(3,SNRim_coil1_sam_shft(RLidx,:),...
    SNRim_coil2_sam_shft(RLidx,:),SNRim_coil3_sam(RLidx,:),...
    SNRim_coil4_sam(RLidx,:),SNRim_coil6_sam_shft(RLidx,:),...
    SNRim_coil7_sam_shft(RLidx,:)),1)).';

figure; 
subplot(121); 
set(gca,'colororder',colorder([1:4,6:7],:),'NextPlot','ReplaceChildren'); 
plot(AP,SNRprofile_sam_AP,'linewidth',1.5)
axis([-inf inf 0 250]); ylabel 'SNR'; xlabel 'A-P [mm]';
set(gca,'fontsize',imfontsize)
subplot(122); 
set(gca,'colororder',colorder([1:4,6:7],:),'NextPlot','ReplaceChildren');
plot(RL,SNRprofile_sam_RL,'linewidth',1.5)
axis([-inf inf 0 250]); ylabel 'SNR'; xlabel 'R-L [mm]'
legend(profilelegend([1:4,6:7]),'location','northwest')
set(gca,'fontsize',imfontsize)

%% noise correlation matrices for array coils

figure; 
subplot(1,3,1); imagesc(abs(ncorr_coil4_cyl)); axis square; title 'Coil #4'
hold on; text(4.5,9.5,sprintf('max=%.1f %%; min=%.1f %%; mean=%.1f %%',...
    max(vec(abs(ncorr_coil4_cyl)-eye(8)))*100,min(vec(abs(...
    ncorr_coil4_cyl)))*100,mean(vec(abs(ncorr_coil4_cyl)-eye(8)))*100),...
    'fontsize',imfontsize,'horizontalalignment','center','fontweight','bold')
set(gca,'fontsize',imfontsize)
subplot(1,3,2); imagesc(abs(ncorr_coil5_cyl)); axis square; title 'Coil #5'
hold on; text(8.5,18.5,sprintf('max=%.1f %%; min=%.1f %%; mean=%.1f %%',...
    max(vec(abs(ncorr_coil5_cyl)-eye(16)))*100,min(vec(abs(...
    ncorr_coil5_cyl)))*100,mean(vec(abs(ncorr_coil5_cyl)-eye(16)))*100),...
    'fontsize',imfontsize,'horizontalalignment','center','fontweight','bold')
set(gca,'fontsize',imfontsize)
subplot(1,3,3); imagesc(abs(ncorr_coil6_cyl)); axis square; title 'Coil #6'
hold on; text(7.5,16.25,sprintf('max=%.1f %%; min=%.1f %%; mean=%.1f %%',...
    max(vec(abs(ncorr_coil6_cyl)-eye(14)))*100,min(vec(abs(...
    ncorr_coil6_cyl)))*100,mean(vec(abs(ncorr_coil6_cyl)-eye(14)))*100),...
    'fontsize',imfontsize,'horizontalalignment','center','fontweight','bold')
set(gca,'fontsize',imfontsize)

%% gmap simulations
% for Cartesian undersampling factor R = 2 and R = 4

% coil #4, cylindrical phantom
csi_path = 'data/coil4_cyl_8ch.mat'; noise_path = []; noiseBW = 0.845;
[~,~,~,~,gmapR2_coil4_cyl] = SNReval_carbon13coils(...
    csi_path,noise_path,noiseBW,QAeval,est_T2s,smap,2,false);
[~,~,~,~,gmapR4_coil4_cyl] = SNReval_carbon13coils(...
    csi_path,noise_path,noiseBW,QAeval,est_T2s,smap,4,false);
% coil #5, cylindrical phantom
csi_path = 'data/coil5_cyl.mat'; noise_path = []; noiseBW = 0.845; 
[~,~,~,~,gmapR2_coil5_cyl] = SNReval_carbon13coils(...
    csi_path,noise_path,noiseBW,QAeval,est_T2s,smap,2,true);
[~,~,~,~,gmapR4_coil5_cyl] = SNReval_carbon13coils(...
    csi_path,noise_path,noiseBW,QAeval,est_T2s,smap,4,true);
% coil #6, cylindrical phantom
csi_path = 'data/coil6_cyl.mat'; noise_path = []; noiseBW = 0.845; 
[~,~,~,~,gmapR2_coil6_cyl] = SNReval_carbon13coils(...
    csi_path,noise_path,noiseBW,QAeval,est_T2s,smap,2,false);
[~,~,~,~,gmapR4_coil6_cyl] = SNReval_carbon13coils(...
    csi_path,noise_path,noiseBW,QAeval,est_T2s,smap,4,false);

figure; subplot(211)
imagesc_row(cat(3,gmapR2_coil4_cyl,gmapR2_coil5_cyl,gmapR2_coil6_cyl));
colorbar; caxis([1 3]); colormap(brighten(gray,0.5))
set(gca,'fontsize',imfontsize);
hold on; text(nx/2+(0:nx:(nx*3-1)),ones(1,3)*(ny-1),...
    {sprintf('max=%.2f; mean=%.2f',max(gmapR2_coil4_cyl(:)),nanmean(gmapR2_coil4_cyl(:))),...
    sprintf('max=%.2f; mean=%.2f',max(gmapR2_coil5_cyl(:)),nanmean(gmapR2_coil5_cyl(:))),...
    sprintf('max=%.2f; mean=%.2f',max(gmapR2_coil6_cyl(:)),nanmean(gmapR2_coil6_cyl(:)))},...
    'fontsize',imfontsize,'color',[1 1 1],'horizontalalignment','center');
subplot(212)
imagesc_row(cat(3,gmapR4_coil4_cyl,gmapR4_coil5_cyl,gmapR4_coil6_cyl));
colorbar; caxis([1 3]); colormap(brighten(gray,0.5))
set(gca,'fontsize',imfontsize)
hold on; text(nx/2+(0:nx:(nx*3-1)),ones(1,3)*(ny-1),...
    {sprintf('max=%.2f; mean=%.2f',max(gmapR4_coil4_cyl(:)),nanmean(gmapR4_coil4_cyl(:))),...
    sprintf('max=%.2f; mean=%.2f',max(gmapR4_coil5_cyl(:)),nanmean(gmapR4_coil5_cyl(:))),...
    sprintf('max=%.2f; mean=%.2f',max(gmapR4_coil6_cyl(:)),nanmean(gmapR4_coil6_cyl(:)))},...
    'fontsize',imfontsize,'color',[1 1 1],'horizontalalignment','center');

