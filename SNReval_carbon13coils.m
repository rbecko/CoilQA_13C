function [SNRim,ncov,ncorr,T2s,gmap] = SNReval_carbon13coils(csi_path,...
    noise_path,noiseBW,QAeval,est_T2s,smap,R,rot,decorr)
%% SNREVAL_CARBON13COILS reconstructs EG-CSI data in SNR units
% EG = ethylene glycol
% The function accepts GE raw p-files (.7 data), Siemens twix-files (.dat 
% data), and mat-files with the following data: 
%           d = [fid points x nx x ny x nz x coils] for the raw csi data
%       noise = [noise sample points x 1] for the raw noise data
%   sw_signal = [1x1] for the spectral bandwidth in Hz
%       txfrq = [1x1] for the transmit frequency in MHz
%          te = [1x1] for the echo time in s
%     dataext = 'string' with either '.7' or '.dat' indicating original
%                  raw data extension
%
% The function currently only accepts multi-channel data from .7 data or
% mat-files
%
% example: [SNRim,ncov,ncorr,T2s,gmap] = SNReval_carbon13coils(csi_path,...
%               noise_path,noiseBW,QAeval,est_T2s,smap,R,rot,decorr)
%
% input:   csi_path  path to CSI raw data
%        noise_path  path to noise raw data
%          noise_BW  noise equivalent bw
%            QAeval  true = plot figures for QA evalution,
%                       false = do not plot figures (opt,default=true)
%           est_T2s  true = estimate T2* based on FWHM, false = skip
%                       (optional, default=false)
%              smap  true = estimate sensitivity maps for coil combination,
%                       false = do sum-of-squares combination
%                       (optional, default=false)
%                 R  acceleration factor >1 for generation of g-factor map
%                       for L/R acceleration (optional, default=1)
%               rot  if true, rotate acceleration pattern by 90 to generate
%                       g-factor map for A/P acceleration instead of L/R
%                       (optional, default=false)
%            decorr  include noise correlations in coil combination for
%                       multi-channel coils
%                       (optional, default=true)
%
% output:     SNRim  SNR image
%              ncov  noise covariance matrix
%             ncorr  noise correlation matrix
%               T2s  T2* in seconds
%              gmap  g-factor map
%
% last edited: October 2019, Rie Beck Olin

%% add paths to required functions
addpath(genpath('utils'))

%% input variables
if ~exist('noise_path','var'),  noise_path = []; end
if ~exist('QAeval','var'),      QAeval = true; end
if ~exist('est_T2s','var'),     est_T2s = false; end
if ~exist('smap','var'),        smap = false; end
if ~exist('R','var'),           R = []; end
if isempty(R),                  R = 1; end
if ~exist('rot','var'),         rot = false; end
if ~exist('decorr','var'),      decorr = true; end

%% read all data
disp('loading data')

[~,~,ext] = fileparts(csi_path);

if strcmp(ext,'.dat')
    raw = io_loadspec_twix(csi_path);
    raw_noise = io_loadspec_twix(noise_path);
    d = raw.fids; noise = vec(double(raw_noise.fids));
    % remove zero-filled fids
    d = double(d(:,any(sum(d),2),any(sum(d),3)));
    % extract data dimensions
    ns = raw.sz(1); nx = size(d,2); ny = size(d,3); nz = 1; nc = 1;
    sw_signal = raw.spectralwidth;
    txfrq = raw.txfrq*1e-6; % Tx frequency in MHz
    te = raw.te*1e-3; % echo time in s
elseif strcmp(ext,'.7')
    [raw,hdr] = read_p(csi_path);
    [raw_noise,~] = read_p(noise_path);
    noise = squeeze(double(raw_noise));
    % extract data dimensions
    ns = hdr.rdb_hdr.da_xres; nx = hdr.rdb_hdr.xcsi;
    ny = hdr.rdb_hdr.ycsi; nz = hdr.rdb_hdr.zcsi;
    nc = (hdr.rdb_hdr.dab(2)-hdr.rdb_hdr.dab(1))+1; % Rx coil elements
    if nc == 1; noise = vec(noise); end
    sw_signal = hdr.rdb_hdr.spectral_width;
    txfrq = hdr.rdb_hdr.ps_mps_freq*1e-7; % Tx frequency in MHz
    te = hdr.rdb_hdr.user25*1e-6; % echo time in s
    % unchop data
    if mod(hdr.rdb_hdr.data_collect_type,2) == 0
        raw(2:2:end,:) = -raw(2:2:end,:);
    end
    % reorder data
    d = zeros(ns,nx,ny,nz,nc);
    for l1 = 1:nx
        for l2 = 1:ny
            for l3 = 1:nz
                ll1 = l1+(l2-1)*nx+(l3-1)*ny*nx;
                for l4 = 1:nc
                    d(:,l1,l2,l3,l4) = raw(ll1,:,1,1,1,l4).';
                end
            end
        end
    end
elseif strcmp(ext,'.mat')
    load(csi_path)
    % check if all data were loaded
    if ~exist('d','var')||~exist('noise','var')||~exist('sw_signal',...
            'var')||~exist('te','var')||~exist('txfrq','var')...
            ||~exist('dataext','var')
        warning('mat-file does not contain all required variables')
        return
    end
    % extract data dimensions
    [ns,nx,ny,nz,nc] = size(d); 
else
    warning('data extension is not .7, .dat, or .mat')
    return
end

% frequency axis
hz = -(-ns/2:ns/2-1)/ns*sw_signal;
ppm = -hz/txfrq;
% time axis
t = (0:size(d,1)-1)/sw_signal;
% check if 2D or 3D spatial data
if nz > 1; dim = 3; else dim = 2; end

fprintf('TE = %.3g ms\n',te*1e3)

%% data pre-processing and signal extraction
% box filter
box_fil = (ppm > -10 & ppm < 10).';
box_fil = box_fil*1/sqrt((1/length(box_fil))*sum(box_fil.^2));

% spatial reconstruction
for n = 2:(dim+1)
    d = sqrt(size(d,n))*ifft(fftshift(d,n),[],n);
    d = ifftshift(d,n);
end

% spectral reconstruction
spec = 1/sqrt(ns)*ifftshift(fft(d,[],1),1);
spec_fil = spec.*repmat(box_fil,1,nx,ny,1,nc); % filtering
d_fil = sqrt(ns)*ifft(fftshift(spec_fil,1),[],1);

% signal extraction in the time domain based on J-coupling and echo time
jc = 1/142; % J-coupling time constant for EG
if jc > te; temp_shift = jc-te; else temp_shift = 2*jc-te; end
tmp = reshape(abs(d_fil),ns,nx*ny*nc); [tx,tidx] = max(max(tmp,[],1),[],2);
[~,lc] = findpeaks(tmp(:,tidx),t,'NPeaks',10,'MinPeakHeight',tx/2,...
    'MinPeakDistance',3*temp_shift/4);

jc = mean(diff(lc));
fprintf('J-coupling = %.4g Hz\n',1/jc)
if jc > te; temp_shift = jc-te; else temp_shift = 2*jc-te; end

sig_idx = round(temp_shift*sw_signal+1);
sig = squeeze(d_fil(sig_idx,:,:,:,:));

% image orientation correction
if strcmp(ext,'.dat')||strcmp(dataext,'.dat')
    sig = permute(sig,[2,1]);
elseif strcmp(ext,'.7')||strcmp(dataext,'.7')
    tmp = zeros(size(sig));
    for l1 = 1:nx
        for l2 = 1:ny
            ll1 = nx-l1+1; ll2 = l2;
            tmp(ll1,ll2,:) = sig(l1,l2,:);
        end
    end
    sig = tmp;
end

if QAeval == 1 && nc > 1
    figure; imagesc_row(abs(sig),[],[],true)
    title 'magnitude signal image(s)'
end

%% noise evaluation
% noise covariance matrix
ncov = cov(noise);
% noise correlation matrix
ncorr = corrcoef(noise);

if nc > 1 && QAeval == 1
    figure; subplot(121)
    imagesc(abs(ncov)), colorbar
    axis square, caxis([0 max(abs(ncov(:)))])
    title 'noise covariance matrix'
    subplot(122); imagesc(abs(ncorr)), colorbar
    axis square, caxis([0 1])
    title 'noise correlation matrix'
    
    % From: Tunnicliffe, E., Graves, M. J., & Robson, M. D. (2011).
    % Use of the noise covariance matrix in array coil quality assurance.
    % In Proc. Intl. Soc. Mag. Reson. Med. (p. 4548).
    % "An appreciable (>15%) drop in SNR was found in images acquired with
    % arrays meeting at least one of the following conditions:
    % 1) The variation of on-diagonal real elements in the NCM is >35%;
    % 2) The largest off-diagonal element is >25% of the largest 
    %    on-diagonal element"
    maxvar = max(abs(diag(real(ncov))-mean(diag(real(ncov))))...
        /mean(diag(real(ncov)))*100);
    maxoff = max(abs(ncov(eye(size(ncov))<1)))...
        /max(abs(ncov(eye(size(ncov))==1)))*100;
    
    disp('ARRAY COIL QA BASED ON NOISE COVARIANCE MATRIX')
    fprintf(['\nGreatest variation from mean of on-diagonal real ',...
        'elements: %.2f %% (should be < 35 %%)\nLargest off-diagonal ',...
        'element relative to largest on-diagonal element: %.2f %% ',...
        '(should be < 25 %%)\n\n'],maxvar,maxoff)
    disp('press any key to continue')
    pause
end

% scale noise with respect to noise equivalent BW factor
noise_cov_scaled = ncov/noiseBW;

%% estimate SNR
% normalize signal data with respect to the noise variance
sig_norm = sig./repmat(reshape(sqrt(diag(noise_cov_scaled)),...
    [1,1,nc]),nx,ny,1);

if nc > 1 && smap == 1
    % data-driven sensitivity-based coil-combination 
    % Perona-Malik anisotropic difussion filtering
    b1_pmfilt = zeros(size(sig_norm));
    for n = 1:nc
        maxsnr = max(abs(vec(sig_norm(:,:,n))));
        b1_pmfilt(:,:,n) = perona_malik(sig_norm(:,:,n),.05*maxsnr,...
            750/maxsnr);
    end
    if decorr == 1
        % include noise correlations
        Rinv = pinv(noise_cov_scaled);
    else
        Rinv = pinv(eye(size(noise_cov_scaled)).*noise_cov_scaled);
    end
    sigvec = reshape(sig,nx*ny,nc).';
    b1vec = reshape(b1_pmfilt,nx*ny,nc).';
    SNRcomplex_vec = zeros(nx*ny,1);
    for n = 1:(nx*ny)
        SNRcomplex_vec(n) = sigvec(:,n).'*Rinv*conj(b1vec(:,n))./...
            sqrt(b1vec(:,n).'*Rinv*conj(b1vec(:,n)));
    end
    SNRcomplex = reshape(SNRcomplex_vec.',nx,ny);
elseif nc > 1 && smap == 0
    SNRcomplex = sos(sig_norm);
else
    SNRcomplex = sig_norm;
end

SNRim = sqrt(2)*abs(SNRcomplex);

if QAeval == 1
    if nc > 1
        figure; imagesc_row(sqrt(2)*abs(sig_norm),[],[],true)
        title(['SNR images for individual coil channels (' ...
            num2str(nc) ' channels)']); 
        colorbar
    end
    
    figure; imagesc_row(SNRim)
    title 'SNR image'; colorbar
end

%% extra: T2* estimation
if est_T2s == 1
    disp('estimating T2* based on FWHM')
    % spectral reconstruction with zero-filling
    spec_0lb = 1/sqrt(ns)*ifftshift(fft(d,ns*4,1),1);
    % image orientation correction
    if strcmp(ext,'.dat')
        spec_0lb = permute(spec_0lb,[1,3,2]);
    elseif strcmp(ext,'.7')
        tmp = zeros(size(spec_0lb));
        for l1 = 1:nx
            for l2 = 1:ny
                ll1 = nx-l1+1; ll2 = l2;
                tmp(:,ll1,ll2,:,:,:) = spec_0lb(:,l1,l2,:,:,:);
            end
        end
        spec_0lb = tmp;
    end

    max_spec_0lb = squeeze(max(abs(spec_0lb),[],1));
    hz_int = interp(hz,4);
    T2s_all = zeros(nx,ny);
    for n = 1:nx
        for m = 1:ny
            for l = 1:nc
                [~,~,widths] = findpeaks(abs(spec_0lb(:,n,m,:,l)),...
                    -hz_int,'SortStr','descend','MinPeakHeight',...
                    .5*max_spec_0lb(n,m,l));
                T2s_all(n,m,l) = 1/widths(1);
            end
        end
    end
    
    % image orientation correction
    if strcmp(ext,'.dat')||strcmp(dataext,'.dat')
        T2s_all = permute(T2s_all,[2,1]);
    elseif strcmp(ext,'.7')||strcmp(dataext,'.7')
        tmp = zeros(size(T2s_all));
        for l1 = 1:nx
            for l2 = 1:ny
                ll1 = nx-l1+1; ll2 = l2;
                tmp(ll1,ll2,:) = T2s_all(l1,l2,:);
            end
        end
        T2s_all = tmp;
    end
    
    if nc > 1
        % remove T2 star estimate based on individual coil SNR values
        tmp = T2s_all; tmp(sig_norm<2) = nan; 
        tmp(repmat(SNRim<10,1,1,nc)) = nan;
        % average across coil elements for multi-channel coils
        T2s = nanmean(tmp,3);
    else
        tmp = T2s_all; tmp(SNRim<10) = nan;
        T2s = tmp;
    end
    if QAeval == 1
        figure; imagesc_row(T2s*1e3);
        h = colorbar; ylabel(h,'T_2* [ms]'); caxis([0 150])
        title(sprintf('mean T2* within mask = %.3g ms',...
            nanmean(T2s(:)*1e3)))
    end
else
    T2s = [];
end

%% extra: gmap estimation for multi-channel coils
if smap == 1 && ~isempty(R) && R>1
    disp(['estimating g-factors for R=' num2str(R)])
    mask = SNRim>10;
    if rot == true
        b1_pmfilt = imrotate(b1_pmfilt,90,'crop');
        mask = imrotate(mask,90,'crop');
    end
    [~,gmap] = ismrm_calculate_sense_unmixing(R, b1_pmfilt);
    % mask gmap
    gmap = gmap.*mask;
    gmap(gmap == 0) = nan;
    if rot == true
        gmap = imrotate(gmap,-90,'crop');
    end
    if QAeval == 1
        figure; imagesc_row(gmap);
        colorbar; caxis([1 max(gmap(:))])
        title(sprintf('mean g-factor = %.3g, max g-factor = %.3g',...
            nanmean(gmap(:)),max(gmap(:))))
    end
else
    gmap = [];
end




