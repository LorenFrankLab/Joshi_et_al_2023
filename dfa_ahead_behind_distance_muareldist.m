function out = dfa_ahead_behind_distance_muareldist(index, excludeperiods, eeg, eeggnd, ~, posdlc, varargin)

    % dfa_ahead_behind_distance_muareldist
    % this code helps you select the forelimb stances that are accompanied with
    % clean decoding estimates. 

%   Any epochs that have movement decoding done, compute:
%	- ahead/behind distance of the peak posterior density to the current body position 
% 	- computes if the ahead-behind distance is closer to the current position of the animal than other times 
% 	- loads posterior, posdlc, linpos etc. WTRACK

    % define defaults 
    appendindex = 0;

    % process varargin if present and overwrite default values
    if (~isempty(varargin))
        assign(varargin{:});
    else
        warning('thresholds empty - cant compute!')
    end
    
    % load the day and epoch number 
    d = index(1);
    e = index(2);
    reftet=index(3);
    
    min_peak_dist=0.1; % for longsweep detection 
    nan_cutoff=10; % for exclusion of a longsweep if data around a plant consists of more than 5 nans 
    nan_eval_window=0.06; % for excusion of a longsweep if data around a pant for nan_cutoff
    buffer=0.06; % make the cycles of rel distance from these peak times 
    interpnan=5; %interpolate over NaNs < 10ms long
    
    % include periods that have velocity above theset threshold
    for i=1:length(excludeperiods)-1
        run_periods(i,1)=excludeperiods(i,2);
        run_periods(i,2)=excludeperiods(i+1,1);
    end
    
    % Load LFP 
    theta_data=eeggnd{1,d}{1,e}{1,reftet}.data;
    theta_time=(eeggnd{1,d}{1,e}{1,reftet}.starttime):(1/eeg{1,d}{1,e}{1,reftet}.samprate):(eeg{1,d}{1,e}{1,reftet}.endtime);

    if ~(size(theta_time,2)==size(theta_data,1))
        theta_data=theta_data(1:size(theta_time,2));
    else 
    end
   
    % Filter LFP
    load('thetafilter.mat');
    theta_filtered_lfp=filtfilt(thetafilter.tf.num,thetafilter.tf.den,double(int16(theta_data)));
    HT = hilbert(theta_filtered_lfp);
    amplitude = sqrt(real(HT).^2 + imag(HT).^2);
    phase_hpc = angle(HT)+pi;

%% load posteriors
    post_path = ''; 
    postfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_cv_classifier_clusterless_vel_0_nose_alltime5x_results.nc',post_path,animal,animal,d,e);  %
   
    if exist(postfile)
        posteriorts_wtrack = ncread(postfile,'time'); % time bins 
        postposbins_wtrack=1+ncread(postfile,'position'); % position on w
        posterior_causal=ncread(postfile,'causal_posterior'); % causal posterior laod 
        posterior_acausal=ncread(postfile,'acausal_posterior'); % acausal posterior load 

        linposfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_linearised_position_nose.nc', post_path, animal, animal, d, e); %
        linpos_nose_wtrack = ncread(linposfile,'linear_position'); % load linpos

        rel_dist_causal=ncread(postfile,'rel_distance_from_animal_position_causal');
        rel_dist_acausal=ncread(postfile,'rel_distance_from_animal_position_acausal');

        abs_dist_causal=abs(rel_dist_causal);
        abs_dist_acausal=abs(rel_dist_acausal);

        credible_interval_95_causal=ncread(postfile,'credible_interval_95_causal');
        credible_interval_50_causal=ncread(postfile,'credible_interval_50_causal');

        credible_interval_95_acausal=ncread(postfile,'credible_interval_95_acausal');
        credible_interval_50_acausal=ncread(postfile,'credible_interval_50_acausal');

    %% Load dlc positions and times 
    dlc_results=posdlc{1,d}{1,e}; 

    %load ptp adjusted timestamps 
    cam_rt_fit=dlc_results.data(:,1);

    %load nose
    dlc_nose_x=dlc_results.data(:,2);
    dlc_nose_y=dlc_results.data(:,3);
    dlc_nose_vel=dlc_results.data(:,4);

    %load tail 
    dlc_tail_x=dlc_results.data(:,5);
    dlc_tail_y=dlc_results.data(:,6);
    dlc_tail_vel=dlc_results.data(:,7);

    %forepawL
    dlc_forepawL_x=dlc_results.data(:,14);
    dlc_forepawL_y=dlc_results.data(:,15);

    %forepawR
    dlc_forepawR_x= dlc_results.data(:,17);
    dlc_forepawR_y= dlc_results.data(:,18);

    %hindpawL
    dlc_hindpawL_x=dlc_results.data(:,20);
    dlc_hindpawL_y=dlc_results.data(:,21);

    %hindpawR
    dlc_hindpawR_x=dlc_results.data(:,23);
    dlc_hindpawR_y=dlc_results.data(:,24);

    dlc_n_records = size(dlc_hindpawL_x,1);
    fprintf('# of DeepLabCut timestamps: %d \n', dlc_n_records)

    % estimated framerate based on camera time 
    est_framerate=median(1./diff(cam_rt_fit));

    % forepawL
    [cam_rt_fit_forepawL_results, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_nose_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, run_periods);

    % HindpawR
    [cam_rt_fit_hindpawR_results, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_nose_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, run_periods);
 
    % ForepawR 
    [cam_rt_fit_forepawR_results, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_nose_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, run_periods);
  
    % HindpawL
    [cam_rt_fit_hindpawL_results, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_nose_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, run_periods);

    %% collect the forelimb plant times in a particular window 
    forelimb_plant=[cam_rt_fit_forepawR_results.midstance; cam_rt_fit_forepawL_results.midstance];
    hindlimb_plant=[cam_rt_fit_hindpawR_results.midstance; cam_rt_fit_hindpawL_results.midstance];

    %sort and remove duplicates 
    forelimb_plant=sort(forelimb_plant); idx_remove=find(diff(forelimb_plant)==0);forelimb_plant(idx_remove)=[]; clear idx_remove;
    hindlimb_plant=sort(hindlimb_plant); idx_remove=find(diff(hindlimb_plant)==0);hindlimb_plant(idx_remove)=[]; clear idx_remove;

    linpos_min_max=[out_min out_max]; % collect stance times only within a window 
    outboundvec = (linpos_nose_wtrack > linpos_min_max(:,1)) & (linpos_nose_wtrack < linpos_min_max (:,2));
    outbound_periods_pos = vec2list(outboundvec,posteriorts_wtrack);

    forelimb_plants_outbound={};
    outbound_ind={};

    outbound_center=[];
    inbound_center=[];
    diff_linpos=[];

    for i=1:length(outbound_periods_pos)
         ind_linpos_start_end=find(posteriorts_wtrack>outbound_periods_pos(i,1)&posteriorts_wtrack<outbound_periods_pos(i,2));

    if ~isempty(ind_linpos_start_end)
         diff_linpos=linpos_nose_wtrack(ind_linpos_start_end(1)) - linpos_nose_wtrack(ind_linpos_start_end(end));

    if  ((abs(diff_linpos))<15)
        outbound_center(i,1)= NaN;
        outbound_center(i,2)=NaN;
        inbound_center(i,1)=NaN;
        inbound_center(i,2)=NaN;  
    else

    if  (diff_linpos < 0)
        outbound_center(i,1)= outbound_periods_pos(i,1);
        outbound_center(i,2)=outbound_periods_pos(i,2);
        inbound_center(i,1)=NaN;
        inbound_center(i,2)=NaN; 
    else
        outbound_center(i,1)= NaN;
        outbound_center(i,2)=NaN;
        inbound_center(i,1)=outbound_periods_pos(i,1);
        inbound_center(i,2)=outbound_periods_pos(i,2);      
    end

    end
    else 
        outbound_center(i,1)= NaN;
        outbound_center(i,2)=NaN;
        inbound_center(i,1)=NaN;
        inbound_center(i,2)=NaN;  
        continue 
    end
    clear diff_linpos
    end

    outbound_center(any(isnan(outbound_center), 2), :) = [];
    inbound_center(any(isnan(inbound_center), 2), :) = [];

    %% Collect forelimb and hindlimb stance times for outbound center and inbound center 
    % forelimb
    for i=1:length(outbound_center)
      ind=find(forelimb_plant>(outbound_center(i,1))&(forelimb_plant<(outbound_center(i,2))));
      forelimb_plant_outbound(i)={forelimb_plant(ind)};
      clear ind
    end
    forelimb_plant_outbound=cell2mat(forelimb_plant_outbound(:));

    for i=1:length(inbound_center)
      ind=find(forelimb_plant>(inbound_center(i,1))&(forelimb_plant<(inbound_center(i,2))));
      forelimb_plant_inbound(i)={forelimb_plant(ind)};
      clear ind
    end
    forelimb_plant_inbound=cell2mat(forelimb_plant_inbound(:));

    % hindlimb
    for i=1:length(outbound_center)
      ind=find(hindlimb_plant>(outbound_center(i,1))&(hindlimb_plant<(outbound_center(i,2))));
      hindlimb_plant_outbound(i)={hindlimb_plant(ind)};
      clear ind
    end
    hindlimb_plant_outbound=cell2mat(hindlimb_plant_outbound(:));

    for i=1:length(inbound_center)
      ind=find(hindlimb_plant>(inbound_center(i,1))&(hindlimb_plant<(inbound_center(i,2))));
      hindlimb_plant_inbound(i)={hindlimb_plant(ind)};
      clear ind
    end
    hindlimb_plant_inbound=cell2mat(hindlimb_plant_inbound(:));
    
     %% collect all info that may have to be plotted later on
     cdat_ahead_behind_distance=imcont('data', rel_dist_causal, 'timestamp', posteriorts_wtrack);
     cdat_abs_distance=imcont('data', abs_dist_causal, 'timestamp', posteriorts_wtrack);
     cdat_ci=imcont('data', credible_interval_50_causal, 'timestamp', posteriorts_wtrack);
     cdat_eeggnd=imcont('data', theta_data,'timestamp', theta_time);
     cdat_theta_filt=imcont('data', theta_filtered_lfp, 'timestamp', theta_time);
     
     %% %create a lowpass filter: [12 25] to remove the wiggles in rel dist and ci 

     filtopt_12_25 = mkfiltopt('name', 'aj_filt', 'filttype', 'lowpass', 'F', [12 25]);
     filtopt_30_50 = mkfiltopt('name', 'aj_filt', 'filttype', 'lowpass', 'F', [30 50]);
     
     cdat_ci_smooth=contfilt(cdat_ci,'filtopt', filtopt_12_25);
     cdat_ahead_behind_distance_smooth=contfilt(cdat_ahead_behind_distance, 'filtopt', filtopt_30_50);
     cdat_abs_distance_smooth=contfilt(cdat_abs_distance, 'filtopt', filtopt_30_50);
     
    % spotchecked that a threshold of 50cm seems reasonable 
    for i=1:size(posteriorts_wtrack)
        if credible_interval_50_causal(i)<ci_thresh
            cdat_ahead_behind_distance_smooth.data(i)=cdat_ahead_behind_distance_smooth.data(i);
            cdat_abs_distance_smooth.data(i)=cdat_abs_distance_smooth.data(i);
        else
            cdat_ahead_behind_distance_smooth.data(i)=NaN;
            cdat_abs_distance_smooth.data(i)=NaN;
        end
    end
    
%% now collect forelimb stances in this window 
% select times where the relative distances are bigger then 5 cms 
% keyboard
    for select = 1:size(outbound_center,1)
        ind_start_end=find(posteriorts_wtrack>outbound_center(select,1)&posteriorts_wtrack<outbound_center(select,2));
        ahead_behind_dist_include=cdat_ahead_behind_distance_smooth.data(ind_start_end);
        ci_50_include=credible_interval_50_causal(ind_start_end);
        time_select=posteriorts_wtrack(ind_start_end);
        
        if ((nanmean(ci_50_include)>ci_mean) || (abs(nanmean(ahead_behind_dist_include))>reldist_mean))
            cdat_ahead_behind_distance_smooth.data(ind_start_end)=NaN;
            cdat_abs_distance_smooth.data(ind_start_end)=NaN;
        else 
            cdat_ahead_behind_distance_smooth.data(ind_start_end)=cdat_ahead_behind_distance_smooth.data(ind_start_end);
        end
        
        clear ahead_behind_dist_include ind_start_end time_select ci_50_include
    end
    
    for select = 1:size(inbound_center,1)
        ind_start_end=find(posteriorts_wtrack>inbound_center(select,1)&posteriorts_wtrack<inbound_center(select,2));
        ahead_behind_dist_include=cdat_ahead_behind_distance_smooth.data(ind_start_end);
        ci_50_include=credible_interval_50_causal(ind_start_end);
        time_select=posteriorts_wtrack(ind_start_end);
        
        if ((nanmean(ci_50_include)>ci_mean) || (abs(nanmean(ahead_behind_dist_include))>reldist_mean))
            cdat_ahead_behind_distance_smooth.data(ind_start_end)=NaN;
        else 
            cdat_ahead_behind_distance_smooth.data(ind_start_end)=cdat_ahead_behind_distance_smooth.data(ind_start_end);
        end
        
        clear ahead_behind_dist_include ind_start_end time_select ci_50_include
    end
     
    forelimb_plant_inbound_mua=forelimb_plant_inbound;
    forelimb_plant_outbound_mua=forelimb_plant_outbound;
    
    %% Interpolate through NaNs in the ahead behind distance - try on July 7 - this works well. interpnan =5 corresponds to regions with NaNs <10ms 
    % find regions with NaNs <10ms big  i.e 5samples at 500 fps resolution

    % Find consecutive rows of NaN that exceed threshold number allowed
    hasNan = all(isnan(cdat_ahead_behind_distance_smooth.data),2);   
    dIdx = find(diff([0;hasNan;0]==1));     %rows that change 1/0
    s1 = dIdx(1:2:end-1);                   %start indices of 1s
    s2 = dIdx(2:2:end);                     %stop indices of 1s
    keepNan = (s2-s1)>=interpnan;                  %which segments have too many consecutive nans
    hasNan(cell2mat(arrayfun(@(x1,x2)x1:x2, s1(~keepNan),s2(~keepNan)-1,'UniformOutput',false)')) = false; 
    % Interp all missing values
    % TT2intrp = retime(cdat_ahead_behind_distance_smooth.data,'regular','linear','TimeStep',seconds(16));
    cdat_ahead_behind_distance_smooth.data = fillmissing(cdat_ahead_behind_distance_smooth.data,'linear');  % This gives you the same results as your line above
    % Replace the NaN values for >10 consecutive rows of missing data
    cdat_ahead_behind_distance_smooth.data(hasNan,:) = NaN; 

    % Find consecutive rows of NaN that exceed threshold number allowed
    hasNan = all(isnan(cdat_abs_distance_smooth.data),2);   
    dIdx = find(diff([0;hasNan;0]==1));     %rows that change 1/0
    s1 = dIdx(1:2:end-1);                   %start indices of 1s
    s2 = dIdx(2:2:end);                     %stop indices of 1s
    keepNan = (s2-s1)>=interpnan;                  %which segments have too many consecutive nans
    hasNan(cell2mat(arrayfun(@(x1,x2)x1:x2, s1(~keepNan),s2(~keepNan)-1,'UniformOutput',false)')) = false; 
    % Interp all missing values
    % TT2intrp = retime(cdat_ahead_behind_distance_smooth.data,'regular','linear','TimeStep',seconds(16));
    cdat_abs_distance_smooth.data = fillmissing(cdat_abs_distance_smooth.data,'linear');  % This gives you the same results as your line above
    % Replace the NaN values for >10 consecutive rows of missing data
    cdat_abs_distance_smooth.data(hasNan,:) = NaN; 

    %     keyboard
    %     % just sanity check plot for the data 
    %     close all; figure(1); hold on; clims=[0 0.3]
    %     imagesc(posteriorts_wtrack,postposbins_wtrack,squeeze(posterior_causal(:,1,:)), clims); colormap gray; cmap=colormap; cmap=flipud(cmap); colormap(cmap); set(gca,'YDir','normal');
    %     plot(posteriorts_wtrack,linpos_nose_wtrack,'y');plot(posteriorts_wtrack,credible_interval_50_causal,'r');
    %     plot(posteriorts_wtrack,cdat_ahead_behind_distance_smooth.data,'m')
    %     plot([outbound_center(:,1) outbound_center(:,1)], [-40 120], 'y');
    %     plot([outbound_center(:,2) outbound_center(:,2)], [-40 120], 'r');
    %     plot([forelimb_plant_outbound(:,1) forelimb_plant_outbound(:,1)], [-40 400], 'm');
    %     t_pre=0.1; t_post=0.1;
    %     periseg('cdat', contchans(cdat_ahead_behind_distance_smooth, 'chans', 1), 'segs', forelimb_plant_outbound, 't_pre',t_pre, 't_post',t_post);

%% only include those forelimb plants that do not have NaN in a window around them. NaNs correspond to data with high credible interval >50cm 
% check for multiple NaNs aorund a plant

% all plants outbound
forelimb_plant_outbound_window=[];
forelimb_plant_outbound_new=[];
forelimb_plant_outbound_window(:,1)=[forelimb_plant_outbound-nan_eval_window];
forelimb_plant_outbound_window(:,2)=[forelimb_plant_outbound+nan_eval_window];
for i=1:size(forelimb_plant_outbound_window,1)
    ind=find(posteriorts_wtrack>forelimb_plant_outbound_window(i,1) & posteriorts_wtrack<forelimb_plant_outbound_window(i,2));
    reldist_nan=sum(isnan(cdat_ahead_behind_distance_smooth.data(ind)));
    
    if reldist_nan>nan_cutoff % this is being very restrictive - but good for the initial test 
        forelimb_plant_outbound_new(i)=NaN;
    else
        forelimb_plant_outbound_new(i)=forelimb_plant_outbound(i);
    end

clear reldist_nan ind
end
clear forelimb_plant_outbound_window forelimb_plant_outbound
idx=isnan(forelimb_plant_outbound_new);forelimb_plant_outbound_new(idx)=[];clear idx;
forelimb_plant_outbound_reldist=forelimb_plant_outbound_new;

% all plants inbound
forelimb_plant_inbound_window=[];
forelimb_plant_inbound_new=[];
forelimb_plant_inbound_window(:,1)=[forelimb_plant_inbound-nan_eval_window];
forelimb_plant_inbound_window(:,2)=[forelimb_plant_inbound+nan_eval_window];
for i=1:size(forelimb_plant_inbound_window,1)
    ind=find(posteriorts_wtrack>forelimb_plant_inbound_window(i,1) & posteriorts_wtrack<forelimb_plant_inbound_window(i,2));
    reldist_nan=sum(isnan(cdat_ahead_behind_distance_smooth.data(ind)));
    
    if reldist_nan>nan_cutoff % this is being very restrictive - but good for the initial test 
        forelimb_plant_inbound_new(i)=NaN;
    else
        forelimb_plant_inbound_new(i)=forelimb_plant_inbound(i);
    end

clear reldist_nan ind
end
clear forelimb_plant_inbound_window forelimb_plant_inbound
idx=isnan(forelimb_plant_inbound_new);forelimb_plant_inbound_new(idx)=[];clear idx;
forelimb_plant_inbound_reldist=forelimb_plant_inbound_new;
forelimb_plant_reldist=sort([forelimb_plant_outbound_reldist';forelimb_plant_inbound_reldist']); % just include those forelimb plants that are in the run periods defined 

%% Now the bit of code that help us analyse the longest rel distance cycles 
cdat_ahead_behind_distance_smooth;

[~, pks_t]=findpeaks(cdat_ahead_behind_distance_smooth.data, posteriorts_wtrack,'MinPeakDistance', min_peak_dist,'MinPeakHeight', rel_dist_cutoff); 
[~, trs_t]=findpeaks(-cdat_ahead_behind_distance_smooth.data, posteriorts_wtrack, 'MinPeakDistance',min_peak_dist,'MinPeakHeight', -rel_dist_cutoff);

% close all;
% plot(posteriorts_wtrack,cdat_ahead_behind_distance_smooth.data); hold on
% plot([pks_t pks_t],[0 20], 'k');
% plot([trs_t trs_t], [-20 0], 'r');
   
for i=1:size(pks_t,1)
    trough_right_ind=find(trs_t>pks_t(i));
    trough_left_ind=find(trs_t<pks_t(i));
    
    if ~(isempty(trough_left_ind)||isempty(trough_right_ind))
        
    max_dist_cycles(i,1)=trs_t(trough_left_ind(end))-buffer;
    max_dist_cycles(i,2)=trs_t(trough_right_ind(1))-buffer;
    else
        continue
    end
    clear trough_right_ind trough_left_ind
end

for i=1:length(max_dist_cycles)
  ind=find(forelimb_plant_reldist>(max_dist_cycles(i,1))&(forelimb_plant_reldist<(max_dist_cycles(i,2))));
  forelimb_plant_max_dist(i)={forelimb_plant_reldist(ind)};
  clear ind
end
forelimb_plant_max_dist=cell2mat(forelimb_plant_max_dist(:));
forelimb_plant_max_dist(forelimb_plant_max_dist==0)=[];

%% collect forelimb plant during outbound or inbound runs

    for i=1:length(outbound_center)
      ind=find(forelimb_plant_max_dist>(outbound_center(i,1))&(forelimb_plant_max_dist<(outbound_center(i,2))));
      forelimb_plant_max_dist_outbound(i)={forelimb_plant_max_dist(ind)};
        
      ind_peaks=find(pks_t>(outbound_center(i,1))&(pks_t<(outbound_center(i,2))));
      ind_trs=find(trs_t>(outbound_center(i,1))&(trs_t<(outbound_center(i,2))));
      peaks_outbound_longsweep(i)={pks_t(ind_peaks)};
      troughs_outbound_longsweep(i)={trs_t(ind_trs)};
      clear ind ind_peaks ind_trs
    end
    forelimb_plant_max_dist_outbound=cell2mat(forelimb_plant_max_dist_outbound(:));
    peaks_outbound_longsweep=cell2mat(peaks_outbound_longsweep(:));
    troughs_outbound_longsweep=cell2mat(troughs_outbound_longsweep(:));
    
    for i=1:length(inbound_center)
      ind=find(forelimb_plant_max_dist>(inbound_center(i,1))&(forelimb_plant_max_dist<(inbound_center(i,2))));
      forelimb_plant_max_dist_inbound(i)={forelimb_plant_max_dist(ind)};
      
      ind_peaks=find(pks_t>(inbound_center(i,1))&(pks_t<(inbound_center(i,2))));
      ind_trs=find(trs_t>(inbound_center(i,1))&(trs_t<(inbound_center(i,2))));
      peaks_inbound_longsweep(i)={pks_t(ind_peaks)};
      troughs_inbound_longsweep(i)={trs_t(ind_trs)};
      clear ind ind_peaks ind_trs
    end
    forelimb_plant_max_dist_inbound=cell2mat(forelimb_plant_max_dist_inbound(:));
    peaks_inbound_longsweep=cell2mat(peaks_inbound_longsweep(:));
    troughs_inbound_longsweep=cell2mat(troughs_inbound_longsweep(:));
    
%% now the main analysis 
% compute the distance distribution of 
if ~isempty(forelimb_plant_max_dist_outbound)   
for i=1:length(forelimb_plant_max_dist_outbound)
    ind=find(posteriorts_wtrack>forelimb_plant_max_dist_outbound(i));
    forelimb_select_outbound_reldist(i)=cdat_ahead_behind_distance_smooth.data(ind(1));
    clear ind
end
else
    forelimb_select_outbound_reldist=nan;
end
    
if ~isempty(forelimb_plant_max_dist_inbound)
for i=1:length(forelimb_plant_max_dist_inbound)
    ind=find(posteriorts_wtrack>forelimb_plant_max_dist_inbound(i));
    forelimb_select_inbound_reldist(i)=cdat_ahead_behind_distance_smooth.data(ind(1));
    clear ind
end
else
    forelimb_select_inbound_reldist=nan;
end

if ~isempty(peaks_outbound_longsweep)
for i=1:length(peaks_outbound_longsweep)
    ind=find(posteriorts_wtrack>peaks_outbound_longsweep(i));
    peaks_select_outbound_reldist(i)=cdat_ahead_behind_distance_smooth.data(ind(1));
    clear ind
end
else
    peaks_select_outbound_reldist=nan;
end

if ~isempty(peaks_inbound_longsweep)   
for i=1:length(peaks_inbound_longsweep)
    ind=find(posteriorts_wtrack>peaks_inbound_longsweep(i));
    peaks_select_inbound_reldist(i)=cdat_ahead_behind_distance_smooth.data(ind(1));
    clear ind
end
else
    peaks_select_inbound_reldist=nan;
end

if ~isempty(troughs_outbound_longsweep)
for i=1:length(troughs_outbound_longsweep)
    ind=find(posteriorts_wtrack>troughs_outbound_longsweep(i));
    troughs_select_outbound_reldist(i)=cdat_ahead_behind_distance_smooth.data(ind(1));
    clear ind
end
else
    trough_select_outbound_reldist=nan;
end

if ~isempty(troughs_inbound_longsweep)   
for i=1:length(troughs_inbound_longsweep)
    ind=find(posteriorts_wtrack>troughs_inbound_longsweep(i));
    troughs_select_inbound_reldist(i)=cdat_ahead_behind_distance_smooth.data(ind(1));
    clear ind
end
else
    troughs_select_inbound_reldist=nan;
end

%% now the main analysis 
% compute the distance distribution of 
if ~isempty(forelimb_plant_max_dist_outbound)
    for i=1:length(forelimb_plant_max_dist_outbound)
        ind=find(posteriorts_wtrack>forelimb_plant_max_dist_outbound(i));
        forelimb_select_outbound_absdist(i)=cdat_abs_distance_smooth.data(ind(1));
        clear ind
    end
    else
        forelimb_select_outbound_absdist=nan;
end

if ~isempty(forelimb_plant_max_dist_inbound)
    for i=1:length(forelimb_plant_max_dist_inbound)
        ind=find(posteriorts_wtrack>forelimb_plant_max_dist_inbound(i));
        forelimb_select_inbound_absdist(i)=cdat_abs_distance_smooth.data(ind(1));
        clear ind
    end
    else
        forelimb_select_inbound_absdist=nan;
end

if ~isempty(peaks_outbound_longsweep)
    for i=1:length(peaks_outbound_longsweep)
        ind=find(posteriorts_wtrack>peaks_outbound_longsweep(i));
        peaks_select_outbound_absdist(i)=cdat_abs_distance_smooth.data(ind(1));
        clear ind
    end
    else
        peaks_select_outbound_absdist=nan;
end

if ~isempty(peaks_inbound_longsweep)
    for i=1:length(peaks_inbound_longsweep)
        ind=find(posteriorts_wtrack>peaks_inbound_longsweep(i));
        peaks_select_inbound_absdist(i)=cdat_abs_distance_smooth.data(ind(1));
        clear ind
    end
    else
        peaks_select_inbound_absdist=nan;
end
    
% close all; plot(posteriorts_wtrack,cdat_abs_distance_smooth.data, 'b')
% hold on 
% plot([forelimb_plant_max_dist_inbound forelimb_plant_max_dist_inbound], [-100 100],'r')
% plot([forelimb_plant_max_dist_outbound forelimb_plant_max_dist_outbound], [-100 100],'k')
% plot(posteriorts_wtrack,cdat_ahead_behind_distance_smooth.data, 'y')
% plot(cam_rt_fit,dlc_forepawR_x); plot(cam_rt_fit,dlc_forepawL_x);
% plot([outbound_center(:,1) outbound_center(:,1)], [-40 400], 'k');
% plot([outbound_center(:,2) outbound_center(:,2)], [-40 400], 'r');

%%  collected stance times for forelimbs and hindlimbs 
    out.index=index; % collection of day epoch 
    out.posteriorts=posteriorts_wtrack;
    
    out.forelimb_inbound_mua=forelimb_plant_inbound_mua;
    out.forelimb_outbound_mua=forelimb_plant_outbound_mua;

    out.forelimb_inbound_reldist=forelimb_plant_inbound_reldist;
    out.forelimb_outbound_reldist=forelimb_plant_outbound_reldist;
    
    out.hindlimb_inbound=hindlimb_plant_inbound;
    out.hindlimb_outbound=hindlimb_plant_outbound;

    out.ahead_behind_smooth=cdat_ahead_behind_distance_smooth;
    out.absolute_dist_smooth= cdat_abs_distance_smooth;
    out.ci_smooth=cdat_ci_smooth;
    
    out.eeggnd=cdat_eeggnd;
    out.theta_filt_5_11=cdat_theta_filt;
    
    out.forelimb_select_outbound_reldist=forelimb_select_outbound_reldist;
    out.forelimb_select_inbound_reldist=forelimb_select_inbound_reldist;
    out.peaks_select_outbound_reldist=peaks_select_outbound_reldist;
    out.peaks_select_inbound_reldist=peaks_select_inbound_reldist;
    out.troughs_select_outbound_reldist=troughs_select_outbound_reldist;
    out.troughs_select_inbound_reldist=troughs_select_inbound_reldist;
    
    out.forelimb_select_outbound_absdist=forelimb_select_outbound_absdist;
    out.forelimb_select_inbound_absdist=forelimb_select_inbound_absdist;
    out.peaks_select_outbound_absdist=peaks_select_outbound_absdist;
    out.peaks_select_inbound_absdist=peaks_select_inbound_absdist;
    
    out.forelimb_inbound_longsweep=forelimb_plant_max_dist_inbound;
    out.forelimb_outbound_longsweep=forelimb_plant_max_dist_outbound;
    out.peaks_outbound_longsweep=peaks_outbound_longsweep;
    out.peaks_inbound_longsweep=peaks_inbound_longsweep;
    
    out.outbound_run=outbound_center;
    out.inbound_run=inbound_center;
    
    else
    out.index=index; % collection of day epoch 
    out.posteriorts=nan;
    
    out.forelimb_inbound_mua=nan;
    out.forelimb_outbound_mua=nan;

    out.forelimb_inbound_reldist=nan;
    out.forelimb_outbound_reldist=nan;
    
    out.hindlimb_inbound=nan;
    out.hindlimb_outbound=nan;

    out.ahead_behind_smooth=nan;
    out.absolute_dist_smooth= nan;
    out.ci_smooth=nan;
    
    out.eeggnd=nan;
    out.theta_filt_5_11=nan;  
    
    out.forelimb_select_outbound_reldist=nan;
    out.forelimb_select_inbound_reldist=nan;
    out.peaks_select_outbound_reldist=nan;
    out.peaks_select_inbound_reldist=nan;
    out.troughs_select_outbound_reldist=nan;
    out.troughs_select_inbound_reldist=nan;
    
    out.forelimb_select_outbound_absdist=nan;
    out.forelimb_select_inbound_absdist=nan;
    out.peaks_select_outbound_absdist=nan;
    out.peaks_select_inbound_absdist=nan;
    
    out.forelimb_inbound_longsweep=nan;
    out.forelimb_outbound_longsweep=nan;
    out.peaks_outbound_longsweep=nan;
    out.peaks_inbound_longsweep=nan;
    
    out.outbound_run=nan;
    out.inbound_run=nan;
    end
end