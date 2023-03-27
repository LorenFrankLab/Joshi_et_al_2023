function out = dfa_speed_steps_theta_forelimb(index, excludeperiods, eeg, eeggnd, ~, posdlc, varargin)

%  dfa_speed_steps_theta_forelimb
%   For any epochs that have movement decoding done, compute:
%	- ahead/behind distance of the peak posterior density to the current body position 
% 	- computes if the ahead-behind distance is closer to the current position of the animal than other times 
% 	- loads posterior, posdlc, linpos etc. WTRACK

% velocity thershold actually does not work as both outbound and inbound
% periods are selected bsaed on their pos rather than intersection with run
% periods

    % define defaults 
    appendindex = 0;

    % process varargin if present and overwrite default values
    if (~isempty(varargin))
        assign(varargin{:});
    else
%          out_min=70;
%          out_max=100;
%          ci_thresh=40;
%          ci_mean=20;
%          reldist_mean=10;
%          time_step=0.60;
    end
    
    % load the day and epoch number 
    d = index(1);
    e = index(2);
    reftet=index(3);
        
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
    
    Fs=eeggnd{1,d}{1,e}{1,reftet}.samprate;
  
    % Filter LFP
    load('thetafilter.mat');
    theta_filtered_lfp=filtfilt(thetafilter.tf.num,thetafilter.tf.den,double(int16(theta_data)));
    HT = hilbert(theta_filtered_lfp);
    amplitude = sqrt(real(HT).^2 + imag(HT).^2);
    phase_hpc = angle(HT)+pi;
    instfreq_theta = Fs/(2*pi)*diff(unwrap(angle(HT)+pi));
    
%% load posteriors
    post_path = ''; 
    postfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_cv_classifier_clusterless_vel_0_nose_alltime5x_results.nc',post_path,animal,animal,d,e);  %
   
    if exist(postfile)
    posteriorts_wtrack = ncread(postfile,'time'); % time bins 
    postposbins_wtrack=1+ncread(postfile,'position'); % position on w

    linposfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_linearised_position_nose.nc', post_path, animal, animal, d, e); %
    linpos_nose_wtrack = ncread(linposfile,'linear_position'); % load linpos

    %% Load dlc positions and times 
    dlc_results=posdlc{1,d}{1,e}; 

    %load ptp adjusted timestamps 
    cam_rt_fit=dlc_results.data(:,1);

    %load nose
    dlc_nose_x=dlc_results.data(:,2);
    dlc_nose_y=dlc_results.data(:,3);

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
    [cam_rt_fit_forepawL_results, cdat_forepawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_tail_vel, est_framerate, run_periods);

    % HindpawR
    [cam_rt_fit_hindpawR_results, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_tail_vel, est_framerate, run_periods);
    
    % ForepawR 
    [cam_rt_fit_forepawR_results, cdat_forepawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_tail_vel, est_framerate, run_periods);
    
    % HindpawL
    [cam_rt_fit_hindpawL_results, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_tail_vel, est_framerate, run_periods);

    % velcoity and inst velocity
    dlc_nose_wtrack_vel=dlc_results.data(:,4); % this value is obtained from 
    
    % acceleration and inst acc 
    nose_acc=diff(dlc_nose_wtrack_vel); % this is the acceleration per camera frame. Multiply this value by the time unit by the est framerate 125
    %to get the acceleration in cm/s2
    
    %forelimbR steps
    filtered_steps_forelimbR=(cdat_forepawR_filt_6_8.data);
    filtered_steps_forelimbR(isnan(filtered_steps_forelimbR))=0;
    cdat_filtered_data_forelimbR=imcont('data', filtered_steps_forelimbR, 'timestamp', cam_rt_fit(1:end-1));
    filtopt_LPF_6_8 = mkfiltopt('name', 'aj_step_LPF', 'filttype', 'bandpass', 'F', [0.5 1  6 7]);
    cdat_step_filt_6_8_forelimbR = contfilt(cdat_filtered_data_forelimbR, 'filtopt', filtopt_LPF_6_8);
    
    HT_steps_forelimbR=hilbert(cdat_step_filt_6_8_forelimbR.data);
    Fs_steps=round(cdat_forepawR_filt_6_8.samplerate);
    instfreq_steps_forelimbR = Fs_steps/(2*pi)*diff(unwrap(angle(HT_steps_forelimbR)));
    
    %forelimbL steps
    filtered_steps_forelimbL=(cdat_forepawL_filt_6_8.data);
    filtered_steps_forelimbL(isnan(filtered_steps_forelimbL))=0;
    cdat_filtered_data_forelimbL=imcont('data', filtered_steps_forelimbL, 'timestamp', cam_rt_fit(1:end-1));
    filtopt_LPF_6_8 = mkfiltopt('name', 'aj_step_LPF', 'filttype', 'bandpass', 'F', [0.5 1  6 7]);
    cdat_step_filt_6_8_forelimbL = contfilt(cdat_filtered_data_forelimbL, 'filtopt', filtopt_LPF_6_8);
    
    HT_steps_forelimbL=hilbert(cdat_step_filt_6_8_forelimbL.data);
    Fs_steps=round(cdat_forepawL_filt_6_8.samplerate);
    instfreq_steps_forelimbL = Fs_steps/(2*pi)*diff(unwrap(angle(HT_steps_forelimbL)));
    
    instfreq_steps=[instfreq_steps_forelimbL+instfreq_steps_forelimbR];
    
    %% collect the inst freq of theta and velocity during the same periods of time 
    
%     ts=cam_rt_fit(1); te=cam_rt_fit(end);
%     n=round((te-ts)/time_step);
%     
%     clear time_windows
%     for i=1:n
%         te_win=ts+time_step*i;
%         ts_win=te_win-time_step;
%         time_windows(i,1)=ts_win;
%         time_windows(i,2)=te_win;
%         %keyboard
%         clear ts_win te_win
%     end
%         
%     for i=2:size(time_windows,1)-1
%         ind_theta=find(theta_time>(time_windows(i,1))&(theta_time<(time_windows(i,2))));
%         ind_vel=find(cam_rt_fit>(time_windows(i,1))&(cam_rt_fit<(time_windows(i,2))));
%         
%         inst_theta(i)=mean(instfreq_theta(ind_theta));
%         inst_vel(i)=mean(dlc_nose_wtrack_vel(ind_vel));
%         inst_acc(i)=mean(nose_acc(ind_vel));
%         inst_steps(i)=mean(instfreq_steps(ind_vel));
%         clear ind_vel ind_theta
%     end
        
%% collect the forelimb plant times in a particular window 
    forelimb_plant=[cam_rt_fit_forepawR_results.midstance; cam_rt_fit_forepawL_results.midstance];
    hindlimb_plant=[cam_rt_fit_hindpawR_results.midstance; cam_rt_fit_hindpawL_results.midstance];

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

    %% Collect forelimb and hindlimb plant times for outbound center and inbound center 
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
    
    %% OUTBOUND Calculations 
    [~,theta_trough_locs]=findpeaks(double(theta_filtered_lfp),theta_time,'MinPeakDistance', 0.08, 'MinPeakHeight', 50);

    theta_trough_locs_run_outbound={}; theta_trough_locs_run_inbound={};  
    for k=1:size(outbound_center, 1)
        locs_collect_idx=find(theta_trough_locs>outbound_center(k,1)&theta_trough_locs<outbound_center(k,2));
        if (~isempty (locs_collect_idx)) && (size(locs_collect_idx,2)>4) && (size(locs_collect_idx,2)<50)
            theta_trough_locs_run_outbound(k)={theta_trough_locs(locs_collect_idx)};
        else 
            theta_trough_locs_run_outbound(k)={nan};
        end
        clear locs_collect_idx;
    end 
    
    for k=1:size(inbound_center, 1)
        locs_collect_idx=find(theta_trough_locs>inbound_center(k,1)&theta_trough_locs<inbound_center(k,2));
        if (~isempty (locs_collect_idx)) && (size(locs_collect_idx,2)>4) && (size(locs_collect_idx,2)<50)
            theta_trough_locs_run_inbound(k)={theta_trough_locs(locs_collect_idx)};
        else 
            theta_trough_locs_run_inbound(k)={nan};
        end
        clear locs_collect_idx;
    end 
        
    theta_trough_locs_run_outbound=cell2mat(theta_trough_locs_run_outbound)';
    theta_trough_locs_run_outbound(any(isnan(theta_trough_locs_run_outbound), 2), :) = [];theta_trough_locs_run_outbound(:,2)=theta_trough_locs_run_outbound; clear idx
    
    theta_trough_locs_run_inbound=cell2mat(theta_trough_locs_run_inbound)';
    theta_trough_locs_run_inbound(any(isnan(theta_trough_locs_run_inbound), 2), :) = [];theta_trough_locs_run_inbound(:,2)=theta_trough_locs_run_inbound; clear idx
    
    %% Calculate speed 
    run_duration_outbound=[]; distance_travelled_outbound=[]; speed_outbound=[]; forelimb_steps_outbound_count=[]; hindlimb_steps_outbound_count=[]; theta_troughs_outbound_count =[]; 
    run_duration_inbound=[]; distance_travelled_inbound=[]; speed_inbound=[]; forelimb_steps_inbound_count=[]; hindlimb_steps_inbound_count=[]; theta_troughs_inbound_count =[];
    inst_vel_inbound=[]; inst_vel_outbound=[]; 
    
    for n=1:size(outbound_center)
        run_duration_outbound(n)=outbound_center(n,2)-outbound_center(n,1);
        idx=find(posteriorts_wtrack>outbound_center(n,1) & posteriorts_wtrack<outbound_center(n,2));

        if size(idx,1)>100
             distance_travelled_outbound(n)=linpos_nose_wtrack(idx(end))-linpos_nose_wtrack(idx(1));
             speed_outbound(n)=abs(distance_travelled_outbound(n))/run_duration_outbound(n);

             % now calculate the number of forelimb steps 
             ind_forelimb_plant=find(forelimb_plant(:,1)>(outbound_center(n,1))&(forelimb_plant(:,1)<(outbound_center(n,2))));
             ind_hindlimb_plant=find(hindlimb_plant(:,1)>(outbound_center(n,1))&(hindlimb_plant(:,1)<(outbound_center(n,2))));
             ind_theta_troughs=find(theta_trough_locs_run_outbound(:,1)>(outbound_center(n,1))&(theta_trough_locs_run_outbound(:,1)<(outbound_center(n,2))));
             
             if ~isempty(ind_theta_troughs)
                 theta_troughs_outbound_count(n)=size(ind_theta_troughs,1);
                 forelimb_steps_outbound_count(n)=size(ind_forelimb_plant,1);
                 hindlimb_steps_outbound_count(n)=size(ind_hindlimb_plant,1);
             else 
                distance_travelled_outbound(n)=nan;
                speed_outbound(n)=NaN;
                forelimb_steps_outbound_count(n)=nan;
                forelimb_steps_outbound_count(n)=nan;
                theta_troughs_outbound_count(n)=nan;
             end
        else
                distance_travelled_outbound(n)=nan;
                speed_outbound(n)=NaN;
                forelimb_steps_outbound_count(n)=nan;
                forelimb_steps_outbound_count(n)=nan;
                theta_troughs_outbound_count(n)=nan;
        end

        clear idx ind_*
    end
    
    for n=1:size(inbound_center)
        run_duration_inbound(n)=inbound_center(n,2)-inbound_center(n,1);
        idx=find(posteriorts_wtrack>inbound_center(n,1) & posteriorts_wtrack<inbound_center(n,2));

        if size(idx,1)>100
             distance_travelled_inbound(n)=linpos_nose_wtrack(idx(end))-linpos_nose_wtrack(idx(1));
             speed_inbound(n)=abs(distance_travelled_inbound(n))/run_duration_inbound(n);

             % now calculate the number of forelimb steps 
             ind_forelimb_plant=find(forelimb_plant(:,1)>(inbound_center(n,1))&(forelimb_plant(:,1)<(inbound_center(n,2))));
             ind_hindlimb_plant=find(hindlimb_plant(:,1)>(inbound_center(n,1))&(hindlimb_plant(:,1)<(inbound_center(n,2))));
             ind_theta_troughs=find(theta_trough_locs_run_inbound(:,1)>(inbound_center(n,1))&(theta_trough_locs_run_inbound(:,1)<(inbound_center(n,2))));
             
             if ~isempty(ind_theta_troughs)
                 theta_troughs_inbound_count(n)=size(ind_theta_troughs,1);
                 forelimb_steps_inbound_count(n)=size(ind_forelimb_plant,1);
                 hindlimb_steps_inbound_count(n)=size(ind_hindlimb_plant,1);
             else 
                distance_travelled_inbound(n)=nan;
                speed_inbound(n)=NaN;
                forelimb_steps_inbound_count(n)=nan;
                forelimb_steps_inbound_count(n)=nan;
                theta_troughs_inbound_count(n)=nan;
             end
        else
                distance_travelled_inbound(n)=nan;
                speed_inbound(n)=NaN;
                forelimb_steps_inbound_count(n)=nan;
                forelimb_steps_inbound_count(n)=nan;
                theta_troughs_inbound_count(n)=nan;
        end
        clear idx ind_*
     end
       
    %% Instantaneous 
    
    % windows in outbound runs 
    for n=1:size(outbound_center)
    ts=outbound_center(n,1); te=outbound_center(n,end);
    num=round((te-ts)/time_step);
    
    if ~(num<1)
    clear time_windows
    for i=1:num
        te_win=ts+time_step*i;
        ts_win=te_win-time_step;
        time_windows(i,1)=ts_win;
        time_windows(i,2)=te_win;
        %keyboard
        clear ts_win te_win
    end
        
    for i=2:size(time_windows,1)-1
        ind_theta=find(theta_time>(time_windows(i,1))&(theta_time<(time_windows(i,2))));
        ind_vel=find(cam_rt_fit>(time_windows(i,1))&(cam_rt_fit<(time_windows(i,2))));
        
        inst_theta_outbound(n,i)=mean(instfreq_theta(ind_theta));
        inst_vel_outbound(n,i)=mean(dlc_nose_wtrack_vel(ind_vel));
        inst_acc_outbound(n,i)=mean(nose_acc(ind_vel));
        inst_steps_outbound(n,i)=mean(instfreq_steps(ind_vel));
        clear ind_vel ind_theta
    end
    
    else
        inst_theta_outbound(n,i)=nan;
        inst_vel_outbound(n,i)=nan;
        inst_acc_outbound(n,i)=nan;
        inst_steps_outbound(n,i)=nan;
    end
    end
    
        % windows in outbound runs 
    for n=1:size(inbound_center)
    ts=inbound_center(n,1); te=inbound_center(n,end);
    num=round((te-ts)/time_step); 
       
    if ~(num<1)
    clear time_windows
    for i=1:num
        te_win=ts+time_step*i;
        ts_win=te_win-time_step;
        time_windows(i,1)=ts_win;
        time_windows(i,2)=te_win;
        %keyboard
        clear ts_win te_win
    end
        
    for i=2:size(time_windows,1)-1
        ind_theta=find(theta_time>(time_windows(i,1))&(theta_time<(time_windows(i,2))));
        ind_vel=find(cam_rt_fit>(time_windows(i,1))&(cam_rt_fit<(time_windows(i,2))));
        
        inst_theta_inbound(n,i)=mean(instfreq_theta(ind_theta));
        inst_vel_inbound(n,i)=mean(dlc_nose_wtrack_vel(ind_vel));
        inst_acc_inbound(n,i)=mean(nose_acc(ind_vel));
        inst_steps_inbound(n,i)=mean(instfreq_steps(ind_vel));
        clear ind_vel ind_theta
    end
    
    else 
        inst_theta_inbound(n,i)=nan;
        inst_vel_inbound(n,i)=nan;
        inst_acc_inbound(n,i)=nan;
        inst_steps_inbound(n,i)=nan;
        
    end
    end

    %%  collected stance times for forelimbs and hindlimbs 
    out.index=index; % collection of day epoch 
    
    %outbound
    out.speed_outbound=speed_outbound;
    out.distance_outbound=distance_travelled_outbound;
    out.forelimb_outbound_steps=forelimb_steps_outbound_count;
    out.hindlimb_outbound_steps=hindlimb_steps_outbound_count;
    out.run_duration_outbound=run_duration_outbound;
    out.theta_count_outbound=theta_troughs_outbound_count;
    
    out.inst_vel_outbound=inst_vel_outbound(:);
    out.inst_theta_outbound=inst_theta_outbound(:);
    out.inst_acc_outbound=inst_acc_outbound(:);
    out.inst_steps_outbound=inst_steps_outbound(:);
    
    %inbound
    out.speed_inbound=speed_inbound;
    out.distance_inbound=distance_travelled_inbound;
    out.forelimb_inbound_steps=forelimb_steps_inbound_count;
    out.hindlimb_inbound_steps=hindlimb_steps_inbound_count;
    out.run_duration_inbound=run_duration_inbound;
    out.theta_count_inbound=theta_troughs_inbound_count;
    
    out.inst_vel_inbound=inst_vel_inbound(:);
    out.inst_theta_inbound=inst_theta_inbound(:);
    out.inst_acc_inbound=inst_acc_inbound(:);
    out.inst_steps_inbound=inst_steps_inbound(:);
    
    else
    out.index=index; % collection of day epoch 
    
    %outbound
    out.speed_outbound=nan;
    out.distance_outbound=nan;
    out.forelimb_outbound_steps=nan;
    out.hindlimb_outbound_steps=nan;
    out.run_duration_outbound=nan;
    out.theta_count_outbound=nan;
    
    out.inst_vel_outbound=nan;
    out.inst_theta_outbound=nan;
    out.inst_acc_outbound=nan;
    out.inst_steps_outbound=nan;
    
    %inbound
    out.speed_inbound=nan;
    out.distance_inbound=nan;
    out.forelimb_inbound_steps=nan;
    out.hindlimb_inbound_steps=nan;
    out.run_duration_inbound=nan;
    out.theta_count_inbound=nan;
    
    out.inst_vel_inbound=nan;
    out.inst_theta_inbound=nan;
    out.inst_acc_inbound=nan;
    out.inst_steps_inbound=nan;
end
end
