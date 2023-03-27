function [out] = dfa_steptheta_lineartrack_forelimb(index,excludeperiods,eeg,eeggnd,~,varargin)

%   dfa_steptheta_lineartrack_forelimb 
%   I calculate the phase preference circular mean direction, mean
%   vector length, and Rayleigh statistical significance for each limb steps
%   during inbound and outbound runs in particular portions on the track
%   Forelimb and Hindlimb steps 

% Detailed explanation goes here
    

    % process varargin if present and overwrite default values
    if (~isempty(varargin))
        assign(varargin{:});
    else
         out_min=70;
         out_max=100;
    end
    
    %% load the relevant files, information 
    d = index(1); e = index(2); reftet=index(3);
  
    % Include periods that have velocity above the set threshold
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

    % Linearised position 
    post_path = '/opt/stelmo/abhilasha/animals/'; 
    linposfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_linearised_position_nose_5x.nc', post_path, animal, animal, d, e); %
    
    if ((~iscell(varargin{1,1}{1,d})) | (isempty(run_periods)) | (isempty(varargin{1,1}{1,d}{1,e})) | (~exist(linposfile)))
    
    %% save computed output
    out.index=index;
    
    out.forelimb_plant_outbound=nan;
    out.hindlimb_plant_outbound=nan;
    
    out.forelimb_lift_outbound=nan;
    out.hindlimb_lift_outbound=nan;
    
    out.forelimb_trough_outbound=nan;
    out.hindlimb_trough_outbound=nan;
    
    out.forelimb_peak_outbound=nan;
    out.hindlimb_peak_outbound=nan;
    
    out.forelimb_phase_plant_outbound=nan;
    out.hindlimb_phase_plant_outbound=nan;
    
    out.forelimb_amp_plant_outbound=nan;
    out.hindlimb_amp_plant_outbound=nan;

    out.forelimb_phase_trough_outbound=nan;
    out.hindlimb_phase_trough_outbound=nan;
    
    out.forelimb_phase_lift_outbound=nan;
    out.hindlimb_phase_lift_outbound=nan;
    
    out.forelimb_amp_lift_outbound=nan;
    out.hindlimb_amp_lift_outbound=nan;

    out.forelimb_phase_peak_outbound=nan;
    out.hindlimb_phase_peak_outbound=nan;
    
    out.forelimb_plant_inbound=nan;
    out.hindlimb_plant_inbound=nan;
    
    out.forelimb_lift_inbound=nan;
    out.hindlimb_lift_inbound=nan;
    
    out.forelimb_trough_inbound=nan;
    out.hindlimb_trough_inbound=nan;
    
    out.forelimb_peak_inbound=nan;
    out.hindlimb_peak_inbound=nan;
    
    out.forelimb_phase_plant_inbound=nan;
    out.hindlimb_phase_plant_inbound=nan;
    
    out.forelimb_amp_plant_inbound=nan;
    out.hindlimb_amp_plant_inbound=nan;

    out.forelimb_phase_trough_inbound=nan;
    out.hindlimb_phase_trough_inbound=nan;
    
    out.forelimb_phase_lift_inbound=nan;
    out.hindlimb_phase_lift_inbound=nan;
    
    out.forelimb_amp_lift_inbound=nan;
    out.hindlimb_amp_lift_inbound=nan;

    out.forelimb_phase_peak_inbound=nan;
    out.hindlimb_phase_peak_inbound=nan;

    out.varargin=nan;
        
    else
    
    linpos_nose_wtrack = ncread(linposfile,'linear_position'); % load linpos
    posteriorts_wtrack = ncread(linposfile,'time');

    % load dlc results from rawpos. These are after Dlc2FFPOS
    dlc_results=varargin{1,1}{1,d}{1,e}; 

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
    dlc_forepawL_x=dlc_results.data(:,8);
    dlc_forepawL_y=dlc_results.data(:,9);

    %forepawR
    dlc_forepawR_x= dlc_results.data(:,11);
    dlc_forepawR_y= dlc_results.data(:,12);

    %hindpawL
    dlc_hindpawL_x=dlc_results.data(:,14);
    dlc_hindpawL_y=dlc_results.data(:,15);

    %hindpawR
    dlc_hindpawR_x=dlc_results.data(:,17);
    dlc_hindpawR_y=dlc_results.data(:,18);

    dlc_n_records = size(dlc_hindpawL_x,1);
    fprintf('# of DeepLabCut timestamps: %d \n', dlc_n_records)

    % estimated framerate based on camera time 
    est_framerate=median(1./diff(cam_rt_fit));
    
    %% create inbound and outbound times 
    linpos_min_max=[out_min out_max] % collect stance times only within a window 
    outboundvec = (linpos_nose_wtrack > linpos_min_max(:,1)) & (linpos_nose_wtrack < linpos_min_max (:,2));
    outbound_periods_pos = vec2list(outboundvec,posteriorts_wtrack);

    outbound_center=[];
    inbound_center=[];
    diff_linpos=[];

    for i=1:length(outbound_periods_pos)

    ind_linpos_start_end=find(posteriorts_wtrack>outbound_periods_pos(i,1)&posteriorts_wtrack<outbound_periods_pos(i,2));

    if ~isempty(ind_linpos_start_end)

    diff_linpos=linpos_nose_wtrack(ind_linpos_start_end(1)) - linpos_nose_wtrack(ind_linpos_start_end(end));

    if  ((abs(diff_linpos))<15)
        outbound_center(i,1)= nan;
        outbound_center(i,2)=nan;
        inbound_center(i,1)=nan;
        inbound_center(i,2)=nan;  
    else

    if  (diff_linpos < 0)
        outbound_center(i,1)= outbound_periods_pos(i,1);
        outbound_center(i,2)=outbound_periods_pos(i,2);
        inbound_center(i,1)=nan;
        inbound_center(i,2)=nan; 
    else
        outbound_center(i,1)= nan;
        outbound_center(i,2)=nan;
        inbound_center(i,1)=outbound_periods_pos(i,1);
        inbound_center(i,2)=outbound_periods_pos(i,2);      
    end

    end
    else 
        outbound_center(i,1)= nan;
        outbound_center(i,2)=nan;
        inbound_center(i,1)=nan;
        inbound_center(i,2)=nan;  
        continue 
    end
    clear diff_linpos
    end

    outbound_center(any(isnan(outbound_center), 2), :) = [];
    inbound_center(any(isnan(inbound_center), 2), :) = [];
    
    %% compute forepaw phase results during outbound portion of the wtrack 

    % Forelimb outbound
    [cam_rt_fit_forepawL_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);
    [cam_rt_fit_forepawR_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);

     cam_rt_fit_forelimb_outbound=structcat(1,cam_rt_fit_forepawL_results_outbound,cam_rt_fit_forepawR_results_outbound);
    [cam_rt_fit_forelimb_outbound]=sort_struct(cam_rt_fit_forelimb_outbound);
    % Forelimb combined phase results  
    [forelimb_phase_results_outbound,results_forelimb_outbound]=thetaphase_assign(cam_rt_fit_forelimb_outbound,theta_time, phase_hpc, amplitude)

    % Hindlimb outbound
    [cam_rt_fit_hindpawR_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);
    [cam_rt_fit_hindpawL_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);
        
    cam_rt_fit_hindlimb_outbound=structcat(1,cam_rt_fit_hindpawL_results_outbound,cam_rt_fit_hindpawR_results_outbound);
    [cam_rt_fit_hindlimb_outbound]=sort_struct(cam_rt_fit_hindlimb_outbound);
    % hindlimb combined phase results  
    [hindlimb_phase_results_outbound,results_hindlimb_outbound]=thetaphase_assign(cam_rt_fit_hindlimb_outbound,theta_time, phase_hpc, amplitude);
    
    %%   
    % Forelimb outbound
    [cam_rt_fit_forepawL_results_inbound, cdat_forepawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);
    [cam_rt_fit_forepawR_results_inbound, cdat_forepawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);

     cam_rt_fit_forelimb_inbound=structcat(1,cam_rt_fit_forepawL_results_inbound,cam_rt_fit_forepawR_results_inbound);
    [cam_rt_fit_forelimb_inbound]=sort_struct(cam_rt_fit_forelimb_inbound);
    % Forelimb combined phase results  
    [forelimb_phase_results_inbound,results_forelimb_inbound]=thetaphase_assign(cam_rt_fit_forelimb_inbound,theta_time, phase_hpc, amplitude)

    % Hindlimb outbound
    [cam_rt_fit_hindpawR_results_inbound, cdat_hindpawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);
    [cam_rt_fit_hindpawL_results_inbound, cdat_hindpawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_v11_cm(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center);
        
    cam_rt_fit_hindlimb_inbound=structcat(1,cam_rt_fit_hindpawL_results_inbound,cam_rt_fit_hindpawR_results_inbound);
    [cam_rt_fit_hindlimb_inbound]=sort_struct(cam_rt_fit_hindlimb_inbound);
    % hindlimb combined phase results  
    [hindlimb_phase_results_inbound,results_hindlimb_inbound]=thetaphase_assign(cam_rt_fit_hindlimb_inbound,theta_time, phase_hpc, amplitude)
    
    %% save computed output
    out.index=index;
    
    out.forelimb_plant_outbound=results_forelimb_outbound(1,:);
    out.hindlimb_plant_outbound=results_hindlimb_outbound(1,:);
    
    out.forelimb_lift_outbound=results_forelimb_outbound(2,:);
    out.hindlimb_lift_outbound=results_hindlimb_outbound(2,:);
    
    out.forelimb_trough_outbound=results_forelimb_outbound(3,:);
    out.hindlimb_trough_outbound=results_hindlimb_outbound(3,:);
    
    out.forelimb_peak_outbound=results_forelimb_outbound(4,:);
    out.hindlimb_peak_outbound=results_hindlimb_outbound(4,:);
    
    out.forelimb_phase_plant_outbound=forelimb_phase_results_outbound.plant;
    out.hindlimb_phase_plant_outbound=hindlimb_phase_results_outbound.plant;
    
    out.forelimb_amp_plant_outbound=forelimb_phase_results_outbound.amp_plant;
    out.hindlimb_amp_plant_outbound=hindlimb_phase_results_outbound.amp_plant;

    out.forelimb_phase_trough_outbound=forelimb_phase_results_outbound.trough;
    out.hindlimb_phase_trough_outbound=hindlimb_phase_results_outbound.trough;
    
    out.forelimb_phase_lift_outbound=forelimb_phase_results_outbound.lift;
    out.hindlimb_phase_lift_outbound=hindlimb_phase_results_outbound.lift;
    
    out.forelimb_amp_lift_outbound=forelimb_phase_results_outbound.amp_lift;
    out.hindlimb_amp_lift_outbound=hindlimb_phase_results_outbound.amp_lift;

    out.forelimb_phase_peak_outbound=forelimb_phase_results_outbound.peak;
    out.hindlimb_phase_peak_outbound=hindlimb_phase_results_outbound.peak;
    
    out.forelimb_plant_inbound=results_forelimb_inbound(1,:);
    out.hindlimb_plant_inbound=results_hindlimb_inbound(1,:);
    
    out.forelimb_lift_inbound=results_forelimb_inbound(2,:);
    out.hindlimb_lift_inbound=results_hindlimb_inbound(2,:);
    
    out.forelimb_trough_inbound=results_forelimb_inbound(3,:);
    out.hindlimb_trough_inbound=results_hindlimb_inbound(3,:);
    
    out.forelimb_peak_inbound=results_forelimb_inbound(4,:);
    out.hindlimb_peak_inbound=results_hindlimb_inbound(4,:);
    
    out.forelimb_phase_plant_inbound=forelimb_phase_results_inbound.plant;
    out.hindlimb_phase_plant_inbound=hindlimb_phase_results_inbound.plant;
    
    out.forelimb_amp_plant_inbound=forelimb_phase_results_inbound.amp_plant;
    out.hindlimb_amp_plant_inbound=hindlimb_phase_results_inbound.amp_plant;

    out.forelimb_phase_trough_inbound=forelimb_phase_results_inbound.trough;
    out.hindlimb_phase_trough_inbound=hindlimb_phase_results_inbound.trough;
    
    out.forelimb_phase_lift_inbound=forelimb_phase_results_inbound.lift;
    out.hindlimb_phase_lift_inbound=hindlimb_phase_results_inbound.lift;
    
    out.forelimb_amp_lift_inbound=forelimb_phase_results_inbound.amp_lift;
    out.hindlimb_amp_lift_inbound=hindlimb_phase_results_inbound.amp_lift;

    out.forelimb_phase_peak_inbound=forelimb_phase_results_inbound.peak;
    out.hindlimb_phase_peak_inbound=hindlimb_phase_results_inbound.peak;

    out.varargin=varargin;
   
    end
end