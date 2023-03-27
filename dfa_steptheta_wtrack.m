function [out] = dfa_steptheta_wtrack(index,excludeperiods,eeg,eeggnd,~,varargin)

%   dfa_steptheta_wtrack 
%   I calculate the phase preference circular mean direction, mean
%   vector length, and Rayleigh statistical significance for each limb steps
%   during inbound and outbound runs in particular portions on the track
%   Forelimb and Hindlimb steps 

    % process varargin if present and overwrite default values
    if (~isempty(varargin))
        assign(varargin{:});
    else
%          out_min=70;
%          out_max=100;
    warning('thresholds empty - cant compute!')
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
    %keyboard
    % Filter LFP
    load('thetafilter.mat');
    theta_filtered_lfp=filtfilt(thetafilter.tf.num,thetafilter.tf.den,double(int16(theta_data)));
    HT = hilbert(theta_filtered_lfp);
    amplitude = sqrt(real(HT).^2 + imag(HT).^2);
    phase_hpc = angle(HT)+pi;

    % Linearised position 
    post_path = ''; 
    postfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_cv_classifier_clusterless_vel_0_nose_alltime5x_results.nc',post_path,animal,animal,d,e);  %
    
    if ((~iscell(varargin{1,1}{1,d})) || (isempty(run_periods)) || (isempty(varargin{1,1}{1,d}{1,e})) || (~exist(postfile)))
    
    out.index=index;
    
    out.forepawR_plant_outbound=nan;
    out.forepawL_plant_outbound=nan;
    out.hindpawR_plant_outbound=nan;
    out.hindpawL_plant_outbound=nan;
    
    out.forepawR_lift_outbound=nan;
    out.forepawL_lift_outbound=nan;
    out.hindpawR_lift_outbound=nan;
    out.hindpawL_lift_outbound=nan;
    
    out.forepawR_trough_outbound=nan;
    out.forepawL_trough_outbound=nan;
    out.hindpawR_trough_outbound=nan;
    out.hindpawL_trough_outbound=nan;
    
    out.forepawR_peak=nan;
    out.forepawL_peak=nan;
    out.hindpawR_peak=nan;
    out.hindpawL_peak=nan;
    
    out.forepawR_phase_plant_outbound=nan;
    out.forepawL_phase_plant_outbound=nan;
    out.hindpawR_phase_plant=nan;
    out.hindpawL_phase_plant=nan;
    
    out.forepawR_amp_plant=nan;
    out.forepawL_amp_plant=nan;
    out.hindpawR_amp_plant=nan;
    out.hindpawL_amp_plant=nan;

    out.forepawR_phase_trough_outbound=nan;
    out.forepawL_phase_trough_outbound=nan;
    out.hindpawR_phase_trough_outbound=nan;
    out.hindpawL_phase_trough_outbound=nan;
    
    out.forepawR_phase_lift_outbound=nan;
    out.forepawL_phase_lift_outbound=nan;
    out.hindpawR_phase_lift_outbound=nan;
    out.hindpawL_phase_lift_outbound=nan;
    
    out.forepawR_amp_lift_outbound=nan;
    out.forepawL_amp_lift_outbound=nan;
    out.hindpawR_amp_lift_outbound=nan;
    out.hindpawL_amp_lift_outbound=nan;

    out.forepawR_phase_peak_outbound=nan;
    out.forepawL_phase_peak_outbound=nan;
    out.hindpawR_phase_peak_outbound=nan;
    out.hindpawL_phase_peak_outbound=nan;
    
    out.forepawR_plant_inbound=nan;
    out.forepawL_plant_inbound=nan;
    out.hindpawR_plant_inbound=nan;
    out.hindpawL_plant_inbound=nan;
    
    out.forepawR_lift_inbound=nan;
    out.forepawL_lift_inbound=nan;
    out.hindpawR_lift_inbound=nan;
    out.hindpawL_lift_inbound=nan;
    
    out.forepawR_trough_inbound=nan;
    out.forepawL_trough_inbound=nan;
    out.hindpawR_trough_inbound=nan;
    out.hindpawL_trough_inbound=nan;
    
    out.forepawR_peak_inbound=nan;
    out.forepawL_peak_inbound=nan;
    out.hindpawR_peak_inbound=nan;
    out.hindpawL_peak_inbound=nan;
    
    out.forepawR_phase_plant_inbound=nan;
    out.forepawL_phase_plant_inbound=nan;
    out.hindpawR_phase_plant_inbound=nan;
    out.hindpawL_phase_plant_inbound=nan;
    
    out.forepawR_amp_plant_inbound=nan;
    out.forepawL_amp_plant_inbound=nan;
    out.hindpawR_amp_plant_inbound=nan;
    out.hindpawL_amp_plant_inbound=nan;

    out.forepawR_phase_trough_inbound=nan;
    out.forepawL_phase_trough_inbound=nan;
    out.hindpawR_phase_trough_inbound=nan;
    out.hindpawL_phase_trough_inbound=nan;
    
    out.forepawR_phase_lift_inbound=nan;
    out.forepawL_phase_lift_inbound=nan;
    out.hindpawR_phase_lift_inbound=nan;
    out.hindpawL_phase_lift_inbound=nan;
    
    out.forepawR_amp_lift_inbound=nan;
    out.forepawL_amp_lift_inbound=nan;
    out.hindpawR_amp_lift_inbound=nan;
    out.hindpawL_amp_lift_inbound=nan;

    out.forepawR_phase_peak_inbound=nan;
    out.forepawL_phase_peak_inbound=nan;
    out.hindpawR_phase_peak_inbound=nan;
    out.hindpawL_phase_peak_inbound=nan;

    out.forepawR_filt_6_8=nan;
    out.forepawL_filt_6_8=nan;
    out.hindpawR_filt_6_8=nan;
    out.hindpawL_filt_6_8=nan;
    
    out.theta_data=nan;
    out.theta_time=nan;
    out.theta_filt_5_11=nan;
    out.theta_amp=nan;
    out.theta_phase=nan;
    out.run_periods=nan;
    out.outbound_periods=nan;
    out.inbound_periods=nan;
    out.varargin=nan;
        
    else
    
    posteriorts_wtrack = ncread(postfile,'time'); % time bins 
    postposbins_wtrack=1+ncread(postfile,'position'); % position on w

    linposfile = sprintf('%s%s/filterframework/decoding_clusterless/%s_%d_%d_linearised_position_nose.nc', post_path, animal, animal, d, e); %
    linpos_nose_wtrack = ncread(linposfile,'linear_position'); % load linpos
    
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
  
    
    %% create inbound and outbound times 
    linpos_min_max=[out_min out_max]; % collect stance times only within a window 
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

    % ForepawL
    [cam_rt_fit_forepawL_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);
    [forepawL_phase_results_outbound,results_forepawL_outbound]=thetaphase_assign(cam_rt_fit_forepawL_results_outbound,theta_time, phase_hpc, amplitude);

    % HindpawR
    [cam_rt_fit_hindpawR_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);
    [hindpawR_phase_results_outbound,results_hindpawR_outbound]=thetaphase_assign(cam_rt_fit_hindpawR_results_outbound,theta_time, phase_hpc, amplitude);
    
    % ForepawR 
    [cam_rt_fit_forepawR_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);
    [forepawR_phase_results_outbound,results_forepawR_outbound]=thetaphase_assign(cam_rt_fit_forepawR_results_outbound,theta_time, phase_hpc, amplitude);
    
    % HindpawL
    [cam_rt_fit_hindpawL_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);
    [hindpawL_phase_results_outbound,results_hindpawL_outbound]=thetaphase_assign(cam_rt_fit_hindpawL_results_outbound,theta_time, phase_hpc, amplitude);

    %%   
    [cam_rt_fit_forepawL_results_inbound, cdat_forepawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);
    [forepawL_phase_results_inbound,results_forepawL_inbound]=thetaphase_assign(cam_rt_fit_forepawL_results_inbound,theta_time, phase_hpc, amplitude);

    % HindpawR
    [cam_rt_fit_hindpawR_results_inbound, cdat_hindpawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);
    [hindpawR_phase_results_inbound,results_hindpawR_inbound]=thetaphase_assign(cam_rt_fit_hindpawR_results_inbound,theta_time, phase_hpc, amplitude);
    
    % ForepawR 
    [cam_rt_fit_forepawR_results_inbound, cdat_forepawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);
    [forepawR_phase_results_inbound,results_forepawR_inbound]=thetaphase_assign(cam_rt_fit_forepawR_results_inbound,theta_time, phase_hpc, amplitude);
 
    % HindpawL
    [cam_rt_fit_hindpawL_results_inbound, cdat_hindpawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);
    [hindpawL_phase_results_inbound,results_hindpawL_inbound]=thetaphase_assign(cam_rt_fit_hindpawL_results_inbound,theta_time, phase_hpc, amplitude);

    %% save computed output
    out.index=index;
    
    out.forepawR_plant_outbound=results_forepawR_outbound(1,:);
    out.forepawL_plant_outbound=results_forepawL_outbound(1,:);
    out.hindpawR_plant_outbound=results_hindpawR_outbound(1,:);
    out.hindpawL_plant_outbound=results_hindpawL_outbound(1,:);
    
    out.forepawR_lift_outbound=results_forepawR_outbound(2,:);
    out.forepawL_lift_outbound=results_forepawL_outbound(2,:);
    out.hindpawR_lift_outbound=results_hindpawR_outbound(2,:);
    out.hindpawL_lift_outbound=results_hindpawL_outbound(2,:);
    
    out.forepawR_trough_outbound=results_forepawR_outbound(3,:);
    out.forepawL_trough_outbound=results_forepawL_outbound(3,:);
    out.hindpawR_trough_outbound=results_hindpawR_outbound(3,:);
    out.hindpawL_trough_outbound=results_hindpawL_outbound(3,:);
    
    out.forepawR_peak=results_forepawR_outbound(4,:);
    out.forepawL_peak=results_forepawL_outbound(4,:);
    out.hindpawR_peak=results_hindpawR_outbound(4,:);
    out.hindpawL_peak=results_hindpawL_outbound(4,:);
    
    out.forepawR_phase_plant_outbound=forepawR_phase_results_outbound.plant;
    out.forepawL_phase_plant_outbound=forepawL_phase_results_outbound.plant;
    out.hindpawR_phase_plant=hindpawR_phase_results_outbound.plant;
    out.hindpawL_phase_plant=hindpawL_phase_results_outbound.plant;
    
    out.forepawR_amp_plant=forepawR_phase_results_outbound.amp_plant;
    out.forepawL_amp_plant=forepawL_phase_results_outbound.amp_plant;
    out.hindpawR_amp_plant=hindpawR_phase_results_outbound.amp_plant;
    out.hindpawL_amp_plant=hindpawL_phase_results_outbound.amp_plant;

    out.forepawR_phase_trough_outbound=forepawR_phase_results_outbound.trough;
    out.forepawL_phase_trough_outbound=forepawL_phase_results_outbound.trough;
    out.hindpawR_phase_trough_outbound=hindpawR_phase_results_outbound.trough;
    out.hindpawL_phase_trough_outbound=hindpawL_phase_results_outbound.trough;
    
    out.forepawR_phase_lift_outbound=forepawR_phase_results_outbound.lift;
    out.forepawL_phase_lift_outbound=forepawL_phase_results_outbound.lift;
    out.hindpawR_phase_lift_outbound=hindpawR_phase_results_outbound.lift;
    out.hindpawL_phase_lift_outbound=hindpawL_phase_results_outbound.lift;
    
    out.forepawR_amp_lift_outbound=forepawR_phase_results_outbound.amp_lift;
    out.forepawL_amp_lift_outbound=forepawL_phase_results_outbound.amp_lift;
    out.hindpawR_amp_lift_outbound=hindpawR_phase_results_outbound.amp_lift;
    out.hindpawL_amp_lift_outbound=hindpawL_phase_results_outbound.amp_lift;

    out.forepawR_phase_peak_outbound=forepawR_phase_results_outbound.peak;
    out.forepawL_phase_peak_outbound=forepawL_phase_results_outbound.peak;
    out.hindpawR_phase_peak_outbound=hindpawR_phase_results_outbound.peak;
    out.hindpawL_phase_peak_outbound=hindpawL_phase_results_outbound.peak;
    
    out.forepawR_plant_inbound=results_forepawR_inbound(1,:);
    out.forepawL_plant_inbound=results_forepawL_inbound(1,:);
    out.hindpawR_plant_inbound=results_hindpawR_inbound(1,:);
    out.hindpawL_plant_inbound=results_hindpawL_inbound(1,:);
    
    out.forepawR_lift_inbound=results_forepawR_inbound(2,:);
    out.forepawL_lift_inbound=results_forepawL_inbound(2,:);
    out.hindpawR_lift_inbound=results_hindpawR_inbound(2,:);
    out.hindpawL_lift_inbound=results_hindpawL_inbound(2,:);
    
    out.forepawR_trough_inbound=results_forepawR_inbound(3,:);
    out.forepawL_trough_inbound=results_forepawL_inbound(3,:);
    out.hindpawR_trough_inbound=results_hindpawR_inbound(3,:);
    out.hindpawL_trough_inbound=results_hindpawL_inbound(3,:);
    
    out.forepawR_peak_inbound=results_forepawR_inbound(4,:);
    out.forepawL_peak_inbound=results_forepawL_inbound(4,:);
    out.hindpawR_peak_inbound=results_hindpawR_inbound(4,:);
    out.hindpawL_peak_inbound=results_hindpawL_inbound(4,:);
    
    out.forepawR_phase_plant_inbound=forepawR_phase_results_inbound.plant;
    out.forepawL_phase_plant_inbound=forepawL_phase_results_inbound.plant;
    out.hindpawR_phase_plant_inbound=hindpawR_phase_results_inbound.plant;
    out.hindpawL_phase_plant_inbound=hindpawL_phase_results_inbound.plant;
    
    out.forepawR_amp_plant_inbound=forepawR_phase_results_inbound.amp_plant;
    out.forepawL_amp_plant_inbound=forepawL_phase_results_inbound.amp_plant;
    out.hindpawR_amp_plant_inbound=hindpawR_phase_results_inbound.amp_plant;
    out.hindpawL_amp_plant_inbound=hindpawL_phase_results_inbound.amp_plant;

    out.forepawR_phase_trough_inbound=forepawR_phase_results_inbound.trough;
    out.forepawL_phase_trough_inbound=forepawL_phase_results_inbound.trough;
    out.hindpawR_phase_trough_inbound=hindpawR_phase_results_inbound.trough;
    out.hindpawL_phase_trough_inbound=hindpawL_phase_results_inbound.trough;
    
    out.forepawR_phase_lift_inbound=forepawR_phase_results_inbound.lift;
    out.forepawL_phase_lift_inbound=forepawL_phase_results_inbound.lift;
    out.hindpawR_phase_lift_inbound=hindpawR_phase_results_inbound.lift;
    out.hindpawL_phase_lift_inbound=hindpawL_phase_results_inbound.lift;
    
    out.forepawR_amp_lift_inbound=forepawR_phase_results_inbound.amp_lift;
    out.forepawL_amp_lift_inbound=forepawL_phase_results_inbound.amp_lift;
    out.hindpawR_amp_lift_inbound=hindpawR_phase_results_inbound.amp_lift;
    out.hindpawL_amp_lift_inbound=hindpawL_phase_results_inbound.amp_lift;

    out.forepawR_phase_peak_inbound=forepawR_phase_results_inbound.peak;
    out.forepawL_phase_peak_inbound=forepawL_phase_results_inbound.peak;
    out.hindpawR_phase_peak_inbound=hindpawR_phase_results_inbound.peak;
    out.hindpawL_phase_peak_inbound=hindpawL_phase_results_inbound.peak;

    out.forepawR_filt_6_8=cdat_forepawR_filt_6_8;
    out.forepawL_filt_6_8=cdat_forepawL_filt_6_8;
    out.hindpawR_filt_6_8=cdat_hindpawR_filt_6_8;
    out.hindpawL_filt_6_8=cdat_hindpawL_filt_6_8;
    
    out.theta_data=theta_data;
    out.theta_time=theta_time;
    out.theta_filt_5_11=theta_filtered_lfp;
    out.theta_amp=amplitude;
    out.theta_phase=phase_hpc;
    out.run_periods=run_periods;
    out.outbound_periods=outbound_center;
    out.inbound_periods=inbound_center;
    out.varargin=varargin;
   
    end
end