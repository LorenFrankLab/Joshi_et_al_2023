function [out] = dfa_steptheta_wtrack_gait(index,excludeperiods,~,~,~,varargin)

%   dfa_steptheta_wtrack_gait
%   I save the forelimb filtered and raw traces, forelimb times, nose
%   position and times to do the cross correlation analysis

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
    
    out.forepawR_plant_inbound=nan;
    out.forepawL_plant_inbound=nan;
    out.hindpawR_plant_inbound=nan;
    out.hindpawL_plant_inbound=nan;
    
    out.forepawR_lift_inbound=nan;
    out.forepawL_lift_inbound=nan;
    out.hindpawR_lift_inbound=nan;
    out.hindpawL_lift_inbound=nan;

    out.forepawR_filt_6_8=nan;
    out.forepawL_filt_6_8=nan;
    out.hindpawR_filt_6_8=nan;
    out.hindpawL_filt_6_8=nan;
    
    out.run_periods=nan;
    out.outbound_periods=nan;
    out.inbound_periods=nan;
     out.xcorr_forepawR_forepawL=nan;
    out.xcorr_forepawR_forepawL_plants=nan;
    out.cam_rt_fit=nan;
    out.speed_val=nan;
    out.est_framerate=nan;
        
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
    linpos_min_max=[out_min out_max]; % collect plant times only within a window 
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

    % HindpawR
    [cam_rt_fit_hindpawR_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);
    
    % ForepawR 
    [cam_rt_fit_forepawR_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);
    
    % HindpawL
    [cam_rt_fit_hindpawL_results_outbound, ~,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, outbound_center, index);

    %%  compute forepaw phase results during inbound portion of the wtrack 
    [cam_rt_fit_forepawL_results_inbound, cdat_forepawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawL_x+dlc_forepawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);

    % HindpawR
    [cam_rt_fit_hindpawR_results_inbound, cdat_hindpawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawR_x+dlc_hindpawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);
    
    % ForepawR 
    [cam_rt_fit_forepawR_results_inbound, cdat_forepawR_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_forepawR_x+dlc_forepawR_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);
 
    % HindpawL
    [cam_rt_fit_hindpawL_results_inbound, cdat_hindpawL_filt_6_8,~, ~]=smooth_dlc_stepcycles_wtrack(dlc_tail_x,(dlc_hindpawL_x+dlc_hindpawL_y), cam_rt_fit, dlc_nose_vel, est_framerate, inbound_center, index);

    
    %% calculate cross correlation all run periods
    maxlag_steps=40; 
    nbouts=size(run_periods,1);

    for k=1:nbouts-1
        % duration of each run
        run_period_duration=run_periods(k,2)-run_periods(k,1);

        % distance traveled in each run 
        idx_step=find((cam_rt_fit>=run_periods(k,1)) & (cam_rt_fit<=run_periods(k,2)));
        speed_val(:,k)=abs((dlc_nose_x(idx_step(end))-dlc_nose_x(idx_step(1))))/(run_period_duration); % average speed in each run

        forepawL_select=cdat_forepawL_filt_6_8.data(idx_step);
        forepawR_select=cdat_forepawR_filt_6_8.data(idx_step);

        % forepawR-forepawL
        [xc,~]=xcorr(forepawR_select,forepawL_select,maxlag_steps, 'coeff');
        xcorr_forepawR_forepawL(k,1:size(xc))=xc;
        clear xc lags forepawL_select forepaw_select;
    end
    
    %% calculate the cross correlation all run periods plant times 
    bin=0.01;
    tmax=0.3;
    
    forelimbR_plant=sort([cam_rt_fit_forepawR_results_outbound.midstance;cam_rt_fit_forepawR_results_inbound.midstance]);
    forelimbL_plant=sort([cam_rt_fit_forepawL_results_outbound.midstance;cam_rt_fit_forepawL_results_inbound.midstance]);
    
    [xcorr_forepawR_forepawL_plants] = spikexcorr(forelimbR_plant, forelimbL_plant, bin,tmax);
    
    %% save computed output
    out.index=index;
    
    out.forepawR_plant_outbound=cam_rt_fit_forepawR_results_outbound.midstance;
    out.forepawL_plant_outbound=cam_rt_fit_forepawL_results_outbound.midstance;
    out.hindpawR_plant_outbound=cam_rt_fit_hindpawR_results_outbound.midstance;
    out.hindpawL_plant_outbound=cam_rt_fit_hindpawL_results_outbound.midstance;
    
    out.forepawR_lift_outbound=cam_rt_fit_forepawR_results_outbound.midswing;
    out.forepawL_lift_outbound=cam_rt_fit_forepawL_results_outbound.midswing;
    out.hindpawR_lift_outbound=cam_rt_fit_hindpawR_results_outbound.midswing;
    out.hindpawL_lift_outbound=cam_rt_fit_hindpawL_results_outbound.midswing;
    
    out.forepawR_plant_inbound=cam_rt_fit_forepawR_results_inbound.midstance;
    out.forepawL_plant_inbound=cam_rt_fit_forepawL_results_inbound.midstance;
    out.hindpawR_plant_inbound=cam_rt_fit_hindpawR_results_inbound.midstance;
    out.hindpawL_plant_inbound=cam_rt_fit_hindpawL_results_inbound.midstance;
    
    out.forepawR_lift_inbound=cam_rt_fit_forepawR_results_inbound.midswing;
    out.forepawL_lift_inbound=cam_rt_fit_forepawL_results_inbound.midswing;
    out.hindpawR_lift_inbound=cam_rt_fit_hindpawR_results_inbound.midswing;
    out.hindpawL_lift_inbound=cam_rt_fit_hindpawL_results_inbound.midswing;

    out.forepawR_filt_6_8=cdat_forepawR_filt_6_8;
    out.forepawL_filt_6_8=cdat_forepawL_filt_6_8;
    out.hindpawR_filt_6_8=cdat_hindpawR_filt_6_8;
    out.hindpawL_filt_6_8=cdat_hindpawL_filt_6_8;
    
    out.run_periods=run_periods;
    out.outbound_periods=outbound_center;
    out.inbound_periods=inbound_center;
    out.xcorr_forepawR_forepawL=xcorr_forepawR_forepawL;
    out.xcorr_forepawR_forepawL_plants=xcorr_forepawR_forepawL_plants;
    out.cam_rt_fit=cam_rt_fit;
    out.speed_val=speed_val;
    out.est_framerate=est_framerate;
   
    end
end
