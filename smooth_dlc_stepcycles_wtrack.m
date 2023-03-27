function [cam_rt_fit_results,cdat_step_filt_6_8,cdat_step_filt_12_40, run_periods]=smooth_dlc_stepcycles_wtrack(~,dlc_results, cam_rt_fit,cdat_vel, ~, run_periods, ~)

cam_rt_fit_results=[];

% converting into cdat format
cdat_vel=imcont('data', cdat_vel,'timestamp', cam_rt_fit, 'chanlabels',{'velocity'});

%% Method1: 
% Smooth with a moving avg of 5

yy=abs(diff(dlc_results));
xx=cam_rt_fit(1:end-1);
cam_rt_fit_smoothed=xx;

% % find lift or plant times 
cdat_step_filt_v1=imcont('data', yy, 'timestamp', xx);

% % plot results
% figure(7); hold on;
%plot(cam_rt_fit(1:end-1),step_cycle_v1, '.r');
% plot(cam_rt_fit_smoothed,step_cycle_smoothed,'b');
% ylim([-20 40]);

%% Method2:  filter using cont functions
% import data into cont format
cdat_step=imcont('data', yy, 'timestamp', xx, 'chanlabels', {'step_cycle'});

%create a lowpass filter: [6 8] to get midstance and midswing
filtopt_LPF_6_8 = mkfiltopt('name', 'aj_step_LPF', 'filttype', 'lowpass', 'F', [6 8]);

% filtering step cycle
cdat_step_filt_6_8 = contfilt(cdat_step, 'filtopt', filtopt_LPF_6_8);

%plot quickly
% close all; hold on;
% plot(cam_rt_fit_smoothed,cdat_step.data,'r'); hold on;
% plot(cam_rt_fit_smoothed,cdat_step_filt_6_8.data, 'b');
% % plot([pks_t{1,1} pks_t{1,1}],[-5 50],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% % plot([trs_t{1,1} trs_t{1,1}],[-5 50],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);
% plot(cam_rt_fit,dlc_results./20);
% plot([pks_ind pks_ind],[-5 50],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([trs_ind trs_ind],[-5 50],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);
% %need to improve on the troughs :(

%% Method3: 

filtopt_LPF_12_40 = mkfiltopt('name', 'aj_step_LPF', 'filttype', 'lowpass', 'F', [12 40]);
cdat_step_filt_12_40 = contfilt(cdat_step, 'filtopt', filtopt_LPF_12_40);

% plot quickly
% close all;
% plot(cam_rt_fit_smoothed,cdat_step.data,'r'); hold on;
% plot(cam_rt_fit_smoothed,cdat_step_filt_12_40.data, 'b');
% plot([segs_min05_swing(:,2) segs_min05_swing(:,2)],[-5 50],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([segs_min05_stance(:,1) segs_min05_stance(:,1)],[-5 50],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);

%% Find the plant and lift times from diff2 | this method works best to identify plant and lift times 
acceleration_paw=diff(cdat_step_filt_6_8.data);

[~,step_lift_locs]=findpeaks(double(acceleration_paw),cam_rt_fit_smoothed(1:end-1),'MinPeakDistance', 0.125, 'MinPeakHeight', 0.02);
[~,step_plant_locs]=findpeaks(double(-acceleration_paw),cam_rt_fit_smoothed(1:end-1),'MinPeakDistance', 0.125, 'MinPeakHeight', 0.02);

% plot to check if lift corresponds to lift and plant correposponds to plant - correct aj april2021
% figure(1); hold on; x
% plot(cam_rt_fit,(dlc_results), '.k');
% plot(cam_rt_fit_smoothed(1:end-1),acceleration_paw.*100,'.');
% plot([step_lift_locs step_lift_locs], [0 200], 'm');
% plot([step_plant_locs step_plant_locs], [0 200], 'k');

% select only those plants and lifts in between run periods defined by
% tails velocity, now find only those troughs in run_periods as defined by the stepping cycle 

%% run periods are now defined as tracktime

step_lift_locs_run=[];

for k=1:size(run_periods, 1)
    locs_collect_idx=find(step_lift_locs>run_periods(k,1)&step_lift_locs<run_periods(k,2));
    step_lift_locs_run=[step_lift_locs_run;step_lift_locs(locs_collect_idx)];
    clear locs_collect_idx;
end 

% colect plant times within run periods
% step_plant_locs=step_plant_locs';
step_plant_locs_run=[];

for k=1:size(run_periods, 1)
    locs_collect_idx=find(step_plant_locs>run_periods(k,1)&step_plant_locs<run_periods(k,2));
    step_plant_locs_run=[step_plant_locs_run;step_plant_locs(locs_collect_idx)];
    clear locs_collect_idx;
end 

% plot to check if lift corresponds to lift and plant correposponds to plant - correct aj april2021
% close all;
% figure(1); hold on;
% plot(cam_rt_fit,dlc_results, '.k');
% plot(cam_rt_fit_smoothed(1:end-1),acceleration_paw.*100,'.');
% plot([step_lift_locs_run step_lift_locs_run], [0 2000], 'm');
% plot([step_plant_locs_run step_plant_locs_run], [0 2000], 'k');

%% make stance windows and swing windows 
% find first plant

for k=1:size(step_plant_locs_run, 1)-1
    find_lift=find(step_lift_locs_run>step_plant_locs_run(k));
    
    if isempty (find_lift)
       plant_lift_windows(k,:)=nan;
    else 
    plant_lift_windows(k,1)=step_plant_locs_run(k);
    plant_lift_windows(k,2)=step_lift_locs_run(find_lift(1));
    clear find_lift;
    end
end

idx_select=~isnan(plant_lift_windows);
plant_lift_period=plant_lift_windows(idx_select(:,1),:); clear idx_select;


lift_plant_windows=[];
% find first lift
for k=1:size(step_lift_locs_run, 1)-1
    find_plant=find(step_plant_locs_run>step_lift_locs_run(k));
    
    if isempty (find_plant)
       lift_plant_windows(k,:)=nan;
    else 
    lift_plant_windows(k,1)=step_lift_locs_run(k);
    lift_plant_windows(k,2)=step_plant_locs_run(find_plant(1));
    clear find_plant;
    end
end

idx_select=~isnan(lift_plant_windows);
lift_plant_period=lift_plant_windows(idx_select(:,1),:); clear idx_select;

% plot to check that the plant-lift and lift-plant periods correspond to
% the correct time - aj checked april 2021 - P-L(:,1)=plant times
% L-P(:,1)=Lift times 

% plot(cam_rt_fit,dlc_results,'r'); hold on;
% plot(cam_rt_fit_smoothed,cdat_step_filt_12_40.data, 'b');
% plot([lift_plant_period(:,2) lift_plant_period(:,2)],[-5 2000],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([lift_plant_period(:,1) lift_plant_period(:,1)],[-5 2000],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);

%% select only those plant and lifts that are inside run epochs

plant_lift_period_include=find(((plant_lift_period(:,2)-plant_lift_period(:,1))<0.5)&((plant_lift_period(:,2)-plant_lift_period(:,1))>-0.5));
lift_plant_period_include=find(((lift_plant_period(:,2)-lift_plant_period(:,1))<0.5)&((lift_plant_period(:,2)-lift_plant_period(:,1))>-0.5));

plant_lift_period=plant_lift_period(plant_lift_period_include, :);
lift_plant_period=lift_plant_period(lift_plant_period_include,:);

% plot(cam_rt_fit,dlc_results,'r'); hold on;
% plot(cam_rt_fit_smoothed,cdat_step_filt_12_40.data, 'b');
% plot([plant_lift_period(:,2) plant_lift_period(:,2)],[-5 2000],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([plant_lift_period(:,1) plant_lift_period(:,1)],[-5 2000],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);

% calculate the swing and stance phase times in cam_rt_fit windows
% calculate the phase of theta during swing phase 
% collect siwng index and lift time

%%
idx_swing_all=[];
idx_lift=nan(size(lift_plant_period,1),1);
idx_midswing=nan(size(lift_plant_period,1),1);

for i=1:size(lift_plant_period)
    idx_swing=find(cam_rt_fit>lift_plant_period(i,1)&cam_rt_fit<lift_plant_period(i,2));
    
    if isempty(idx_swing)
        idx_swing_all(i,1)=nan;
    else
        
    idx_select=nonzeros(idx_swing(round(1+0.1.*(size(idx_swing,1))):round(1+0.3.*(size(idx_swing,1)))));
    % aj testing terminal swing 
    idx_swing_all(i,round(1+0.1.*(size(idx_swing,1))):round(1+0.3.*size(idx_swing,1)))=nonzeros(idx_swing(round(1+0.1.*(size(idx_swing,1))):round(1+0.3.*(size(idx_swing,1)))));
    idx_lift(i)=idx_select(1);
    idx_midswing(i)=round((idx_select(1)+idx_select(end))/2);
    clear idx_swing idx_select;
    end
end

idx_swing_all=idx_swing_all(:);
idx_swing_all=nonzeros(idx_swing_all);
% remove nans
idx_swing_all(isnan(idx_swing_all))=[];
idx_lift(isnan(idx_lift))=[];
idx_midswing=nonzeros(idx_midswing);
idx_midswing(isnan(idx_midswing))=[];

cam_rt_fit_pks=pks_t;

%Same for stance. Collect stance index
segs_min05_stance=nan(size(lift_plant_period));
for j=1:size(lift_plant_period)-2
    segs_min05_stance(j,1)=lift_plant_period(j,2);
    segs_min05_stance(j,2)=lift_plant_period(j+1,1);
end

a=segs_min05_stance(:,2)-segs_min05_stance(:,1);
idx_stance_select=find(a>0.05&a<0.5);
run_stance_cycle_1=segs_min05_stance(idx_stance_select,1);
run_stance_cycle_2=segs_min05_stance(idx_stance_select,2);
run_stance_cycles(:,1)=run_stance_cycle_1;
run_stance_cycles(:,2)=run_stance_cycle_2;

idx_stance_all=[];
idx_plant=nan(size(run_stance_cycles,1),1);
idx_midstance=nan(size(run_stance_cycles,1),1);

for k=1:size(run_stance_cycles)
    idx_stance=find(cam_rt_fit>run_stance_cycles(k,1)&cam_rt_fit<run_stance_cycles(k,2));
        
    if isempty(idx_stance)
        idx_stance_all(k,1)=nan;
    else
    
    idx_select=idx_stance(round(1+0.1.*(size(idx_stance,1))):round(1+0.3.*(size(idx_stance,1))));
    % aj testing terminal stance 
    idx_stance_all(k,round(1+0.1.*(size(idx_stance,1))):round(0.3.*size(idx_stance,1)))=idx_stance(round(1+0.1.*(size(idx_stance,1))):round(0.3.*(size(idx_stance,1))));    
    %idx_plant(k)=idx_stance_all(k);
    idx_plant(k)=idx_select(1);
    
    idx_midstance(k)=round((idx_select(1)+idx_select(end))/2);
    clear idx_stance idx_select;
    end
end

idx_stance_all=idx_stance_all(:);
idx_stance_all=nonzeros(idx_stance_all);

idx_stance_all=idx_stance_all(:);
idx_stance_all=nonzeros(idx_stance_all);
% remove nans
idx_stance_all(isnan(idx_stance_all))=[];
idx_plant(isnan(idx_plant))=[];
idx_midstance=nonzeros(idx_midstance);
idx_midstance(isnan(idx_midstance))=[];

cam_rt_fit_stance=cam_rt_fit(idx_stance_all); 
cam_rt_fit_swing=cam_rt_fit(idx_swing_all); 
cam_rt_fit_midstance=cam_rt_fit(idx_midstance); 
cam_rt_fit_midswing=cam_rt_fit(idx_midswing); 

% % % % plot quickly
% close all;
% %raw data
% plot(cam_rt_fit,dlc_results./100, 'r');hold on;
% % method1
% plot(cam_rt_fit_smoothed,step_cycle_smoothed,'b');
% %method2
% plot(cam_rt_fit_smoothed,cdat_step_filt_6_8.data, 'k');
% %method3
% plot(cam_rt_fit_smoothed,cdat_step_filt_12_40.data, 'm');
% %lift and plant times as detremined from method 3
% plot([pks_t{1,1} pks_t{1,1}],[-5 50],'b');
% run_swingcyceles(:,1 is the plant time and (:,2 is the lift time))
% plot([lift_plant_period(:,2) lift_plant_period(:,2)],[-5 200],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([lift_plant_period(:,1) lift_plant_period(:,1)],[-5 200],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);
% plot([cam_rt_fit_stance cam_rt_fit_stance],[0 1], 'k')
% plot([cam_rt_fit_midstance cam_rt_fit_midstance], [0 5], 'k')
% plot([cam_rt_fit_swing cam_rt_fit_swing],[0 1], 'm')
% plot([cam_rt_fit_midswing cam_rt_fit_midswing], [0 5], 'm')
% %xlim([2210 2230])
% ylim([-5 15])

cam_rt_fit_results.plant=lift_plant_period(:,2);
cam_rt_fit_results.lift=lift_plant_period(:,1);
cam_rt_fit_results.swing=cam_rt_fit_swing;
cam_rt_fit_results.stance=cam_rt_fit_stance;
cam_rt_fit_results.midswing=cam_rt_fit_midswing;

% PLANT TIMES used in Joshi et al 2023. This is the value that corresponds 
% to the limb of the rat firmly touching the surface of the track 
% ED Figure 1. 
cam_rt_fit_results.midstance=cam_rt_fit_midstance; 

cam_rt_fit_results.peak=pks_t;
cam_rt_fit_results.trough=trs_t;

end