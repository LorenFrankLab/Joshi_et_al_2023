function [cam_rt_fit_results,cdat_step_filt_6_8,cdat_step_filt_12_50, run_periods]=smooth_dlc_stepcycles_lineartrack(~,dlc_results, cam_rt_fit,cdat_vel, est_framerate, run_periods)

cam_rt_fit_results=[];

% converting into cdat format
cdat_vel=imcont('data', cdat_vel,'timestamp', cam_rt_fit, 'chanlabels',{'velocity'});

%% Method1: 

yy=abs(diff(dlc_results));
xx=cam_rt_fit(1:end-1);
cam_rt_fit_smoothed=xx;

% % find lift or plant times | 
cdat_step_filt_v1=imcont('data', yy, 'timestamp', xx);

% % plot results
% figure(7); hold on;
%plot(cam_rt_fit(1:end-1),step_cycle_v1, '.r');
% plot(cam_rt_fit_smoothed,step_cycle_smoothed,'b');
% ylim([-20 40]);

%% Method2:
% import data into cont format
cdat_step=imcont('data', yy, 'timestamp', xx, 'chanlabels', {'step_cycle'});

%create a lowpass filter: [6 8] to get midstance and midswing
filtopt_LPF_6_8 = mkfiltopt('name', 'aj_step_LPF', 'filttype', 'lowpass', 'F', [6 8]);

% filtering step cycle
cdat_step_filt_6_8 = contfilt(cdat_step, 'filtopt', filtopt_LPF_6_8); % 

% select peaks and troughs crossings | these are the midstance and midswing
% times 
[pks_t, ~]= contpeaks(cdat_step_filt_6_8, 'type', 'peaks', 'thresh', 0.4, 'segs', run_periods);
[trs_t, ~]=contpeaks(cdat_step_filt_6_8,'type', 'valleys', 'thresh', 0.1, 'segs', run_periods);

% %plot quickly
% close all;
% plot(cam_rt_fit_smoothed,cdat_step.data,'r'); hold on;
% plot(cam_rt_fit_smoothed,cdat_step_filt_6_8.data, 'b');
% plot([pks_t{1,1} pks_t{1,1}],[-5 50],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([trs_t{1,1} trs_t{1,1}],[-5 50],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);
% plot(cam_rt_fit,dlc_results./20);

%% Method3:

filtopt_LPF_12_50 = mkfiltopt('name', 'aj_step_LPF', 'filttype', 'lowpass', 'F', [12 40]);
cdat_step_filt_12_50 = contfilt(cdat_step, 'filtopt', filtopt_LPF_12_50);

% plot quickly
close all;
plot(cam_rt_fit_smoothed,cdat_step.data,'r'); hold on;
plot(cam_rt_fit_smoothed,cdat_step_filt_12_50.data, 'b');
% plot([segs_min05_swing(:,2) segs_min05_swing(:,2)],[-5 50],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
% plot([segs_min05_stance(:,1) segs_min05_stance(:,1)],[-5 50],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);

%% Find the plant and lift times from diff2 | this method works best to identfy plant and lift times
acceleration_paw=diff(cdat_step_filt_6_8.data);

[~,step_lift_locs]=findpeaks(double(acceleration_paw),cam_rt_fit_smoothed(1:end-1),'MinPeakDistance', 0.100, 'MinPeakHeight', 0.03);
[~,step_plant_locs]=findpeaks(double(-acceleration_paw),cam_rt_fit_smoothed(1:end-1),'MinPeakDistance', 0.100, 'MinPeakHeight', 0.03);

% close all; figure(1); hold on;
% plot(cam_rt_fit,(dlc_results)./10, '.k');
% plot(cam_rt_fit(1:end-2),acceleration_paw.*20-5,'-');
% plot([step_lift_locs step_lift_locs], [-5 15], 'm');
% plot([step_plant_locs step_plant_locs], [-5 15], 'k');
% ylim([-7 12])
% xlim([2448 2451])
% select only those plants and lifts in between run periods defined by
% tails velocity, now find only those troughs in run_periods as defined by the stepping cycle 

%% run periods are now defined as tracktime

% step_lift_locs=step_lift_locs';
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

if length(step_plant_locs_run)>5
    
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

    %select only those plant and lifts that are inside run epochs

    % % calculate the swing and stance phase times in cam_rt_fit windows
    % calculate the phase of theta during swing phase 
    % collect siwng index and lift time

    plant_lift_period_include=find(((plant_lift_period(:,2)-plant_lift_period(:,1))<0.5)&((plant_lift_period(:,2)-plant_lift_period(:,1))>-0.5));
    lift_plant_period_include=find(((lift_plant_period(:,2)-lift_plant_period(:,1))<0.5)&((lift_plant_period(:,2)-lift_plant_period(:,1))>-0.5));

    plant_lift_period=plant_lift_period(plant_lift_period_include, :);
    lift_plant_period=lift_plant_period(lift_plant_period_include,:);

    idx_swing_all=[];
    idx_lift=nan(size(lift_plant_period,1),1);
    idx_midswing=nan(size(lift_plant_period,1),1);

    for i=1:size(lift_plant_period)
        idx_swing=find(cam_rt_fit>lift_plant_period(i,1)&cam_rt_fit<lift_plant_period(i,2));

        if isempty(idx_swing)
            idx_swing_all(i,1)=nan;
        else

        idx_select=idx_swing(round(1+0.1.*(size(idx_swing,1))):round(1+0.3.*(size(idx_swing,1))));
        idx_swing_all(i,round(1+0.1.*(size(idx_swing,1))):round(1+0.3.*size(idx_swing,1)))=idx_swing(round(1+0.1.*(size(idx_swing,1))):round(1+0.3.*(size(idx_swing,1))));

        %idx_lift(i)=idx_swing(1);
        idx_lift(i)=idx_select(1);
        idx_midswing(i)=round((idx_select(1)+idx_select(end))/2);
        %idx_midswing(i,1:size(idx_swing)/5)=idx_swing(1:idx_swing/5);
        clear idx_swing; idx_select;
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
        idx_stance_all(k,round(1+0.1.*(size(idx_stance,1))):round(1+0.30.*size(idx_stance,1)))=idx_stance(round(1+0.1.*(size(idx_stance,1))):round(1+0.3.*(size(idx_stance,1))));    
        idx_plant(k)=idx_select(1);
        idx_midstance(k)=round((idx_select(1)+idx_select(end))/2);
        clear idx_stance;clear idx_select;
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

    % % % plot quickly
    % close all;
    % %raw data
    % plot(cam_rt_fit,dlc_results./100, '.r');hold on;
    % % % method1
    % % plot(cam_rt_fit_smoothed,step_cycle_smoothed,'b');
    % % %method2
    % % plot(cam_rt_fit_smoothed,cdat_step_filt_6_8.data, 'k');
    % % %method3
    % % plot(cam_rt_fit_smoothed,cdat_step_filt_12_50.data, 'm');
    % % %lift and plant times as detremined from method 3
    % % plot([pks_t{1,1} pks_t{1,1}],[-5 50],'b');
    % % run_swingcyceles(:,1 is the plant time and (:,2 is the lift time))
    % plot([lift_plant_period(:,2) lift_plant_period(:,2)],[-5 2000],'color',[0.3686 0.2353 0.6000], 'LineWidth', 0.5);
    % plot([lift_plant_period(:,1) lift_plant_period(:,1)],[-5 2000],'color',[0.4 0.6353 0.4000], 'LineWidth', 0.5);
    % plot([cam_rt_fit_stance cam_rt_fit_stance],[0 1], 'k')
    % plot([cam_rt_fit_swing cam_rt_fit_swing],[0 1], 'm')
    % %xlim([2210 2230])
    % ylim([-5 15])

    cam_rt_fit_results.plant=lift_plant_period(:,2);
    cam_rt_fit_results.lift=lift_plant_period(:,1);
    cam_rt_fit_results.swing=cam_rt_fit_swing;
    cam_rt_fit_results.stance=cam_rt_fit_stance;
    cam_rt_fit_results.midswing=cam_rt_fit_midswing;
    cam_rt_fit_results.peak=pks_t{1,1};
    cam_rt_fit_results.trough=trs_t{1,1};

    % This is the value that corresponds to the limb of the rat firmly touching
    % the surface of the track (ED Figure 1).   
    cam_rt_fit_results.midstance=cam_rt_fit_midstance;

    else

    cam_rt_fit_results.plant=NaN;
    cam_rt_fit_results.lift=NaN;
    cam_rt_fit_results.swing=NaN;
    cam_rt_fit_results.stance=NaN;

end
end