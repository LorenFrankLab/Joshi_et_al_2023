% dfs_speed_steps_theta.m 
% -> dfa_speed_steps_theta_forelimb
% -> Figure1EF_EDFigure2.mat

 %% Example to create task phase specific data from proprocessed files  
% % set filters
% clear all; close all; clc; 
% 
% %animal filters | enter animal names to anlyse 
% animals={'Jaq' 'Roqui', 'Lotus', 'Monty', 'Peanut'};
% 
% for a = 1:length(animals)
% 
%     % epoch filters | select wtrack epochs
%     epochfilter{1} = 'isequal($type, ''run'') && isequal($environment, ''wtrack'')';
% 
%     % iterators | these are goimg to define your base unit of analysis, such as per cell or per epoch
%     iterator= 'epocheeganal'; 
% 
%     % goes into tetinfo and searches for the correct tetrodes to be included 
%     tetfilter = 'isequal($area,''ca1Rref'')';
% 
%     % timefilter is a filter that excludes a selected time window
%     timefilter{1} = {'aj_get2dstate_dlc', '($mobility == 1)','mobility_velocity',[04 200], 'mobility_buffer',0.25, 'min_durations',1, 'vel2choose', 'nose'};
% 
%     f(a) = createfilter('animal',animals{a},'epochs',epochfilter,'excludetime', timefilter, 'eegtetrodes',tetfilter, 'iterator', iterator);
%     f(a) = setfilterfunction(f(a), 'dfa_speed_steps_theta_forelimb', {'eeg', 'eeggnd', 'theta' ,'posdlc'},'animal',animals{a},'out_min',30, 'out_max',100, 'time_step',0.125, 'a',a);
%     f(a)= runfilter(f(a));
% end
% 
% fname ='';
% save(fname,'f', '-v7.3')
  
%% Example Outbound runs task phase
close all; clear all; 
results_filename='';
load([results_filename, '.mat'])
clearvars -except f;

savefigs=0; savefig_scatter=0; 
% parameters
ci_alpha=0.01; % corresponds to 99% CI 
destdir='';

animals={'Jaq' 'Roqui', 'Lotus', 'Monty', 'Peanut'};

close all
for a=[1,2,3,4,5]
    
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    
    ind_del=[];
    for i=1:size(ind,1)
        
        if ((isnan(f(a).output{1, 1}(i).inst_acc_outbound(1)))|(isnan(f(a).output{1, 1}(i).speed_outbound(1))))
             ind_del=[ind_del; i];
        end
    end
    f(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    
    %keyboard
    count=size(ind_select,1);
    
    %keyboard
    speed_outbound(a) = {cell2mat(arrayfun(@(x) x.speed_outbound,f(a).output{1}(ind_select),'Un',0))'};
    distance_outbound(a) = {cell2mat(arrayfun(@(x) x.distance_outbound,f(a).output{1}(ind_select),'Un',0))'};
    forelimb_steps_outbound(a) = {cell2mat(arrayfun(@(x) x.forelimb_outbound_steps,f(a).output{1}(ind_select),'Un',0))'};
    theta_count_outbound(a) = {cell2mat(arrayfun(@(x) x.theta_count_outbound,f(a).output{1}(ind_select),'Un',0))'};
    run_duration_outbound(a)={cell2mat(arrayfun(@(x) x.run_duration_outbound,f(a).output{1}(ind_select),'Un',0))'};
    
    inst_vel_outbound(a)= {cell2mat((arrayfun(@(x) x.inst_vel_outbound,f(a).output{1}(ind_select),'Un',0))')};
    inst_acc_outbound(a)= {cell2mat((arrayfun(@(x) x.inst_acc_outbound,f(a).output{1}(ind_select),'Un',0))')};
    inst_freq_outbound(a)={cell2mat((arrayfun(@(x) x.inst_theta_outbound,f(a).output{1}(ind_select),'Un',0))')};
    inst_steps_outbound(a)={cell2mat((arrayfun(@(x) x.inst_steps_outbound,f(a).output{1}(ind_select),'Un',0))')};
    
    %plot data per animal
    %keyboard
    speed_all=(speed_outbound{a}(:));
    forelimb_steps_all=(forelimb_steps_outbound{a}(:));
    theta_count_all=(theta_count_outbound{a}(:));
    run_duration_all=(run_duration_outbound{a}(:));
    distance_all=(distance_outbound{a}(:));
    idx=find(run_duration_all(:)<4);
    
    %plot data per animal 
    inst_vel_all{a}=(inst_vel_outbound{a}(:));
    inst_acc_all{a}=(inst_acc_outbound{a}(:));
    inst_theta_all{a}=(inst_freq_outbound{a}(:));
    inst_steps_all{a}=(inst_steps_outbound{a}(:));
    
    idx_include{a}=find((inst_vel_all{a}>0) & (inst_steps_all{a}>0) & (inst_theta_all{a}>0));
    
    counter(a)=count;
    clear count speed_all forelimb_steps_all theta_count_all run_duration_all distance_all idx ind_select
end

idx_include_all=cell2mat(idx_include(:));
inst_steps_all=cell2mat(inst_steps_all(:));
inst_vel_all=cell2mat(inst_vel_all(:));
inst_theta_all=cell2mat(inst_theta_all(:));
inst_acc_all=cell2mat(inst_acc_all(:));

idx_include=find((inst_vel_all>vel_thresh_min) & (inst_steps_all>steps_thresh_min) & (inst_steps_all<steps_thresh_max) & (inst_theta_all>theta_thresh_min));
idx_include_vel=find((inst_vel_all>vel_thresh_min)&(inst_vel_all<vel_thresh_max));
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);

    %% SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'
    clear all; load('Figure1EF_EDFigure2.mat');

    % step-theta 
    close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
    scatter(center_outbound.inst_steps_all_plot,center_outbound.inst_theta_all_plot, 'w');
    binscatter(center_outbound.inst_steps_all_plot,center_outbound.inst_theta_all_plot, 100);lsline;
    [r3, p3]=corrcoef(center_outbound.inst_steps_all_plot,center_outbound.inst_theta_all_plot, 'rows', 'complete');
    xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2)]), 'Color', 'r', 'FontSize', 12);

%% or continue on from here if doing your own analysis :)

total_epochs=sum(counter); cd(destdir);

close all; figure(1); hold on; xlim([0 120]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r1, p1]=corrcoef(inst_vel_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r1(2) p1(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstVelTheta_30_100_outbound_nose.fig');gcf;
end

close all; figure(2); hold on; xlim([-1 1]); ylim([5 12]); alpha(0.5); 
scatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r2, p2]=corrcoef(inst_acc_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete'); % note this accelartion need to be multiplied by the estimated framerate
xlabel('Instantaneous Acceleration'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r2(2) p2(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstAccTheta_30_100_outbound_nose.fig');gcf;
end

close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r3, p3]=corrcoef(inst_steps_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstStepsTheta_30_100_outbound_nose.fig');gcf;
end

close all; figure(4); hold on; xlim([0 120]); ylim([0 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'w');
binscatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 100);lsline;
[r4, p4]=corrcoef(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Forelimb'); text(4,4.7,num2str([r4(2) p4(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1 
    savefig('InstVelSteps_30_100_outbound_nose.fig');gcf;
end

close all; figure(5); hold on; xlim([0 120]); alpha(0.5); %lsline;
hist(inst_vel_all(idx_include_vel),100);
xlabel('Instantaneous Velocity'); ylabel('Count'); 
if savefig_scatter==1 
    savefig('InstVelHist_30_100_outbound_nose.fig');gcf;
end

%% save variables for regression analysis
inst_steps_all_outbound=inst_steps_all;
inst_theta_all_outbound=inst_theta_all;
inst_acc_all_outbound=inst_acc_all;
inst_vel_all_outbound=inst_vel_all;
idx_include_outbound=idx_include;

rp_outbound(:,1)=[r1(2); r2(2); r3(2); r4(2)];
rp_outbound(:,2)=[p1(2); p2(2); p3(2); p4(2)];

%% plot the scatter plot for inbound and outbound 
inst_steps_all_outbound_z=zscore(inst_steps_all_outbound(idx_include_outbound));
inst_theta_all_outbound_z=zscore(inst_theta_all_outbound(idx_include_outbound));
inst_acc_all_outbound_z=zscore(inst_acc_all_outbound(idx_include_outbound));
inst_vel_all_outbound_z=zscore(inst_vel_all_outbound(idx_include_outbound));

lm1 = fitlm(inst_vel_all_outbound_z,inst_theta_all_outbound_z);
ci1 = coefCI(lm1, ci_alpha); 

lm2 = fitlm(inst_acc_all_outbound_z,inst_theta_all_outbound_z);
ci2 = coefCI(lm2, ci_alpha); 

lm3 = fitlm(inst_steps_all_outbound_z,inst_theta_all_outbound_z);
ci3 = coefCI(lm3, ci_alpha); 

lm4 = fitlm(inst_vel_all_outbound_z,inst_steps_all_outbound_z);
ci4 = coefCI(lm4, ci_alpha); 

rp_outbound(:,3)=[ci1(2); ci2(2); ci3(2); ci4(2)]; % lower 95%ci 
rp_outbound(:,4)=[ci1(4); ci2(4); ci3(4); ci4(4)]; % upper 95%ci
%% INBOUND
clearvars -except f savefig_scatter destdir inst_steps_all_outbound inst_theta_all_outbound idx_include_outbound speed_outbound rp_outbound ci_alpha
animals={'Jaq' 'Roqui', 'Lotus', 'Monty', 'Peanut'};

for a=[1,2,3,4,5]
    
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    
    ind_del=[];
    for i=1:size(ind,1)
        
        if ((isnan(f(a).output{1, 1}(i).inst_acc_inbound(1)))||(isnan(f(a).output{1, 1}(i).speed_inbound(1))))
             ind_del=[ind_del; i];
        end
    end
    f(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    
    count=size(ind_select,1);

    %keyboard
    speed_inbound(a) = {cell2mat(arrayfun(@(x) x.speed_inbound,f(a).output{1}(ind_select),'Un',0))'};
    distance_inbound(a) = {cell2mat(arrayfun(@(x) x.distance_inbound,f(a).output{1}(ind_select),'Un',0))'};
    forelimb_steps_inbound(a) = {cell2mat(arrayfun(@(x) x.forelimb_inbound_steps,f(a).output{1}(ind_select),'Un',0))'};
    theta_count_inbound(a) = {cell2mat(arrayfun(@(x) x.theta_count_inbound,f(a).output{1}(ind_select),'Un',0))'};
    run_duration_inbound(a)={cell2mat(arrayfun(@(x) x.run_duration_inbound,f(a).output{1}(ind_select),'Un',0))'};
    
    inst_vel_inbound(a)= {cell2mat((arrayfun(@(x) x.inst_vel_inbound,f(a).output{1}(ind_select),'Un',0))')};
    inst_acc_inbound(a)= {cell2mat((arrayfun(@(x) x.inst_acc_inbound,f(a).output{1}(ind_select),'Un',0))')};
    inst_freq_inbound(a)={cell2mat((arrayfun(@(x) x.inst_theta_inbound,f(a).output{1}(ind_select),'Un',0))')};
    inst_steps_inbound(a)={cell2mat((arrayfun(@(x) x.inst_steps_inbound,f(a).output{1}(ind_select),'Un',0))')};
    
    %plot data per animal
    %keyboard
    speed_all=(speed_inbound{a}(:));
    forelimb_steps_all=(forelimb_steps_inbound{a}(:));
    theta_count_all=(theta_count_inbound{a}(:));
    run_duration_all=(run_duration_inbound{a}(:));
    distance_all=(distance_inbound{a}(:));
    idx=find(run_duration_all(:)<4);
    
    %plot data per animal 
    inst_vel_all{a}=(inst_vel_inbound{a}(:));
    inst_acc_all{a}=(inst_acc_inbound{a}(:));
    inst_theta_all{a}=(inst_freq_inbound{a}(:));
    inst_steps_all{a}=(inst_steps_inbound{a}(:));
    
    idx_include{a}=find((inst_vel_all{a}>0) & (inst_steps_all{a}>0) & (inst_theta_all{a}>0));
     
    counter(a)=count;
     clear count speed_all forelimb_steps_all theta_count_all run_duration_all distance_all idx ind_select

end

idx_include_all=cell2mat(idx_include(:));
inst_steps_all=cell2mat(inst_steps_all(:));
inst_vel_all=cell2mat(inst_vel_all(:));
inst_theta_all=cell2mat(inst_theta_all(:));
inst_acc_all=cell2mat(inst_acc_all(:));

idx_include=find((inst_vel_all>vel_thresh_min) & (inst_steps_all>steps_thresh_min) & (inst_steps_all<steps_thresh_max) & (inst_theta_all>theta_thresh_min));
idx_include_vel=find((inst_vel_all>vel_thresh_min)&(inst_vel_all<vel_thresh_max));
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);

    %% SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'
    clear all; load('Figure1EF_EDFigure2.mat');

    % step-theta 
    close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
    scatter(center_inbound.inst_steps_all_plot,center_inbound.inst_theta_all_plot, 'w');
    binscatter(center_inbound.inst_steps_all_plot,center_inbound.inst_theta_all_plot, 100);lsline;
    [r3, p3]=corrcoef(center_inbound.inst_steps_all_plot,center_inbound.inst_theta_all_plot, 'rows', 'complete');
    xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2)]), 'Color', 'r', 'FontSize', 12);

%% or continue on from here if doing your own analysis :)
total_epochs=sum(counter); cd(destdir);

close all; figure(1); hold on; xlim([0 120]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r1, p1]=corrcoef(inst_vel_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r1(2) p1(2) total_epochs]), 'Color', 'r', 'FontSize', 12);

if savefig_scatter==1
    savefig('InstVelTheta_30_100_inbound_nose.fig');gcf;
end

close all; figure(2); hold on; xlim([-1 1]); ylim([5 12]); alpha(0.5); 
scatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r2, p2]=corrcoef(inst_acc_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Acceleration'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r2(2) p2(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstAccTheta_30_100_inbound_nose.fig');gcf;
end

close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r3, p3]=corrcoef(inst_steps_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstStepsTheta_30_100_inbound_nose.fig');gcf;
end

close all; figure(4); hold on; xlim([0 120]); ylim([0 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'w');
binscatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 100);lsline;
[r4, p4]=corrcoef(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Forelimb'); text(4,4.7,num2str([r4(2) p4(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1 
    savefig('InstVelSteps_30_100_inbound_nose.fig');gcf;
end

close all; figure(5); hold on; xlim([0 120]); alpha(0.5); %lsline;
hist(inst_vel_all(idx_include_vel),100);
xlabel('Instantaneous Velocity'); ylabel('Count'); 
if savefig_scatter==1 
    savefig('InstVelHist_30_100_inbound_nose.fig');gcf;
end

%% save variables for regression analysis
inst_steps_all_inbound=inst_steps_all;
inst_theta_all_inbound=inst_theta_all;
inst_acc_all_inbound=inst_acc_all;
inst_vel_all_inbound=inst_vel_all;
idx_include_inbound=idx_include;

rp_inbound(:,1)=[r1(2); r2(2); r3(2); r4(2)];
rp_inbound(:,2)=[p1(2) p2(2) p3(2) p4(2)];

%% CI fr inbound
inst_steps_all_inbound_z=zscore(inst_steps_all_inbound(idx_include_inbound));
inst_theta_all_inbound_z=zscore(inst_theta_all_inbound(idx_include_inbound));
inst_acc_all_inbound_z=zscore(inst_acc_all_inbound(idx_include_inbound));
inst_vel_all_inbound_z=zscore(inst_vel_all_inbound(idx_include_inbound));

lm1 = fitlm(inst_vel_all_inbound_z,inst_theta_all_inbound_z);
ci1 = coefCI(lm1, ci_alpha); 

lm2 = fitlm(inst_acc_all_inbound_z,inst_theta_all_inbound_z);
ci2 = coefCI(lm2, ci_alpha); 

lm3 = fitlm(inst_steps_all_inbound_z,inst_theta_all_inbound_z);
ci3 = coefCI(lm3, ci_alpha); 

lm4 = fitlm(inst_vel_all_inbound_z,inst_steps_all_inbound_z);
ci4 = coefCI(lm4, ci_alpha); 

rp_inbound(:,3)=[ci1(2); ci2(2); ci3(2); ci4(2)]; % lower 95%ci 
rp_inbound(:,4)=[ci1(4); ci2(4); ci3(4); ci4(4)]; % upper 95%ci

%% Outer arms combined - Inbound
clearvars -except f savefig_scatter destdir rp_inbound rp_outbound ci_alpha

results_dir=''; cd(results_dir);
results_filename=[''];
load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};
rel_dist_forelimb_outbound=[];rel_dist_forlimb_inbound=[];rel_dist_hindlimb_outbound=[];rel_dist_hinlimb_inbound=[];
f=f1;

for a=[1,2,3,4,5]
    
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    
    ind_del=[]; ind_del_f1=[]; ind_del_f2=[];
    for i=1:size(ind,1)
        
        if ((isnan(f1(a).output{1, 1}(i).inst_acc_inbound(1)))||(isnan(f1(a).output{1, 1}(i).speed_inbound(1))))
             ind_del_f1=[ind_del_f1; i];
        end
    end

    for i=1:size(ind,1)
        
        if ((isnan(f2(a).output{1, 1}(i).inst_acc_inbound(1)))||(isnan(f2(a).output{1, 1}(i).speed_inbound(1))))
             ind_del_f2=[ind_del_f2; i];
        end
    end
    
    ind_del=unique([ind_del_f1;ind_del_f2], 'rows');
    
    f1(a).output{1, 1}(ind_del)=[];
    f2(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    count=size(ind_select,1);
    
    % keyboard 
    inst_vel_inbound_f1= cell2mat((arrayfun(@(x) x.inst_vel_inbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_acc_inbound_f1= cell2mat((arrayfun(@(x) x.inst_acc_inbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_freq_inbound_f1=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_steps_inbound_f1=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f1(a).output{1}(ind_select),'Un',0))');
    
    inst_vel_inbound_f2=cell2mat((arrayfun(@(x) x.inst_vel_inbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_acc_inbound_f2=cell2mat((arrayfun(@(x) x.inst_acc_inbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_freq_inbound_f2=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_steps_inbound_f2=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f2(a).output{1}(ind_select),'Un',0))');
    
    % plot data per animal 
    inst_vel_all{a}=([inst_vel_inbound_f1(:); inst_vel_inbound_f2(:)]);
    inst_acc_all{a}=([inst_acc_inbound_f1(:); inst_acc_inbound_f2(:)]);
    inst_theta_all{a}=([inst_freq_inbound_f1(:); inst_freq_inbound_f2(:)]);
    inst_steps_all{a}=([inst_steps_inbound_f1(:); inst_steps_inbound_f2(:)]);
    
    idx_include{a}=find((inst_vel_all{a}>0) & (inst_steps_all{a}>0) & (inst_theta_all{a}>0));
     
    counter(a)=count;
    clear count speed_all forelimb_steps_all theta_count_all run_duration_all distance_all idx ind_select inst_vel_inbound_f1 inst_acc_inbound_f1 inst_freq_inbound_f1 ...
        inst_steps_inbound_f1 inst_vel_inbound_f2 inst_acc_inbound_f2 inst_freq_inbound_f2 inst_steps_inbound_f2

end

idx_include_all=cell2mat(idx_include(:));
inst_steps_all=cell2mat(inst_steps_all(:));
inst_vel_all=cell2mat(inst_vel_all(:));
inst_theta_all=cell2mat(inst_theta_all(:));
inst_acc_all=cell2mat(inst_acc_all(:));

idx_include=find((inst_vel_all>vel_thresh_min) & (inst_steps_all>steps_thresh_min) & (inst_steps_all<steps_thresh_max) & (inst_theta_all>theta_thresh_min));
idx_include_vel=find((inst_vel_all>vel_thresh_min)&(inst_vel_all<vel_thresh_max));
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);

    %% SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'
    clear all; load('Figure1EF_EDFigure2.mat');

    % step-theta 
    close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
    scatter(outer_inbound.inst_steps_all_plot,outer_inbound.inst_theta_all_plot, 'w');
    binscatter(outer_inbound.inst_steps_all_plot,outer_inbound.inst_theta_all_plot, 100);lsline;
    [r3, p3]=corrcoef(outer_inbound.inst_steps_all_plot, outer_inbound.inst_theta_all_plot, 'rows', 'complete');
    xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2)]), 'Color', 'r', 'FontSize', 12);

%% or continue on from here if doing your own analysis :)
% total number of epochs  
total_epochs=sum(counter); cd(destdir);

close all; figure(1); hold on; xlim([0 120]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r1, p1]=corrcoef(inst_vel_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r1(2) p1(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstVelTheta_outerInbound_nose.fig');gcf;
end

close all; figure(2); hold on; xlim([-1 1]); ylim([5 12]); alpha(0.5); 
scatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r2, p2]=corrcoef(inst_acc_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Acceleration'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r2(2) p2(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstAccTheta_outerInbound_nose.fig');gcf;
end

close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r3, p3]=corrcoef(inst_steps_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstStepsTheta_outerInbound_nose.fig');gcf;
end

close all; figure(4); hold on; xlim([0 120]); ylim([0 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'w');
binscatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 100);lsline;
[r4, p4]=corrcoef(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Forelimb'); text(4,4.7,num2str([r4(2) p4(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1 
    savefig('InstVelSteps_outerInbound_nose.fig');gcf;
end

close all; figure(5); hold on; xlim([0 120]); alpha(0.5); %lsline;
hist(inst_vel_all(idx_include_vel),100);
xlabel('Instantaneous Velocity'); ylabel('Count'); 
if savefig_scatter==1 
    savefig('InstVelHist_outerInbound_nose.fig');gcf;
end

%% save variables for regression analysis
inst_steps_all_inbound=inst_steps_all;
inst_theta_all_inbound=inst_theta_all;
inst_acc_all_inbound=inst_acc_all;
inst_vel_all_inbound=inst_vel_all;
idx_include_inbound=idx_include;

rp_outer_inbound(:,1)=[r1(2); r2(2); r3(2); r4(2)];
rp_outer_inbound(:,2)=[p1(2) p2(2) p3(2) p4(2)];

% CI for inbound
inst_steps_all_inbound_z=zscore(inst_steps_all_inbound(idx_include_inbound));
inst_theta_all_inbound_z=zscore(inst_theta_all_inbound(idx_include_inbound));
inst_acc_all_inbound_z=zscore(inst_acc_all_inbound(idx_include_inbound));
inst_vel_all_inbound_z=zscore(inst_vel_all_inbound(idx_include_inbound));

lm1 = fitlm(inst_vel_all_inbound_z,inst_theta_all_inbound_z);
ci1 = coefCI(lm1, ci_alpha); 

lm2 = fitlm(inst_acc_all_inbound_z,inst_theta_all_inbound_z);
ci2 = coefCI(lm2, ci_alpha); 

lm3 = fitlm(inst_steps_all_inbound_z,inst_theta_all_inbound_z);
ci3 = coefCI(lm3, ci_alpha); 

lm4 = fitlm(inst_vel_all_inbound_z,inst_steps_all_inbound_z);
ci4 = coefCI(lm4, ci_alpha); 

rp_outer_inbound(:,3)=[ci1(2); ci2(2); ci3(2); ci4(2)]; % lower 95%ci 
rp_outer_inbound(:,4)=[ci1(4); ci2(4); ci3(4); ci4(4)]; % upper 95%ci

%% Outer arms combined - Outbound
clearvars -except f f1 f2 savefig_scatter destdir animals rp_outbound ci_alpha rp_inbound rp_outer_inbound

for a=[1,2,3,4,5]
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    
    ind_del=[]; ind_del_f1=[]; ind_del_f2=[];
    for i=1:size(ind,1)
        if ((isnan(f1(a).output{1, 1}(i).inst_acc_outbound (1)))||(isnan(f1(a).output{1, 1}(i).speed_outbound (1))))
             ind_del_f1=[ind_del_f1; i];
        end
    end

    for i=1:size(ind,1)
        if ((isnan(f2(a).output{1, 1}(i).inst_acc_outbound (1)))||(isnan(f2(a).output{1, 1}(i).speed_outbound (1))))
             ind_del_f2=[ind_del_f2; i];
        end
    end
    
    ind_del=unique([ind_del_f1;ind_del_f2], 'rows');
    
    f1(a).output{1, 1}(ind_del)=[];
    f2(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    count=size(ind_select,1);

    % keyboard 
    inst_vel_outbound_f1= cell2mat((arrayfun(@(x) x.inst_vel_outbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_acc_outbound_f1= cell2mat((arrayfun(@(x) x.inst_acc_outbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_freq_outbound_f1=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_steps_outbound_f1=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f1(a).output{1}(ind_select),'Un',0))');
    
    inst_vel_outbound_f2=cell2mat((arrayfun(@(x) x.inst_vel_outbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_acc_outbound_f2=cell2mat((arrayfun(@(x) x.inst_acc_outbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_freq_outbound_f2=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_steps_outbound_f2=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f2(a).output{1}(ind_select),'Un',0))');
    
    % plot data per animal 
    inst_vel_all{a}=([inst_vel_outbound_f1(:); inst_vel_outbound_f2(:)]);
    inst_acc_all{a}=([inst_acc_outbound_f1(:); inst_acc_outbound_f2(:)]);
    inst_theta_all{a}=([inst_freq_outbound_f1(:); inst_freq_outbound_f2(:)]);
    inst_steps_all{a}=([inst_steps_outbound_f1(:); inst_steps_outbound_f2(:)]);
    
    idx_include{a}=find((inst_vel_all{a}>0) & (inst_steps_all{a}>0) & (inst_theta_all{a}>0));
     
    counter(a)=count;
    clear count speed_all forelimb_steps_all theta_count_all run_duration_all distance_all idx ind_select inst_vel_outbound_f1 inst_acc_outbound_f1 inst_freq_outbound_f1 ...
        inst_steps_outbound_f1 inst_vel_outbound_f2 inst_acc_outbound_f2 inst_freq_outbound_f2 inst_steps_outbound_f2

end

idx_include_all=cell2mat(idx_include(:));
inst_steps_all=cell2mat(inst_steps_all(:));
inst_vel_all=cell2mat(inst_vel_all(:));
inst_theta_all=cell2mat(inst_theta_all(:));
inst_acc_all=cell2mat(inst_acc_all(:));

idx_include=find((inst_vel_all>vel_thresh_min) & (inst_steps_all>steps_thresh_min) & (inst_steps_all<steps_thresh_max) & (inst_theta_all>theta_thresh_min));
idx_include_vel=find((inst_vel_all>vel_thresh_min)&(inst_vel_all<vel_thresh_max));
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);

    %% SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'
    clear all; load('Figure1EF_EDFigure2.mat');

    % step-theta 
    close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
    scatter(outer_outbound.inst_steps_all_plot,outer_outbound.inst_theta_all_plot, 'w');
    binscatter(outer_outbound.inst_steps_all_plot,outer_outbound.inst_theta_all_plot, 100);lsline;
    [r3, p3]=corrcoef(outer_outbound.inst_steps_all_plot,outer_outbound.inst_theta_all_plot, 'rows', 'complete');
    xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2)]), 'Color', 'r', 'FontSize', 12);

%% or continue on from here if doing your own analysis :)
% total number of epochs  
total_epochs=sum(counter); cd(destdir);

close all; figure(1); hold on; xlim([0 120]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r1, p1]=corrcoef(inst_vel_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r1(2) p1(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstVelTheta_outeroutbound_nose.fig');gcf;
end

close all; figure(2); hold on; xlim([-1 1]); ylim([5 12]); alpha(0.5); 
scatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r2, p2]=corrcoef(inst_acc_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Acceleration'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r2(2) p2(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstAccTheta_outeroutbound_nose.fig');gcf;
end

close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r3, p3]=corrcoef(inst_steps_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstStepsTheta_outeroutbound_nose.fig');gcf;
end

close all; figure(4); hold on; xlim([0 120]); ylim([0 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'w');
binscatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 100);lsline;
[r4, p4]=corrcoef(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Forelimb'); text(4,4.7,num2str([r4(2) p4(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1 
    savefig('InstVelSteps_outeroutbound_nose.fig');gcf;
end

close all; figure(5); hold on; xlim([0 120]); alpha(0.5); %lsline;
hist(inst_vel_all(idx_include_vel),100);
xlabel('Instantaneous Velocity'); ylabel('Count'); 
if savefig_scatter==1 
    savefig('InstVelHist_outeroutbound_nose.fig');gcf;
end

%% save variables for regression analysis
inst_steps_all_outbound=inst_steps_all;
inst_theta_all_outbound=inst_theta_all;
inst_acc_all_outbound=inst_acc_all;
inst_vel_all_outbound=inst_vel_all;
idx_include_outbound=idx_include;

rp_outer_outbound(:,1)=[r1(2); r2(2); r3(2); r4(2)];
rp_outer_outbound(:,2)=[p1(2) p2(2) p3(2) p4(2)];

% plot the scatter plot for inbound and outbound 
inst_steps_all_outbound_z=zscore(inst_steps_all_outbound(idx_include_outbound));
inst_theta_all_outbound_z=zscore(inst_theta_all_outbound(idx_include_outbound));
inst_acc_all_outbound_z=zscore(inst_acc_all_outbound(idx_include_outbound));
inst_vel_all_outbound_z=zscore(inst_vel_all_outbound(idx_include_outbound));

lm1 = fitlm(inst_vel_all_outbound_z,inst_theta_all_outbound_z);
ci1 = coefCI(lm1, ci_alpha); 

lm2 = fitlm(inst_acc_all_outbound_z,inst_theta_all_outbound_z);
ci2 = coefCI(lm2, ci_alpha); 

lm3 = fitlm(inst_steps_all_outbound_z,inst_theta_all_outbound_z);
ci3 = coefCI(lm3, ci_alpha); 

lm4 = fitlm(inst_vel_all_outbound_z,inst_steps_all_outbound_z);
ci4 = coefCI(lm4, ci_alpha); 

rp_outer_outbound(:,3)=[ci1(2); ci2(2); ci3(2); ci4(2)]; % lower 95%ci 
rp_outer_outbound(:,4)=[ci1(4); ci2(4); ci3(4); ci4(4)]; % upper 95%ci

%% T Junction - Outbound
clearvars -except savefig_scatter destdir rp_outbound ci_alpha rp_inbound rp_outer_inbound rp_outer_outbound

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};
rel_dist_forelimb_outbound=[];rel_dist_forlimb_inbound=[];rel_dist_hindlimb_outbound=[];rel_dist_hinlimb_inbound=[];

f=f1;

for a=[1,2,3,4,5]
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    
    ind_del=[]; ind_del_f1=[]; ind_del_f2=[];
    for i=1:size(ind,1)
        if ((isnan(f1(a).output{1, 1}(i).inst_acc_outbound (1)))||(isnan(f1(a).output{1, 1}(i).speed_outbound (1))))
             ind_del_f1=[ind_del_f1; i];
        end
    end

    for i=1:size(ind,1)
        if ((isnan(f2(a).output{1, 1}(i).inst_acc_outbound (1)))||(isnan(f2(a).output{1, 1}(i).speed_outbound (1))))
             ind_del_f2=[ind_del_f2; i];
        end
    end
    
    ind_del=unique([ind_del_f1;ind_del_f2], 'rows');
    
    f1(a).output{1, 1}(ind_del)=[];
    f2(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    count=size(ind_select,1);

    % keyboard 
    inst_vel_outbound_f1= cell2mat((arrayfun(@(x) x.inst_vel_outbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_acc_outbound_f1= cell2mat((arrayfun(@(x) x.inst_acc_outbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_freq_outbound_f1=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_steps_outbound_f1=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f1(a).output{1}(ind_select),'Un',0))');
    
    inst_vel_outbound_f2=cell2mat((arrayfun(@(x) x.inst_vel_outbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_acc_outbound_f2=cell2mat((arrayfun(@(x) x.inst_acc_outbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_freq_outbound_f2=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_steps_outbound_f2=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f2(a).output{1}(ind_select),'Un',0))');
    
    % plot data per animal 
    inst_vel_all{a}=([inst_vel_outbound_f1(:); inst_vel_outbound_f2(:)]);
    inst_acc_all{a}=([inst_acc_outbound_f1(:); inst_acc_outbound_f2(:)]);
    inst_theta_all{a}=([inst_freq_outbound_f1(:); inst_freq_outbound_f2(:)]);
    inst_steps_all{a}=([inst_steps_outbound_f1(:); inst_steps_outbound_f2(:)]);
    
    idx_include{a}=find((inst_vel_all{a}>0) & (inst_steps_all{a}>0) & (inst_theta_all{a}>0));
     
    counter(a)=count;
    clear count speed_all forelimb_steps_all theta_count_all run_duration_all distance_all idx ind_select inst_vel_outbound_f1 inst_acc_outbound_f1 inst_freq_outbound_f1 ...
        inst_steps_outbound_f1 inst_vel_outbound_f2 inst_acc_outbound_f2 inst_freq_outbound_f2 inst_steps_outbound_f2

end

idx_include_all=cell2mat(idx_include(:));
inst_steps_all=cell2mat(inst_steps_all(:));
inst_vel_all=cell2mat(inst_vel_all(:));
inst_theta_all=cell2mat(inst_theta_all(:));
inst_acc_all=cell2mat(inst_acc_all(:));

idx_include=find((inst_vel_all>vel_thresh_min) & (inst_steps_all>steps_thresh_min) & (inst_steps_all<steps_thresh_max) & (inst_theta_all>theta_thresh_min));
idx_include_vel=find((inst_vel_all>vel_thresh_min)&(inst_vel_all<vel_thresh_max));
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);

    %% SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'
    clear all; load('Figure1EF_EDFigure2.mat');

    % step-theta 
    close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
    scatter(tjunc_outbound.inst_steps_all_plot,tjunc_outbound.inst_theta_all_plot, 'w');
    binscatter(tjunc_outbound.inst_steps_all_plot,tjunc_outbound.inst_theta_all_plot, 100);lsline;
    [r3, p3]=corrcoef(tjunc_outbound.inst_steps_all_plot,tjunc_outbound.inst_theta_all_plot, 'rows', 'complete');
    xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2)]), 'Color', 'r', 'FontSize', 12);

%% or continue on from here if doing your own analysis :)
% total number of epochs  
total_epochs=sum(counter); cd(destdir);

close all; figure(1); hold on; xlim([0 120]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r1, p1]=corrcoef(inst_vel_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r1(2) p1(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstVelTheta_conOutbound_nose.fig');gcf;
end

close all; figure(2); hold on; xlim([-1 1]); ylim([5 12]); alpha(0.5); 
scatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r2, p2]=corrcoef(inst_acc_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Acceleration'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r2(2) p2(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstAccTheta_conOutbound_nose.fig');gcf;
end

close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r3, p3]=corrcoef(inst_steps_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstStepsTheta_conOutbound_nose.fig');gcf;
end

close all; figure(4); hold on; xlim([0 120]); ylim([0 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'w');
binscatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 100);lsline;
[r4, p4]=corrcoef(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Forelimb'); text(4,4.7,num2str([r4(2) p4(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1 
    savefig('InstVelSteps_conOutbound_nose.fig');gcf;
end

close all; figure(5); hold on; xlim([0 120]); alpha(0.5); %lsline;
hist(inst_vel_all(idx_include_vel),100);
xlabel('Instantaneous Velocity'); ylabel('Count'); 
if savefig_scatter==1 
    savefig('InstVelHist_conOutbound_nose.fig');gcf;
end

%% save variables for regression analysis
inst_steps_all_outbound=inst_steps_all;
inst_theta_all_outbound=inst_theta_all;
inst_acc_all_outbound=inst_acc_all;
inst_vel_all_outbound=inst_vel_all;
idx_include_outbound=idx_include;

rp_tjunc_outbound(:,1)=[r1(2); r2(2); r3(2); r4(2)];
rp_tjunc_outbound(:,2)=[p1(2) p2(2) p3(2) p4(2)];

% plot the scatter plot for inbound and outbound 
inst_steps_all_outbound_z=zscore(inst_steps_all_outbound(idx_include_outbound));
inst_theta_all_outbound_z=zscore(inst_theta_all_outbound(idx_include_outbound));
inst_acc_all_outbound_z=zscore(inst_acc_all_outbound(idx_include_outbound));
inst_vel_all_outbound_z=zscore(inst_vel_all_outbound(idx_include_outbound));

lm1 = fitlm(inst_vel_all_outbound_z,inst_theta_all_outbound_z);
ci1 = coefCI(lm1, ci_alpha); 

lm2 = fitlm(inst_acc_all_outbound_z,inst_theta_all_outbound_z);
ci2 = coefCI(lm2, ci_alpha); 

lm3 = fitlm(inst_steps_all_outbound_z,inst_theta_all_outbound_z);
ci3 = coefCI(lm3, ci_alpha); 

lm4 = fitlm(inst_vel_all_outbound_z,inst_steps_all_outbound_z);
ci4 = coefCI(lm4, ci_alpha); 

rp_tjunc_outbound(:,3)=[ci1(2); ci2(2); ci3(2); ci4(2)]; % lower 95%ci 
rp_tjunc_outbound(:,4)=[ci1(4); ci2(4); ci3(4); ci4(4)]; % upper 95%ci

%% T Junction - Inbound
clearvars -except f f1 f2 savefig_scatter destdir animals rp_outbound ci_alpha rp_inbound rp_outer_inbound rp_outer_outbound rp_tjunc_outbound

for a=[1,2,3,4,5]
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    
    ind_del=[]; ind_del_f1=[]; ind_del_f2=[];
    for i=1:size(ind,1)
        if ((isnan(f1(a).output{1, 1}(i).inst_acc_inbound(1)))||(isnan(f1(a).output{1, 1}(i).speed_inbound (1))))
             ind_del_f1=[ind_del_f1; i];
        end
    end

    for i=1:size(ind,1)
        if ((isnan(f2(a).output{1, 1}(i).inst_acc_inbound (1)))||(isnan(f2(a).output{1, 1}(i).speed_inbound (1))))
             ind_del_f2=[ind_del_f2; i];
        end
    end
    
    ind_del=unique([ind_del_f1;ind_del_f2], 'rows');
    
    f1(a).output{1, 1}(ind_del)=[];
    f2(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    count=size(ind_select,1);

    % keyboard 
    inst_vel_inbound_f1= cell2mat((arrayfun(@(x) x.inst_vel_inbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_acc_inbound_f1= cell2mat((arrayfun(@(x) x.inst_acc_inbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_freq_inbound_f1=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f1(a).output{1}(ind_select),'Un',0))');
    inst_steps_inbound_f1=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f1(a).output{1}(ind_select),'Un',0))');
    
    inst_vel_inbound_f2=cell2mat((arrayfun(@(x) x.inst_vel_inbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_acc_inbound_f2=cell2mat((arrayfun(@(x) x.inst_acc_inbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_freq_inbound_f2=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f2(a).output{1}(ind_select),'Un',0))');
    inst_steps_inbound_f2=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f2(a).output{1}(ind_select),'Un',0))');
    
    % plot data per animal 
    inst_vel_all{a}=([inst_vel_inbound_f1(:); inst_vel_inbound_f2(:)]);
    inst_acc_all{a}=([inst_acc_inbound_f1(:); inst_acc_inbound_f2(:)]);
    inst_theta_all{a}=([inst_freq_inbound_f1(:); inst_freq_inbound_f2(:)]);
    inst_steps_all{a}=([inst_steps_inbound_f1(:); inst_steps_inbound_f2(:)]);
    
    idx_include{a}=find((inst_vel_all{a}>0) & (inst_steps_all{a}>0) & (inst_theta_all{a}>0));
     
    counter(a)=count;
    clear count speed_all forelimb_steps_all theta_count_all run_duration_all distance_all idx ind_select inst_vel_inbound_f1 inst_acc_inbound_f1 inst_freq_inbound_f1 ...
        inst_steps_inbound_f1 inst_vel_inbound_f2 inst_acc_inbound_f2 inst_freq_inbound_f2 inst_steps_inbound_f2

end

idx_include_all=cell2mat(idx_include(:));
inst_steps_all=cell2mat(inst_steps_all(:));
inst_vel_all=cell2mat(inst_vel_all(:));
inst_theta_all=cell2mat(inst_theta_all(:));
inst_acc_all=cell2mat(inst_acc_all(:));

idx_include=find((inst_vel_all>vel_thresh_min) & (inst_steps_all>steps_thresh_min) & (inst_steps_all<steps_thresh_max) & (inst_theta_all>theta_thresh_min));
idx_include_vel=find((inst_vel_all>vel_thresh_min)&(inst_vel_all<vel_thresh_max));
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);
 
    %% SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'
    clear all; load('Figure1EF_EDFigure2.mat');

    % step-theta 
    close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
    scatter(tjunc_inbound.inst_steps_all_plot,tjunc_inbound.inst_theta_all_plot, 'w');
    binscatter(tjunc_inbound.inst_steps_all_plot,tjunc_inbound.inst_theta_all_plot, 100);lsline;
    [r3, p3]=corrcoef(tjunc_inbound.inst_steps_all_plot,tjunc_inbound.inst_theta_all_plot, 'rows', 'complete');
    xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2)]), 'Color', 'r', 'FontSize', 12);

%% or continue on from here if doing your own analysis :)
total_epochs=sum(counter); cd(destdir);

close all; figure(1); hold on; xlim([0 120]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_vel_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r1, p1]=corrcoef(inst_vel_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r1(2) p1(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstVelTheta_conInbound_nose.fig');gcf;
end

close all; figure(2); hold on; xlim([-1 1]); ylim([5 12]); alpha(0.5); 
scatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_acc_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r2, p2]=corrcoef(inst_acc_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Acceleration'); ylabel('Instantaneous Frequency Theta'); text(2,4.7,num2str([r2(2) p2(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstAccTheta_conInbound_nose.fig');gcf;
end

close all; figure(3); hold on; xlim([3 12]); ylim([5 12]); alpha(0.5); %lsline;
scatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 'w');
binscatter(inst_steps_all(idx_include),inst_theta_all(idx_include), 100);lsline;
[r3, p3]=corrcoef(inst_steps_all(idx_include),inst_theta_all(idx_include), 'rows', 'complete');
xlabel('Instantaneous Frequency Forelimb'); ylabel('Instantaneous Frequency Theta'); text(4,4.7,num2str([r3(2) p3(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1
    savefig('InstStepsTheta_conInbound_nose.fig');gcf;
end

close all; figure(4); hold on; xlim([0 120]); ylim([0 12]); alpha(0.5); %lsline;
scatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'w');
binscatter(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 100);lsline;
[r4, p4]=corrcoef(inst_vel_all(idx_include_vel),inst_steps_all(idx_include_vel), 'rows', 'complete');
xlabel('Instantaneous Velocity'); ylabel('Instantaneous Frequency Forelimb'); text(4,4.7,num2str([r4(2) p4(2) total_epochs]), 'Color', 'r', 'FontSize', 12);
if savefig_scatter==1 
    savefig('InstVelSteps_conInbound_nose.fig');gcf;
end

close all; figure(5); hold on; xlim([0 120]); alpha(0.5); %lsline;
hist(inst_vel_all(idx_include_vel),100);
xlabel('Instantaneous Velocity'); ylabel('Count'); 
if savefig_scatter==1 
    savefig('InstVelHist_conInbound_nose.fig');gcf;
end

%% save variables for regression analysis
inst_steps_all_inbound=inst_steps_all;
inst_theta_all_inbound=inst_theta_all;
inst_acc_all_inbound=inst_acc_all;
inst_vel_all_inbound=inst_vel_all;
idx_include_inbound=idx_include;

rp_tjunc_inbound(:,1)=[r1(2); r2(2); r3(2); r4(2)];
rp_tjunc_inbound(:,2)=[p1(2) p2(2) p3(2) p4(2)];

%% plot the scatter plot for inbound and inbound 
inst_steps_all_inbound_z=zscore(inst_steps_all_inbound(idx_include_inbound));
inst_theta_all_inbound_z=zscore(inst_theta_all_inbound(idx_include_inbound));
inst_acc_all_inbound_z=zscore(inst_acc_all_inbound(idx_include_inbound));
inst_vel_all_inbound_z=zscore(inst_vel_all_inbound(idx_include_inbound));

lm1 = fitlm(inst_vel_all_inbound_z,inst_theta_all_inbound_z);
ci1 = coefCI(lm1, ci_alpha); 

lm2 = fitlm(inst_acc_all_inbound_z,inst_theta_all_inbound_z);
ci2 = coefCI(lm2, ci_alpha); 

lm3 = fitlm(inst_steps_all_inbound_z,inst_theta_all_inbound_z);
ci3 = coefCI(lm3, ci_alpha); 

lm4 = fitlm(inst_vel_all_inbound_z,inst_steps_all_inbound_z);
ci4 = coefCI(lm4, ci_alpha); 

rp_tjunc_inbound(:,3)=[ci1(2); ci2(2); ci3(2); ci4(2)]; % lower 95%ci 
rp_tjunc_inbound(:,4)=[ci1(4); ci2(4); ci3(4); ci4(4)]; % upper 95%ci