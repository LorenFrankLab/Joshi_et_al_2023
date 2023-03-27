%% dfs_gait_motifs_analysis.m
% -> EDFigure3.mat
% -> dfa_steptheta_wtrack_gait.m

% This script will plot the cross correlations between rigth and left
% forelimb steps during every run and help us estimate if these are out of
% phase in our data. 

%% set filters
% clear all; close all; clc; 
% 
% %animal filters | enter animal names to anlyse 
% animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};
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
%     tetfilter = 'isequal($area,''ca1Rref'') '; %isequal($area,''ca1Lref'') |
% 
%     % timefilter is a filter that excludes a selected time window
%     timefilter{1} = {'aj_get2dstate_dlc', '($mobility == 1)','mobility_velocity',[10 200], 'mobility_buffer',0.25, 'min_durations',2.5, 'vel2choose', 'nose'};
%     f(a) = createfilter('animal',animals{a},'epochs',epochfilter,'excludetime', timefilter, 'eegtetrodes',tetfilter, 'iterator', iterator);
%     f(a) = setfilterfunction(f(a), 'dfa_steptheta_wtrack_gait', {'eeg', 'eeggnd', 'theta' ,'posdlc'},'animal',animals{a},'out_min',30, 'out_max',400, 'a',a);
%     f(a)= runfilter(f(a));
% end
% 
% fname ='';
% save(fname,'f', '-v7.3')

%% SAVE AND PLOT 
% variables saved in EDFigure3.mat
% run from here to plot

% EDFigure1BC
clear all; close all; load('EDFigure3.mat');

%animal filters | enter animal names to anlyse 
animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};

maxlag_steps=40;
est_framerate=125;
lags_steps=-maxlag_steps/est_framerate:1/est_framerate:maxlag_steps/est_framerate;

for a = [1 2 3 4 5]
    count=0; animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)|(ind(:,1)==6));
   
    xcorr_forepawR_forepawL_all{a}=(arrayfun(@(x) x.xcorr_forepawR_forepawL,f(a).output{1}(ind_select),'Un',0))';
    speed_val_all{a}=(arrayfun(@(x) x.speed_val,f(a).output{1}(ind_select),'Un',0))';
    
    for i=1:size(xcorr_forepawR_forepawL_all{a},1)
        if size(xcorr_forepawR_forepawL_all{a}{i},1)<2
            xcorr_forepawR_forepawL_all{a}{i}=[];
            speed_val_all{a}{i}=[];
            idx(i)=0;
        else
            idx(i)=1;
        end
    end
    
    % per animal computation
    xcorr_to_plot=cell2mat(xcorr_forepawR_forepawL_all{a}(idx>0));
    speed_val_to_plot=(speed_val_all{a}(idx>0));
    speed_val_to_plot=cell2mat(speed_val_to_plot(:)');
    clim=[(mean(min(xcorr_to_plot))) (mean(max(xcorr_to_plot)))];
    
    b=1:size(xcorr_to_plot,1);
    %[~,b]=sort(speed_val_to_plot(speed_val_to_plot>0));
   
    subplot(2,5,a)
    imagesc(lags_steps,[1 b],xcorr_to_plot(b,:), clim);
    xlabel('lags secs'); ylabel('run periods'); title(['forepawR-forepawL', animal]);
    vline(0)

    % Metric - for each run, calculate the min and max index and value 
    for i=1:size(xcorr_to_plot)
        [min_val(i), min_ind(i)]=min(xcorr_to_plot(i,:));
        [max_val(i), max_ind(i)]=max(xcorr_to_plot(i,:));
    end

    subplot(2,5,a+5)
    histogram(min_val);
    [hist1,c1]=histogram(min_val);
    hold on 
    histogram(max_val);
    [hist2,c2]=histogram(max_val);
    % calculate the overlaps between two normal distributions
    % [overlap2] = calc_overlap_twonormal(s1,s2,mu1,mu2,xstart,xend,xinterval)
    
    bothHistograms = [hist1'; hist2'];
    minCounts = min(bothHistograms, [], 1);
    maxCounts = max(bothHistograms, [], 1);

%     pd1=fitdist(hist1','Normal');
%     pd2=fitdist(hist2','Normal');
%     [overlap] = calc_overlap_twonormal(pd1.sigma,pd2.sigma,pd1.mu,pd2.mu,0.01);
%      
    ratios(a) = (minCounts ./ maxCounts).*100;
    clear ind* xcorr_to_plot animal idx b speed_val_to_plot clim minCounts maxCounts bothHistograms min_val max_val
end

% EDFigure1A
close all; 
bin=0.01; tmax=0.3;

for a = [1,2,3,4,5]
    count=0; animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)|(ind(:,1)==6));
   
    xcorr_forepawR_forepawL_all{a}=(arrayfun(@(x) x.xcorr_forepawR_forepawL_plants,f(a).output{1}(ind_select),'Un',0))';
    speed_val_all{a}=(arrayfun(@(x) x.speed_val,f(a).output{1}(ind_select),'Un',0))';
    
end

for a=[1 2 3 4 5]
    animal= animals(a);
        
    % per animal computation
    for i=1:size(xcorr_forepawR_forepawL_all{a},1)
        left_right_forelimb_collect_count{i}=xcorr_forepawR_forepawL_all{a}{i}.c1vsc2;
        left_right_forelimb_collect_time{i}=xcorr_forepawR_forepawL_all{a}{i}.time;
    end
    
    xcorr_to_plot=cell2mat(left_right_forelimb_collect_count(:));
    time_to_plot=(left_right_forelimb_collect_time{1});

    subplot(1,5,a)
    plot(time_to_plot,mean(xcorr_to_plot(:,:))./ max(mean(xcorr_to_plot(:,:)))); hold on 
    xlim([-0.3 0.3]);
    xlabel('Lags secs'); ylabel('Normalized Count'); title(['Left-Right Plant xcorr', animal]);
    vline(0)
%     
%     subplot(2,5,a+5)
%     hist(mean(xcorr_to_plot(:,:)), 20); hold on 
%     xlabel('Lags secs'); ylabel('Count'); title(['forepawR-forepawL', animal]);
%     vline(0)
%     keyboard
    clear ind* xcorr_to_plot animal idx b speed_val_to_plot clim minCounts maxCounts bothHistograms min_val max_val time_to_plot left_right_forelimb_collect_count left_right_forelimb_collect_time
end
