%% dfs_plotpower_outbound_inbound_peranimal.m 
% -> dfa_steptheta_wtrack.m
% -> smooth_dlc_stepcycles_wtrack.m
% -> Figure1ACD.mat

%% set filters
% clear all; close all; clc; 

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
%     timefilter{1} = {'aj_get2dstate_dlc', '($mobility == 1)','mobility_velocity',[04 200], 'mobility_buffer',0.25, 'min_durations',1, 'vel2choose', 'nose'};
%     f(a) = createfilter('animal',animals{a},'epochs',epochfilter,'excludetime', timefilter, 'eegtetrodes',tetfilter, 'iterator', iterator);
%     f(a) = setfilterfunction(f(a), 'dfa_steptheta_wtrack', {'eeg', 'eeggnd', 'theta' ,'posdlc'},'animal',animals{a},'out_min',30, 'out_max',100, 'a',a);
%     f(a)= runfilter(f(a));
% end

% fname ='';
% save(fname,'f', '-v7.3')

%% OUTBOUND Forelimb RIGHT 
% set parameters

clear all; close all; clc;

% load results file
load(''); clearvars -except f t_window minpeakheight noise_thresh;
animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};

for a = [1 2 3 4 5]
     
    count=0;
        
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)|(ind(:,1)==6));
    %ind_select=((ind(:,1)==1))
    
    cdat_forepawR{a}=(arrayfun(@(x) x.forepawR_filt_6_8',f(a).output{1}(ind_select),'Un',0))';
    theta_data_all{a}=(arrayfun(@(x) x.theta_filt_5_11',f(a).output{1}(ind_select),'Un',0))';
    theta_time_all{a}=(arrayfun(@(x) x.theta_time',f(a).output{1}(ind_select),'Un',0))';
    run_periods_all{a}=(arrayfun(@(x) x.outbound_periods',f(a).output{1}(ind_select),'Un',0))';
    
    % calculate the correlations coeffiecient for each epoch 
    for e=1:size(cdat_forepawR{1,a},1)
        
        cdat_forepawR_filt=cdat_forepawR{1,a}{e,1};
        theta_data=theta_data_all{1,a}{e,1};
        theta_time=theta_time_all{1,a}{e,1};
        run_periods=run_periods_all{1,a}{e,1}';
        run_periods=run_periods(1:end-1,:);  
        
        if isstruct(cdat_forepawR_filt)
         
        count=count+1;
        cdat_forepawR_filt.data(isnan(cdat_forepawR_filt.data))=0; % contpsd doesnt deal wwell with nan 
        subplot(1,2,1); hold on;
        [p1,f1]=contpsd(cdat_forepawR_filt, 'window_t',t_window , 'segs', run_periods);
        p1_norm=p1/max(p1); hold on 
        plot(f1,p1_norm, 'color', [(0.5/a) (0.5/a) (0.5/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        subplot(1,2,2); hold on;
        cdat_lfp=imcont('data', theta_data','timestamp', theta_time'); hold on;
        [p2,f2]=contpsd(cdat_lfp, 'window_t', 1, 'segs', run_periods);
        p2_norm=p2/max(p2);
        plot(f2,p2_norm, 'color',[(0.6/a) (0.6/a) (0.6/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        psd_res_steps{a,e}=[p1_norm]';
        psd_res_theta{a,e}=[p2_norm]';
        
        if ~isempty(f1)
            [~,idx_limbpeak]= findpeaks(p1_norm, 'MinPeakHeight', minpeakheight);
            psd_peak_forelimbR{a,e}=max((f1(idx_limbpeak)));
        else
            psd_peak_forelimbR{a,e}=[];
        end

        mydata_steps_f1=f1;
        mydata_theta_f2=f2;
        
        clear p1 p1_norm run_periods run_select p2 f2 p2_norm theta_time theta_data cdat_lfp idx_limbpeak cdat_forepawR_filt
        else 
            
        clear run_periods theta_time theta_data cdat_lfp cdat_forepawR_filt
        end
        counter(a)=count;
    end   
    
    forelimbR_outbound_peak{a}=(psd_peak_forelimbR(a,:));
    clear ind* cdat_*
end
forelimbR_outbound_peak_all=cell2mat(psd_peak_forelimbR(~cellfun('isempty',psd_peak_forelimbR)));
mydata_steps=psd_res_steps(:);
mydata_theta=psd_res_theta(:);

mydata_steps = mydata_steps(~cellfun('isempty', mydata_steps'));
mydata_theta = mydata_theta(~cellfun('isempty', mydata_theta'));

mydata_steps = cell2mat(mydata_steps);
mydata_theta = cell2mat(mydata_theta);

    %% SAVE AND PLOT 
    % variables above saved in Figure1ACD.mat
    % run from here to plot 
    close all; clear all; load('Figure1ACD.mat')

    mydata_steps_f1=mydata_steps_f1_out_fr;
    mydata_steps=mydata_steps_out_fr;
    
    plot(mydata_steps_f1,mean(mydata_steps,1), 'DisplayName', 'Stepping Rhythm', 'Linewidth', 2, 'color', [0.597656250000000	0.195312500000000	0.796875000000000]);hold on % Identify Grp1 by giving it a tag number
    [~,ind]=findpeaks(mean(mydata_steps,1));
    % step_peaks=mydata_steps_f1(ind); step_peaks(step_peaks>noise_thresh)=[];

    % plot(mydata_theta_f2,mean(mydata_theta,1), 'DisplayName', 'Hippocampal EEG' , 'Linewidth', 2); % Identify Grp2 by giving it a tag number
    % [val,ind]=findpeaks(mean(mydata_theta,1));
    % theta_peaks=mydata_theta_f2(ind); theta_peaks(theta_peaks>noise_thresh)=[];

    x_steps = [mydata_steps_f1; flip(mydata_steps_f1); mydata_steps_f1(1)];
    s = plot(mydata_steps_f1,mean(mydata_steps,1), 'color', [0.5 0.5 0.5], 'Linewidth', 2);
    y1 = mean(mydata_steps,1)-std(mydata_steps,[],1);
    y2 = mean(mydata_steps,1)+std(mydata_steps,[],1);
    patch_in_steps = [ y1, fliplr(y2), y1(1)];
    hold on; p = patch(x_steps,patch_in_steps',repmat(size(mydata_steps,1),size(mydata_steps,2)*2+1,1),'facealpha', 0.25,'linestyle','none');
    clear y1 y2 

    % x_theta = [mydata_theta_f2; flip(mydata_theta_f2); mydata_theta_f2(1)]
    % s = plot(mydata_theta_f2,mean(mydata_theta,1));
    % y1 = mean(mydata_theta,1)-std(mydata_theta,[],1);
    % y2 = mean(mydata_theta,1)+std(mydata_theta,[],1);
    % patch_in_theta = [ y1, fliplr(y2), y1(1)];
    % hold on; p = patch(x_theta',patch_in_theta',repmat(size(mydata_theta,1),size(mydata_theta,2)*2+1,1),'facealpha', 0.5,'linestyle','none' )
    % clear y1 y2

    xlim([0 15]); ylim([0 1.2]);
    title(['W Track - OUT fr'], 'FontName', 'Arial', 'FontSize', 18);
    xlabel('Frequency Hz', 'FontName', 'Arial', 'FontSize', 18);
    ylabel('Power (a.u.)', 'FontName', 'Arial', 'FontSize', 18);
    legend({'Stepping Rhythm'});
    %'Hippocampal EEG'})

%savefig('StepFreq_forelimbR_wtrack_30_100_outbound.fig');

%% OUTBOUND Forelimb LEFT 
clearvars -except f t_window minpeakheight noise_thresh forelimbR_outbound_peak*
animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};

for a = [1 2 3 4 5]
     
    count=0;
        
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)|(ind(:,1)==6));
    %ind_select=((ind(:,1)==1))
    
    cdat_forepawL{a}=(arrayfun(@(x) x.forepawL_filt_6_8',f(a).output{1}(ind_select),'Un',0))';
    theta_data_all{a}=(arrayfun(@(x) x.theta_filt_5_11',f(a).output{1}(ind_select),'Un',0))';
    theta_time_all{a}=(arrayfun(@(x) x.theta_time',f(a).output{1}(ind_select),'Un',0))';
    run_periods_all{a}=(arrayfun(@(x) x.outbound_periods',f(a).output{1}(ind_select),'Un',0))';

    for e=1:length(cdat_forepawL{1,a})
        
        cdat_forepawL_filt=cdat_forepawL{1,a}{e,1};
        theta_data=theta_data_all{1,a}{e,1};
        theta_time=theta_time_all{1,a}{e,1};
        run_periods=run_periods_all{1,a}{e,1}';
        run_periods=run_periods(1:end-1,:);  
        
        if isstruct(cdat_forepawL_filt)
         
        count=count+1;
        cdat_forepawL_filt.data(isnan(cdat_forepawL_filt.data))=0; % contpsd doesnt deal wwell with nan 
%        keyboard
        subplot(1,2,1); hold on;
        [p1,f1]=contpsd(cdat_forepawL_filt, 'window_t',t_window , 'segs', run_periods);
        p1_norm=p1/max(p1); hold on 
        plot(f1,p1_norm, 'color', [(0.5/a) (0.5/a) (0.5/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        subplot(1,2,2); hold on;
        cdat_lfp=imcont('data', theta_data','timestamp', theta_time'); hold on;
        [p2,f2]=contpsd(cdat_lfp, 'window_t', 1, 'segs', run_periods);
        p2_norm=p2/max(p2);
        plot(f2,p2_norm, 'color',[(0.6/a) (0.6/a) (0.6/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        psd_res_steps{a,e}=[p1_norm]';
        psd_res_theta{a,e}=[p2_norm]';
 %      keyboard
 
         if ~isempty(f1)
            [~,idx_limbpeak]= findpeaks(p1_norm, 'MinPeakHeight', minpeakheight);
            psd_peak_forelimbL{a,e}=max((f1(idx_limbpeak)));
        else
            psd_peak_forelimbL{a,e}=[];
        end
        
        mydata_steps_f1=f1;
        mydata_theta_f2=f2;
        
                clear p1 p1_norm run_periods run_select p2 f2 p2_norm theta_time theta_data cdat_lfp idx_limbpeak cdat_forepawR_filt
        else 
            
        clear run_periods theta_time theta_data cdat_lfp
        end
        counter(a)=count;
    end   
    
    forelimbL_outbound_peak{a}=(psd_peak_forelimbL(a,:));
    clear ind* cdat_* 
end
forelimbL_outbound_peak_all=cell2mat(psd_peak_forelimbL(~cellfun('isempty',psd_peak_forelimbL)));
mydata_steps=psd_res_steps(:);
mydata_theta=psd_res_theta(:);

mydata_steps = mydata_steps(~cellfun('isempty', mydata_steps'));
mydata_theta = mydata_theta(~cellfun('isempty', mydata_theta'));

mydata_steps = cell2mat(mydata_steps);
mydata_theta = cell2mat(mydata_theta);

    %% SAVE AND PLOT 
    % variables above saved in Figure1ACD.mat
    % run from here to plot 
    close all; clear all; load('Figure1ACD.mat')

    mydata_steps_f1=mydata_steps_f1_out_fl;
    mydata_steps=mydata_steps_out_fl;
    
    plot(mydata_steps_f1,mean(mydata_steps,1), 'DisplayName', 'Stepping Rhythm', 'Linewidth', 2, 'color', [0.597656250000000	0.195312500000000	0.796875000000000]);hold on % Identify Grp1 by giving it a tag number
    [~,ind]=findpeaks(mean(mydata_steps,1));
    % step_peaks=mydata_steps_f1(ind); step_peaks(step_peaks>noise_thresh)=[];

    % plot(mydata_theta_f2,mean(mydata_theta,1), 'DisplayName', 'Hippocampal EEG' , 'Linewidth', 2); % Identify Grp2 by giving it a tag number
    % [val,ind]=findpeaks(mean(mydata_theta,1));
    % theta_peaks=mydata_theta_f2(ind); theta_peaks(theta_peaks>noise_thresh)=[];

    x_steps = [mydata_steps_f1; flip(mydata_steps_f1); mydata_steps_f1(1)];
    s = plot(mydata_steps_f1,mean(mydata_steps,1), 'color', [0.5 0.5 0.5], 'Linewidth', 2);
    y1 = mean(mydata_steps,1)-std(mydata_steps,[],1);
    y2 = mean(mydata_steps,1)+std(mydata_steps,[],1);
    patch_in_steps = [ y1, fliplr(y2), y1(1)];
    hold on; p = patch(x_steps,patch_in_steps',repmat(size(mydata_steps,1),size(mydata_steps,2)*2+1,1),'facealpha', 0.25,'linestyle','none');
    clear y1 y2 

    % x_theta = [mydata_theta_f2; flip(mydata_theta_f2); mydata_theta_f2(1)]
    % s = plot(mydata_theta_f2,mean(mydata_theta,1));
    % y1 = mean(mydata_theta,1)-std(mydata_theta,[],1);
    % y2 = mean(mydata_theta,1)+std(mydata_theta,[],1);
    % patch_in_theta = [ y1, fliplr(y2), y1(1)];
    % hold on; p = patch(x_theta',patch_in_theta',repmat(size(mydata_theta,1),size(mydata_theta,2)*2+1,1),'facealpha', 0.5,'linestyle','none' )
    % clear y1 y2

    xlim([0 15]); ylim([0 1.2]);
    title(['W Track - OUT fl'], 'FontName', 'Arial', 'FontSize', 18);
    xlabel('Frequency Hz', 'FontName', 'Arial', 'FontSize', 18);
    ylabel('Power (a.u.)', 'FontName', 'Arial', 'FontSize', 18);
    legend({'Stepping Rhythm'});

%savefig('StepFreq_forelimbL_wtrack_30_100_outbound.fig');

%% INBOUND FORELIMB LEFT 
clearvars -except f t_window minpeakheight noise_thresh forelimbR_outbound_peak* forelimbL_outbound_peak*
animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};

for a = [1 2 3 4 5]
     
    count=0;
        
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)|(ind(:,1)==6));
    %ind_select=((ind(:,1)==1))
    
    cdat_forepawL{a}=(arrayfun(@(x) x.forepawL_filt_6_8',f(a).output{1}(ind_select),'Un',0))';
    theta_data_all{a}=(arrayfun(@(x) x.theta_filt_5_11',f(a).output{1}(ind_select),'Un',0))';
    theta_time_all{a}=(arrayfun(@(x) x.theta_time',f(a).output{1}(ind_select),'Un',0))';
    run_periods_all{a}=(arrayfun(@(x) x.inbound_periods',f(a).output{1}(ind_select),'Un',0))';

    for e=1:length(cdat_forepawL{1,a})
        
        cdat_forepawL_filt=cdat_forepawL{1,a}{e,1};
        theta_data=theta_data_all{1,a}{e,1};
        theta_time=theta_time_all{1,a}{e,1};
        run_periods=run_periods_all{1,a}{e,1}';
        run_periods=run_periods(1:end-1,:);  
        
        if isstruct(cdat_forepawL_filt)
         
        count=count+1;
        cdat_forepawL_filt.data(isnan(cdat_forepawL_filt.data))=0; % contpsd doesnt deal wwell with nan 
%        keyboard
        subplot(1,2,1); hold on;
        [p1,f1]=contpsd(cdat_forepawL_filt, 'window_t',t_window , 'segs', run_periods);
        p1_norm=p1/max(p1); hold on 
        plot(f1,p1_norm, 'color', [(0.5/a) (0.5/a) (0.5/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        subplot(1,2,2); hold on;
        cdat_lfp=imcont('data', theta_data','timestamp', theta_time'); hold on;
        [p2,f2]=contpsd(cdat_lfp, 'window_t', 1, 'segs', run_periods);
        p2_norm=p2/max(p2);
        plot(f2,p2_norm, 'color',[(0.6/a) (0.6/a) (0.6/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        psd_res_steps{a,e}=[p1_norm]';
        psd_res_theta{a,e}=[p2_norm]';
 %      keyboard
 
        if ~isempty(f1)
            [~,idx_limbpeak]= findpeaks(p1_norm, 'MinPeakHeight', minpeakheight);
            psd_peak_forelimbL{a,e}=max((f1(idx_limbpeak)));
        else
            psd_peak_forelimbL{a,e}=[];
        end
        
        mydata_steps_f1=f1;
        mydata_theta_f2=f2;
        
                clear p1 p1_norm run_periods run_select p2 f2 p2_norm theta_time theta_data cdat_lfp idx_limbpeak cdat_forepawR_filt
        else 
            
        clear run_periods theta_time theta_data cdat_lfp
        end
        counter(a)=count;
    end   
    
    forelimbL_inbound_peak{a}=(psd_peak_forelimbL(a,:));
    clear ind* cdat_* 
end
forelimbL_inbound_peak_all=cell2mat(psd_peak_forelimbL(~cellfun('isempty',psd_peak_forelimbL)));
mydata_steps=psd_res_steps(:);
mydata_theta=psd_res_theta(:);

mydata_steps = mydata_steps(~cellfun('isempty', mydata_steps'));
mydata_theta = mydata_theta(~cellfun('isempty', mydata_theta'));

mydata_steps = cell2mat(mydata_steps);
mydata_theta = cell2mat(mydata_theta);

    %% SAVE AND PLOT 
    % variables above saved in Figure1ACD.mat
    % run from here to plot 
    close all; clear all; load('Figure1ACD.mat')

    mydata_steps_f1=mydata_steps_f1_in_fl;
    mydata_steps=mydata_steps_in_fl;

    plot(mydata_steps_f1,mean(mydata_steps,1), 'DisplayName', 'Stepping Rhythm', 'Linewidth', 2, 'color', [0.597656250000000	0.195312500000000	0.796875000000000]);hold on % Identify Grp1 by giving it a tag number
    [~,ind]=findpeaks(mean(mydata_steps,1));
    % step_peaks=mydata_steps_f1(ind); step_peaks(step_peaks>noise_thresh)=[];

    % plot(mydata_theta_f2,mean(mydata_theta,1), 'DisplayName', 'Hippocampal EEG' , 'Linewidth', 2); % Identify Grp2 by giving it a tag number
    % [val,ind]=findpeaks(mean(mydata_theta,1));
    % theta_peaks=mydata_theta_f2(ind); theta_peaks(theta_peaks>noise_thresh)=[];

    x_steps = [mydata_steps_f1; flip(mydata_steps_f1); mydata_steps_f1(1)];
    s = plot(mydata_steps_f1,mean(mydata_steps,1), 'color', [0.5 0.5 0.5], 'Linewidth', 2);
    y1 = mean(mydata_steps,1)-std(mydata_steps,[],1);
    y2 = mean(mydata_steps,1)+std(mydata_steps,[],1);
    patch_in_steps = [ y1, fliplr(y2), y1(1)];
    hold on; p = patch(x_steps,patch_in_steps',repmat(size(mydata_steps,1),size(mydata_steps,2)*2+1,1),'facealpha', 0.25,'linestyle','none');
    clear y1 y2 

    % x_theta = [mydata_theta_f2; flip(mydata_theta_f2); mydata_theta_f2(1)]
    % s = plot(mydata_theta_f2,mean(mydata_theta,1));
    % y1 = mean(mydata_theta,1)-std(mydata_theta,[],1);
    % y2 = mean(mydata_theta,1)+std(mydata_theta,[],1);
    % patch_in_theta = [ y1, fliplr(y2), y1(1)];
    % hold on; p = patch(x_theta',patch_in_theta',repmat(size(mydata_theta,1),size(mydata_theta,2)*2+1,1),'facealpha', 0.5,'linestyle','none' )
    % clear y1 y2

    xlim([0 15]); ylim([0 1.2]);
    title(['W Track - IN fl)'], 'FontName', 'Arial', 'FontSize', 18);
    xlabel('Frequency Hz', 'FontName', 'Arial', 'FontSize', 18);
    ylabel('Power (a.u.)', 'FontName', 'Arial', 'FontSize', 18);
    legend({'Stepping Rhythm'});
    %'Hippocampal EEG'})

%savefig('StepFreq_forelimbL_wtrack_30_100_inbound.fig');

%% INBOUND Forelimb RIGHT
clearvars -except f t_window minpeakheight noise_thresh forelimbR_outbound_peak* forelimbL_outbound_peak* forelimbR_inbound_peak* forelimbL_inbound_peak*
animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};

for a = [1 2 3 4 5]
     
    count=0;
        
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)|(ind(:,1)==6));
    %ind_select=((ind(:,1)==1))
    
    cdat_forepawR{a}=(arrayfun(@(x) x.forepawR_filt_6_8',f(a).output{1}(ind_select),'Un',0))';
    theta_data_all{a}=(arrayfun(@(x) x.theta_filt_5_11',f(a).output{1}(ind_select),'Un',0))';
    theta_time_all{a}=(arrayfun(@(x) x.theta_time',f(a).output{1}(ind_select),'Un',0))';
    run_periods_all{a}=(arrayfun(@(x) x.inbound_periods',f(a).output{1}(ind_select),'Un',0))';

    for e=1:length(cdat_forepawR{1,a})
        
        cdat_forepawR_filt=cdat_forepawR{1,a}{e,1};
        theta_data=theta_data_all{1,a}{e,1};
        theta_time=theta_time_all{1,a}{e,1};
        run_periods=run_periods_all{1,a}{e,1}';
        run_periods=run_periods(1:end-1,:);  
        
        if isstruct(cdat_forepawR_filt)
         
        count=count+1;
        cdat_forepawR_filt.data(isnan(cdat_forepawR_filt.data))=0; % contpsd doesnt deal wwell with nan 
%        keyboard
        subplot(1,2,1); hold on;
        [p1,f1]=contpsd(cdat_forepawR_filt, 'window_t',t_window , 'segs', run_periods);
        p1_norm=p1/max(p1); hold on 
        plot(f1,p1_norm, 'color', [(0.5/a) (0.5/a) (0.5/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        subplot(1,2,2); hold on;
        cdat_lfp=imcont('data', theta_data','timestamp', theta_time'); hold on;
        [p2,f2]=contpsd(cdat_lfp, 'window_t', 1, 'segs', run_periods);
        p2_norm=p2/max(p2);
        plot(f2,p2_norm, 'color',[(0.6/a) (0.6/a) (0.6/a)], 'Linewidth', 0.1);
        xlim([0 15]);
        
        psd_res_steps{a,e}=[p1_norm]';
        psd_res_theta{a,e}=[p2_norm]';
 %      keyboard
 
         if ~isempty(f1)
            [~,idx_limbpeak]= findpeaks(p1_norm, 'MinPeakHeight', minpeakheight);
            psd_peak_forelimbR{a,e}=max((f1(idx_limbpeak)));
        else
            psd_peak_forelimbR{a,e}=[];
        end
        
        mydata_steps_f1=f1;
        mydata_theta_f2=f2;
        
        clear p1 p1_norm run_periods run_select p2 f2 p2_norm theta_time theta_data cdat_lfp
        else 
            
        clear run_periods theta_time theta_data cdat_lfp
        end
        counter(a)=count;
    end
    
    forelimbR_inbound_peak{a}=(psd_peak_forelimbR(a,:));
    clear ind* cdat_* 
end
forelimbR_inbound_peak_all=cell2mat(psd_peak_forelimbR(~cellfun('isempty',psd_peak_forelimbR)));
mydata_steps=psd_res_steps(:);
mydata_theta=psd_res_theta(:);

mydata_steps = mydata_steps(~cellfun('isempty', mydata_steps'));
mydata_theta = mydata_theta(~cellfun('isempty', mydata_theta'));

mydata_steps = cell2mat(mydata_steps);
mydata_theta = cell2mat(mydata_theta);

    %% SAVE AND PLOT 
    % variables above saved in Figure1ACD.mat
    % run from here to plot 
    close all; clear all; load('Figure1ACD.mat')

    mydata_steps_f1=mydata_steps_f1_in_fr;
    mydata_steps=mydata_steps_in_fr;
    
    plot(mydata_steps_f1,mean(mydata_steps,1), 'DisplayName', 'Stepping Rhythm', 'Linewidth', 2, 'color', [0.597656250000000	0.195312500000000	0.796875000000000]);hold on % Identify Grp1 by giving it a tag number
    [val,ind]=findpeaks(mean(mydata_steps,1));
    % step_peaks=mydata_steps_f1(ind); step_peaks(step_peaks>noise_thresh)=[];

    % plot(mydata_theta_f2,mean(mydata_theta,1), 'DisplayName', 'Hippocampal EEG' , 'Linewidth', 2); % Identify Grp2 by giving it a tag number
    % [val,ind]=findpeaks(mean(mydata_theta,1));
    % theta_peaks=mydata_theta_f2(ind); theta_peaks(theta_peaks>noise_thresh)=[];

    x_steps = [mydata_steps_f1; flip(mydata_steps_f1); mydata_steps_f1(1)];
    s = plot(mydata_steps_f1,mean(mydata_steps,1), 'color', [0.5 0.5 0.5], 'Linewidth', 2);
    y1 = mean(mydata_steps,1)-std(mydata_steps,[],1);
    y2 = mean(mydata_steps,1)+std(mydata_steps,[],1);
    patch_in_steps = [ y1, fliplr(y2), y1(1)];
    hold on; p = patch(x_steps,patch_in_steps',repmat(size(mydata_steps,1),size(mydata_steps,2)*2+1,1),'facealpha', 0.25,'linestyle','none');
    clear y1 y2 

    % x_theta = [mydata_theta_f2; flip(mydata_theta_f2); mydata_theta_f2(1)]
    % s = plot(mydata_theta_f2,mean(mydata_theta,1));
    % y1 = mean(mydata_theta,1)-std(mydata_theta,[],1);
    % y2 = mean(mydata_theta,1)+std(mydata_theta,[],1);
    % patch_in_theta = [ y1, fliplr(y2), y1(1)];
    % hold on; p = patch(x_theta',patch_in_theta',repmat(size(mydata_theta,1),size(mydata_theta,2)*2+1,1),'facealpha', 0.5,'linestyle','none' )
    % clear y1 y2

    xlim([0 15]); ylim([0 1.2]);
    title(['W Track - IN fr'], 'FontName', 'Arial', 'FontSize', 18);
    xlabel('Frequency Hz', 'FontName', 'Arial', 'FontSize', 18);
    ylabel('Power (a.u.)', 'FontName', 'Arial', 'FontSize', 18);
    legend({'Stepping Rhythm'});
    %'Hippocampal EEG'})

%savefig('StepFreq_forelimbR_wtrack_30_100_inbound.fig');

%% STATS 

close all
for a=[1,2,3,4,5]
  forelimbR_outbound_peak_plot=cell2mat(forelimbR_outbound_peak{a}(~cellfun('isempty',forelimbR_outbound_peak{a})));
  forelimbR_inbound_peak_plot=cell2mat(forelimbR_inbound_peak{a}(~cellfun('isempty',forelimbR_inbound_peak{a})));
  
  forelimbL_outbound_peak_plot=cell2mat(forelimbL_outbound_peak{a}(~cellfun('isempty',forelimbL_outbound_peak{a})));
  forelimbL_inbound_peak_plot=cell2mat(forelimbL_inbound_peak{a}(~cellfun('isempty',forelimbL_inbound_peak{a})));
  
  data_all=[(forelimbR_outbound_peak_plot+forelimbL_outbound_peak_plot)'; (forelimbR_inbound_peak_plot+forelimbL_inbound_peak_plot)'];
  id_all=[ones(size((forelimbR_outbound_peak_plot+forelimbL_outbound_peak_plot)))'; 1+ones(size((forelimbR_inbound_peak_plot+forelimbL_inbound_peak_plot)))'];

  % PLOT means per animal 
  data=[mean((forelimbR_outbound_peak_plot+forelimbL_outbound_peak_plot))'; mean((forelimbR_inbound_peak_plot+forelimbL_inbound_peak_plot))'];
  id=[1;2];

%   %keyboard
%   ax=figure(1);
%   plot([id],[data], 'k'); hold on;
%   [p1,~]=kruskalwallis(data_all,id_all); % collect the significance per animal 
%   %keyboard
%   p_all{a}=p1;
%   figure(1); text(1.5,a,['Rat', num2str(a), ' p: ' num2str(p_all{a})]);
%   %r_all{a}=r1;
  
  clear data* id* p1 r1 forelimbR_outbound_peak_plot forelimbR_inbound_peak_plot forelimbL_outbound_peak_plot forelimbL_inbound_peak_plot
end

ax=figure(1);
% empty circles for r value per epoch 
data_all=[(forelimbR_outbound_peak_all+forelimbL_outbound_peak_all); (forelimbR_inbound_peak_all+forelimbL_inbound_peak_all)]; %aniamls are mostly trotting 
id_all=[ones(size((forelimbR_outbound_peak_all+forelimbL_outbound_peak_all))); (1+ ones(size((forelimbR_inbound_peak_all+forelimbL_inbound_peak_all))))]; %aniamls are mostly trotting 
% scatter([id_all],[data_all],'k','jitter','on'); hold on 
boxplot([data_all], [id_all],'Notch', 'on'); hold on 

ylim([0 10]);  xlim([0.5 2.5]);
% legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5')
xlabel('Portion on WTrack');
ylabel('Peak Frequency [Hz]');
labels=({'C-OUT' 'C-IN'});
xticklabels(labels)

all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);