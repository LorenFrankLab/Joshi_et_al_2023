%% dfs_mua_analysis_allcombination.m
% -> EDFigure45C_EDFigure7DE_mua_mod_score.mat
% -> dfa_ahead_behind_distance_muareldist.m
% -> dfa_ahead_behind_distance_muareldist_lintrack.m

%% Brief Description 
% Forelimb plants are collected during different parts of the track
% A corresponding plant triggered average is calculated 
% These plant times are then shuffled n times (5000) and an equivalent
% Number of plants as are in that epoch are collected m times (1000) - this
% is my shuffled distribution per epoch
% A modulation score is then computed per epoch for the actual plant times
% and for the shuffles

% modulation score: can be just peak-trough for mua 
% mod_score_stdev: is calculated as the summed square deviations from the mean
% stats on this data will look like a comparsion to the shuffled
% distributions. 

clear all;

%% initialization parameters 
mod_score_window=0.04; % corresponds to 80ms around the plant time 
t_pre=0.07; t_post=0.07; % the usal window around which we compute one cycle of MUA fluctuations ~125 ms is a theta cycle 
OffsetMin=-0.07; OffsetMax=0.07; % for shuffle analysis 
number_of_shuffles=5000; % per epoch 
savefigs=0; %saves the mua trace per wtrack portion for all aniamls 
shuffles=0; n_samps=1000; %does the shuffles and saves them for each portion of the track in results_dir

resdir='';

%% OUTBOUND MUA
results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};   

for a=[1,2,3,4,5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_outbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        forelimb_inbound_mua = f(a).output{1,1}(i).forelimb_inbound_mua;
        forelimb_outbound_mua = f(a).output{1,1}(i).forelimb_outbound_mua;
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
        
        forelimb_outbound_mua(:,2)=forelimb_outbound_mua;
        forelimb_inbound_mua(:,2)=forelimb_inbound_mua;

        limb_select_mua=forelimb_outbound_mua; number_steps=[number_steps; size(forelimb_outbound_mua,1)];
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis

            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_mua_shuff(:,2)=forelimb_outbound_mua_shuff;
                limb_select_mua_shuff=forelimb_outbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
            clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
    
% degree of theta modulation Frank et al 2001
% theta_mod=maximum value -minimum value
% saveworkspace 
if shuffles==1
    cd(results_dir); save('')
end

%% %% INBOUND MUA
clearvars -except f t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles mod_score_window

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};   

for a=[1,2,3,4,5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_inbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        forelimb_inbound_mua = f(a).output{1,1}(i).forelimb_inbound_mua;
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
        
        forelimb_inbound_mua(:,2)=forelimb_inbound_mua;

        limb_select_mua=forelimb_inbound_mua; number_steps=[number_steps; size(forelimb_inbound_mua,1)];
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;


%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_mua_shuff(:,2)=forelimb_inbound_mua_shuff;
                limb_select_mua_shuff=forelimb_inbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
           clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
    
% degree of theta modulation Frank et al 2001
% theta_mod=maximum value -minimum value
% saveworkspace 
if shuffles==1
    cd(results_dir); save('')
end

%% outbound outer arms mua
clearvars -except t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles mod_score_window

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};

f=f1;
for a=[1,2,3,4,5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
 
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_outbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
         
        forelimb_outbound_mua_outer1 = f1(a).output{1,1}(i).forelimb_outbound_mua;
        forelimb_outbound_mua_outer2 = f2(a).output{1,1}(i).forelimb_outbound_mua;
        
        forelimb_outbound_mua=[forelimb_outbound_mua_outer1;forelimb_outbound_mua_outer2];
        forelimb_outbound_mua(:,2)=forelimb_outbound_mua;
        limb_select_mua=forelimb_outbound_mua;
        
        number_steps=[number_steps; size(limb_select_mua,1)];
                
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;


%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_mua_shuff(:,2)=forelimb_outbound_mua_shuff;
                limb_select_mua_shuff=forelimb_outbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
            clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
            end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

if shuffles==1
    cd(results_dir); cd(results_dir); save('')
end

%% INBOUND outer arms mua
clearvars -except t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles mod_score_window

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};

f=f1;
for a=[1,2,3,4,5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
 
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_outbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
         
        forelimb_inbound_mua_outer1 = f1(a).output{1,1}(i).forelimb_inbound_mua;
        forelimb_inbound_mua_outer2 = f2(a).output{1,1}(i).forelimb_inbound_mua;
        
        forelimb_inbound_mua=[forelimb_inbound_mua_outer1;forelimb_inbound_mua_outer2];
        forelimb_inbound_mua(:,2)=forelimb_inbound_mua;
        limb_select_mua=forelimb_inbound_mua;
        
        number_steps=[number_steps; size(limb_select_mua,1)];
                
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;


%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_mua_shuff(:,2)=forelimb_inbound_mua_shuff;
                limb_select_mua_shuff=forelimb_inbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
                  clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
            end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

if shuffles==1
    cd(results_dir); cd(results_dir); save('')
end


%% INBOUND TJUNC mua
clearvars -except t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles mod_score_window

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};

f=f1;
for a=[1,2,3,4,5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
 
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_outbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
         
        forelimb_inbound_mua_outer1 = f1(a).output{1,1}(i).forelimb_inbound_mua;
        forelimb_inbound_mua_outer2 = f2(a).output{1,1}(i).forelimb_inbound_mua;
        
        forelimb_inbound_mua=[forelimb_inbound_mua_outer1;forelimb_inbound_mua_outer2];
        forelimb_inbound_mua(:,2)=forelimb_inbound_mua;
        limb_select_mua=forelimb_inbound_mua;
        
        number_steps=[number_steps; size(limb_select_mua,1)];
                
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;


%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_mua_shuff(:,2)=forelimb_inbound_mua_shuff;
                limb_select_mua_shuff=forelimb_inbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
            clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
            end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

if shuffles==1
    cd(results_dir); cd(results_dir); save('')
end

%% OUTBOUND TJUNC mua
clearvars -except animals f f1 f2 t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles mod_score_window
   f=f1;
   
for a=[1,2,3,4,5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
 
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_outbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
         
        forelimb_outbound_mua_outer1 = f1(a).output{1,1}(i).forelimb_outbound_mua;
        forelimb_outbound_mua_outer2 = f2(a).output{1,1}(i).forelimb_outbound_mua;
        
        forelimb_outbound_mua=[forelimb_outbound_mua_outer1;forelimb_outbound_mua_outer2];
        forelimb_outbound_mua(:,2)=forelimb_outbound_mua;
        limb_select_mua=forelimb_outbound_mua;
        
        number_steps=[number_steps; size(limb_select_mua,1)];
                
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;


%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_mua_shuff(:,2)=forelimb_outbound_mua_shuff;
                limb_select_mua_shuff=forelimb_outbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
             clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
            end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

if shuffles==1
    cd(results_dir); cd(results_dir); save('')
end

%% initialization parameters 
mod_score_window=0.04; % corresponds to 80ms around the plant time 
t_pre=0.07; t_post=0.07; % the usal window around which we compute one cycle of MUA fluctuations ~125 ms is a theta cycle 
OffsetMin=-0.07; OffsetMax=0.07; % for shuffle analysis 
number_of_shuffles=5000; % per epoch 
savefigs=1; %saves the mua trace per wtrack portion for all aniamls 
shuffles=1; n_samps=1000; %does the shuffles and saves them for each portion of the track in results_dir

resdir='';

% LINEARTRACK
results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 

animals={'Jaq', 'Roqui', 'Lotus'};   

for a=[1,2,3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=2:1:size(ind_select) % first epoch is exculded 
        counter=[counter; size(i,1)];
                
        if ~isnan(f(a).output{1,1}(i).forelimb_outbound_mua)
            
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        %keyboard
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);
        
        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        forelimb_inbound_mua_rl = f(a).output{1,1}(i).forelimb_inbound_mua;
        forelimb_outbound_mua_lr = f(a).output{1,1}(i).forelimb_outbound_mua;
        
        cdat_mua=imcont('data', mua', 'timestamp', muatimes');
        forelimb_outbound_mua=sort([forelimb_inbound_mua_rl;forelimb_outbound_mua_lr]);
        forelimb_outbound_mua(:,2)=forelimb_outbound_mua;

        limb_select_mua=forelimb_outbound_mua; number_steps=[number_steps; size(forelimb_outbound_mua,1)];
        trig_plot_mua=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua, 't_pre',t_pre, 't_post',t_post);
        lags_mua=trig_plot_mua.lags;

        % stdev mod score for all data 
        mydata_mua{a,i}=trig_plot_mua.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_mua{a,i}(:)-mean(mydata_mua{a,i})));

        % stdev mod score for all data rescaled 
        mydata_mua_rescale{a,i}=(mydata_mua{a,i})./(max(mydata_mua{a,i}));
        mod_score_perepoch_rescale_stdev{a,i}= sum(abs(mydata_mua_rescale{a,i}(:)-mean(mydata_mua_rescale{a,i})));

        % peak-trough mod score rescaled
        [pks_val,~]=findpeaks(mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        [trs_val,~]=findpeaks(-mydata_mua_rescale{a,i}, lags_mua,'MinPeakDistance', mod_score_window);
        mod_score_perepoch{a,i}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val

%        keyboard
%         plot(mydata_mua{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_mua.peritrig_data,1)
%             if sum(isnan(trig_plot_mua.peritrig_data(n_rows,:)))<2
%                 mydata_mua_nonnan{n_rows}=trig_plot_mua.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_mua_nonnan(:))), 'r'); hold on;


%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_mua_shuff=limb_select_mua(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_mua_shuff(:,2)=forelimb_outbound_mua_shuff;
                limb_select_mua_shuff=forelimb_outbound_mua_shuff;

                trig_plot_mua_shuff=periseg('cdat', contchans(cdat_mua, 'chans', 1), 'segs', limb_select_mua_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_mua_shuff{a,i,i_shuff}=trig_plot_mua_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_mua_shuff trig_plot_mua_shuff;
            end
            mydata_mua_shuff_all=cell2mat(squeeze(((mydata_mua_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_mua_shuff_all,1),size(limb_select_mua(:,1),1));
                    mydata_mua_nulldist{a,i,n_shuff}=nanmean(mydata_mua_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % calculate also for the rescaled version 
                    mydata_mua_nulldist_rescale{a,i,n_shuff}=(nanmean(mydata_mua_shuff_all(randIdcs,:))./(max(nanmean(mydata_mua_shuff_all(randIdcs,:)))));
                    mod_score_perepoch_rescale_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_mua_nulldist{a,i,n_shuff}(:)-mean(mydata_mua_nulldist{a,i,n_shuff})));
                    
                    % peak-trough mod score rescaled
                    [pks_val,~]=findpeaks(mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    [trs_val,~]=findpeaks(-mydata_mua_nulldist_rescale{a,i,n_shuff}, lags_mua,'MinPeakDistance', mod_score_window);
                    mod_score_perepoch_pershuff{a,i,n_shuff}=max(pks_val)-min(abs(trs_val)); clear trs_val pks_val
                    clear randIdcs
            end
        end
            clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all mod_score_window ...
            lags_mua mydata_mua_nulldist mod_score_perepoch_pershuff mod_score_perepoch mod_score_perepoch_rescale_pershuff_stdev ...
            mod_score_perepoch_rescale_stdev 
            end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
    
% degree of theta modulation Frank et al 2001
% theta_mod=maximum value -minimum value
% saveworkspace 
if shuffles==1
    cd(results_dir); save('')
end

%% CENTER ARM OUTBOUND | PER EPOCH ANALSYIS 
clear all;
results_dir=''; cd(results_dir);
load('');
p_95=0.05; p_75=0.25;

for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
    
    % calculate p values here for per animal 
    % comapre z scores 
    results=z_score(a,:); % collect the modulation scores per epoch per animal  
    ind= find(results~=0); % just remove the epochs where modulation score was not compted 
    z_score_animal=z_score(a,ind); clear ind % modulation index z_score for each epoch 

    results_shuff=z_score_shuff(a,:,:); % collect the modulation scores per epoch per animal  
    ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
    z_score_shuff_animal=z_score_shuff(a,ind); clear ind; % modulation index z_score for each epoch 
    [~,p_ttest_animal{a}]=ttest(z_score_animal,mean(z_score_shuff_animal)); % ttest: performs a t-test of the hypothesis that the data in X come from a distribution with mean M.  M must be a scalar.
    % keyboard
    clear z_score_shuff_animal results_shuff z_score_animal results
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff));

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;

bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density

sprintf('OUTBOUND CENTER, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.outbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.outbound=[z_score-mean(z_score_shuff)];

%% CENTER ARM INBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table p_95 p_75 mod_score_z
cd(results_dir)
load('');

for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
%     keyboard
    % calculate p values here for per animal 
    % comapre z scores 
    results=z_score(a,:); % collect the modulation scores per epoch per animal  
    ind= find(results~=0); % just remove the epochs where modulation score was not compted 
    z_score_animal=z_score(a,ind); clear ind % modulation index z_score for each epoch 

    results_shuff=z_score_shuff(a,:,:); % collect the modulation scores per epoch per animal  
    ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
    z_score_shuff_animal=z_score_shuff(a,ind); clear ind; % modulation index z_score for each epoch 
    [~,p_ttest_animal{a}]=ttest(z_score_animal,mean(z_score_shuff_animal)); % ttest: performs a t-test of the hypothesis that the data in X come from a distribution with mean M.  M must be a scalar.
    % keyboard
    clear z_score_shuff_animal results_shuff z_score_animal results
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); 

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;

bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density

sprintf('INBOUND CENTER, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.inbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.inbound=[z_score-mean(z_score_shuff)];

median(mod_score_z.outbound-mod_score_z.inbound) % 1.43

%% OUTER ARM OUTBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table p_95 p_75 mod_score_z
cd(results_dir)
load('')
% calculate a modulation_sig score that just does 1 for yes modulated above
% 95 % CI 2 if not
for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); 

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;

bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density
sprintf('OUTBOUND OUTER, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.outer_outbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.outer_outbound=[z_score-mean(z_score_shuff)];

%% OUTER ARM INBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table p_95 p_75 mod_score_z
cd(results_dir)
load('')

for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); 

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;

bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density
sprintf('INBOUND OUTER, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.outer_inbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.outer_inbound=[z_score-mean(z_score_shuff)];

%% %% TJUNC OUTBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table p_95 p_75 mod_score_z
cd(results_dir)
load('')

for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); 

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;

bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density
sprintf('OUTBOUND TJUNC, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.tjunc_outbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.tjunc_outbound=[z_score-mean(z_score_shuff)];

%% TJUNC INBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 

clearvars -except results_dir  results_table p_95 p_75 mod_score_z
cd(results_dir)
load('')

for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
    
    results=z_score(a,:); % collect the modulation scores per epoch per animal  
    ind= find(results~=0); % just remove the epochs where modulation score was not compted 
    z_score_animal=z_score(a,ind); clear ind % modulation index z_score for each epoch 

    results_shuff=z_score_shuff(a,:,:); % collect the modulation scores per epoch per animal  
    ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
    z_score_shuff_animal=z_score_shuff(a,ind); clear ind; % modulation index z_score for each epoch 
    [~,p_ttest_animal{a}]=ttest(z_score_animal,mean(z_score_shuff_animal));
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); 

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;
 
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density
sprintf('INBOUND TJUNC, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.tjunc_inbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.tjunc_inbound=[z_score-mean(z_score_shuff)];

%% LINEARTRACK| PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table p_95 p_75 mod_score_z
cd(results_dir)
load('');
p_95=0.05; p_75=0.25;

for a=1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
    for i=1:size(mod_score_perepoch_pershuff_stdev,2) % iterate through epochs 
            if ~isempty(cell2mat(mod_score_perepoch_pershuff_stdev(a,i,:))) % only collect mod score from those animal/epochs where data is valid 
                
                mod_score_perepoch_collect=mod_score_perepoch_stdev(a,i); % collect the actual modulation for mean of the epoch
                mod_score_perepoch_pershuff_stdev_collect=mod_score_perepoch_pershuff_stdev(a,i,1:end); % collect the distribution of modulation score per shuffle mean       

                x=squeeze(cell2mat(mod_score_perepoch_pershuff_stdev_collect(:,:,:))); % collecting the distribution of modulation score per epoch
                z_score(a,i)=(cell2mat(mod_score_perepoch_collect)-mean(x))./(std(x));
                z_score_shuff(a,i,:)=((x(:))-mean(x))./(std(x));
                
                CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % percentile method
                CI_95 = CIFcn(x,95); % calculate the 95% CI
                CI_75 = CIFcn(x,75); % calculate the 75% CI
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-mua modulation per epoch 
                
                  % OPTIONAL PLOT TO CHECK THE DISTRIBUTION
%                 hist(x,50); 
%                 arrayfun(@(x)xline(x,'-m','prctile'),CI);
%                 arrayfun(@(y)xline(y,'-r','Mean'),y);

                % collect 1 for outside 95% CI
                if (y<CI_95(1))||(y>CI_95(2))
                    modulation_sig95(a,i)=1;
                else
                    modulation_sig95(a,i)=2;
                end

                % collect 1 for outside 90% CI
                if (y<CI_75(1))||(y>CI_75(2))
                    modulation_sig75(a,i)=1;
                else
                    modulation_sig75(a,i)=2;
                end
                
                clear x y CI_75 CI_95 mod_score_perepoch_collect mod_score_perepoch_pershuff_stdev_collect
            end
    end
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); 

total_epochs=(size(find(modulation_sig75==1 | modulation_sig75==2),1));
epochs_75=size(find(modulation_sig75==1),1); %/total_epochs.*100;;
epochs_95=size(find(modulation_sig95==1),1); %/total_epochs.*100;;

bin_prob_density=binopdf(0:total_epochs,total_epochs,p_95);
p_val_95CI=bin_prob_density(epochs_95+1); clear bin_prob_density
bin_prob_density=binopdf(0:total_epochs,total_epochs,p_75);
p_val_75CI=bin_prob_density(epochs_75+1); clear bin_prob_density
sprintf('LINEARTRACK, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.lintrack=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.lintrack=[z_score-mean(z_score_shuff)];

%% SAVE AND PLOT 
% variables above saved in EDFigure45C_EDFigure7DE_mua_mod_score.mat
% run from here to plot 

clear all; load('EDFigure5C_EDFigure7DE_mua_mod_score.mat');
clearvars -except results_dir  results_table mod_score_z

% plot(results_table.outbound(1,4:6), 'X-'); hold on 
% plot(results_table.inbound(1,4:6), 'X-'); hold on 
% plot(results_table.outer_outbound(1,4:6), 'X-'); hold on 
% plot(results_table.outer_inbound(1,4:6), 'X-'); hold on 
% plot(results_table.tjunc_outbound(1,4:6), 'X-'); hold on 
% plot(results_table.tjunc_inbound(1,4:6), 'X-'); hold on 
% legend('OUTBOUND - CENTER', 'INBOUND - CENTER', 'OUTBOUND - OUTER', 'INBOUND - CENTER', 'OUTBOUND - TJUNC', 'INBOUND - TJUNC')

close all;
subplot(1,2,1);
plot(1+0.2*rand(1), results_table.outbound(1,4), 'X-', 'MarkerSize', 12, 'Color', 'g'); hold on 
plot(1+0.2*rand(1), results_table.inbound(1,4), 'X-', 'MarkerSize', 12, 'Color', 'r'); hold on 
plot(1+0.2*rand(1), results_table.outer_outbound(1,4), 'X-', 'MarkerSize', 12,  'Color', 'k'); hold on 
plot(1+0.2*rand(1), results_table.outer_inbound(1,4), 'X-', 'MarkerSize', 12,  'Color', 'm'); hold on 
plot(1+0.2*rand(1), results_table.tjunc_outbound(1,4), 'X-', 'MarkerSize', 12,  'Color', 'y'); hold on 
plot(1+0.2*rand(1), results_table.tjunc_inbound(1,4), 'X-', 'MarkerSize', 12,  'Color', 'b'); hold on 
plot(1+0.2*rand(1), results_table.lintrack(1,4), 'X-', 'MarkerSize', 12,  'Color', 'c'); hold on 

plot(2+0.2*rand(1), results_table.outbound(1,5), 'X-', 'MarkerSize', 12, 'Color', 'g'); hold on 
plot(2+0.2*rand(1), results_table.inbound(1,5), 'X-', 'MarkerSize', 12,  'Color', 'r'); hold on 
plot(2+0.2*rand(1), results_table.outer_outbound(1,5), 'X-', 'MarkerSize', 12,  'Color', 'k'); hold on 
plot(2+0.2*rand(1), results_table.outer_inbound(1,5), 'X-', 'MarkerSize', 12,  'Color', 'm'); hold on 
plot(2+0.2*rand(1), results_table.tjunc_outbound(1,5), 'X-', 'MarkerSize', 12,  'Color', 'y'); hold on 
plot(2+0.2*rand(1), results_table.tjunc_inbound(1,5), 'X-', 'MarkerSize', 12,  'Color', 'b'); hold on 
plot(2+0.2*rand(1), results_table.lintrack(1,5), 'X-', 'MarkerSize', 12,  'Color', 'c'); hold on 

plot(3+0.2*rand(1), results_table.outbound(1,6), 'X-', 'MarkerSize', 12, 'Color', 'g'); hold on 
plot(3+0.2*rand(1), results_table.inbound(1,6), 'X-', 'MarkerSize', 12,  'Color', 'r'); hold on 
plot(3+0.2*rand(1), results_table.outer_outbound(1,6), 'X-', 'MarkerSize', 12,  'Color', 'k'); hold on 
plot(3+0.2*rand(1), results_table.outer_inbound(1,6), 'X-', 'MarkerSize', 12,  'Color', 'm'); hold on 
plot(3+0.2*rand(1), results_table.tjunc_outbound(1,6), 'X-', 'MarkerSize', 12,  'Color', 'y'); hold on 
plot(3+0.2*rand(1), results_table.tjunc_inbound(1,6), 'X-', 'MarkerSize', 12,  'Color', 'b'); hold on 
plot(3+0.2*rand(1), results_table.lintrack(1,6), 'X-', 'MarkerSize', 12,  'Color', 'c'); hold on 

plot([0 4], [0.05 0.05], '-k');
legend('OUTBOUND - CENTER', 'INBOUND - CENTER', 'OUTBOUND - OUTER', 'INBOUND - OUTER', 'OUTBOUND - TJUNC', 'INBOUND - TJUNC', 'LINTRACK')
xticklabels({'' 'Bin. Prob. 95% CI' 'Bin. Prob. 75% CI' 'T Test Shuff'})
xlabel('Statistical Test')
ylabel('P-Value of Comparison')

subplot(1,2,2);
bar(1, -log(results_table.outbound(1,6))); hold on 
bar(2, -log(results_table.inbound(1,6))); hold on 
bar(3, -log(results_table.outer_outbound(1,6))); hold on 
bar(4, -log(results_table.outer_inbound(1,6))); hold on 
bar(5, -log(results_table.tjunc_outbound(1,6))); hold on 
bar(6, -log(results_table.tjunc_inbound(1,6))); hold on 
bar(7, -log(results_table.lintrack(1,6))); hold on 
plot([0 8], [-log(0.05) -log(0.05)], '-k');
ylabel('-log [P-Value of Comparison]');
legend('OUTBOUND - CENTER', 'INBOUND - CENTER', 'OUTBOUND - OUTER', 'INBOUND - OUTER', 'OUTBOUND - TJUNC', 'INBOUND - TJUNC', 'LINTRACK');
title('T Test - Task Phase')

%keyboard
% save figure in a figures 
resdir='';
cd(resdir)
%savefigs(['Forelimb_center_wtrack_mua_modulation_stats']); close(gcf)

mod_score_z_all=[mod_score_z.outbound;mod_score_z.inbound; mod_score_z.outer_outbound; mod_score_z.outer_inbound; mod_score_z.tjunc_outbound; mod_score_z.tjunc_inbound; mod_score_z.lintrack];
id=[ones(size(mod_score_z.outbound));1+ones(size(mod_score_z.inbound)); 2+ones(size(mod_score_z.outer_outbound)); 3+ones(size(mod_score_z.outer_inbound)); 4+ones(size(mod_score_z.tjunc_outbound)); 5+(ones(size(mod_score_z.tjunc_inbound))); 6+(ones(size(mod_score_z.lintrack)))];
[r,p,stats]=kruskalwallis(mod_score_z_all,id);

close all;
ax=figure(1); hold on 
boxplot(mod_score_z_all, id, 'Notch', 'on');
all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);

xticklabels({'OUTBOUND - CENTER', 'INBOUND - CENTER', 'OUTBOUND - OUTER', 'INBOUND - OUTER', 'OUTBOUND - TJUNC', 'INBOUND - TJUNC', 'LINTRACK'});
ylabel('Modulation Score [z]');

%keyboard
% save figure in a figures 
resdir='';
cd(resdir)
%savefigs(['']); close(gcf)

multcompare(stats)
%keyboard
% save figure in a figures 
resdir='';
cd(resdir)
%savefigs(['']); close(gcf)

%% INBOUND OUTBOUND COMPARISON 

[r,p]=signrank(mod_score_z.outbound',mod_score_z.inbound');

% average difference
mean(mod_score_z.outbound-mod_score_z.inbound); %1.6165
mean(mod_score_z.outbound) % 1.4820
mean(mod_score_z.inbound) % -0.1345

%% OUTBOUND INBOUND-OUTER arm 
[r,p]=signrank(mod_score_z.outbound',mod_score_z.outer_inbound');
% signrank p=0.37, rat1, rat2, rat3 rat4, rat5 p>0.05