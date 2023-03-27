%% dfs_reldist_analysis_allcombination.m
% -> dfa_ahead_behind_distance_muareldist.m
% -> dfa_ahead_behind_distance_muareldist_lintrack.m
% -> aj_get2dstate_dlc.m
% -> Figure23_plant_triggered_avgs.mat
% -> Figure2D_EDFigure7BC_decode_to_animal_dist.mat

%% CENTER ARM 
% clear all; close all; clc; 

% %animal filters | enter animal names to anlyse 
% animals={'Jaq','Roqui', 'Lotus', 'Monty', 'Peanut'};

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
%     timefilter{1} = {'aj_get2dstate_dlc', '($mobility == 1)','mobility_velocity',[10 200], 'mobility_buffer',0, 'min_durations',0.75, 'vel2choose', 'nose'};
% 
%     f(a) = createfilter('animal',animals{a},'epochs',epochfilter,'excludetime', timefilter, 'eegtetrodes',tetfilter, 'iterator', iterator);
%     f(a) = setfilterfunction(f(a), 'dfa_ahead_behind_distance_muareldist', {'eeg', 'eeggnd', 'theta' ,'posdlc'},'animal',animals{a},'out_min',60, 'out_max',100, 'a',a, 'ci_thresh', 50, 'ci_mean', 50, 'reldist_mean', 50, 'rel_dist_cutoff',10);
%     f(a)= runfilter(f(a));
% end
% 
% fname ='/';
% save(fname,'f', '-v7.3')

%% LINEARTRACK 
% % % set filters
% clear all; close all; clc; 
% 
% %animal filters | enter animal names to anlyse 
% animals={'Jaq','Roqui', 'Lotus'};
% %animals={'Roqui'};
% 
% for a = 1:length(animals)
% 
%     % epoch filters | select wtrack epochs
%     epochfilter{1} = 'isequal($type, ''run'') && isequal($environment, ''lineartrack'')';
% 
%     % iterators | these are goimg to define your base unit of analysis, such as per cell or per epoch
%     iterator= 'epocheeganal'; 
% 
%     % goes into tetinfo and searches for the correct tetrodes to be included 
%     tetfilter = 'isequal($area,''ca1Rref'') '; %isequal($area,''ca1Lref'') |
% 
%     % timefilter is a filter that excludes a selected time window
%     timefilter{1} = {'aj_get2dstate_dlc', '($mobility == 1)','mobility_velocity',[10 200], 'mobility_buffer',0.25 'min_durations',0.75, 'vel2choose', 'nose'};
% 
%     f(a) = createfilter('animal',animals{a},'epochs',epochfilter,'excludetime', timefilter, 'eegtetrodes',tetfilter, 'iterator', iterator);
%     f(a) = setfilterfunction(f(a), 'dfa_ahead_behind_distance_muareldist_lintrack', {'eeg', 'eeggnd', 'theta' ,'posdlc'},'animal',animals{a},'out_min',40, 'out_max',80, 'a',a, 'ci_thresh', 50, 'ci_mean', 50, 'reldist_mean', 50, 'rel_dist_cutoff',10);
%     f(a)= runfilter(f(a));
% end
% 
% fname ='';
% save(fname,'f', '-v7.3')

%% initialization parameters 
clear all;
mod_score_window=0.04; % corresponds to 80ms around the plant time 
t_pre=0.07; t_post=0.07; % the usal window around which we compute one cycle of MUA fluctuations ~125 ms is a theta cycle 
OffsetMin=-0.07; OffsetMax=0.07; % for shuffle analysis 
number_of_shuffles=5000; % per epoch 
savefigs=0; %saves the mua trace per wtrack portion for all aniamls 
shuffles=0; n_samps=1000; %does the shuffles and saves them for each portion of the track in results_dir

resdir='';

%% OUTBOUND RELDIST
results_dir='/'; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
%clearvars -except f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};   

for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx % these epochs have good quality decoding

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
    end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[2]
    ind= cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat2_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
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
    cd(results_dir); save(''); 
end

%% INBOUND RELDIST 
clearvars -except f t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};   
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx 

        if size(f(a).output{1,1}(i).forelimb_inbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_reldist_shuff(:,2)=forelimb_inbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_inbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat2_idx

        if size(f(a).output{1,1}(i).forelimb_inbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_reldist_shuff(:,2)=forelimb_inbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_inbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx

        if size(f(a).output{1,1}(i).forelimb_inbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_reldist_shuff(:,2)=forelimb_inbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_inbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx

        if size(f(a).output{1,1}(i).forelimb_inbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_inbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_inbound_reldist_shuff(:,2)=forelimb_inbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_inbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

% saveworkspace
if shuffles ==1
    cd(results_dir); save(''); 
end 

%% INBOUND outer arms RELDIST
clearvars -except t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles

results_dir='/'; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};

f=f1;
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx 

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig; 
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for  i=rat2_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
                    clear randIdcs
            end
        end
             clearvars -except animals results_dir lags_reldist mydata_reldist_nulldist results_dir  a t_pre t_post shuffles number_of_shuffles n_samps OffsetMin OffsetMax f f1 f2 counter number_steps number_steps_correct number_steps_incorrect ...
            mod_score_mua_perepoch_stdev_correct mod_score_mua_perepoch_stdev_incorrect mod_score_mua_perepoch_stdev lags_mua_correct ...
            mod_score_mua_perepoch_pershuff_stdev mod_score_mua_perepoch_pershuff_stdev_incorrect_control ...
            mydata_mua_correct mydata_mua_incorrect mydata_mua ...
            mod_score_theta_perepoch_stdev_correct mod_score_theta_perepoch_stdev_incorrect mod_score_theta_perepoch_stdev lags_theta_correct ...
            mod_score_theta_perepoch_pershuff_stdev mod_score_theta_perepoch_pershuff_stdev_incorrect_control ...
            mydata_theta_correct mydata_theta_incorrect mydata_theta ... 
            mod_score_perepoch_stdev_correct mod_score_perepoch_stdev_incorrect mod_score_perepoch_stdev lags_reldist_correct ...
            mod_score_perepoch_pershuff_stdev mod_score_perepoch_pershuff_stdev_incorrect_control ...
            mydata_reldist_correct mydata_reldist_incorrect mydata_reldist ...
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
% saveworkspace
if shuffles==1
    cd(results_dir); cd(results_dir); save('');  
end

%% OUTBOUND outer arms RELDIST
clearvars -except animals f f1 f2  t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles
   f=f1;
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx 

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for  i=rat2_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        
        end
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        
        end
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
% saveworkspace
if shuffles==1
    cd(results_dir); save('');  
end

%% INBOUND Choice point T RELDIST
clearvars -except t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
f1=f; clear f;

results_filename=[''];
load([results_filename,'.mat']); 
f2=f; clear f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};
rel_dist_forelimb_outbound_mua=[];rel_dist_forlimb_inbound=[];rel_dist_hindlimb_outbound=[];rel_dist_hinlimb_inbound=[];

f=f1;
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx 

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for  i=rat2_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_inbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_inbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_inbound_longsweep;
        
        forelimb_inbound_reldist=[forelimb_inbound_reldist_outer1;forelimb_inbound_reldist_outer2];
        forelimb_inbound_reldist(:,2)=forelimb_inbound_reldist;
        limb_select_reldist=forelimb_inbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

% saveworkspace
if shuffles==1
    cd(results_dir); save('');  
end

%% OUTBOUND Choice point T
clearvars -except animals f f1 f2 t_pre t_post n_samps savefigs resdir OffsetMin OffsetMax number_of_shuffles results_dir shuffles
   f=f1;
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx 

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for  i=rat2_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        
        forelimb_outbound_reldist_outer1 = f1(a).output{1,1}(i).forelimb_outbound_longsweep;
        forelimb_outbound_reldist_outer2 = f2(a).output{1,1}(i).forelimb_outbound_longsweep;
        
        forelimb_outbound_reldist=[forelimb_outbound_reldist_outer1;forelimb_outbound_reldist_outer2];
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

% saveworkspace
if shuffles==1
    cd(results_dir); save('')
end

%% LINEARTRACK
clear all;
% initialization parameters 
mod_score_window=0.04; % corresponds to 80ms around the plant time 
t_pre=0.07; t_post=0.07; % the usal window around which we compute one cycle of MUA fluctuations ~125 ms is a theta cycle 
OffsetMin=-0.07; OffsetMax=0.07; % for shuffle analysis 
number_of_shuffles=5000; % per epoch 
savefigs=1; %saves the mua trace per wtrack portion for all aniamls 
shuffles=1; n_samps=1000; %does the shuffles and saves them for each portion of the track in results_dir

resdir='';

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
%clearvars -except f;

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};   

for a=[1,2,3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=2:1:size(ind_select) % these epochs have good quality decoding

        if size(f(a).output{1,1}(i).forelimb_outbound_reldist,2)>1      
         counter=[counter; size(i,1)];
                
        animal=animals{1,a};
        base_dir='/opt/stelmo/abhilasha/animals';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        
        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);
        [muatimes, mua] = getMUAtrace(marks{d}{e}(tets),eeg{d}{e}{tets(1)});
        
        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist_lr = f(a).output{1,1}(i).forelimb_outbound_longsweep;
  
        forelimb_outbound_reldist_rl = f(a).output{1,1}(i).forelimb_inbound_longsweep;
        forelimb_outbound_reldist=sort([forelimb_outbound_reldist_lr;forelimb_outbound_reldist_rl]);
        forelimb_outbound_reldist(:,2)=forelimb_outbound_reldist;
        limb_select_reldist=forelimb_outbound_reldist;
        number_steps=[number_steps; size(limb_select_reldist,1)];
        
        trig_plot_reldist=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
        lags_reldist=trig_plot_reldist.lags;
        mydata_reldist{a,i}=trig_plot_reldist.peritrig;
        mod_score_perepoch_stdev{a,i}= sum(abs(mydata_reldist{a,i}(:)-mean(mydata_reldist{a,i})));
        
%        keyboard
%         plot(mydata_reldist{a,i}, 'k'); hold on;
%         for n_rows=1:size(trig_plot_reldist.peritrig_data,1)
%             if sum(isnan(trig_plot_reldist.peritrig_data(n_rows,:)))<2
%                 mydata_reldist_nonnan{n_rows}=trig_plot_reldist.peritrig_data(n_rows,:);
%             end
%         end
%         plot(nanmean(cell2mat(mydata_reldist_nonnan(:))), 'r'); hold on;
        % just for fun, remove those row with nans

%        keyboard
        if shuffles == 1
            %% Iteration for shuffling analysis
            %RandomOffsets=OffsetMin:(1/(number_of_shuffles*10)):OffsetMax;
            RandomOffsets=OffsetMin+(OffsetMax-OffsetMin).*rand(number_of_shuffles,1); % in seconds 
            
            % this is my shuffled distribution
            for i_shuff=1:size(RandomOffsets,1)
                %number_of_shuffles
                forelimb_outbound_reldist_shuff=limb_select_reldist(:,1)+RandomOffsets(i_shuff);
                forelimb_outbound_reldist_shuff(:,2)=forelimb_outbound_reldist_shuff;
                limb_select_reldist_shuff=forelimb_outbound_reldist_shuff;

                trig_plot_reldist_shuff=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist_shuff, 't_pre',t_pre, 't_post',t_post);  
                mydata_reldist_shuff{a,i,i_shuff}=trig_plot_reldist_shuff.peritrig_data; % doing peritrig_data instead of peritrig ensures that all shuffled plants get to go in the massive shuffle 
                
            close all; clear limb_select_reldist_shuff trig_plot_reldist_shuff;
            end
            mydata_reldist_shuff_all=cell2mat(squeeze(((mydata_reldist_shuff(a,i,:)))));
           
            % now create a null distribution based on the shuffles i.e from
            % the shuffled relationships, sample the same number of plants
            % that you have in the epoch a 1000 times 
            
            for n_shuff=1:n_samps
                    randIdcs=randperm(size(mydata_reldist_shuff_all,1),size(limb_select_reldist(:,1),1));
                    mydata_reldist_nulldist{a,i,n_shuff}=nanmean(mydata_reldist_shuff_all(randIdcs,:));
                    mod_score_perepoch_pershuff_stdev{a,i,n_shuff}= sum(abs(mydata_reldist_nulldist{a,i,n_shuff}(:)-mean(mydata_reldist_nulldist{a,i,n_shuff})));
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
            number_of_epochs_all number_steps_all number_steps_correct_all number_steps_incorrect_all
        end       
    end
        number_of_epochs_all(a)=sum(counter);         
        number_steps_all(a)=sum(number_steps);
        clear counter number_steps
end

% saveworkspace
if shuffles==1
    cd(results_dir); save('')
end

keyboard
%% CENTER ARM OUTBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clear all;
results_dir=''; cd(results_dir);
load('');
p_95=0.05; p_75=0.25;

% calculate a modulation_sig score that just does 1 for yes modulated above
% 95 % CI 2 if not
for a=[1,2,3,5]
    % 1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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
%     keyboard
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
load('');

% calculate a modulation_sig score that just does 1 for yes modulated above
% 95 % CI 2 if not
for a=[1,2,3,5]
    %1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 

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
%     keyboard
    clear z_score_shuff_animal results_shuff z_score_animal results
end

% comapre z scores 
results=z_score(:); % collect the modulation scores per epoch 
ind= find(results~=0); % just remove the epochs where modulation score was not compted 
z_score=z_score(ind); clear ind % modulation index z_score for each epoch 

results_shuff=z_score_shuff(:); % collect the modulation scores per epoch 
ind = find(results_shuff~=0); % just remove the epochs where modulation score was not compted 
z_score_shuff=z_score_shuff(ind); clear ind; % modulation index z_score for each epoch 
[~,p_ttest]=ttest(z_score,mean(z_score_shuff)); % | ttest: performs a t-test of the hypothesis that the data in X come from a distribution with mean M.  M must be a scalar.

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

%% OUTER ARM OUTBOUND | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table p_95 p_75 mod_score_z
load('wtrack_allanimals_reldist_shuff_outer_outbound_longsweep_sep29.mat')
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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
load('wtrack_allanimals_reldist_shuff_outer_inbound_longsweep_sep29.mat')
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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

% calculates p values per animal| individual aniamls were ns 
for a=[1,2,3,5]
    %1:size(mod_score_perepoch_pershuff_stdev,1) % iterate through animals 
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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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
%     keyboard
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
sprintf('INBOUND TJUNC, TOTAL EPOCHS (%d) | 95 | 75 | Binomial Prob. | Ttest with shuffled distribution = %d %d %d %d %d', total_epochs, round(epochs_95), round(epochs_75), p_val_95CI, p_val_75CI, p_ttest)

% best to save in a table for a final plot 
results_table.tjunc_inbound=[total_epochs, epochs_95, epochs_75, p_val_95CI, p_val_75CI, p_ttest];
mod_score_z.tjunc_inbound=[z_score-mean(z_score_shuff)];

%% LINEARTRACK  | PLOT ALL ANIMALS individually and together | PER EPOCH ANALSYIS 
clearvars -except results_dir  results_table mod_score_z

load('');
p_95=0.05; p_75=0.25;

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
                y=abs(cell2mat(mod_score_perepoch_collect)); % plant-reldist modulation per epoch 
                
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
% variables saved in Figure2D_EDFigure7BC_decode_to_animal_dist.mat
% RUN from here to plot 

clear all; load('Figure2D_EDFigure7BC_decode_to_animal_dist');
% clearvars -except results_dir  results_table mod_score_z

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
%savefig(['Forelimb_center_wtrack_reldist_modulation_stats']); close(gcf)

mod_score_z_all=[mod_score_z.outbound;mod_score_z.inbound; mod_score_z.outer_outbound; mod_score_z.outer_inbound; mod_score_z.tjunc_outbound; mod_score_z.tjunc_inbound; mod_score_z.lintrack];
id=[ones(size(mod_score_z.outbound));1+ones(size(mod_score_z.inbound)); 2+ones(size(mod_score_z.outer_outbound)); 3+ones(size(mod_score_z.outer_inbound)); 4+ones(size(mod_score_z.tjunc_outbound)); 5+(ones(size(mod_score_z.tjunc_inbound))); 6+(ones(size(mod_score_z.lintrack)))];
[~,~,stats]=kruskalwallis(mod_score_z_all,id);

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
%savefig(['']); close(gcf)

multcompare(stats)
%keyboard
% %save figure in a figures 
resdir='';
cd(resdir)
%savefig(['']); close(gcf)

%% INBOUND OUTBOUND COMPARISON 
mod_score_z_all=[mod_score_z.outbound;mod_score_z.inbound];
[r,p]=signrank(mod_score_z.outbound',mod_score_z.inbound');

% average diff 
median(mod_score_z.outbound-mod_score_z.inbound) % 1.64

%% SAVE AND PLOT 
% variables saved in Figure23_plant_triggered_avgs.mat
% OUTBOUND INBOUND PLOT 
close all; clear all; load('Figure23_plant_triggered_avgs.mat')

subplot(2,2,1)
hist(z_score_shuff_outbound,50); 
arrayfun(@(z_score_shuff_outbound)xline(z_score_shuff_outbound,'-m','95 prctile'),CI_95_outbound);
arrayfun(@(y)xline(y_outbound,'-r','Mean'),y_outbound);
xlabel('Modulation Score - Shuffles');
ylabel('#');
title('OUTBOUND CENTER');

subplot(2,2,3)
hist(z_score_shuff_inbound,50); 
arrayfun(@(z_score_shuff_inbound)xline(z_score_shuff_inbound,'-m','95 prctile'),CI_95_outbound);
arrayfun(@(y)xline(y_inbound,'-r','Mean'),y_inbound);
xlabel('Modulation Score - Shuffles');
ylabel('#');
title('INBOUND CENTER');

subplot(2,2,2)
plot(lags_reldist,y_mean_outbound, 'Color', [17/255 17/255 17/255], 'linewidth', 1);
hold on
plot(lags_reldist,ci_upper95_outbound, 'Color', [100/255 100/255 100/255], 'linewidth', 1);
plot(lags_reldist,ci_lower95_outbound, 'Color', [100/255 100/255 100/255], 'linewidth', 1);
hold on 
plot(lags_reldist,(mydata_reldist_outbound));
xlim([lags_reldist(1) lags_reldist(end)]);
xlabel('Time since Plants [s]');
ylabel('Relative Distance [cm]');
title('OUTBOUND CENTER');
ylim([-6 6])

subplot(2,2,4)
plot(lags_reldist,y_mean_inbound, 'Color', [17/255 17/255 17/255], 'linewidth', 1);
hold on
plot(lags_reldist,ci_upper95_inbound, 'Color', [100/255 100/255 100/255], 'linewidth', 1);
plot(lags_reldist,ci_lower95_inbound, 'Color', [100/255 100/255 100/255], 'linewidth', 1);
hold on 
plot(lags_reldist,(mydata_reldist_inbound));
xlim([lags_reldist(1) lags_reldist(end)]);
xlabel('Time since Plants [s]');
ylabel('Relative Distance [cm]');
title('INBOUND CENTER');
ylim([-6 6])