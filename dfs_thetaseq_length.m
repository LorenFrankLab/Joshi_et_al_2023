%% dfs_thetaseq_length.m
% -> EDFigure6.mat

%% 1. collect peak reldist median per epochs 
clear all;
results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};   
load('decode_quality_idx.mat')

% median peak reldist ini outbound and inbound 
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat1_idx % these epochs have good quality decoding
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        d=index(1); e=index(2);

        peak_outbound_reldist{a,i}= median(double(f(a).output{1,1}(i).peaks_select_outbound_reldist));
        peak_inbound_reldist{a,i}= median(double(f(a).output{1,1}(i).peaks_select_inbound_reldist));
    end
end
for a=[2]
    ind= cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat2_idx
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        d=index(1); e=index(2);

        peak_outbound_reldist{a,i} = median(double(f(a).output{1,1}(i).peaks_select_outbound_reldist));
        peak_inbound_reldist{a,i} = median(double(f(a).output{1,1}(i).peaks_select_inbound_reldist));
    end
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat3_idx
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        d=index(1); e=index(2);

        peak_outbound_reldist{a,i} = median(double(f(a).output{1,1}(i).peaks_select_outbound_reldist));
        peak_inbound_reldist{a,i} = median(double(f(a).output{1,1}(i).peaks_select_inbound_reldist));
    end
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; number_steps=0;
    
    for i=rat5_idx
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        d=index(1); e=index(2);

        peak_outbound_reldist{a,i} = median(double(f(a).output{1,1}(i).peaks_select_outbound_reldist));
        peak_inbound_reldist{a,i} = median(double(f(a).output{1,1}(i).peaks_select_inbound_reldist));
    end
end

% peak 
peak_outbound_reldist_all=cell2mat(peak_outbound_reldist(:)');
peak_inbound_reldist_all=cell2mat(peak_inbound_reldist(:)');

peak_outbound_reldist_all=peak_outbound_reldist_all(peak_outbound_reldist_all<50);
peak_inbound_reldist_all=peak_inbound_reldist_all(peak_inbound_reldist_all<50);

data=[peak_outbound_reldist_all';peak_inbound_reldist_all'];
id=[ones(length(peak_outbound_reldist_all),1); 1+(ones(length(peak_inbound_reldist_all),1))];
kruskalwallis(data,id)

% plot values 
ax=figure(1); hold on 
boxplot(data,id, 'Notch', 'on'); 
all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);
xlabel(['OUT', 'IN']);
ylim([0 40])

%% 2. Complimnetary analysis to caluclate the length of theta sequences in a windows 
% check out also the peak of theta sequence upon scaling with theta sequences 
% OUTBOUND ALL PLANTS 
clear all;

% initialization parameters 
mod_score_window=0.04; % corresponds to 80ms around the plant time 
t_pre=0.07; t_post=0.07; % the usal window around which we compute one cycle of MUA fluctuations ~125 ms is a theta cycle 
savefigs=0;

results_dir=''; cd(results_dir);
results_filename=[''];

load([results_filename,'.mat']); 
load('decode_quality_idx.mat')

animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i= rat1_idx
        %
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_outbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_stance_outbound_run;
        theta_troughs_outbound=f(a).output{1,1}(i).theta_troughs_outbound; speed_outbound{a,i}=f(a).output{1,1}(i).speed_outbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_outbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_outbound_reldist{run_i}) && (size(forelimb_outbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_outbound{run_i}) && size(theta_troughs_outbound{run_i},2)>1));
        end
        forelimb_outbound_reldist=forelimb_outbound_reldist(idx); 
        theta_troughs_outbound=theta_troughs_outbound(idx); clear idx;

%        keyboard
        % collect individual run avergaes
        if size(forelimb_outbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
        else
        for run_i=1:size(forelimb_outbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_outbound_reldist_analyse(:,1)=forelimb_outbound_reldist{run_i}; 
            forelimb_outbound_reldist_analyse(:,2)=forelimb_outbound_reldist_analyse;
            limb_select_reldist=forelimb_outbound_reldist_analyse;

            theta_trough_outbound_analyse(:,1)=theta_troughs_outbound{run_i};
            theta_trough_outbound_analyse(:,2)=theta_trough_outbound_analyse;
            trough_select=theta_trough_outbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end

       for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_outbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_outbound_analyse 
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        %keyboard
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_outbound_reldist trig_plot_mua ...
            forelimb_outbound forelimb_outbound hindlimb_outbound hindlimb_outbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean

end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat2_idx
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_outbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_stance_outbound_run;
        theta_troughs_outbound=f(a).output{1,1}(i).theta_troughs_outbound; speed_outbound{a,i}=f(a).output{1,1}(i).speed_outbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_outbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_outbound_reldist{run_i}) && (size(forelimb_outbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_outbound{run_i}) && size(theta_troughs_outbound{run_i},2)>1));
        end
        forelimb_outbound_reldist=forelimb_outbound_reldist(idx); 
        theta_troughs_outbound=theta_troughs_outbound(idx); clear idx;

%        keyboard
        % collect individual run avergaes
 if size(forelimb_outbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_outbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_outbound_reldist_analyse(:,1)=forelimb_outbound_reldist{run_i}; 
            forelimb_outbound_reldist_analyse(:,2)=forelimb_outbound_reldist_analyse;
            limb_select_reldist=forelimb_outbound_reldist_analyse;

            theta_trough_outbound_analyse(:,1)=theta_troughs_outbound{run_i};
            theta_trough_outbound_analyse(:,2)=theta_trough_outbound_analyse;
            trough_select=theta_trough_outbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end


            for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_outbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_outbound_analyse
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_outbound_reldist trig_plot_mua ...
            forelimb_outbound forelimb_outbound hindlimb_outbound hindlimb_outbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat3_idx
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_outbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_stance_outbound_run;
        theta_troughs_outbound=f(a).output{1,1}(i).theta_troughs_outbound; speed_outbound{a,i}=f(a).output{1,1}(i).speed_outbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_outbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_outbound_reldist{run_i}) && (size(forelimb_outbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_outbound{run_i}) && size(theta_troughs_outbound{run_i},2)>1));
        end
        forelimb_outbound_reldist=forelimb_outbound_reldist(idx); 
        theta_troughs_outbound=theta_troughs_outbound(idx); clear idx;

%        keyboard
        % collect individual run avergaes
 if size(forelimb_outbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_outbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_outbound_reldist_analyse(:,1)=forelimb_outbound_reldist{run_i}; 
            forelimb_outbound_reldist_analyse(:,2)=forelimb_outbound_reldist_analyse;
            limb_select_reldist=forelimb_outbound_reldist_analyse;

            theta_trough_outbound_analyse(:,1)=theta_troughs_outbound{run_i};
            theta_trough_outbound_analyse(:,2)=theta_trough_outbound_analyse;
            trough_select=theta_trough_outbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end


            for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_outbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_outbound_analyse
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_outbound_reldist trig_plot_mua ...
            forelimb_outbound forelimb_outbound hindlimb_outbound hindlimb_outbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat5_idx
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_outbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_outbound_reldist = f(a).output{1,1}(i).forelimb_stance_outbound_run;
        theta_troughs_outbound=f(a).output{1,1}(i).theta_troughs_outbound; speed_outbound{a,i}=f(a).output{1,1}(i).speed_outbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_outbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_outbound_reldist{run_i}) && (size(forelimb_outbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_outbound{run_i}) && size(theta_troughs_outbound{run_i},2)>1));
        end
        forelimb_outbound_reldist=forelimb_outbound_reldist(idx); 
        theta_troughs_outbound=theta_troughs_outbound(idx); clear idx;

%        keyboard
        % collect individual run avergaes
 if size(forelimb_outbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_outbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_outbound_reldist_analyse(:,1)=forelimb_outbound_reldist{run_i}; 
            forelimb_outbound_reldist_analyse(:,2)=forelimb_outbound_reldist_analyse;
            limb_select_reldist=forelimb_outbound_reldist_analyse;

            theta_trough_outbound_analyse(:,1)=theta_troughs_outbound{run_i};
            theta_trough_outbound_analyse(:,2)=theta_trough_outbound_analyse;
            trough_select=theta_trough_outbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end


            for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_outbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_outbound_analyse
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_outbound_reldist trig_plot_mua ...
            forelimb_outbound forelimb_outbound hindlimb_outbound hindlimb_outbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean
end

% length of theta sequence in theta phase units 
cdat_interpt_theta=cdat_interpt_theta(:);
cdat_interpt_theta=cdat_interpt_theta(~cellfun('isempty', cdat_interpt_theta'));

% whole animal all troughs
for i=1:size(cdat_interpt_theta,1)        
    if sum(isnan(cdat_interpt_theta{i}))<1
        cdat_interpt_theta_plot(i)={cdat_interpt_theta{i}'};
        cdat_interpt_theta_plot_length_outbound(i)=max(cdat_interpt_theta{i})-min(cdat_interpt_theta{i});
    else 
    end
end
cdat_interpt_theta_plot=cdat_interpt_theta_plot(~cellfun('isempty', cdat_interpt_theta_plot));
cdat_interpt_theta_plot=cell2mat(cdat_interpt_theta_plot');

% whole animal all troughs mean per epoch
cdat_interpt_theta_length_epoch_all=cdat_interpt_theta_length_epoch_all(:);
cdat_interpt_theta_length_epoch_all=cdat_interpt_theta_length_epoch_all(cdat_interpt_theta_length_epoch_all~=0);
cdat_interpt_theta_length_outbound_epoch_all=cdat_interpt_theta_length_epoch_all(~isnan(cdat_interpt_theta_length_epoch_all));

reldist_time= [xi'; flip(xi'); xi(1)];
y1 = nanmean(cdat_interpt_theta_plot,1)-(std(cdat_interpt_theta_plot,[],1)/sqrt(size(cdat_interpt_theta_plot,1)));
y2 = nanmean(cdat_interpt_theta_plot,1)+(std(cdat_interpt_theta_plot,[],1)/sqrt(size(cdat_interpt_theta_plot,1)));
patch_in_dist = [ y1, fliplr(y2), y1(1)];
hold on; p2 = patch(reldist_time,patch_in_dist',repmat(size(cdat_interpt_theta_plot,1),size(cdat_interpt_theta_plot,2)*2+1,1),'facealpha', 0.5,'linestyle','none');
%clear y1 y2
plot(xi,nanmean(cdat_interpt_theta_plot,1));
text(0.5,-0.05,num2str([size(cdat_interpt_theta_plot,1)]), 'Color', 'r', 'FontSize', 12);
xlabel('Time [fraction of time between theta troughs]');
ylabel('Relative Distance [cm]')

%% INBOUND ALL PLANTS 
clearvars -except f t_pre t_post mod_score_window speed_outbound cdat_interpt_theta_plot_length_outbound cdat_interpt_theta_length_outbound_epoch_all
load('decode_quality_idx.mat')
animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};
for a=[1]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat1_idx
        %
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_inbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_stance_inbound_run;
       theta_troughs_inbound=f(a).output{1,1}(i).theta_troughs_inbound; speed_inbound{a,i}=f(a).output{1,1}(i).speed_inbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_inbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_inbound_reldist{run_i}) && (size(forelimb_inbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_inbound{run_i}) && size(theta_troughs_inbound{run_i},2)>1));
        end
        forelimb_inbound_reldist=forelimb_inbound_reldist(idx); 
        theta_troughs_inbound=theta_troughs_inbound(idx); clear idx;

 if size(forelimb_inbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_inbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_inbound_reldist_analyse(:,1)=forelimb_inbound_reldist{run_i}; 
            forelimb_inbound_reldist_analyse(:,2)=forelimb_inbound_reldist_analyse;
            limb_select_reldist=forelimb_inbound_reldist_analyse;

            theta_trough_inbound_analyse(:,1)=theta_troughs_inbound{run_i};
            theta_trough_inbound_analyse(:,2)=theta_trough_inbound_analyse;
            trough_select=theta_trough_inbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end

        for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_inbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_inbound_analyse 
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        %keyboard
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_inbound_reldist trig_plot_mua ...
            forelimb_inbound forelimb_inbound hindlimb_inbound hindlimb_inbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean

end
for a=[2]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat2_idx
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_inbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_stance_inbound_run;
       theta_troughs_inbound=f(a).output{1,1}(i).theta_troughs_inbound; speed_inbound{a,i}=f(a).output{1,1}(i).speed_inbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_inbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_inbound_reldist{run_i}) && (size(forelimb_inbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_inbound{run_i}) && size(theta_troughs_inbound{run_i},2)>1));
        end
        forelimb_inbound_reldist=forelimb_inbound_reldist(idx); 
        theta_troughs_inbound=theta_troughs_inbound(idx); clear idx;

 if size(forelimb_inbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_inbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_inbound_reldist_analyse(:,1)=forelimb_inbound_reldist{run_i}; 
            forelimb_inbound_reldist_analyse(:,2)=forelimb_inbound_reldist_analyse;
            limb_select_reldist=forelimb_inbound_reldist_analyse;

            theta_trough_inbound_analyse(:,1)=theta_troughs_inbound{run_i};
            theta_trough_inbound_analyse(:,2)=theta_trough_inbound_analyse;
            trough_select=theta_trough_inbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end

       for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_inbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_inbound_analyse 
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        %keyboard
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_inbound_reldist trig_plot_mua ...
            forelimb_inbound forelimb_inbound hindlimb_inbound hindlimb_inbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean
end
for a=[3]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat3_idx
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_inbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_stance_inbound_run;
       theta_troughs_inbound=f(a).output{1,1}(i).theta_troughs_inbound; speed_inbound{a,i}=f(a).output{1,1}(i).speed_inbound;

        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_inbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_inbound_reldist{run_i}) && (size(forelimb_inbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_inbound{run_i}) && size(theta_troughs_inbound{run_i},2)>1));
        end
        forelimb_inbound_reldist=forelimb_inbound_reldist(idx); 
        theta_troughs_inbound=theta_troughs_inbound(idx); clear idx;
      
 if size(forelimb_inbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_inbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_inbound_reldist_analyse(:,1)=forelimb_inbound_reldist{run_i}; 
            forelimb_inbound_reldist_analyse(:,2)=forelimb_inbound_reldist_analyse;
            limb_select_reldist=forelimb_inbound_reldist_analyse;

            theta_trough_inbound_analyse(:,1)=theta_troughs_inbound{run_i};
            theta_trough_inbound_analyse(:,2)=theta_trough_inbound_analyse;
            trough_select=theta_trough_inbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end

       for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_inbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_inbound_analyse 
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        %keyboard
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_inbound_reldist trig_plot_mua ...
            forelimb_inbound forelimb_inbound hindlimb_inbound hindlimb_inbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean
end
for a=[5]
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=((ind(:,1)==1)|(ind(:,1)==2)|(ind(:,1)==3)|(ind(:,1)==4)|(ind(:,1)==5)| (ind(:,1)==6)); counter=0; num_runs=0; 
    
    for i=rat5_idx
        counter=[counter; size(i,1)];
        
        if size(f(a).output{1,1}(i).forelimb_stance_inbound_run,2)>1      
        
        animal=animals{1,a};
        base_dir='';
        destdir=[base_dir,'/',animal,'/','filterframework', '/'];
        
        index=f(a).output{1,1}(i).index;close all
        
        d=index(1); e=index(2);
        
        task = loaddatastruct(destdir,animal,'task',[d e]);
        trials = loaddatastruct(destdir,animal,'trials',[d e]);
        tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
        posteriorts = f(a).output{1,1}(i).posteriorts;

        tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
        tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
        tets=vertcat(tets_ca1R,tets_ca1L);

        eeg = loadeegstruct(destdir,animal,'eeggnd',d,e,tets);
        eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
        
        marks = loaddatastruct(destdir, animal,'marks',[d e ]);

        ahead_behind_dist_smooth = f(a).output{1,1}(i).ahead_behind_smooth;
        forelimb_inbound_reldist = f(a).output{1,1}(i).forelimb_stance_inbound_run;
       theta_troughs_inbound=f(a).output{1,1}(i).theta_troughs_inbound; speed_inbound{a,i}=f(a).output{1,1}(i).speed_inbound;
      
        % keyboard
        % remove empty cells from runs
        for run_i=1:size(forelimb_inbound_reldist,2)
            idx(run_i)=(~isempty(forelimb_inbound_reldist{run_i}) && (size(forelimb_inbound_reldist{run_i},1)>1) && (~isempty(theta_troughs_inbound{run_i}) && size(theta_troughs_inbound{run_i},2)>1));
        end
        forelimb_inbound_reldist=forelimb_inbound_reldist(idx); 
        theta_troughs_inbound=theta_troughs_inbound(idx); clear idx;

 if size(forelimb_inbound_reldist,2)<1
%             trig_plot_reldist={};
%             trig_plot_theta={};
%             mod_score_sec={};
%             mod_score_theta={};
 else
        for run_i=1:size(forelimb_inbound_reldist,2)
            num_runs=[num_runs; size(run_i,1)];
            forelimb_inbound_reldist_analyse(:,1)=forelimb_inbound_reldist{run_i}; 
            forelimb_inbound_reldist_analyse(:,2)=forelimb_inbound_reldist_analyse;
            limb_select_reldist=forelimb_inbound_reldist_analyse;

            theta_trough_inbound_analyse(:,1)=theta_troughs_inbound{run_i};
            theta_trough_inbound_analyse(:,2)=theta_trough_inbound_analyse;
            trough_select=theta_trough_inbound_analyse;

            trig_plot_reldist(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', limb_select_reldist, 't_pre',t_pre, 't_post',t_post);
            trig_plot_theta(run_i)=periseg('cdat', contchans(ahead_behind_dist_smooth, 'chans', 1), 'segs', trough_select, 't_pre',t_pre, 't_post',t_post);

            %keyboard
            % caluclate modulation score per run 
            lags_reldist=trig_plot_reldist(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_reldist(run_i).peritrig, lags_reldist,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_sec(run_i)=max(pks_val)-min(trig_plot_reldist(run_i).peritrig(lags_reldist==trs_ind)); clear trs_ind pks_val trs_val;

            lags_theta=trig_plot_theta(run_i).lags;
            [pks_val,~]=findpeaks(trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            [trs_val,trs_ind]=findpeaks(-trig_plot_theta(run_i).peritrig, lags_theta,'MinPeakDistance', mod_score_window);
            trs_ind=trs_ind(find(trs_val==min(trs_val))); 
            mod_score_theta(run_i)=max(pks_val)-min(trig_plot_theta(run_i).peritrig(lags_theta==trs_ind)); clear trs_ind pks_val trs_val;

            xi=0:1/500:1;
            for i_plant=1:size(limb_select_reldist,1)-1
                        idx_i_plant=(find((posteriorts>limb_select_reldist(i_plant,1)) & (posteriorts<limb_select_reldist(i_plant+1,1))));

                        if ~isempty(idx_i_plant)
                            cdat_plant_plant{a,i,run_i,i_plant}.y=ahead_behind_dist_smooth.data(idx_i_plant)';
                            cdat_plant_plant{a,i,run_i,i_plant}.x=rescale(1:size(idx_i_plant,1),0,1);
                            %plot(cdat_plant_plant_x{a,i,run_i,i_plant},cdat_plant_plant_y{a,i,run_i,i_plant}); hold on 
                             
                            cdat_interpt{a,i,run_i,i_plant}=interp1q(cdat_plant_plant{a,i,run_i,i_plant}.x',cdat_plant_plant{a,i,run_i,i_plant}.y',xi');
                            clear idx_i_plant;
                        else
                            cdat_plant_plant{a,i,run_i,i_plant}.x=NaN;
                            cdat_plant_plant{a,i,run_i,i_plant}.y=NaN;
                        end
            end

       for i_trough=1:size(trough_select,1)-1
                        idx_i_troughs=(find((posteriorts>trough_select(i_trough,1)) & (posteriorts<trough_select(i_trough+1,1))));

                        if ~isempty(idx_i_troughs)
                            cdat_trough_trough{a,i,run_i,i_trough}.y=ahead_behind_dist_smooth.data(idx_i_troughs)';
                            cdat_trough_trough{a,i,run_i,i_trough}.x=rescale(1:size(idx_i_troughs,1),0,1);
                            cdat_interpt_theta{a,i,run_i,i_trough}=interp1q(cdat_trough_trough{a,i,run_i,i_trough}.x',cdat_trough_trough{a,i,run_i,i_trough}.y',xi');
                            cdat_interpt_theta_length_epoch(run_i,i_trough)=max(cdat_interpt_theta{a,i,run_i,i_trough})-min(cdat_interpt_theta{a,i,run_i,i_trough});
                            clear idx_i_plant;
                        else
                            cdat_trough_trough{a,i,run_i,i_trough}.x=NaN;
                            cdat_trough_trough{a,i,run_i,i_trough}.y=NaN;
                        end
            end
            clear forelimb_inbound_reldist_analyse limb_select_reldist RandomOffsets trough_select theta_trough_inbound_analyse 
        end

        lags_reldist=trig_plot_reldist(run_i).lags;
        mydata_reldist{a,i}=trig_plot_reldist;
        mydata_theta{a,i}=trig_plot_theta;
        mod_score_sec_all{a,i}=mod_score_sec;
        mod_score_theta_all{a,i}=mod_score_theta;
        cdat_interpt_theta_length_epoch=cdat_interpt_theta_length_epoch(:);
        cdat_interpt_theta_length_epoch_all(a,i)=double(median(cdat_interpt_theta_length_epoch(cdat_interpt_theta_length_epoch~=0)));
       
        %keyboard
        clear mod_score_theta mod_score_sec trig_plot_reldist trig_plot_theta i_plant posteriorts tets tets_ca1R tets_ca1L RandomOffsets ahead_behind_dist_smooth marks eeg eegtimes cdat_eeg cdat_mua muatimes index d e forelimb_inbound_reldist trig_plot_mua ...
            forelimb_inbound forelimb_inbound hindlimb_inbound hindlimb_inbound cdat_mua eeg eegtimes run_i cdat_interpt_theta_length_epoch

        end
        else 
        end 
     end       
        number_of_epochs_all(a)=sum(counter);
        number_runs_all(a)=sum(num_runs);
        clear counter num_runs num_runs_clean
end

% time in theta phase units
cdat_interpt_theta=cdat_interpt_theta(:);
cdat_interpt_theta=cdat_interpt_theta(~cellfun('isempty', cdat_interpt_theta'));

for i=1:size(cdat_interpt_theta,1)        
    if sum(isnan(cdat_interpt_theta{i}))<1
        cdat_interpt_theta_plot(i)={cdat_interpt_theta{i}'};
        cdat_interpt_theta_plot_length_inbound(i)=max(cdat_interpt_theta{i})-min(cdat_interpt_theta{i});
    else 
    end
end
cdat_interpt_theta_plot=cdat_interpt_theta_plot(~cellfun('isempty', cdat_interpt_theta_plot));
cdat_interpt_theta_plot=cell2mat(cdat_interpt_theta_plot');

% whole animal all troughs mean per epoch
cdat_interpt_theta_length_epoch_all=cdat_interpt_theta_length_epoch_all(:);
cdat_interpt_theta_length_epoch_all=cdat_interpt_theta_length_epoch_all(cdat_interpt_theta_length_epoch_all~=0);
cdat_interpt_theta_length_inbound_epoch_all=cdat_interpt_theta_length_epoch_all(~isnan(cdat_interpt_theta_length_epoch_all));

reldist_time= [xi'; flip(xi'); xi(1)];
y1 = nanmean(cdat_interpt_theta_plot,1)-(std(cdat_interpt_theta_plot,[],1)/sqrt(size(cdat_interpt_theta_plot,1)));
y2 = nanmean(cdat_interpt_theta_plot,1)+(std(cdat_interpt_theta_plot,[],1)/sqrt(size(cdat_interpt_theta_plot,1)));
patch_in_dist = [ y1, fliplr(y2), y1(1)];
hold on; p2 = patch(reldist_time,patch_in_dist',repmat(size(cdat_interpt_theta_plot,1),size(cdat_interpt_theta_plot,2)*2+1,1),'facealpha', 0.5,'linestyle','none');
%clear y1 y2
plot(xi,nanmean(cdat_interpt_theta_plot,1));
text(0.5,-0.05,num2str([size(cdat_interpt_theta_plot,1)]), 'Color', 'r', 'FontSize', 12);
xlabel('Time [fraction of time between theta troughs]');
ylabel('Relative Distance [cm]')

% stats 
a= cdat_interpt_theta_plot_length_outbound(cdat_interpt_theta_plot_length_outbound~=0);
b= cdat_interpt_theta_plot_length_inbound(cdat_interpt_theta_plot_length_inbound~=0);

data=[a,b]; id=[ones(length(a),1); 1+(ones(length(b),1))]; clear a b;
kruskalwallis(data,id); clear data id 

% stats 
a= cdat_interpt_theta_length_outbound_epoch_all;
b= cdat_interpt_theta_length_inbound_epoch_all;

data=[a',b']'; id=[ones(length(a),1); 1+(ones(length(b),1))];
kruskalwallis(data,id);

% plot values 
close all
ax=figure(1); hold on 
boxplot(data,id, 'Notch', 'on'); 
all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);
xlabel(['','OUT','IN']);
ylim([0 40])

%% SAVE AND PLOT 
% variables above saved in EDFigure6A_boxplots.mat
% run from here to plot 
clear all; load('EDFigure6.mat')

% plot values 
close all
ax=figure(1); hold on 
boxplot(data_peak_dist_median_perepoch,id_peak_dist_median_perepoch, 'Notch', 'on'); 
all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);
xlabel(['','OUT','IN']);
ylim([0 40])

% plot values 
ax=figure(2); hold on 
boxplot(data_peak_theta_seq,id_peak_theta_seq, 'Notch', 'on'); 
all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);
xlabel(['OUT', 'IN']);
ylim([0 40])

%% SAVE and PLOT 
% variables saved in in EDFigure6.mat
% epoch examples of plant-theta phase relationship

close all;
rose(rat1_epoch10_forelimb_plant_phasedist, 15);