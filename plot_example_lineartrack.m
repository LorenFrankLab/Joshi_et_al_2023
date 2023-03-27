%% plot_example_lineartrack.m 
% -> dfa_clusterfilter
% -> Figure1ACD.mat
% -> assumes data is preprocessed in filterframework format

clear all; close all; clc;
animal=''; %animal ID
d = 2; e =4; ref=13; fps = 125; %epoch ID
basenm=[animal,'_0',num2str(d),'_',num2str(e)];
destdir=[base_dir,'/',animal,'/','filterframework', '/'];
resdir=[base_dir,'/',animal,'/','results', '/'];

%% Load task info
task = loaddatastruct(destdir,animal,'task',[d e]);
trials = loaddatastruct(destdir,animal,'trials',[d e]);
rips = loaddatastruct(destdir,animal,'ca1rippleskons',[d e]);
tetinfo = loaddatastruct(destdir,animal,'tetinfo',[d e]);
tets_ca1R = evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1R'')');
tets_ca1L=evaluatefilter(tetinfo{d}{e},'isequal($area,''ca1L'')');
tets=vertcat(tets_ca1R, tets_ca1L);

eeg = loadeegstruct(destdir,animal,'eeg',d,e,tets);
eegtimes = geteegtimes(eeg{d}{e}{tets(1)});
posdlc = loaddatastruct(destdir,animal,'posdlc',[d e]);
rawpos = loaddatastruct(destdir,animal,'rawpos',[d e]);

%% Load theta 
% Collect LFP unreferenced i.e. grounded to cerebellum from FFeeg folder
LFPgnd_filename=[animal,'eeggnd','0',num2str(d), '-', num2str(e),'-',num2str(ref),'.mat'];
LFPgnd_data=load(LFPgnd_filename);
LFP_rt_dat=[(LFPgnd_data.eeggnd{1, d}{1, e}{1, ref}.starttime):(1/LFPgnd_data.eeggnd{1, d}{1, e}{1, ref}.samprate):(LFPgnd_data.eeggnd{1, d}{1, e}{1, ref}.endtime)]';
LFPgnd_data=LFPgnd_data.eeggnd{1,d}{1,e}{1,ref}.data;

if ~(size(LFP_rt_dat,1)==size(LFPgnd_data,1))
    LFPgnd_data=LFPgnd_data(1:size(LFP_rt_dat,1));
else 
end

%% Load the adjusted timestamps. This step requirs pos data to be in fliterframework format

% load the pos file. This was created during Binary2FFPos with 
% argument, use_ptp=true. 
load_file_cam_rt_fit=[animal,'posdlc0',num2str(d),'.mat'];
cam_rt_fit_pos=load(load_file_cam_rt_fit);

% adjusted camera time. This is the ptp adjusted camera time. 
cam_rt_fit=cam_rt_fit_pos.posdlc{1,d}{1,e}.data(:,1);

% corresponding camera frame. These are the corresponding camera frame
% numbers for the ptp adjusted time. 
cam_rt_fit_ind=cam_rt_fit_pos.posdlc{1,d}{1,e}.data(:,6);

% estimated framerate based on camera time 
est_framerate=median(1./diff(cam_rt_fit));

%% 3. Load dlc_results posdlc

% load dlc results from rawpos. These are after Dlc2FFPOS
dlc_results=cam_rt_fit_pos.posdlc{1,d}{1,e}; 

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
dlc_forepawL_x=dlc_results.data(:,8);
dlc_forepawL_y=dlc_results.data(:,9);

%forepawR
dlc_forepawR_x= dlc_results.data(:,11);
dlc_forepawR_y= dlc_results.data(:,12);

%hindpawL
dlc_hindpawL_x=dlc_results.data(:,14);
dlc_hindpawL_y=dlc_results.data(:,15);

%hindpawR
dlc_hindpawR_x=dlc_results.data(:,17);
dlc_hindpawR_y=dlc_results.data(:,18);

dlc_n_records = size(dlc_hindpawL_x,1);
fprintf('# of DeepLabCut timestamps: %d \n', dlc_n_records)

clearvars dlc_results

%% Load spiking data 
animals={''};

for a = 1:length(animals)

     animal=animals(a);
     animalinfo =animaldef(animal{1});
     animaldir = animalinfo{2};
     animalprefix = animalinfo{3};

     timefilter = {};
    % epoch filters | select wtrack epochs
     epochfilter = '(isequal($type,''run'') && isequal($environment,''lineartrack''))';

    % iterators | these are goimg to define your base unit of analysis, such as per cell or per epoch
     iterator = 'singlecellanal';
     cellfilter_all= ['($accepted== 1)']; %% keep an eye out for format seems important 
     f_all= createfilter('animal', animal,'epochs', epochfilter, 'cells', cellfilter_all, 'excludetime', timefilter, 'iterator', iterator);
     % Set Analysis Function
     f_all= setfilterfunction(f_all, 'dfa_clusterfilter', {'spikes','posdlc','task','cellinfo'},animaldir,animalprefix);
     % Run Analysis
     f_all= runfilter(f_all);

end
% f.output now contains spike times. 
% for each epoch collect all spike across tetrodes 
num_cells_all=size((f_all.output{1,1}),2);
for i=1:num_cells_all
    day_number_all(i)=f_all.output{1,1}(i).index(:,1);
    epoch_number_all(i)=f_all.output{1,1}(i).index(:,2);
    day_epoch_all(i,1)=[day_number_all(i)];
    day_epoch_all(i,2)=[epoch_number_all(i)];
end
num_epochs_all=unique(day_epoch_all,'rows');

%% Collect all spikes from accepted cells | lets collect all spikes from unique epochs into a struct | 

for i=5 % for day 2 epoch 4 
    %i=1:length(num_epochs_all) % if you want to look at all epochs. 
    animal=''; a=1;
    
    d=num_epochs_all(i,1); e=num_epochs_all(i,2);
    ind = cell2mat(arrayfun(@(x) x.index',f_all.output{1},'Un',0))'; idx=zeros(length(ind));
    day_select=(ind(:,1)==d);
    epoch_select=(ind(:,2)==e);
    idx=find((day_select+epoch_select)==2);
    ind_select(idx)=1; 
    results_dir='/results'; cd(results_dir); results_filename=[''];
    
    load(results_filename)
    run_periods=thetamod_f(1).output{1, 1}(9).run_periods; clear f
     
     % ForepawR 
     [cam_rt_fit_forepawR_results, ~,~]=smooth_dlc_stepcycles_lineartrack(dlc_tail_x,dlc_forepawR_x, cam_rt_fit, dlc_tail_vel, est_framerate, run_periods);

     % plotted run priod example
     xleft=6242.5; xright=6246; %start and end of the run period
     
    % Forelimb Plant and Lift times used in Joshi et al 2023. This is the 
    % value that corresponds to the limb of the rat firmly touching the 
    % surface of the track as in ED Figure 1. 
    valid_forelimbR_lift_collect 
    valid_forelimbR_plant_collect
    
    % Load all spikes in the epoch 
    for idx_select=1:(size(idx,1))
                
        k=idx(idx_select);
        spiketimes{idx_select}=cell2mat(arrayfun(@(x) x.spiketimes',f_all.output{1,1}(k),'Un',0))'; %save spiketimes k

        % keyboard
        spiketimes_=spiketimes{idx_select};
        spiketimes_select=find((spiketimes_>xleft) & (spiketimes_<xright));
        num_spikes_run(idx_select)=length(spiketimes_select);
        spiketimes_mean(idx_select)=median(spiketimes_(spiketimes_select));
        clear spiketimes_ spiketimes_select k
    end

%% SAVE AND PLOT 
% Variables saved in Figure1ACD.mat. Run from here if using 
    
    load('Figure1ACD.mat');
    [val_, idx_toplot] = sort(spiketimes_mean',1);

    % load all cell spikes and plot rasters 
    for idx_select=1:(size(idx,1))
                
        k=idx_toplot(idx_select);
        spiketimes_toplot{idx_select}=spiketimes{k};

        % can also plot the spiketimerasters 
        % keyboard
        spiketimes_=spiketimes_toplot{idx_select};
        spiketimes_select=find(spiketimes_>xleft & spiketimes_<xright);
        
        if ~isempty(spiketimes_) && ~isempty(spiketimes_select) && size(spiketimes_select,1)>10
            % close all; 
            figure(1000+e); hold on; 
            ax(1)=subplot(14,1,1:5); %ax(1).YDir = 'reverse';
            plot([spiketimes_(spiketimes_select) spiketimes_(spiketimes_select)], [idx_select-0.5 idx_select+0.5],'k'); hold on 
             
%             plot_order=length(spiketimes_(spiketimes_select));
%             if plot_order<45 % putative pyramidal cells
%                 ax(1)=subplot(14,1,1:5)
%                 plot([spiketimes_(spiketimes_select) spiketimes_(spiketimes_select)], [idx_select-0.5 idx_select+0.5],'k'); hold on 
%             elseif plot_order>=45 % puttaive interneurons
%                 ax(2)=subplot(14,1,6:7)
%                 plot([spiketimes_(spiketimes_select) spiketimes_(spiketimes_select)], [idx_select-0.5 idx_select+0.5],'k'); hold on 
%             end
%             clear plot_order
            
            % keyboard
            clear spiketimes_ spiketimes_select k
        else
            % keyboard
            clear spiketimes_ spiketimes_select k
            continue            
        end
          clear spiketimes_ spiketimes_select k
    end
    
    %ylim([0 size(idx,1)+1]);
    set(gca,'box','off'); set(gca,'xtick',[]);
    xlim([xleft xright]); hold on
    ylabel('Unit ID');
    
    % Plot LFP 
    ax(2)=subplot(14,1,6:7);
    LFP_rt_dat_collect=LFP_rt_dat(LFP_rt_dat>xleft & LFP_rt_dat<xright); % save LFP_rt_dat_collect
    LFPgnd_data_collect=LFPgnd_data(LFP_rt_dat>xleft & LFP_rt_dat<xright); % save LFPgnd_data_collect
    plot(LFP_rt_dat_collect,LFPgnd_data_collect, 'k'); hold on 
    xlim([xleft xright]); set(gca,'box','off'); set(gca,'xtick',[])

    % Plot Pose 
    ax(4)=subplot(14,1,8:14); hold on 
    plot(cam_rt_fit,dlc_nose_x, 'color', 'b'); % save cam_rt_fit_dlc nose_x dlc_tails_x
    plot(cam_rt_fit,(dlc_forepawR_x), 'color', [0.597656250000000	0.195312500000000	0.796875000000000]); % save dlc_forelimbR_x
    plot([valid_forelimbR_plant valid_forelimbR_plant],[0 120],'color','k', 'LineWidth', 1);
    plot([valid_forelimbR_lift_collect valid_forelimbR_lift_collect],[0 120],'color', 'r', 'LineWidth', 1);
    plot(cam_rt_fit,dlc_tail_x,'color', [0.5 0.5 0.5]);
    xlim([xleft xright]); set(gca,'YGrid','on'); 
    set(gca,'box','off'); %set(gca,'xtick',[]); %ax(4).YDir = 'reverse'; 
    linkaxes(ax, 'x');  xlim([xleft xright]); set(gca,'box','off'); xlabel('Time [s]');
    
%%
%     output_all.muatimes_spikes{i}=muatimes_spikes;
%     output_all.muatrace_spikes{i}=muatrace_spikes;
%     output_all.muatimes_marks{i}=muatimes_marks;
%     output_all.muatrace_marks{i}=muatrace_marks;
%     output_all.index{i}=[d e];
%    clear idx_toplot_ spiketimes spiketimes_all muatimes_marks muatrace_marks muatimes_spikes muatrace_spikes ind ind_select day_select epoch_select
%     clear xright xleft
end
