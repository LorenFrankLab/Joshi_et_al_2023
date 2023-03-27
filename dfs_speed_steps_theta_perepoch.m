%% dfs_speed_steps_theta_perepoch
% -> Figure1EF_EDFigure2.mat

%% OUTBOUND CENTER  
close all; clear all; 
results_filename='';
load([results_filename, '.mat'])
clearvars -except f;

savefigs=0; savefig_scatter=0; 
% parameters
destdir='';

animals={'Jaq' 'Roqui', 'Lotus', 'Monty', 'Peanut'};
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
    
    % calculate the correlation coefficient for each epoch
    for e=1:size(ind_select,1)
        
        inst_vel_outbound{a,e}= {cell2mat((arrayfun(@(x) x.inst_vel_outbound,f(a).output{1}(e),'Un',0))')};
        inst_acc_outbound{a,e}= {cell2mat((arrayfun(@(x) x.inst_acc_outbound,f(a).output{1}(e),'Un',0))')};
        inst_freq_outbound{a,e}={cell2mat((arrayfun(@(x) x.inst_theta_outbound,f(a).output{1}(e),'Un',0))')};
        inst_steps_outbound{a,e}={cell2mat((arrayfun(@(x) x.inst_steps_outbound,f(a).output{1}(e),'Un',0))')};

        inst_vel_epoch=cell2mat(inst_vel_outbound{a,e}(:));
        inst_acc_epoch=cell2mat(inst_acc_outbound{a,e}(:));
        inst_theta_epoch=cell2mat(inst_freq_outbound{a,e}(:));
        inst_steps_epoch=cell2mat(inst_steps_outbound{a,e}(:));

        idx_include=find((inst_vel_epoch>vel_thresh_min) & (inst_steps_epoch>steps_thresh_min) & (inst_steps_epoch<steps_thresh_max) & (inst_theta_epoch>theta_thresh_min)); % excluding the noisy data 
        idx_include_vel=find((inst_vel_epoch>vel_thresh_min)&(inst_vel_epoch<vel_thresh_max)); % excluding the noisy data 
        
        %% inst vel versus theta 
        [r1, p1]=corrcoef(inst_vel_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r1) 
            r_veltheta{a,e}=r1(2);
            p_veltheta{a,e}=p1(2);
        else 
            r_veltheta{a,e}=[];
            p_veltheta{a,e}=[];
        end
         
        %% inst acc versus theta 
        [r2, p2]=corrcoef(inst_acc_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r2) 
            r_acctheta{a,e}=r2(2);
            p_acctheta{a,e}=p2(2);
        else 
            r_acctheta{a,e}=[];
            p_acctheta{a,e}=[];
        end
        
         %% inst steps versus inst theta 
        [r3, p3]=corrcoef(inst_steps_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');

        if ~isnan(r3) 
            r_stepstheta{a,e}=r3(2);
            p_stepstheta{a,e}=p3(2);
        else 
            r_stepstheta{a,e}=[];
            p_stepstheta{a,e}=[];
        end
        
        %% inst steps versus inst theta 
         [r4, p4]=corrcoef(inst_vel_epoch(idx_include_vel),inst_steps_epoch(idx_include_vel), 'rows', 'complete');
         
        if ~isnan(r4) 
            r_stepsvel{a,e}=r4(2);
            p_stepsvel{a,e}=p4(2);
        else 
            r_stepsvel{a,e}=[];
            p_stepsvel{a,e}=[];
        end
        %%
        clear r1 r2 r3 r4 p1 p2 p3 p4 inst_vel_epoch inst_acc_epoch inst_theta_epoch inst_steps_epoch idx_include idx_include_vel
    end
    
    counter(a)=count;
    
    %save per animal per epoch results 
    r_veltheta_outbound_center{a}=(r_veltheta(a,:));
    r_acctheta_outbound_center{a}=(r_acctheta(a,:));
    r_stepstheta_outbound_center{a}=(r_stepstheta(a,:));
    r_stepsvel_outbound_center{a}=(r_stepsvel(a,:));

    p_veltheta_outbound_center{a}=(p_veltheta(a,:));
    p_acctheta_outbound_center{a}=(p_acctheta(a,:));
    p_stepstheta_outbound_center{a}=(p_stepstheta(a,:));
    p_stepsvel_outbound_center{a}=(p_stepsvel(a,:));
end

r_veltheta_outbound_center_all=cell2mat(r_veltheta(~cellfun('isempty',r_veltheta)));
r_acctheta_outbound_center_all=cell2mat(r_acctheta(~cellfun('isempty',r_acctheta)));
r_stepstheta_outbound_center_all=cell2mat(r_stepstheta(~cellfun('isempty',r_stepstheta)));
r_stepsvel_outbound_center_all=cell2mat(r_stepsvel(~cellfun('isempty',r_stepsvel)));
    
p_veltheta_outbound_center_all=cell2mat(p_veltheta(~cellfun('isempty',p_veltheta)));
p_acctheta_outbound_center_all=cell2mat(p_acctheta(~cellfun('isempty',p_acctheta)));
p_stepstheta_outbound_center_all=cell2mat(p_stepstheta(~cellfun('isempty',p_stepstheta)));
p_stepsvel_outbound_center_all=cell2mat(p_stepsvel(~cellfun('isempty',p_stepsvel)));

%% INBOUND CENTER 
clearvars -except f savefig_scatter destdir ci_alpha r_* p_*
animals={'Jaq' 'Roqui', 'Lotus', 'Monty', 'Peanut'};

for a=[1,2,3,4,5]
    
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    
    ind_del=[];
    for i=1:size(ind,1)
        
        if ((isnan(f(a).output{1, 1}(i).inst_acc_inbound(1)))|(isnan(f(a).output{1, 1}(i).speed_inbound(1))))
             ind_del=[ind_del; i];
        end
    end
    f(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    
    %keyboard
    count=size(ind_select,1);
    
    % calculate the correlation coefficient for each epoch
    
    for e=1:size(ind_select,1)
        inst_vel_inbound{a,e}= {cell2mat((arrayfun(@(x) x.inst_vel_inbound,f(a).output{1}(e),'Un',0))')};
        inst_acc_inbound{a,e}= {cell2mat((arrayfun(@(x) x.inst_acc_inbound,f(a).output{1}(e),'Un',0))')};
        inst_freq_inbound{a,e}={cell2mat((arrayfun(@(x) x.inst_theta_inbound,f(a).output{1}(e),'Un',0))')};
        inst_steps_inbound{a,e}={cell2mat((arrayfun(@(x) x.inst_steps_inbound,f(a).output{1}(e),'Un',0))')};

        inst_vel_epoch=cell2mat(inst_vel_inbound{a,e}(:));
        inst_acc_epoch=cell2mat(inst_acc_inbound{a,e}(:));
        inst_theta_epoch=cell2mat(inst_freq_inbound{a,e}(:));
        inst_steps_epoch=cell2mat(inst_steps_inbound{a,e}(:));

        idx_include=find((inst_vel_epoch>vel_thresh_min) & (inst_steps_epoch>steps_thresh_min) & (inst_steps_epoch<steps_thresh_max) & (inst_theta_epoch>theta_thresh_min)); % excluding the noisy data 
        idx_include_vel=find((inst_vel_epoch>vel_thresh_min)&(inst_vel_epoch<vel_thresh_max)); % excluding the noisy data 
        
        %% inst vel versus theta 
        [r1, p1]=corrcoef(inst_vel_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r1) 
            r_veltheta{a,e}=r1(2);
            p_veltheta{a,e}=p1(2);
        else 
            r_veltheta{a,e}=[];
            p_veltheta{a,e}=[];
        end
         
        %% inst acc versus theta 
        [r2, p2]=corrcoef(inst_acc_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r2) 
            r_acctheta{a,e}=r2(2);
            p_acctheta{a,e}=p2(2);
        else 
            r_acctheta{a,e}=[];
            p_acctheta{a,e}=[];
        end
        
         %% inst steps versus inst theta 
        [r3, p3]=corrcoef(inst_steps_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');

        if ~isnan(r3) 
            r_stepstheta{a,e}=r3(2);
            p_stepstheta{a,e}=p3(2);
        else 
            r_stepstheta{a,e}=[];
            p_stepstheta{a,e}=[];
        end
        
        %% inst steps versus inst theta 
         [r4, p4]=corrcoef(inst_vel_epoch(idx_include_vel),inst_steps_epoch(idx_include_vel), 'rows', 'complete');
         
        if ~isnan(r4) 
            r_stepsvel{a,e}=r4(2);
            p_stepsvel{a,e}=p4(2);
        else 
            r_stepsvel{a,e}=[];
            p_stepsvel{a,e}=[];
        end
        %%
        clear r1 r2 r3 r4 p1 p2 p3 p4 inst_vel_epoch inst_acc_epoch inst_theta_epoch inst_steps_epoch idx_include idx_include_vel
    end
    
    counter(a)=count;
    
    %save per animal per epoch results 
    r_veltheta_inbound_center{a}=(r_veltheta(a,:));
    r_acctheta_inbound_center{a}=(r_acctheta(a,:));
    r_stepstheta_inbound_center{a}=(r_stepstheta(a,:));
    r_stepsvel_inbound_center{a}=(r_stepsvel(a,:));
    
    p_veltheta_inbound_center{a}=(p_veltheta(a,:));
    p_acctheta_inbound_center{a}=(p_acctheta(a,:));
    p_stepstheta_inbound_center{a}=(p_stepstheta(a,:));
    p_stepsvel_inbound_center{a}=(p_stepsvel(a,:));
end

r_veltheta_inbound_center_all=cell2mat(r_veltheta(~cellfun('isempty',r_veltheta)));
r_acctheta_inbound_center_all=cell2mat(r_acctheta(~cellfun('isempty',r_acctheta)));
r_stepstheta_inbound_center_all=cell2mat(r_stepstheta(~cellfun('isempty',r_stepstheta)));
r_stepsvel_inbound_center_all=cell2mat(r_stepsvel(~cellfun('isempty',r_stepsvel)));
    
p_veltheta_inbound_center_all=cell2mat(p_veltheta(~cellfun('isempty',p_veltheta)));
p_acctheta_inbound_center_all=cell2mat(p_acctheta(~cellfun('isempty',p_acctheta)));
p_stepstheta_inbound_center_all=cell2mat(p_stepstheta(~cellfun('isempty',p_stepstheta)));
p_stepsvel_inbound_center_all=cell2mat(p_stepsvel(~cellfun('isempty',p_stepsvel)));

%% Outer arms combined - Inbound
clearvars -except savefig_scatter destdir ci_alpha r_* p_*

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
    
    %calculate per epoch
    for e=1:size(ind_select,1)
        % keyboard 
        inst_vel_inbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_vel_inbound,f1(a).output{1}(e),'Un',0))');
        inst_acc_inbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_acc_inbound,f1(a).output{1}(e),'Un',0))');
        inst_freq_inbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f1(a).output{1}(e),'Un',0))');
        inst_steps_inbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f1(a).output{1}(e),'Un',0))');

        inst_vel_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_vel_inbound,f2(a).output{1}(e),'Un',0))');
        inst_acc_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_acc_inbound,f2(a).output{1}(e),'Un',0))');
        inst_freq_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f2(a).output{1}(e),'Un',0))');
        inst_steps_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f2(a).output{1}(e),'Un',0))');

        inst_vel_epoch=([(inst_vel_inbound_f1{a,e}(:)); (inst_vel_inbound_f2{a,e}(:))]); 
        inst_acc_epoch=([(inst_acc_inbound_f1{a,e}(:)); (inst_acc_inbound_f2{a,e}(:))]);
        inst_theta_epoch=([(inst_freq_inbound_f1{a,e}(:)); (inst_freq_inbound_f2{a,e}(:))]);
        inst_steps_epoch=([(inst_steps_inbound_f1{a,e}(:)); (inst_steps_inbound_f2{a,e}(:))]);
        
        idx_include=find((inst_vel_epoch>vel_thresh_min) & (inst_steps_epoch>steps_thresh_min) & (inst_steps_epoch<steps_thresh_max) & (inst_theta_epoch>theta_thresh_min)); % excluding the noisy data 
        idx_include_vel=find((inst_vel_epoch>vel_thresh_min)&(inst_vel_epoch<vel_thresh_max)); % excluding the noisy data 
        
        %% inst vel versus theta 
        [r1, p1]=corrcoef(inst_vel_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r1) 
            r_veltheta{a,e}=r1(2);
            p_veltheta{a,e}=p1(2);
        else 
            r_veltheta{a,e}=[];
            p_veltheta{a,e}=[];
        end
         
        %% inst acc versus theta 
        [r2, p2]=corrcoef(inst_acc_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r2) 
            r_acctheta{a,e}=r2(2);
            p_acctheta{a,e}=p2(2);
        else 
            r_acctheta{a,e}=[];
            p_acctheta{a,e}=[];
        end
        
         %% inst steps versus inst theta 
        [r3, p3]=corrcoef(inst_steps_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');

        if ~isnan(r3) 
            r_stepstheta{a,e}=r3(2);
            p_stepstheta{a,e}=p3(2);
        else 
            r_stepstheta{a,e}=[];
            p_stepstheta{a,e}=[];
        end
        
        %% inst steps versus inst theta 
         [r4, p4]=corrcoef(inst_vel_epoch(idx_include_vel),inst_steps_epoch(idx_include_vel), 'rows', 'complete');
         
        if ~isnan(r4) 
            r_stepsvel{a,e}=r4(2);
            p_stepsvel{a,e}=p4(2);
        else 
            r_stepsvel{a,e}=[];
            p_stepsvel{a,e}=[];
        end
        %%
        clear r1 r2 r3 r4 p1 p2 p3 p4 inst_vel_epoch inst_acc_epoch inst_theta_epoch inst_steps_epoch idx_include idx_include_vel
    end
    
    counter(a)=count;
    
    %save per animal per epoch results 
    r_veltheta_inbound_outer{a}=(r_veltheta(a,:));
    r_acctheta_inbound_outer{a}=(r_acctheta(a,:));
    r_stepstheta_inbound_outer{a}=(r_stepstheta(a,:));
    r_stepsvel_inbound_outer{a}=(r_stepsvel(a,:));
    
    p_veltheta_inbound_outer{a}=(p_veltheta(a,:));
    p_acctheta_inbound_outer{a}=(p_acctheta(a,:));
    p_stepstheta_inbound_outer{a}=(p_stepstheta(a,:));
    p_stepsvel_inbound_outer{a}=(p_stepsvel(a,:));
end

r_veltheta_inbound_outer_all=cell2mat(r_veltheta(~cellfun('isempty',r_veltheta)));
r_acctheta_inbound_outer_all=cell2mat(r_acctheta(~cellfun('isempty',r_acctheta)));
r_stepstheta_inbound_outer_all=cell2mat(r_stepstheta(~cellfun('isempty',r_stepstheta)));
r_stepsvel_inbound_outer_all=cell2mat(r_stepsvel(~cellfun('isempty',r_stepsvel)));
    
p_veltheta_inbound_outer_all=cell2mat(p_veltheta(~cellfun('isempty',p_veltheta)));
p_acctheta_inbound_outer_all=cell2mat(p_acctheta(~cellfun('isempty',p_acctheta)));
p_stepstheta_inbound_outer_all=cell2mat(p_stepstheta(~cellfun('isempty',p_stepstheta)));
p_stepsvel_inbound_outer_all=cell2mat(p_stepsvel(~cellfun('isempty',p_stepsvel)));

%% Outer arms combined - Outbound
clearvars -except f1 f2 savefig_scatter destdir ci_alpha r_* p_*
animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};

for a=[1,2,3,4,5]
    
    count=0;
    animal= animals(a);
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    
    ind_del=[]; ind_del_f1=[]; ind_del_f2=[];
    for i=1:size(ind,1)
        
        if ((isnan(f1(a).output{1, 1}(i).inst_acc_outbound(1)))||(isnan(f1(a).output{1, 1}(i).speed_outbound(1))))
             ind_del_f1=[ind_del_f1; i];
        end
    end

    for i=1:size(ind,1)
        
        if ((isnan(f2(a).output{1, 1}(i).inst_acc_outbound(1)))||(isnan(f2(a).output{1, 1}(i).speed_outbound(1))))
             ind_del_f2=[ind_del_f2; i];
        end
    end
    
    ind_del=unique([ind_del_f1;ind_del_f2], 'rows');
    
    f1(a).output{1, 1}(ind_del)=[];
    f2(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    count=size(ind_select,1);
    
    for e=1:size(ind_select,1)
        % keyboard 
        inst_vel_outbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_vel_outbound,f1(a).output{1}(e),'Un',0))');
        inst_acc_outbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_acc_outbound,f1(a).output{1}(e),'Un',0))');
        inst_freq_outbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f1(a).output{1}(e),'Un',0))');
        inst_steps_outbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f1(a).output{1}(e),'Un',0))');

        inst_vel_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_vel_outbound,f2(a).output{1}(e),'Un',0))');
        inst_acc_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_acc_outbound,f2(a).output{1}(e),'Un',0))');
        inst_freq_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f2(a).output{1}(e),'Un',0))');
        inst_steps_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f2(a).output{1}(e),'Un',0))');

        inst_vel_epoch=([(inst_vel_outbound_f1{a,e}(:)); (inst_vel_outbound_f2{a,e}(:))]); 
        inst_acc_epoch=([(inst_acc_outbound_f1{a,e}(:)); (inst_acc_outbound_f2{a,e}(:))]);
        inst_theta_epoch=([(inst_freq_outbound_f1{a,e}(:)); (inst_freq_outbound_f2{a,e}(:))]);
        inst_steps_epoch=([(inst_steps_outbound_f1{a,e}(:)); (inst_steps_outbound_f2{a,e}(:))]);
        
        idx_include=find((inst_vel_epoch>vel_thresh_min) & (inst_steps_epoch>steps_thresh_min) & (inst_steps_epoch<steps_thresh_max) & (inst_theta_epoch>theta_thresh_min)); % excluding the noisy data 
        idx_include_vel=find((inst_vel_epoch>vel_thresh_min)&(inst_vel_epoch<vel_thresh_max)); % excluding the noisy data 
        
        %% inst vel versus theta 
        [r1, p1]=corrcoef(inst_vel_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r1) 
            r_veltheta{a,e}=r1(2);
            p_veltheta{a,e}=p1(2);
        else 
            r_veltheta{a,e}=[];
            p_veltheta{a,e}=[];
        end
         
        %% inst acc versus theta 
        [r2, p2]=corrcoef(inst_acc_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r2) 
            r_acctheta{a,e}=r2(2);
            p_acctheta{a,e}=p2(2);
        else 
            r_acctheta{a,e}=[];
            p_acctheta{a,e}=[];
        end
        
         %% inst steps versus inst theta 
        [r3, p3]=corrcoef(inst_steps_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');

        if ~isnan(r3) 
            r_stepstheta{a,e}=r3(2);
            p_stepstheta{a,e}=p3(2);
        else 
            r_stepstheta{a,e}=[];
            p_stepstheta{a,e}=[];
        end
        
        %% inst steps versus inst theta 
         [r4, p4]=corrcoef(inst_vel_epoch(idx_include_vel),inst_steps_epoch(idx_include_vel), 'rows', 'complete');
         
        if ~isnan(r4) 
            r_stepsvel{a,e}=r4(2);
            p_stepsvel{a,e}=p4(2);
        else 
            r_stepsvel{a,e}=[];
            p_stepsvel{a,e}=[];
        end
        %%
        clear r1 r2 r3 r4 p1 p2 p3 p4 inst_vel_epoch inst_acc_epoch inst_theta_epoch inst_steps_epoch idx_include idx_include_vel
    end
    
    counter(a)=count;
    
    %save per animal per epoch results 
    r_veltheta_outbound_outer{a}=(r_veltheta(a,:));
    r_acctheta_outbound_outer{a}=(r_acctheta(a,:));
    r_stepstheta_outbound_outer{a}=(r_stepstheta(a,:));
    r_stepsvel_outbound_outer{a}=(r_stepsvel(a,:));
    
    p_veltheta_outbound_outer{a}=(p_veltheta(a,:));
    p_acctheta_outbound_outer{a}=(p_acctheta(a,:));
    p_stepstheta_outbound_outer{a}=(p_stepstheta(a,:));
    p_stepsvel_outbound_outer{a}=(p_stepsvel(a,:));
end

r_veltheta_outbound_outer_all=cell2mat(r_veltheta(~cellfun('isempty',r_veltheta)));
r_acctheta_outbound_outer_all=cell2mat(r_acctheta(~cellfun('isempty',r_acctheta)));
r_stepstheta_outbound_outer_all=cell2mat(r_stepstheta(~cellfun('isempty',r_stepstheta)));
r_stepsvel_outbound_outer_all=cell2mat(r_stepsvel(~cellfun('isempty',r_stepsvel)));
    
p_veltheta_outbound_outer_all=cell2mat(p_veltheta(~cellfun('isempty',p_veltheta)));
p_acctheta_outbound_outer_all=cell2mat(p_acctheta(~cellfun('isempty',p_acctheta)));
p_stepstheta_outbound_outer_all=cell2mat(p_stepstheta(~cellfun('isempty',p_stepstheta)));
p_stepsvel_outbound_outer_all=cell2mat(p_stepsvel(~cellfun('isempty',p_stepsvel)));

%% T Junction - Outbound
clearvars -except savefig_scatter destdir ci_alpha r_* p_*

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
        
        if ((isnan(f1(a).output{1, 1}(i).inst_acc_outbound(1)))||(isnan(f1(a).output{1, 1}(i).speed_outbound(1))))
             ind_del_f1=[ind_del_f1; i];
        end
    end

    for i=1:size(ind,1)
        
        if ((isnan(f2(a).output{1, 1}(i).inst_acc_outbound(1)))||(isnan(f2(a).output{1, 1}(i).speed_outbound(1))))
             ind_del_f2=[ind_del_f2; i];
        end
    end
    
    ind_del=unique([ind_del_f1;ind_del_f2], 'rows');
    
    f1(a).output{1, 1}(ind_del)=[];
    f2(a).output{1, 1}(ind_del)=[];
    
    ind = cell2mat(arrayfun(@(x) x.index',f1(a).output{1},'Un',0))';
    ind_select=(ind(:,1)==1|ind(:,2)==2| ind(:,1)==3|ind(:,1)==4 |ind(:,1)==5 |ind(:,1)==6); %wtrack
    count=size(ind_select,1);
    
    for e=1:size(ind_select,1)
        % keyboard 
        inst_vel_outbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_vel_outbound,f1(a).output{1}(e),'Un',0))');
        inst_acc_outbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_acc_outbound,f1(a).output{1}(e),'Un',0))');
        inst_freq_outbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f1(a).output{1}(e),'Un',0))');
        inst_steps_outbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f1(a).output{1}(e),'Un',0))');

        inst_vel_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_vel_outbound,f2(a).output{1}(e),'Un',0))');
        inst_acc_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_acc_outbound,f2(a).output{1}(e),'Un',0))');
        inst_freq_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_outbound,f2(a).output{1}(e),'Un',0))');
        inst_steps_outbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_outbound,f2(a).output{1}(e),'Un',0))');

        inst_vel_epoch=([(inst_vel_outbound_f1{a,e}(:)); (inst_vel_outbound_f2{a,e}(:))]); 
        inst_acc_epoch=([(inst_acc_outbound_f1{a,e}(:)); (inst_acc_outbound_f2{a,e}(:))]);
        inst_theta_epoch=([(inst_freq_outbound_f1{a,e}(:)); (inst_freq_outbound_f2{a,e}(:))]);
        inst_steps_epoch=([(inst_steps_outbound_f1{a,e}(:)); (inst_steps_outbound_f2{a,e}(:))]);
        
        idx_include=find((inst_vel_epoch>vel_thresh_min) & (inst_steps_epoch>steps_thresh_min) & (inst_steps_epoch<steps_thresh_max) & (inst_theta_epoch>theta_thresh_min)); % excluding the noisy data 
        idx_include_vel=find((inst_vel_epoch>vel_thresh_min)&(inst_vel_epoch<vel_thresh_max)); % excluding the noisy data 
        
        %% inst vel versus theta 
        [r1, p1]=corrcoef(inst_vel_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r1) 
            r_veltheta{a,e}=r1(2);
            p_veltheta{a,e}=p1(2);
        else 
            r_veltheta{a,e}=[];
            p_veltheta{a,e}=[];
        end
         
        %% inst acc versus theta 
        [r2, p2]=corrcoef(inst_acc_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r2) 
            r_acctheta{a,e}=r2(2);
            p_acctheta{a,e}=p2(2);
        else 
            r_acctheta{a,e}=[];
            p_acctheta{a,e}=[];
        end
        
         %% inst steps versus inst theta 
        [r3, p3]=corrcoef(inst_steps_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');

        if ~isnan(r3) 
            r_stepstheta{a,e}=r3(2);
            p_stepstheta{a,e}=p3(2);
        else 
            r_stepstheta{a,e}=[];
            p_stepstheta{a,e}=[];
        end
        
        %% inst steps versus inst theta 
         [r4, p4]=corrcoef(inst_vel_epoch(idx_include_vel),inst_steps_epoch(idx_include_vel), 'rows', 'complete');
         
        if ~isnan(r4) 
            r_stepsvel{a,e}=r4(2);
            p_stepsvel{a,e}=p4(2);
        else 
            r_stepsvel{a,e}=[];
            p_stepsvel{a,e}=[];
        end
        %%
        clear r1 r2 r3 r4 p1 p2 p3 p4 inst_vel_epoch inst_acc_epoch inst_theta_epoch inst_steps_epoch idx_include idx_include_vel
    end
    
    counter(a)=count;
    
    %save per animal per epoch results 
    r_veltheta_outbound_tjunc{a}=(r_veltheta(a,:));
    r_acctheta_outbound_tjunc{a}=(r_acctheta(a,:));
    r_stepstheta_outbound_tjunc{a}=(r_stepstheta(a,:));
    r_stepsvel_outbound_tjunc{a}=(r_stepsvel(a,:));
    
    p_veltheta_outbound_tjunc{a}=(p_veltheta(a,:));
    p_acctheta_outbound_tjunc{a}=(p_acctheta(a,:));
    p_stepstheta_outbound_tjunc{a}=(p_stepstheta(a,:));
    p_stepsvel_outbound_tjunc{a}=(p_stepsvel(a,:));
end

r_veltheta_outbound_tjunc_all=cell2mat(r_veltheta(~cellfun('isempty',r_veltheta)));
r_acctheta_outbound_tjunc_all=cell2mat(r_acctheta(~cellfun('isempty',r_acctheta)));
r_stepstheta_outbound_tjunc_all=cell2mat(r_stepstheta(~cellfun('isempty',r_stepstheta)));
r_stepsvel_outbound_tjunc_all=cell2mat(r_stepsvel(~cellfun('isempty',r_stepsvel)));
    
p_veltheta_outbound_tjunc_all=cell2mat(p_veltheta(~cellfun('isempty',p_veltheta)));
p_acctheta_outbound_tjunc_all=cell2mat(p_acctheta(~cellfun('isempty',p_acctheta)));
p_stepstheta_outbound_tjunc_all=cell2mat(p_stepstheta(~cellfun('isempty',p_stepstheta)));
p_stepsvel_outbound_tjunc_all=cell2mat(p_stepsvel(~cellfun('isempty',p_stepsvel)));

%% T Junction - Inbound
clearvars -except f1 f2 savefig_scatter destdir ci_alpha r_* p_*
animals={'Jaq', 'Roqui', 'Lotus', 'Monty', 'Peanut'};

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
    
    for e=1:size(ind_select,1)
        % keyboard 
        inst_vel_inbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_vel_inbound,f1(a).output{1}(e),'Un',0))');
        inst_acc_inbound_f1{a,e}= cell2mat((arrayfun(@(x) x.inst_acc_inbound,f1(a).output{1}(e),'Un',0))');
        inst_freq_inbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f1(a).output{1}(e),'Un',0))');
        inst_steps_inbound_f1{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f1(a).output{1}(e),'Un',0))');

        inst_vel_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_vel_inbound,f2(a).output{1}(e),'Un',0))');
        inst_acc_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_acc_inbound,f2(a).output{1}(e),'Un',0))');
        inst_freq_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_theta_inbound,f2(a).output{1}(e),'Un',0))');
        inst_steps_inbound_f2{a,e}=cell2mat((arrayfun(@(x) x.inst_steps_inbound,f2(a).output{1}(e),'Un',0))');

        inst_vel_epoch=([(inst_vel_inbound_f1{a,e}(:)); (inst_vel_inbound_f2{a,e}(:))]); 
        inst_acc_epoch=([(inst_acc_inbound_f1{a,e}(:)); (inst_acc_inbound_f2{a,e}(:))]);
        inst_theta_epoch=([(inst_freq_inbound_f1{a,e}(:)); (inst_freq_inbound_f2{a,e}(:))]);
        inst_steps_epoch=([(inst_steps_inbound_f1{a,e}(:)); (inst_steps_inbound_f2{a,e}(:))]);
        
        idx_include=find((inst_vel_epoch>vel_thresh_min) & (inst_steps_epoch>steps_thresh_min) & (inst_steps_epoch<steps_thresh_max) & (inst_theta_epoch>theta_thresh_min)); % excluding the noisy data 
        idx_include_vel=find((inst_vel_epoch>vel_thresh_min)&(inst_vel_epoch<vel_thresh_max)); % excluding the noisy data 
        
        %% inst vel versus theta 
        [r1, p1]=corrcoef(inst_vel_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r1) 
            r_veltheta{a,e}=r1(2);
            p_veltheta{a,e}=p1(2);
        else 
            r_veltheta{a,e}=[];
            p_veltheta{a,e}=[];
        end
         
        %% inst acc versus theta 
        [r2, p2]=corrcoef(inst_acc_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');
        
        if ~isnan(r2) 
            r_acctheta{a,e}=r2(2);
            p_acctheta{a,e}=p2(2);
        else 
            r_acctheta{a,e}=[];
            p_acctheta{a,e}=[];
        end
        
         %% inst steps versus inst theta 
        [r3, p3]=corrcoef(inst_steps_epoch(idx_include),inst_theta_epoch(idx_include), 'rows', 'complete');

        if ~isnan(r3) 
            r_stepstheta{a,e}=r3(2);
            p_stepstheta{a,e}=p3(2);
        else 
            r_stepstheta{a,e}=[];
            p_stepstheta{a,e}=[];
        end
        
        %% inst steps versus inst theta 
         [r4, p4]=corrcoef(inst_vel_epoch(idx_include_vel),inst_steps_epoch(idx_include_vel), 'rows', 'complete');
         
        if ~isnan(r4) 
            r_stepsvel{a,e}=r4(2);
            p_stepsvel{a,e}=p4(2);
        else 
            r_stepsvel{a,e}=[];
            p_stepsvel{a,e}=[];
        end
        %%
        clear r1 r2 r3 r4 p1 p2 p3 p4 inst_vel_epoch inst_acc_epoch inst_theta_epoch inst_steps_epoch idx_include idx_include_vel
    end
    
    counter(a)=count;
    
    %save per animal per epoch results 
    r_veltheta_inbound_tjunc{a}=(r_veltheta(a,:));
    r_acctheta_inbound_tjunc{a}=(r_acctheta(a,:));
    r_stepstheta_inbound_tjunc{a}=(r_stepstheta(a,:));
    r_stepsvel_inbound_tjunc{a}=(r_stepsvel(a,:));
    
    p_veltheta_inbound_tjunc{a}=(p_veltheta(a,:));
    p_acctheta_inbound_tjunc{a}=(p_acctheta(a,:));
    p_stepstheta_inbound_tjunc{a}=(p_stepstheta(a,:));
    p_stepsvel_inbound_tjunc{a}=(p_stepsvel(a,:));
end

r_veltheta_inbound_tjunc_all=cell2mat(r_veltheta(~cellfun('isempty',r_veltheta)));
r_acctheta_inbound_tjunc_all=cell2mat(r_acctheta(~cellfun('isempty',r_acctheta)));
r_stepstheta_inbound_tjunc_all=cell2mat(r_stepstheta(~cellfun('isempty',r_stepstheta)));
r_stepsvel_inbound_tjunc_all=cell2mat(r_stepsvel(~cellfun('isempty',r_stepsvel)));
    
p_veltheta_inbound_tjunc_all=cell2mat(p_veltheta(~cellfun('isempty',p_veltheta)));
p_acctheta_inbound_tjunc_all=cell2mat(p_acctheta(~cellfun('isempty',p_acctheta)));
p_stepstheta_inbound_tjunc_all=cell2mat(p_stepstheta(~cellfun('isempty',p_stepstheta)));
p_stepsvel_inbound_tjunc_all=cell2mat(p_stepsvel(~cellfun('isempty',p_stepsvel)));

%% COMPARISON ACROSS TASK PORTIONS 
% per epoch - save these

% VEL THETA 
ax(1)=subplot(3,1,1);
data=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all; r_stepstheta_outbound_outer_all; r_stepstheta_inbound_outer_all; r_stepstheta_outbound_tjunc_all; r_stepstheta_inbound_tjunc_all];
id=[ones(size(r_stepstheta_outbound_center_all)); 1+ ones(size(r_stepstheta_inbound_center_all)); 2+ ones(size(r_stepstheta_outbound_outer_all)); 3+ones(size(r_stepstheta_inbound_outer_all)); 4+ ones(size(r_stepstheta_outbound_tjunc_all)); 5+ones(size(r_stepstheta_inbound_tjunc_all))];
boxplot(data,id, 'Notch', 'on'); clear data id 
title('Forelimb Vs Theta')

ax(2)=subplot(3,1,2);
data=[r_veltheta_outbound_center_all; r_veltheta_inbound_center_all; r_veltheta_outbound_outer_all; r_veltheta_inbound_outer_all; r_veltheta_outbound_tjunc_all; r_veltheta_inbound_tjunc_all];
id=[ones(size(r_veltheta_outbound_center_all)); 1+ ones(size(r_veltheta_inbound_center_all)); 2+ ones(size(r_veltheta_outbound_outer_all)); 3+ones(size(r_veltheta_inbound_outer_all)); 4+ ones(size(r_veltheta_outbound_tjunc_all)); 5+ones(size(r_veltheta_inbound_tjunc_all))];
boxplot(data,id, 'Notch', 'on'); clear data id 
title('Vel Vs Theta')

ax(3)=subplot(3,1,3);
data=[r_acctheta_outbound_center_all; r_acctheta_inbound_center_all; r_acctheta_outbound_outer_all; r_acctheta_inbound_outer_all; r_acctheta_outbound_tjunc_all; r_acctheta_inbound_tjunc_all];
id=[ones(size(r_acctheta_outbound_center_all)); 1+ ones(size(r_acctheta_inbound_center_all)); 2+ ones(size(r_acctheta_outbound_outer_all)); 3+ones(size(r_acctheta_inbound_outer_all)); 4+ ones(size(r_acctheta_outbound_tjunc_all)); 5+ones(size(r_acctheta_inbound_tjunc_all))];
boxplot(data,id, 'Notch', 'on'); clear data id 
title('Acc Vs Theta')

xlabel('Portion on WTrack');
ylabel('Correlation Coefficient');
ylim([-0.07 0.6]); xlim([0 7]);
labels=({'C-OUT' 'C-IN' 'O-OUT' 'O-IN' 'T-OUT' 'T-IN' ''});
xticklabels(labels)
linkaxes(ax)
text()

all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);

if savefig_scatter==1
    cd('')
    savefig('.fig'); gcf;
end

%% MULTCOMPARE FOR STATS of comparsion between groups 

% STEPS THETA 
ax(1)=subplot(3,1,1);
data=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all; r_stepstheta_outbound_outer_all; r_stepstheta_inbound_outer_all; r_stepstheta_outbound_tjunc_all; r_stepstheta_inbound_tjunc_all];
id=[ones(size(r_stepstheta_outbound_center_all)); 1+ ones(size(r_stepstheta_inbound_center_all)); 2+ ones(size(r_stepstheta_outbound_outer_all)); 3+ones(size(r_stepstheta_inbound_outer_all)); 4+ ones(size(r_stepstheta_outbound_tjunc_all)); 5+ones(size(r_stepstheta_inbound_tjunc_all))];
[~,~,stats]=kruskalwallis(data,id); clear data id 
multcompare(stats)
if savefig_scatter==1
    cd('')
    savefig(''); gcf;
end

ax(2)=subplot(3,1,2);
data=[r_veltheta_outbound_center_all; r_veltheta_inbound_center_all; r_veltheta_outbound_outer_all; r_veltheta_inbound_outer_all; r_veltheta_outbound_tjunc_all; r_veltheta_inbound_tjunc_all];
id=[ones(size(r_veltheta_outbound_center_all)); 1+ ones(size(r_veltheta_inbound_center_all)); 2+ ones(size(r_veltheta_outbound_outer_all)); 3+ones(size(r_veltheta_inbound_outer_all)); 4+ ones(size(r_veltheta_outbound_tjunc_all)); 5+ones(size(r_veltheta_inbound_tjunc_all))];
[~,~,stats]=kruskalwallis(data,id); clear data id 
multcompare(stats)
if savefig_scatter==1
    cd('')
    savefig(''); gcf;
end

ax(3)=subplot(3,1,3);
data=[r_acctheta_outbound_center_all; r_acctheta_inbound_center_all; r_acctheta_outbound_outer_all; r_acctheta_inbound_outer_all; r_acctheta_outbound_tjunc_all; r_acctheta_inbound_tjunc_all];
id=[ones(size(r_acctheta_outbound_center_all)); 1+ ones(size(r_acctheta_inbound_center_all)); 2+ ones(size(r_acctheta_outbound_outer_all)); 3+ones(size(r_acctheta_inbound_outer_all)); 4+ ones(size(r_acctheta_outbound_tjunc_all)); 5+ones(size(r_acctheta_inbound_tjunc_all))];
[~,~,stats]=kruskalwallis(data,id); clear data id 
multcompare(stats)
if savefig_scatter==1
    cd(')
    savefig(''); gcf;
end

%% PLOT correlation coefficients for OUTBOUND AND INBOUND PORTIONS on the center arm
% Scatter PLOT
% this will plot all the data and lines per animal
close all
for a=[1,2,3,4,5]
  r_stepstheta_outbound_center_plot=cell2mat(r_stepstheta_outbound_center{a}(~cellfun('isempty',r_stepstheta_outbound_center{a})));
  r_stepstheta_inbound_center_plot=cell2mat(r_stepstheta_inbound_center{a}(~cellfun('isempty',r_stepstheta_inbound_center{a})));
  
   data_all=[r_stepstheta_outbound_center_plot'; r_stepstheta_inbound_center_plot'];
   id_all=[(ones(size(r_stepstheta_outbound_center_plot))'); ((1+ ones(size(r_stepstheta_inbound_center_plot)))')];

  % PLOT means per animal 
  data=[mean(r_stepstheta_outbound_center_plot)'; mean(r_stepstheta_inbound_center_plot)'];
  id=[1;2];

  figure(1);
  plot([id],[data], 'k'); hold on 
  [p{a},~]=kruskalwallis(data_all,id_all); % collect the significance per animal 
  figure(1); text(1.5,0.05*a,['Rat', num2str(a), ' p: ' num2str(p{a})]);
  % keyboard
  clear data* id* r_stepstheta_outbound_center_plot r_stepstheta_inbound_center_plot
end

figure(1)
% now scatter all the dots per epoch
% collect correlations with p<0.05
p_select_outbound=find(p_stepstheta_outbound_center_all<0.05);
p_select_inbound=find(p_stepstheta_inbound_center_all<0.05);

% empty circles for r value per epoch 
data_all=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all];
id_all=[ones(size(r_stepstheta_outbound_center_all)); (1+ ones(size(r_stepstheta_inbound_center_all)))];
scatter([id_all],[data_all],'k','jitter','on'); hold on 

% filled circles for significant r values 
data_all_p_select=[r_stepstheta_outbound_center_all(p_select_outbound); r_stepstheta_inbound_center_all(p_select_inbound)];
id_all_p_select=[ones(size(r_stepstheta_outbound_center_all(p_select_outbound))); (1+ ones(size(r_stepstheta_inbound_center_all(p_select_inbound))))];
scatter([id_all_p_select],[data_all_p_select],'k','filled','jitter','on'); hold on 

legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5')
xlabel('Portion on WTrack');
ylabel('Correlation Coefficient');
ylim([-0.25 0.45]); xlim([0 3]);
labels=({'' 'C-OUT' 'C-IN'});
xticklabels(labels)

if savefig_scatter==1
    cd('')
    savefig(''); gcf;
end

%% BOXPLOT STEP THETA

for a=[1,2,3,4,5]
  r_stepstheta_outbound_center_plot=cell2mat(r_stepstheta_outbound_center{a}(~cellfun('isempty',r_stepstheta_outbound_center{a})));
  r_stepstheta_inbound_center_plot=cell2mat(r_stepstheta_inbound_center{a}(~cellfun('isempty',r_stepstheta_inbound_center{a})));
  
   data_all=[r_stepstheta_outbound_center_plot'; r_stepstheta_inbound_center_plot'];
   id_all=[(ones(size(r_stepstheta_outbound_center_plot))'); ((1+ ones(size(r_stepstheta_inbound_center_plot)))')];

  % PLOT means per animal 
  data=[mean(r_stepstheta_outbound_center_plot)'; mean(r_stepstheta_inbound_center_plot)'];
  id=[1;2];

  ax=figure(1);
  plot([id],[data], 'k'); hold on 
  [p{a},~]=kruskalwallis(data_all,id_all); % collect the significance per animal 
  [p_steptheta{a},~]=signrank(r_stepstheta_outbound_center_plot,r_stepstheta_inbound_center_plot); % collect the significance per animal 
  
  figure(1); text(1.5,0.05*a,['Rat', num2str(a), ' p: ' num2str(p{a})]);

  clear data* id* r_stepstheta_outbound_center_plot r_stepstheta_inbound_center_plot
end

    % SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'

    close all; clear all; load('Figure1EF_EDFigure2.mat')

    ax=figure(1);
    % now scatter all the dots per epoch
    data_all=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all];
    id_all=[ones(size(r_stepstheta_outbound_center_all)); (1+ ones(size(r_stepstheta_inbound_center_all)))];
    boxplot([data_all],[id_all],'Notch','on'); hold on; 
    %legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5');
    xlabel('Portion on WTrack');
    ylabel('Correlation Coefficient');
    %ylim([-0.25 0.45]); xlim([0 3]);
    labels=({'C-OUT' 'C-IN'});
    xticklabels(labels)

    all_lines = findobj(ax,'Type','Line');
    arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);

    if savefig_scatter==1
        cd('')
        savefig(''); gcf;
    end

%% BOXPLOT & Scatter
close all
for a=[1,2,3,4,5]
  r_stepstheta_outbound_center_plot=cell2mat(r_stepstheta_outbound_center{a}(~cellfun('isempty',r_stepstheta_outbound_center{a})));
  r_stepstheta_inbound_center_plot=cell2mat(r_stepstheta_inbound_center{a}(~cellfun('isempty',r_stepstheta_inbound_center{a})));
  
   data_all=[r_stepstheta_outbound_center_plot'; r_stepstheta_inbound_center_plot'];
   id_all=[(ones(size(r_stepstheta_outbound_center_plot))'); ((1+ ones(size(r_stepstheta_inbound_center_plot)))')];

  % PLOT means per animal 
  data=[mean(r_stepstheta_outbound_center_plot)'; mean(r_stepstheta_inbound_center_plot)'];
  id=[1;2];

  ax=figure(1);
  plot([id],[data], 'k'); hold on 
  [p{a},~]=kruskalwallis(data_all,id_all); % collect the significance per animal 
  figure(1); text(1.5,0.05*a,['Rat', num2str(a), ' p: ' num2str(p{a})]);

  clear data* id* r_stepstheta_outbound_center_plot r_stepstheta_inbound_center_plot
end

ax=figure(1);
% now scatter all the dots per epoch
data_all=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all];
id_all=[ones(size(r_stepstheta_outbound_center_all)); (1+ ones(size(r_stepstheta_inbound_center_all)))];
boxplot([data_all],[id_all],'Notch','on'); hold on 

% empty circles for r value per epoch 
data_all=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all];
id_all=[ones(size(r_stepstheta_outbound_center_all)); (1+ ones(size(r_stepstheta_inbound_center_all)))];
scatter([id_all],[data_all],'k','jitter','on'); hold on 

% filled circles for significant r values 
data_all_p_select=[r_stepstheta_outbound_center_all(p_select_outbound); r_stepstheta_inbound_center_all(p_select_inbound)];
id_all_p_select=[ones(size(r_stepstheta_outbound_center_all(p_select_outbound))); (1+ ones(size(r_stepstheta_inbound_center_all(p_select_inbound))))];
scatter([id_all_p_select],[data_all_p_select],'k','filled','jitter','on'); hold on 

legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5')
xlabel('Portion on WTrack');
ylabel('Correlation Coefficient');
ylim([-0.25 0.45]); xlim([0 3]);
labels=({'C-OUT' 'C-IN'});
xticklabels(labels)

all_lines = findobj(ax,'Type','Line');
arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);

if savefig_scatter==1
    cd('')
    savefig(''); gcf;
end

keyboard

%% BOXPLOT VEL THETA
close all
for a=[1,2,3,4,5]
  r_veltheta_outbound_center_plot=cell2mat(r_veltheta_outbound_center{a}(~cellfun('isempty',r_veltheta_outbound_center{a})));
  r_veltheta_inbound_center_plot=cell2mat(r_veltheta_inbound_center{a}(~cellfun('isempty',r_veltheta_inbound_center{a})));
  
  data_all=[r_veltheta_outbound_center_plot'; r_veltheta_inbound_center_plot'];
  id_all=[(ones(size(r_veltheta_outbound_center_plot))'); ((1+ ones(size(r_veltheta_inbound_center_plot)))')];

  % PLOT means per animal 
  data=[mean(r_veltheta_outbound_center_plot)'; mean(r_veltheta_inbound_center_plot)'];
  id=[1;2];

  ax=figure(1);
  plot([id],[data], 'k'); hold on 
   [p{a},~]=kruskalwallis(data_all,id_all); % collect the significance per animal 
   [p_vel_ttest_out{a},~]=ttest(r_veltheta_outbound_center_plot,0); % collect the significance per animal 
   [p_vel_ttest_in{a},~]=ttest(r_veltheta_inbound_center_plot,0); % collect the significance per animal 
   [p_vel{a},~]=signrank(r_veltheta_outbound_center_plot,r_veltheta_inbound_center_plot); % collect the significance per animal 
  figure(1); text(1.5,0.05*a,['Rat', num2str(a), ' p: ' num2str(p{a})]);

  clear data* id* r_veltheta_outbound_center_plot r_veltheta_inbound_center_plot
end

    % SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'

    close all; clear all; load('Figure1EF_EDFigure2.mat')

    ax=figure(1);
    % now scatter all the dots per epoch
    data_all=[r_veltheta_outbound_center_all; r_veltheta_inbound_center_all];
    id_all=[ones(size(r_veltheta_outbound_center_all)); (1+ ones(size(r_veltheta_inbound_center_all)))];
    boxplot([data_all],[id_all],'Notch','on'); hold on;
    %legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5');
    xlabel('Portion on WTrack');
    ylabel('Correlation Coefficient');
    %ylim([-0.25 0.45]); xlim([0 3]);
    labels=({'C-OUT' 'C-IN'});
    xticklabels(labels);

    all_lines = findobj(ax,'Type','Line');
    arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);

    if savefig_scatter==1
        cd('')
        savefig(''); gcf;
    end

%% BOXPLOT ACC THETA
close all
for a=[1,2,3,4,5]
  r_acctheta_outbound_center_plot=cell2mat(r_acctheta_outbound_center{a}(~cellfun('isempty',r_acctheta_outbound_center{a})));
  r_acctheta_inbound_center_plot=cell2mat(r_acctheta_inbound_center{a}(~cellfun('isempty',r_acctheta_inbound_center{a})));
  
   data_all=[r_acctheta_outbound_center_plot'; r_acctheta_inbound_center_plot'];
   id_all=[(ones(size(r_acctheta_outbound_center_plot))'); ((1+ ones(size(r_acctheta_inbound_center_plot)))')];

  % PLOT means per animal 
  data=[mean(r_acctheta_outbound_center_plot)'; mean(r_acctheta_inbound_center_plot)'];
  id=[1;2];

  ax=figure(1);
  plot([id],[data], 'k'); hold on 
  [p{a},~]=kruskalwallis(data_all,id_all); % collect the significance per animal 
%   [p_acc_ttest_out{a},~]=ttest(r_acctheta_outbound_center_plot,0); % collect the significance per animal 
%   [p_acc_ttest_in{a},~]=ttest(r_acctheta_inbound_center_plot,0); % collect the significance per animal 
%  [p_acc{a},~]=signrank(r_acctheta_outbound_center_plot,r_acctheta_inbound_center_plot); % collect the significance per animal 
  figure(1); text(1.5,0.05*a,['Rat', num2str(a), ' p: ' num2str(p{a})]);
keyboard
  clear data* id* r_acctheta_outbound_center_plot r_acctheta_inbound_center_plot
end

    % SAVE AND PLOT 
    % The variables above are saved in Figure1EF_EDFigure2.mat
    % Run this section if using 'Figure1EF_EDFigure2.mat'

    close all; clear all; load('Figure1EF_EDFigure2.mat')

    ax=figure(1);
    % now scatter all the dots per epoch
    data_all=[r_acctheta_outbound_center_all; r_acctheta_inbound_center_all];
    id_all=[ones(size(r_acctheta_outbound_center_all)); (1+ ones(size(r_acctheta_inbound_center_all)))];
    boxplot([data_all],[id_all],'Notch','on'); hold on 
    %legend('Rat1', 'Rat2', 'Rat3', 'Rat4', 'Rat5')
    xlabel('Portion on WTrack');
    ylabel('Correlation Coefficient');
    %ylim([-0.25 0.45]); xlim([0 3]);
    labels=({'C-OUT' 'C-IN'});
    xticklabels(labels);

    all_lines = findobj(ax,'Type','Line');
    arrayfun(@(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines);

    if savefig_scatter==1
        cd('')
        savefig(''); gcf;
    end

%%
% steptheta 
clear data
data=[r_stepstheta_outbound_center_all; r_stepstheta_inbound_center_all; r_stepstheta_outbound_outer_all; r_stepstheta_inbound_outer_all; r_stepstheta_outbound_tjunc_all; r_stepstheta_inbound_tjunc_all];
[~,p_steptheta_t{1}]=ttest(r_stepstheta_outbound_center_all,0);
[~,p_steptheta_t{2}]=ttest(r_stepstheta_inbound_center_all,0);
[~,p_steptheta_t{3}]=ttest(r_stepstheta_outbound_outer_all,0);
[~,p_steptheta_t{4}]=ttest(r_stepstheta_inbound_outer_all,0);
[~,p_steptheta_t{5}]=ttest(r_stepstheta_outbound_tjunc_all,0);
[~,p_steptheta_t{6}]=ttest(r_stepstheta_inbound_tjunc_all,0);
clear data 

% veltheta
data=[r_veltheta_outbound_center_all; r_veltheta_inbound_center_all; r_veltheta_outbound_outer_all; r_veltheta_inbound_outer_all; r_veltheta_outbound_tjunc_all; r_veltheta_inbound_tjunc_all];
[~,p_veltheta_t{1}]=ttest(r_veltheta_outbound_center_all,0);
[~,p_veltheta_t{2}]=ttest(r_veltheta_inbound_center_all,0);
[~,p_veltheta_t{3}]=ttest(r_veltheta_outbound_outer_all,0);
[~,p_veltheta_t{4}]=ttest(r_veltheta_inbound_outer_all,0);
[~,p_veltheta_t{5}]=ttest(r_veltheta_outbound_tjunc_all,0);
[~,p_veltheta_t{6}]=ttest(r_veltheta_inbound_tjunc_all,0);
clear data 

% acctheta
data=[r_acctheta_outbound_center_all; r_acctheta_inbound_center_all; r_acctheta_outbound_outer_all; r_acctheta_inbound_outer_all; r_acctheta_outbound_tjunc_all; r_acctheta_inbound_tjunc_all];
[~,p_acctheta_t{1}]=ttest(r_acctheta_outbound_center_all,0);
[~,p_acctheta_t{2}]=ttest(r_acctheta_inbound_center_all,0);
[~,p_acctheta_t{3}]=ttest(r_acctheta_outbound_outer_all,0);
[~,p_acctheta_t{4}]=ttest(r_acctheta_inbound_outer_all,0);
[~,p_acctheta_t{5}]=ttest(r_acctheta_outbound_tjunc_all,0);
[~,p_acctheta_t{6}]=ttest(r_acctheta_inbound_tjunc_all,0);
clear data 
