%% validate_plant_times.m
% -> EDFigure1.mat

clear all; close all; clc;
load('EDFigure1.mat')
 
hand_splay=hand_splay_left_limb; % these are human annotated "hand splaying" events of the left forelimb. 
% Collected these blind by manually going through one video of
% lineartrack data and "grabbing" the frame in DLC that corresponded to the first
% frame where all fingers of the rat firmly spalyed, and the left forelimb does not
% have any forward motion. n=114 plants selected throughout the video. 

mirror_plants=mirror_plants_left_limb; % as a second line of evidence, for the same video we also installed 
% mirror at roughly 45 deg. Collected these frames by looking at the
% mirrow view manually going though all frames and selecting the first frames 
% that corresponded to the rats left forelimb showing touchdown and no
% forward movement. n=66 plants selected throuout the video. 

%% calculate the offset between human annotated frames and algorithm defined frames 

% for every frame identified as a hand_splay 
% collect the closest plant frame 
for n=1:length(hand_splay)
    [diff,index] = min(abs((n_frames_left_limb)-hand_splay(n)));
    if diff<50
        offset_hand_splay_frame(n)=(n_frames_left_limb(index)-hand_splay(n));
        offset_hand_splay_time(n)=(n_frames_left_limb(index)-hand_splay(n)).*(1/est_framerate);
        clear diff index
    else 
        clear diff index
    end
end
    
% for every frame identified as a m
% collect the closest plant frame 
for n=1:length(mirror_plants)
    [diff,index] = min(abs((n_frames_left_limb)-mirror_plants(n)));
    
    if diff<50
        offset_mirror_plants_frame(n)=(n_frames_left_limb(index)-mirror_plants(n));
        offset_mirror_plants_time(n)=(n_frames_left_limb(index)-mirror_plants(n)).*(1/est_framerate);
        clear diff index
    else 
        clear diff index
    end
    
end

% %% plot results frame
% figure(1);
% subplot(1,2,1)
% hist(offset_hand_splay_frame,5);
% xlabel('Frames from Plant - Hand Splay');
% ylabel('Count');
% 
% subplot(1,2,2)
% hist(offset_mirror_plants_frame,5);
% xlabel('Frames from Plant - Mirror Plants');
% ylabel('Count');
% 
% results.median_framediff_handsplay=median(offset_hand_splay_frame); %median offset of handsplay frames from the plant 
% results.iqr_framediff_handsplay=iqr(offset_hand_splay_frame);
% results.median_framediff_mirrorplants=median(offset_mirror_plants_frame); %median offset of mirrorplant frames from the plant 
% results.iqr_framediff_mirrorplants=iqr(offset_mirror_plants_frame);

%% plot results time 
figure(2);
subplot(1,2,1)
hist(offset_hand_splay_time,5);hold on;
plot([0 0],[0 50], '-r'); xlim([-0.04 +0.04]);
xlabel('Time from Plant - Hand Splay');
ylabel('Count');

subplot(1,2,2)
hist(offset_mirror_plants_time,5); hold on;
plot([0 0],[0 50], '-r'); xlim([-0.04 +0.04]);
xlabel('Time from Plant - Mirror Plants');
ylabel('Count');

results.median_timediff_handsplay=median(offset_hand_splay_time); %median offset of handsplay frames from the plant 
results.iqr_timediff_handsplay=iqr(offset_hand_splay_time);
results.median_timediff_mirrorplants=median(offset_mirror_plants_time); %median offset of mirrorplant frames from the plant 
results.iqr_timediff_mirrorplants=iqr(offset_hand_splay_time);