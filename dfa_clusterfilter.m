function [out] = dfa_clusterfilter(index, ~, spikes, ~, ~, ~,~,~, varargin)
   
% dfa_clusterfilter
% retrieve spikes from accepted cells 

if (index(3)<size(spikes{index(1)}{index(2)},2)) && (index(4)<size(spikes{index(1)}{index(2)}{index(3)},2))
    
if (~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})) && (~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data))
    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
else
    spiketimes = [];
end

else 
    spiketimes = [];
end

out.index=index;
out.spiketimes=spiketimes;

end