function [spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, ...
    pcFeatures, pcFeatureIdx, channelPositions] = loadEphysData(ephys_path, savePath, datasetidx)
% JF, Load ephys data (1-indexed)
% ------
% Inputs
% ------
% ephys_path: character array defining the path to your kilosorted output files
% datasetidx: 1 x 1 double vector, only use if you have chronically
% recorded stitched datasets.
% ------
% Outputs
% ------
% spikeTimes_samples: nSpikes × 1 uint64 vector giving each spike time in samples (*not* seconds)
% spikeTemplates: nSpikes × 1 uint32 vector giving the identity of each
%   spike's matched template
% templateWaveforms: nTemplates × nTimePoints × nChannels single matrix of
%   template waveforms for each template and channel
% templateAmplitudes: nSpikes × 1 double vector of the amplitude scaling factor
%   that was applied to the template when extracting that spike
% pcFeatures: nSpikes × nFeaturesPerChannel × nPCFeatures  single
%   matrix giving the PC values for each spike
% pcFeatureIdx: nTemplates × nPCFeatures uint32  matrix specifying which
%   channels contribute to each entry in dim 3 of the pc_features matrix
% channelPositions: nChannels x 2 double matrix, each row gives the x and y
%   coordinates of each channel
%

% load spike templates (= waveforms)
spike_templates_0idx = readNPY([ephys_path, filesep, 'spike_templates.npy']);
spikeTemplates = spike_templates_0idx + 1;


% load spike times
if exist(fullfile(ephys_path, 'spike_times_corrected.npy')) % When running pyKS stitched you need the 'aligned / corrected' spike times
    spikeTimes_samples = double(readNPY([ephys_path, filesep, 'spike_times_corrected.npy']));
    spikeTimes_datasets = double(readNPY([ephys_path, filesep, 'spike_datasets.npy'])) + 1; %  which dataset? (zero-indexed so +1)
else
    spikeTimes_samples = double(readNPY([ephys_path, filesep, 'spike_times.npy']));
    spikeTimes_datasets = ones(size(spikeTimes_samples));
end

templateAmplitudes = double(readNPY([ephys_path, filesep, 'amplitudes.npy'])); % ensure double (KS4 saves as single)


% Load and unwhiten templates
templateWaveforms_whitened = readNPY([ephys_path, filesep, 'templates.npy']);
winv = readNPY([ephys_path, filesep, 'whitening_mat_inv.npy']);
templateWaveforms = zeros(size(templateWaveforms_whitened));
for t = 1:size(templateWaveforms, 1)
    templateWaveforms(t, :, :) = squeeze(templateWaveforms_whitened(t, :, :)) * winv;
end

if exist(fullfile([ephys_path, filesep, 'pc_features.npy']))
    pcFeatures = readNPY([ephys_path, filesep, 'pc_features.npy']);
    pcFeatureIdx = readNPY([ephys_path, filesep, 'pc_feature_ind.npy']) + 1;
else % not computed in early kilosort3 version - the distance and drift metrics (which are based on the PCs) will not be calculated
    pcFeatures = NaN;
    pcFeatureIdx = NaN;
end
channelPositions = readNPY([ephys_path, filesep, 'channel_positions.npy']);
%goodChannels = readNPY([ephys_path filesep  'channel_map.npy']) + 1;


% in any merging/splitting has been done in phy, create the corresponding
% templates & pc_features
if exist(fullfile([ephys_path, filesep, 'spike_clusters.npy']))
    spike_clusters_0idx = readNPY([ephys_path, filesep, 'spike_clusters.npy']); % already manually-curated
    spikeClusters = int32(spike_clusters_0idx) + 1;


    newTemplates = unique(spikeClusters(~ismember(spikeClusters, int32(spikeTemplates))));


    if ~isempty(newTemplates)
        % initialize templates and pc features
        try
        templateWaveforms = [templateWaveforms; zeros(max(newTemplates)-size(templateWaveforms, 1), size(templateWaveforms, 2), size(templateWaveforms, 3))];
        catch
            keyboard
        end
        pcFeatureIdx = [pcFeatureIdx; zeros(max(newTemplates)-size(pcFeatureIdx, 1), size(pcFeatureIdx, 2))];
        for iNewTemplate = newTemplates'
            % find corresponding pre merge/split templates and PCs
            oldTemplates = spikeTemplates(spikeClusters == iNewTemplate);
            if length(unique(oldTemplates)) > 1 % average if merge
                newWaveform = mean(templateWaveforms(unique(oldTemplates), :, :), 1);
                newPcFeatureIdx = mean(pcFeatureIdx(unique(oldTemplates), :), 1);
            else % just take value if split
                newWaveform = templateWaveforms(unique(oldTemplates), :, :);
                newPcFeatureIdx = pcFeatureIdx(unique(oldTemplates), :);
            end
            templateWaveforms(iNewTemplate, :, :) = newWaveform;
            pcFeatureIdx(iNewTemplate, :) = newPcFeatureIdx;
        end
        % check raw waveforms 
        bc.load.checkAndConvertRawWaveforms(savePath, spikeTemplates, spikeClusters)


    else
        spikeClusters = spikeTemplates;
         % check raw waveforms 
        bc.load.checkAndConvertRawWaveforms(savePath, spikeTemplates, spikeClusters)

    end
else
    spikeClusters = spikeTemplates;
     % check raw waveforms 
     bc.load.checkAndConvertRawWaveforms(savePath, spikeTemplates, spikeClusters)
end

%% Only use data set of interest - for unit match
if nargin > 3 && ~isempty(datasetidx) %- for unit match

    spikeTimes_samples = spikeTimes_samples(spikeTimes_datasets == datasetidx);
    spikeTemplates = spikeTemplates(spikeTimes_datasets == datasetidx);
    templateAmplitudes = templateAmplitudes(spikeTimes_datasets == datasetidx);
    pcFeatures = pcFeatures(spikeTimes_datasets == datasetidx, :, :);
end

end
