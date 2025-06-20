function param_complete = checkParameterFields(param)
% JF, Check input structure has all necessary fields + add them with
% default values if not. This is to ensure backcompatibility when any new
% paramaters are introduced. By default, any parameters not already present
% will be set so that the quality metrics are calculated in the same way as
% they were before these new parameters were introduced. (i.e. this does
% not change anything for the user!).
% ------
% Inputs
% ------
% - param 
% ------
% Outputs
% ------
% - param_complete 

%% Name changes
[~, param] = bc.qm.prettify_names([], param); % some names were changed for added clarity

%% Default values for fields
% duplicate spikes
defaultValues.removeDuplicateSpikes = 0;
defaultValues.duplicateSpikeWindow_s = 0.0001;
defaultValues.saveSpikes_withoutDuplicates = 1;
defaultValues.recomputeDuplicateSpikes = 0;

% raw waveforms 
defaultValues.detrendWaveform = 0;
defaultValues.extractRaw = 1;

defaultValues.computeSpatialDecay = 1;

% amplitude 
defaultValues.gain_to_uV = NaN;

% phy saving 
defaultValues.saveAsTSV = 0;
defaultValues.unitType_for_phy = 0;

% separate good from mua in non-somatic
defaultValues.splitGoodAndMua_NonSomatic = 0;

% refactory period violations
defaultValues.hillOrLlobetMethod = 1;

% waveform - noise
defaultValues.normalizeSpDecay = 0;% whether to normalize spatial decay points relative to 
% maximum - this makes the spatrial decay slop calculation more invariant to the 
% spike-sorting algorithm used
defaultValues.spDecayLinFit = 1;
defaultValues.minSpatialDecaySlopeExp = 0.01; % in a.u./um
defaultValues.maxSpatialDecaySlopeExp = 0.1; % in a.u./um
defaultValues.maxScndPeakToTroughRatio_noise = 0.8; % peak must be less than this x the trough 

% waveform - non-somatic
defaultValues.maxMainPeakToTroughRatio_nonSomatic = 0.8; % peak must be less than this x the trough 
defaultValues.minWidthFirstPeak_nonSomatic = 0; % in samples 
defaultValues.minWidthMainTrough_nonSomatic = 0; % in samples
defaultValues.minTroughToPeak2Ratio_nonSomatic = 0; % trough should be min 5 x bigger than 1rst peak to count as non-somatic 
defaultValues.maxPeak1ToPeak2Ratio_nonSomatic = Inf; % if units have an initial peak before the trough,
    % it must be at least firstPeakRatio times larger than the peak after
    % the trough to qualify as a non-somatic unit. 0 means this value is
    % not used.

%% Check for missing fields and add them with default value
[param_complete, missingFields] = bc.qm.addMissingFieldsWithDefault(param, defaultValues);

%% Display result
if ~isempty(missingFields)
    disp('Missing param fields filled in with default values');
    disp(missingFields);
end

end
