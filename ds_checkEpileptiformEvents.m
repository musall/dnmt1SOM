%% define experiments to use
recLabels = {'Animal','ExperimentDate','Location','Expertise','Folder','Probe','Background','Path','Ramppower','useRec','hasGamma'};
recInfo{1} = {'SOM2554','26/01/2022','S1','Naive','SOM2554_20221026','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2554','27/01/2022','S1','Naive','SOM2554_20221027','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2554','28/01/2022','S1','Naive','SOM2554_20221028','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2554','29/01/2022','S1','Naive','SOM2554_20221029','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','26/01/2022','S1','Naive','SOM2563_20221026','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','27/01/2022','S1','Naive','SOM2563_20221027','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','28/01/2022','S1','Naive','SOM2563_20221028','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','29/01/2022','S1','Naive','SOM2563_20221029','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1'};
recInfo{2} = {'SOM2626','08/02/2023','S1','Naive','SOM2626_20230208','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2626','09/02/2023','S1','Naive','SOM2626_20230209','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2626','10/02/2023','S1','Naive','SOM2626_20230210','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2626','11/02/2023','S1','Naive','SOM2626_20230211','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','24/01/2023','S1','Naive','SOM2627_20230124','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','25/01/2023','S1','Naive','SOM2627_20230125','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','26/01/2023','S1','Naive','SOM2627_20230126','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','27/01/2023','S1','Naive','SOM2627_20230127','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1'};

%% basic variables
clear opts
opts.Location = {'V1'};
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = false;
opts.groups = {'SOM-Cre', 'SOM-DNMT1'};
opts.baseDur = 0.5;
opts.postStim = 1;
opts.stepSize = 4;
opts.stimType = 1; %1 for vision, 2 for audio, 3 for tactile
opts.brainRange = [-200, 1500]; %depth range for plotting
opts.verbose = false; %flag to supress some of the text outputs
opts.brainThresh = 1;
opts.groupColors = {[0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]};
opts.optoGenetics = false; %flag to isolate sensory responses during optogenetics
opts.showPlots = false;
opts.verbose = true;
opts.winSize = [0 15];
opts.reload = false;
opts.savePath = 'F:\DNMT1_project\';
opts.eventThresh = 1.5;
opts.eventFreqRange = [1, 3];
opts.loadSC = false;

groupColors = {...
[0.0, 0.5, 0.5],
[1.0, 0.498, 0.314],
[0.42, 0.35, 0.8]
};

%% run over groups
nrGroups = length(opts.groups);
allOut = cell(1, length(opts.groups));
usedRecs = cell(1, length(opts.groups));
failedRecs = cell(1, length(opts.groups));
for iGroups = 1 : nrGroups

    allOut{iGroups} = struct();
    recIdx = strcmpi(recLabels, 'Folder');
    sessionIdx = strcmpi(recLabels, 'sessionCnt');
    animalIdx = strcmpi(recLabels, 'Animal');
    allInfo{iGroups} = recInfo{iGroups};
    
    for iRecs = 1 : size(recInfo{iGroups}, 1)
        try
            out = pC_findOscillatoryEvents_DNMT1(recInfo{iGroups}(iRecs,:), recLabels, opts);
            if ~isfield(out, 'peakMeanPower')
                out.peakMeanPower = [];
                out.peakVals = {};
            end
            animalCnt = find(ismember(recInfo{iGroups}(:, animalIdx), recInfo{iGroups}{iRecs, animalIdx}), 1);
            out.animalEventCnt = repmat(animalCnt, 1 ,length(out.useEvents));
                        
            % Loop over fields in the output struct and append them to allOut
            fields = fieldnames(out);
            for iField = 1:numel(fields)
                field = fields{iField};
                if ~isfield(allOut{iGroups}, field)
                    % Initialize the field as a cell array if it doesn't exist
                    allOut{iGroups}.(field) = cell(size(recInfo{iGroups}, 1), 1);
                end
                % Store the result in the corresponding cell
                allOut{iGroups}.(field){iRecs} = out.(field);
            end
            
            usedRecs{iGroups} = [usedRecs{iGroups}, {recInfo{iGroups}{iRecs,recIdx}}];
        catch ME
            disp(['Aborted: ' ME.message]);
            failedRecs{iGroups} = [failedRecs{iGroups}, {recInfo{iGroups}{iRecs,recIdx}}];
        end
    end
end

%% check peak values
clear temp peakVals meanPower sessionPower sessionEvent peakCounts sessionCnt sessionIDs animalIDs
figure;
minEventDur = 4;
for iGroups = 1 : length(allOut)
    %
    subplot(2,2,1);
    sessionPower{iGroups} = cellfun(@(x) nanmean(single(x)), allOut{iGroups}.peakMeanPower)';
    useIdx = ~isnan(sessionPower{iGroups});
    sessionPower{iGroups} = sessionPower{iGroups}(useIdx);
    peakCounts{iGroups} = cellfun(@length, allOut{iGroups}.peakVals);
    peakCounts{iGroups} = peakCounts{iGroups}(useIdx)';
    
    a = Violin(sessionPower{iGroups}, iGroups, 'ViolinColor', groupColors{iGroups}, 'ViolinAlpha', 0.2, 'bandwidth', 0.6);
    a.EdgeColor = groupColors{iGroups};
    a.BoxColor = 'k';
    a.ScatterPlot.MarkerFaceAlpha = 0.8;
        
    axis square
    xlim([0.25 2.75])
    ylim([0 14])
    ylabel('Mean spectral power (dBmV)')
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabel = opts.groups;
    a.ViolinPlot.LineWidth = 4;
    
    if iGroups == 2
        niceFigure
        
    disp('=====================')
    fprintf('Mean WT Session Power = %.2f +- %.2f\n', mean(sessionPower{1}), sem(sessionPower{1}));
    fprintf('Mean RD10 Session Power = %.2f +- %.2f\n', mean(sessionPower{2}), sem(sessionPower{2}));
    end
    
%% show mean power between 1-3 Hz in all events
    subplot(2,2,2);
    meanPower{iGroups} = cat(2, allOut{iGroups}.peakMeanPower{:});
    eventDurations{iGroups} = cat(1, allOut{iGroups}.eventDurations{:})';
    meanPower{iGroups}(eventDurations{iGroups}<minEventDur) = nan;
    useIdx = ~isnan(meanPower{iGroups});
    meanPower{iGroups} = meanPower{iGroups}(useIdx);
    animalIDs{iGroups} = cat(2, allOut{iGroups}.animalEventCnt{:});
    animalIDs{iGroups} = animalIDs{iGroups}(useIdx);

    a = Violin(meanPower{iGroups}, iGroups, 'ViolinColor', groupColors{iGroups}, 'bandwidth', 0.7);
    a.EdgeColor = groupColors{iGroups};
    a.BoxColor = 'k';
    axis square
    niceFigure
    xlim([0.25 2.75])
    ylim([0 14])
    ylabel('Mean session spectral power (dBmV)')
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabel = opts.groups;
    
    if iGroups == 2
        fprintf('Mean WT Event Power = %.2f +- %.2f in %i events\n', mean(meanPower{1}), sem(meanPower{1}), length(meanPower{1}));
        fprintf('Mean RD10 Event Power = %.2f +- %.2f in %i events\n', mean(meanPower{2}), sem(meanPower{2}), length(meanPower{2}));

    end
        
    %% get peak valuies
    for iRecs = 1 : length(allOut{iGroups}.peakVals)
        if ~isempty(allOut{iGroups}.peakVals{iRecs})
            try
                temp{iGroups}{iRecs} = cat(1, allOut{iGroups}.peakVals{iRecs}{:});
                tempAnimal{iGroups}{iRecs} = repmat(allOut{iGroups}.animalEventCnt{iRecs}(1), length(temp{iGroups}{iRecs}), 1);                
            catch
                temp{iGroups}{iRecs} = cat(2, allOut{iGroups}.peakVals{iRecs}{:})';
                tempAnimal{iGroups}{iRecs} = repmat(allOut{iGroups}.animalEventCnt{iRecs}(1), 1, length(temp{iGroups}{iRecs}));                
            end
        else
            temp{iGroups}{iRecs} = [];
            tempSession{iGroups}{iRecs} = [];
            tempAnimal{iGroups}{iRecs} = [];
        end
    end
    
    %% show duration individual events
    subplot(2,2,3);
    eventDurations{iGroups} = cat(1, allOut{iGroups}.eventDurations{:})';
    eventDurations{iGroups} = eventDurations{iGroups}(useIdx);

    a = Violin(eventDurations{iGroups}, iGroups, 'ViolinColor', groupColors{iGroups}, 'bandwidth', 0.7);
    a.EdgeColor = groupColors{iGroups};
    a.BoxColor = 'k';
    axis square
    niceFigure
    xlim([0.25 2.75])
%     ylim([0 14])
    ylabel('Event durations (s)')
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabel = opts.groups;
    
    if iGroups == 2
        fprintf('Median WT Event duration = %.2f, (%.2f - %.2f)\n', median(eventDurations{1}), min(eventDurations{1}), max(eventDurations{1}));
        fprintf('Median RD10 Event duration = %.2f, (%.2f - %.2f)\n', median(eventDurations{2}), min(eventDurations{2}), max(eventDurations{2}));
    end
    
    %% show peak power of individual events
    subplot(2,2,4);
    [peakVals{iGroups}, rejIdx] = rmoutliers(cat(1, temp{iGroups}{:}), 'percentile', [0 98]);
    peakVals{iGroups} = peakVals{iGroups}(~isnan(peakVals{iGroups}));

    a = Violin(peakVals{iGroups}, iGroups, 'ViolinColor', groupColors{iGroups});
    a.EdgeColor = groupColors{iGroups};
    a.BoxColor = 'k';
    axis square
    niceFigure
    xlim([0.25 2.75])
    ylabel('LFP voltage (mV)')
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabel = opts.groups;
    if iGroups == 2
        nhline(0, '--k');
        fprintf('Mean WT Peak Magnitude = %.2f +- %.2f\n', mean(peakVals{1}), sem(peakVals{1}));
        fprintf('Mean RD10 Peak Magnitude = %.2f +- %.2f\n', mean(peakVals{2}), sem(peakVals{2}));
        disp('=====================')

    end
    ylim([-4 0.5])
end

%% same as before but only within spont period
checkPath = [opts.savePath filesep 'spont_periods'];
allFiles = dir([checkPath filesep '*.mat']);
allFiles = {allFiles(:).name};
eventPower = cell(1,2);
fs = 2500;
for iGroups = 1 : length(allOut)
    
    for iRecs = 1 : length(usedRecs{iGroups})
        
        cIdx = contains(allFiles, usedRecs{iGroups}{iRecs});
        if sum(cIdx) > 0
        cFile = allFiles{cIdx};
        
        load(fullfile(checkPath, cFile));
        spontOn = find(spontMask,1, 'first') / fs;
        spontOff = find(spontMask,1, 'last') / fs;
        
        if sum(spontMask) > 0
        useIdx = allOut{iGroups}.useEvents{iRecs} > spontOn & allOut{iGroups}.useEvents{iRecs} < spontOff;
        useIdx = useIdx & allOut{iGroups}.eventDurations{iRecs} > 5.5;

        if ~isempty(useIdx)
            eventPower{iGroups} = [eventPower{iGroups} allOut{iGroups}.peakMeanPower{iRecs}(useIdx)];
        end
        end
        end
    end
end
