function EEG = c_TMSEEG_Preprocess_AARATEPPipeline(varargin)
% c_TMSEEG_Preprocess_AARATEPPipeline - Preprocess raw TMS-EEG data
% 
% Inputs:
%   EEG: input dataset in EEGLab struct format, assumed to not yet have been epoched
%   'pulseEvent': string indicating event type marking each pulse 
%   'outputDir': path to a folder for saving output files. Will be created automatically. If already
%                 exists, the existing dir will be moved to a <outputDir>_old# folder to avoid
%                 overwriting old results
%   'epochTimespan': timespan, in seconds, around which to epoch. Duration should be less than the
%                     inter-trial interval. E.g. [-1 2] epochs 1 s before to 2 s after each pulse
%
% Optional inputs;
%   'outputFilePrefix': string to include as a prefix in every output filename
%   'artifactTimespan': timespan, in seconds, within which to replace data with interpolated values
%   'baselineTimespan': timespan, in seconds, to use for per-epoch baseline subtraction
%   'doEarlyArtifactInterpolation': whether to interpolate the artifact timespan before early processing
%   'doPostICAArtifactInterpolation': whether to interpolate the artifact timespan after ICA cleanup
%   'doBaselineCorrection': whether to apply baseline subtraction after downsampling
%   'doHighPassFiltering': whether to apply the high-pass filtering stage
%   'doBasicChannelRejection': whether to run the early bad-channel detection stage
%   'doInitialEyeBlinkRemoval': whether to run the initial eye-blink ICA rejection stage
%   'downsampleTo': frequency, in Hz, to downsample to before most preprocessing stages. Note that
%                    it is critical to pass in data before downsampling to do some artifact 
%                    interpolation prior to this downsampling stage.
%   'bandpassFreqSpan': frequency span, in Hz, for lower and upper bandpass filter cutoffs.
%   'badChannelDetectionMethod': string, method for bad channel detection
%   'badChannelThreshold': scalar, threshold at which to classify a channel as bad or good. Exact
%                           meaning of this threshold depends on the badChannelDetectionMethod
%   'SOUNDlambda': lambda regularization parameter used in SOUND
%   'initialEyeComponentThreshold': threshold for rejecting eye-related ICs during early eye-blink 
%						rejection stage. Set to 1 or greater to skip early eye IC rejection stage.
%   'leadFieldPath': path to lead field file or lead field matrix itself to be used by SOUND. If 
%                     empty, SOUND will use a template lead field
%   'doDecayRemoval': whether to perform decay artifact fitting/removal stage
%   'doDecayRemovalPerTrial': whether to do decay fitting and removal partially on a per-trial basis
%                              (true) or to subtract the exact same decay timecourse from every 
%                              trial (false)
%   'postNotchCleanupMethod': method to run after notch filtering. One of {'ICA','SSPSIR','BOTH'}
%   'doSSPSIRBeforeICA': legacy flag, retained for backward compatibility
%   'SSPSIRTimeRange': time range, in ms, used by SSP-SIR for artifact PCA (e.g. [0 20])
%   'SSPSIRPC': SSP-SIR PC setting. Can be a scalar (fixed number of PCs) or e.g. {'data', 90}
%   'ICAType': string, method of ICA to be used, as referenced in c_EEG_ICA
%   'stimMuscleComponentThreshold': scalar, threshold over which to classify ICs as TMS-induced 
%                                    muscle artifact
%   'onOverRejection': string, what to do if too many components are rejected.
%   'doDebug': bool, whether to keep intermediate results and generate additional plots
%   'doPlotFinalTimtopo': bool, whether to generate final timtopo plot after preprocessing
%   'maximizePlotsToMonitor': which monitor index on which to show debug plots
%   'plotXLim': timespan, in seconds, for final and debug timtopo plots
%   'plotTPOIs': time-points of interest, in seconds, for final and debug timtopo plots (e.g. 
%                 typical TEP latencies)
%   'plotChans': channels of interest, either as a cell list of labels or vector of channel indices,
%                 to be included in debug erpimage plots.

%% add dependencies to path
persistent pathModified
if isempty(pathModified)
	mfilepath = fileparts(which(mfilename));
	addpath(fullfile(mfilepath, './Common'));
	addpath(fullfile(mfilepath, './Common/EEGAnalysisCode'));
	addpath(fullfile(mfilepath, './Common/GUI'));
	c_EEG_openEEGLabIfNeeded();
	pathModified = true;
end

%% parse inputs
%defaultTPOIs = [15 25 45 60 95 115 180 280]*1e-3;
defaultTPOIs = [25 33 43 55 90 135 250]*1e-3;

p = inputParser();
p.addRequired('EEG', @isstruct);
p.addParameter('pulseEvents', '', @(x) ischar(x) || iscellstr(x));
p.addParameter('outputDir', '', @ischar);
p.addParameter('outputFilePrefix', 'PreprocessedResults', @ischar);
p.addParameter('epochTimespan', [], @c_isSpan);
p.addParameter('artifactTimespan', [-0.002, 0.009], @c_isSpan); %Upper bound was originally 12 ms
p.addParameter('baselineTimespan', [-0.5 -0.01], @c_isSpan);
p.addParameter('doEarlyArtifactInterpolation', true, @islogical);
p.addParameter('doBaselineCorrection', true, @islogical);
p.addParameter('doHighPassFiltering', true, @islogical);
p.addParameter('doBasicChannelRejection', true, @islogical);
p.addParameter('doInitialEyeBlinkRemoval', true, @islogical);
p.addParameter('filterPrePostExtrapolationDurations', [0 0], @(x) isnumeric(x) && length(x)==2);
p.addParameter('downsampleTo', 1000, @isscalar);
p.addParameter('bandpassFreqSpan', [1 200], @c_isSpan);
p.addParameter('badChannelDetectionMethod', {'PREP_deviation', 'TESA_DDWiener_PerTrial'}, @(x) ischar(x) || iscellstr(x));
p.addParameter('badChannelThreshold', 10, @isscalar);
p.addParameter('initialEyeComponentThreshold', 0.9, @isscalar);
p.addParameter('SOUNDlambda', 10^-1.5, @isscalar);
p.addParameter('leadFieldPath', '', @(x) ischar(x) || ismatrix(x));  % used by SOUND
p.addParameter('doDecayRemoval', true, @islogical);
p.addParameter('doDecayRemovalPerTrial', true, @islogical);
p.addParameter('lineNoiseFreq', 60, @isscalar);  % in Hz
p.addParameter('lineNoiseNumHarmonics', 1, @isscalar);
p.addParameter('postNotchCleanupMethod', 'ICA', @(x) ischar(x) || isstring(x));
p.addParameter('doSSPSIRBeforeICA', true, @islogical);
p.addParameter('SSPSIRTimeRange', [0 20], @(x) isnumeric(x) && length(x)==2);
p.addParameter('SSPSIRPC', {'data', 90}, @(x) isnumeric(x) || iscell(x));
p.addParameter('ICAType', 'fastica', @ischar);
p.addParameter('brainComponentThreshold', 0.3, @isscalar);
p.addParameter('stimMuscleComponentThreshold', 8, @isscalar);
p.addParameter('onOverRejection', 'pause', @ischar);
p.addParameter('doPostICAArtifactInterpolation', false, @islogical);
p.addParameter('doDebug', false, @islogical);
p.addParameter('doPlotFinalTimtopo', true, @islogical);
p.addParameter('maximizePlotsToMonitor', 2, @isschar);
p.addParameter('plotXLim', [-0.1 0.35], @c_isSpan);
p.addParameter('plotTPOIs', defaultTPOIs, @isnumeric);
p.addParameter('plotChans', {'Fz', 'Pz'}, @(x) iscellstr(x) || isnumeric(x));
p.parse(varargin{:});
s = p.Results;
EEG = s.EEG;

s.postNotchCleanupMethod = upper(strtrim(char(s.postNotchCleanupMethod)));
if strcmp(s.postNotchCleanupMethod, 'SSPSIR+ICA') || strcmp(s.postNotchCleanupMethod, 'SSPSIR_ICA') || strcmp(s.postNotchCleanupMethod, 'SSPSIRANDICA') || strcmp(s.postNotchCleanupMethod, 'BOTH')
	s.postNotchCleanupMethod = 'BOTH';
end
if ~ismember(s.postNotchCleanupMethod, {'ICA', 'SSPSIR', 'BOTH'})
	error('postNotchCleanupMethod must be one of {''ICA'',''SSPSIR'',''BOTH''}');
end

% Backward compatibility with legacy flag if new parameter was not explicitly passed.
if ismember('postNotchCleanupMethod', p.UsingDefaults) && s.doSSPSIRBeforeICA
	s.postNotchCleanupMethod = 'SSPSIR';
end

md = struct();
md.pipelineVersion = '2.1.1';
md.postNotchCleanupMethod = s.postNotchCleanupMethod;

% check for required named arguments
assert(~isempty(s.pulseEvents), 'Pulse event(s) must be specified');
assert(~isempty(s.outputDir),'Output directory must be specified');
assert(~isempty(s.epochTimespan), 'Epoch timespan must be specified');

if s.doDebug
    intermediateEEGs = {};
    intermediateLabels = {};
end

if ischar(s.pulseEvents)
	s.pulseEvents = {s.pulseEvents};
end

%% make sure no pulses are unexpectedly close together
EEG = c_TMSEEG_handleBurstEvents(EEG,...
	'pulseEvents', s.pulseEvents,...
	'method', 'error',...
	'burstMaxIPI', max(abs(s.epochTimespan)));

%% make output directory
assert(~isempty(s.outputDir),'Output directory must be specified');

if c_exist(s.outputDir,'dir')
	% if outputDir already exists, assume it has old results that we should
	%  move somewhere else rather than deleting / ovewriting
	listing = dir([s.outputDir '_old*']);
	existingDirs = arrayfun(@(listItem) fullfile(listItem.folder, listItem.name), listing, 'UniformOutput',false);
	movePrevOutputTo = c_str_makeUnique(existingDirs, [s.outputDir '_old']);
	movefile(s.outputDir, movePrevOutputTo);
	pause(1); % give time for move to process before making dir again below
	c_saySingle('Moved previous output to %s',...
		c_path_convert(movePrevOutputTo, 'makeRelativeTo', fileparts(s.outputDir)));
end

if ~c_exist(s.outputDir,'dir')
	c_saySingle('Making output directory at %s', s.outputDir);
	mkdir(s.outputDir);
end

%% epoch
c_say('Epoching');
fn = @() pop_epoch(EEG, s.pulseEvents,s.epochTimespan);
if true
	[~,EEG] = evalc('fn()');
else
	EEG = fn();
end
c_sayDone('Extracted %d epochs', EEG.trials);

%     if s.doDebug
%         intermediateEEGs{end+1} = EEG;
%         intermediateLabels{end+1} = 'Epoched';
%     end

%% interpolate artifact timespan
if s.doEarlyArtifactInterpolation
	c_say('Interpolating artifact timespan');
	EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
		'timespanToReplace', s.artifactTimespan,...
		'method', 'ARExtrapolation',...
		'prePostFitDurations', [20 20]*1e-3);
	c_sayDone();
else
	c_saySingle('Skipping artifact interpolation stage.');
end

%% downsample
if EEG.srate > s.downsampleTo
	c_say('Downsampling to %.1f Hz', s.downsampleTo);
	EEG = pop_resample(EEG, s.downsampleTo);
	c_sayDone();

%         if s.doDebug
%             intermediateEEGs{end+1} = EEG;
%             intermediateLabels{end+1} = 'Downsampled';
%         end
end

%% baseline correct
if s.doBaselineCorrection
	baselineTimespan = s.baselineTimespan;
	baselineTimespan(2) = min(baselineTimespan(2), s.artifactTimespan(1));
	tIndices = EEG.times >= baselineTimespan(1)*1e3 & EEG.times < baselineTimespan(2)*1e3;
	toSubtract = nanmean(EEG.data(:, tIndices, :), 2);
	EEG.data = EEG.data - toSubtract;

	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		intermediateLabels{end+1} = 'Baseline subtracted';
	end
else
	c_saySingle('Skipping baseline correction stage.');
end

%% High-pass filtering
if s.doHighPassFiltering
	c_say('Highpass filtering')
	timeToExtend = 0.5;
	maxTimeToExtend = min(abs(s.epochTimespan - s.artifactTimespan*3));
	timeToExtend = min(timeToExtend, maxTimeToExtend);
	EEG = c_TMSEEG_applyModifiedBandpassFilter(EEG,...
		'piecewiseTimeToExtend', timeToExtend,...
		'lowCutoff', s.bandpassFreqSpan(1),...
		'artifactTimespan', s.artifactTimespan*3,...
		'prePostExtrapolationDurations', s.filterPrePostExtrapolationDurations);
	c_sayDone();
	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		intermediateLabels{end+1} = 'High-pass filtered';
	end
else
	c_saySingle('Skipping high-pass filtering stage.');
end

%% basic channel rejection
if s.doBasicChannelRejection
	doReplaceBadChanImmediately = true;

	[EEG, misc] = c_TMSEEG_detectBadChannels(EEG,...
		'detectionMethod', s.badChannelDetectionMethod,...
		'artifactTimespan', s.artifactTimespan*2,...
		'threshold', s.badChannelThreshold,...
		'replaceMethod', c_if(doReplaceBadChanImmediately, 'interpolate', 'none'),...
		'doPlot', true);

	badChannels = misc.badChannelIndices;
	rejectedChannelsLegend = getRejectedChannelsLegendText(EEG, badChannels);
	addRejectedChannelsTextNote(misc.hf, rejectedChannelsLegend, 'center');
	figure(misc.hf);
	drawnow;

	c_FigurePrinter.copyToFile(fullfile(s.outputDir, [s.outputFilePrefix '_QC_EarlyChannelRejection.png']),...
		'magnification', 2, 'resolution', 450, 'doTransparent', false, 'hf', misc.hf);
	close(misc.hf);

	md.earlyRejectedChannels = badChannels;

	if s.doDebug
		if ~isempty(badChannels)
			if doReplaceBadChanImmediately
				intermediateEEGs{end+1} = EEG;
			else
				intermediateEEGs{end+1} = pop_select(EEG, 'nochannel', badChannels);
			end
			intermediateLabels{end+1} = 'After channel rejection';
		else
			intermediateEEGs{end+1} = EEG;
			intermediateLabels{end+1} = '(After channel rejection)';
		end
		tmpEEG = [];
	end
else
	c_saySingle('Skipping basic channel rejection stage.');
	badChannels = [];
	md.earlyRejectedChannels = badChannels;
end

%% Initial eye-blink removal
% Note: do this prior to SOUND, since SOUND tends to make eye blink topographies look more "brain-like" in a way
%  that makes them more difficult to reject later.

if s.doInitialEyeBlinkRemoval && s.initialEyeComponentThreshold < 1
	c_say('Rereferencing');
	EEG = pop_reref(EEG, []);
	c_sayDone();
	
	c_say('Running early ICA for eye artifacts');
	EEG = c_EEG_ICA(EEG, 'method', 'fastica');
	c_sayDone();
	
	c_say('Labeling ICs using ICLabel');
	[EEG, misc] = c_TMSEEG_runICLabel(EEG,...
		'eyeComponentThreshold', s.initialEyeComponentThreshold,...
		'muscleComponentThreshold', nan,...
		'brainComponentThreshold', nan,...
		'otherComponentThreshold', nan,...
		'doPlot', true,...
		'doRejection', false);
	c_sayDone();

	md.eyeICA_numComp = length(misc.rejectComponents);
	md.eyeICA_numRejComp = sum(misc.rejectComponents);

	if true
		if isfield(misc, 'hf') && ~isempty(misc.hf) && isgraphics(misc.hf, 'figure')
			hf = misc.hf;
		else
			hf = figure;
			misc.hf = hf;
		end
		figure(hf);
		adjustICAComponentFigureLayout(hf);
		c_FigurePrinter.copyToFile(fullfile(s.outputDir, [s.outputFilePrefix, '_QC_EyeICA_ClassifiedComponents.png']),...
			'magnification', 2,...
			'resolution', 350,...
			'doTransparent', false);
		close(hf);
	end

	if  md.eyeICA_numRejComp > 0
		c_say('Rejecting %d/%d components', md.eyeICA_numRejComp, md.eyeICA_numComp);
		if false
			EEG = pop_subcomp(EEG, find(misc.rejectComponents));
		else
			compproj_toRemove = EEG.icawinv(:, misc.rejectComponents) * eeg_getdatact(EEG, 'component', find(misc.rejectComponents), 'reshape', '2d');
			compproj_toRemove = reshape(compproj_toRemove, size(compproj_toRemove, 1), EEG.pnts, EEG.trials);
			if false
				EEG_a = pop_subcomp(EEG, find(misc.rejectComponents));
				EEG_b = EEG;
				EEG_b.data = EEG_b.data - compproj_toRemove;
				c_EEG_plotRawComparison({EEG_a, EEG_b}, 'descriptors', {'pop_subcomp', 'custom subcomp'});
				EEG_c = EEG;
				EEG_c.data = EEG_b.data - EEG_a.data;
				pop_eegplot(EEG_c, 1, 1, 0);
			end
			EEG.data = EEG.data - compproj_toRemove;
		end
		c_sayDone();
	end
	
	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		intermediateLabels{end+1} = sprintf('Eye components%s removed', c_if(md.eyeICA_numRejComp > 0, '', ' (not)'));
	end
else
	c_saySingle('Skipping early eye component rejection stage.')
	md.eyeICA_numComp = 0;
	md.eyeICA_numRejComp = 0;
end


if true
	% save EEG prior to SOUND for some analyses
	outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_preSOUND.mat']);
	c_say('Prior to ICA rejection, saving results to %s', outputPath);
	save(outputPath, 'EEG', 'md');
	c_sayDone();
end

%% SOUND

c_say('Running SOUND');
EEG = c_TMSEEG_runSOUND(EEG,...
	'doRereferenceBeforeSOUND', false,...
	'replaceChannels', badChannels,...
	'lambda', s.SOUNDlambda,...
	'leadFieldPath', s.leadFieldPath);
c_sayDone();

if s.doDebug
	intermediateEEGs{end+1} = EEG;
	intermediateLabels{end+1} = 'After SOUND';
end

if true
    % save EEG prior to component rejection to allow different component rejection choices later
    outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_preDecayRemoval.mat']);
    c_say('Prior to ICA rejection, saving results to %s', outputPath);
    c_save(outputPath, 'EEG', 'md');
    c_sayDone();
end

%% decay fitting and removal
if s.doDecayRemoval
	[EEG, misc] = c_TMSEEG_fitAndRemoveDecayArtifact(EEG,...
		'artifactTimespan', s.artifactTimespan,...
		'trialAggMethod_timeCourseRemoval', c_if(s.doDecayRemovalPerTrial, 'none', 'mean'),...
		'aggTrimPercent', 10,...
		'doPlot', true);

	figure(misc.hf);

	% Make the decay-fit line plot text a bit smaller for readability.
	allAxes = findall(misc.hf, 'type', 'axes');
	if ~isempty(allAxes)
		set(allAxes, 'FontSize', 10, 'TitleFontSizeMultiplier', 0.9, 'LabelFontSizeMultiplier', 0.9);
	end
	allText = findall(misc.hf, 'type', 'text');
	if ~isempty(allText)
		set(allText, 'FontSize', 10);
	end
	allLegends = findall(misc.hf, 'type', 'legend');
	if ~isempty(allLegends)
		set(allLegends, 'FontSize', 10);
	end

	drawnow;
	c_FigurePrinter.copyToFile(fullfile(s.outputDir, [s.outputFilePrefix '_QC_DecayFitAndRemoval.png']), 'magnification', 2, 'resolution', 450, 'doTransparent', false);
	close(misc.hf);

	md.didRemoveDecay = misc.didRemoveDecay;

	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		intermediateLabels{end+1} = 'After decay removal';
	end
else
	c_saySingle('Skipping decay fitting/removal stage.');
	md.didRemoveDecay = false;

	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		intermediateLabels{end+1} = '(Decay removal skipped)';
	end
end

%% interpolate artifact timespan (to remove any large decay residuals in artifactTimespan)
if s.doPostICAArtifactInterpolation
	c_say('Interpolating artifact timespan');
	EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
		'timespanToReplace', s.artifactTimespan,...
		'method', 'ARExtrapolation',...
		'prePostFitDurations', [20 20]*1e-3);
	c_sayDone();
else
	c_saySingle('Skipping artifact interpolation stage.');
end

if s.doDebug
	intermediateEEGs{end+1} = EEG;
	intermediateLabels{end+1} = 'Artifact interpolated';
end

%% notch
c_say('Notch filtering');
if false
	EEG = pop_eegfiltnew(EEG, s.lineNoiseFreq-2, s.lineNoiseFreq+2, 2*EEG.srate, 1);
else
	EEG = c_EEG_filter_butterworth(EEG, s.lineNoiseFreq + [-2, 2], 'type', 'bandstop');
	if true
		% also (optionally) filter higher harmonics
		% even though they will be low-passed later, can still interfere with ICA
		for iHarmonic = 2:s.lineNoiseNumHarmonics
			EEG = c_EEG_filter_butterworth(EEG, [-2, 2] + s.lineNoiseFreq*iHarmonic, 'type', 'bandstop');
		end
	end
end
c_sayDone();

%% SSP-SIR (TESA)
if strcmp(s.postNotchCleanupMethod, 'SSPSIR') || strcmp(s.postNotchCleanupMethod, 'BOTH')
	assert(exist('pop_tesa_sspsir', 'file') == 2,...
		'TESA function pop_tesa_sspsir was not found on the MATLAB path');

	if true
		% save EEG prior to SSP-SIR for alternate cleanup analyses
		outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_preSSPSIR.mat']);
		c_say('Prior to SSP-SIR, saving results to %s', outputPath);
		c_save(outputPath, 'EEG', 'md');
		c_sayDone();
	end

	EEG_beforeSSPSIR = EEG;
	allFigsBeforeSSPSIR = findall(0, 'type', 'figure');

	c_say('Running SSP-SIR');
	EEG = pop_tesa_sspsir(EEG,...
		'artScale', 'manual',...
		'timeRange', s.SSPSIRTimeRange,...
		'PC', s.SSPSIRPC);
	c_sayDone();

	drawnow;
	allFigsAfterSSPSIR = findall(0, 'type', 'figure');
	newSSPSIRFigs = setdiff(allFigsAfterSSPSIR, allFigsBeforeSSPSIR);
	for iFig = 1:length(newSSPSIRFigs)
		if isgraphics(newSSPSIRFigs(iFig), 'figure')
			close(newSSPSIRFigs(iFig));
		end
	end

	hf = figure;
	set(hf, 'Position', [50 50 1700 1100]);
	ht = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
	title(ht, 'SSP-SIR QC: all-channel averages + GMFA before vs after');

	beforeMean = squeeze(mean(EEG_beforeSSPSIR.data, 3));
	afterMean = squeeze(mean(EEG.data, 3));
	numCh = size(beforeMean, 1);
	gmfaBefore = sqrt(mean((beforeMean - mean(beforeMean, 1)).^2, 1));
	gmfaAfter = sqrt(mean((afterMean - mean(afterMean, 1)).^2, 1));

	xCoords = nan(numCh, 1);
	yCoords = nan(numCh, 1);
	for iCh = 1:numCh
		if iCh <= length(EEG.chanlocs) && isfield(EEG.chanlocs, 'X') && isfield(EEG.chanlocs, 'Y') && ...
				~isempty(EEG.chanlocs(iCh).X) && ~isempty(EEG.chanlocs(iCh).Y)
			xCoords(iCh) = EEG.chanlocs(iCh).X;
			yCoords(iCh) = EEG.chanlocs(iCh).Y;
		end
	end

	if any(isnan(xCoords) | isnan(yCoords))
		for iCh = 1:numCh
			if iCh <= length(EEG.chanlocs) && isfield(EEG.chanlocs, 'theta') && isfield(EEG.chanlocs, 'radius') && ...
					~isempty(EEG.chanlocs(iCh).theta) && ~isempty(EEG.chanlocs(iCh).radius)
				% EEGLab convention: x = sin(theta)*radius (left-right), y = cos(theta)*radius (front-back)
				xCoords(iCh) = EEG.chanlocs(iCh).radius * sin(EEG.chanlocs(iCh).theta * pi/180);
				yCoords(iCh) = EEG.chanlocs(iCh).radius * cos(EEG.chanlocs(iCh).theta * pi/180);
			end
		end
	end

	if any(isnan(xCoords) | isnan(yCoords))
		tmpTheta = linspace(0, 2*pi, numCh+1);
		tmpTheta = tmpTheta(1:end-1);
		xCoords = 0.48*cos(tmpTheta(:));
		yCoords = 0.48*sin(tmpTheta(:));
	else
		r = sqrt(xCoords.^2 + yCoords.^2);
		rMax = max(r);
		if rMax > 0
			xCoords = xCoords / rMax * 0.18;
			yCoords = yCoords / rMax * 0.18;
		end
	end

	channelColors = zeros(numCh, 3);
	xMin = min(xCoords);
	xMax = max(xCoords);
	if xMax > xMin
		xNorm = (xCoords - xMin) / (xMax - xMin);
	else
		xNorm = 0.5*ones(numCh, 1);
	end
	yMin = min(yCoords);
	yMax = max(yCoords);
	if yMax > yMin
		yNorm = (yCoords - yMin) / (yMax - yMin);
	else
		yNorm = 0.5*ones(numCh, 1);
	end
	for iCh = 1:numCh
		t = xNorm(iCh);
		if t < 0.5
			u = t / 0.5;
			baseColor = [0, u, 1-u]; % blue -> green
		else
			u = (t - 0.5) / 0.5;
			baseColor = [u, 1-u, 0]; % green -> red
		end
		shadeScale = 1 - 0.35*(1 - yNorm(iCh)); % darker front (high yNorm = front), lighter back
		channelColors(iCh, :) = max(0, min(1, baseColor * shadeScale));
	end

	axBefore = nexttile(1);
	hold(axBefore, 'on');
	for iCh = 1:numCh
		plot(axBefore, EEG_beforeSSPSIR.times, beforeMean(iCh, :), 'LineWidth', 0.6, 'Color', channelColors(iCh, :));
	end
	hGMFABefore = plot(axBefore, EEG_beforeSSPSIR.times, gmfaBefore, 'k', 'LineWidth', 2.2);
	xline(axBefore, 0, 'k--');
	xlabel('Time (ms)');
	ylabel('Amplitude (\muV)');
	title(sprintf('Before SSP-SIR (%d channels)', numCh));
	grid on;

	axAfter = nexttile(2);
	hold(axAfter, 'on');
	for iCh = 1:numCh
		plot(axAfter, EEG.times, afterMean(iCh, :), 'LineWidth', 0.6, 'Color', channelColors(iCh, :));
	end
	hGMFAAfter = plot(axAfter, EEG.times, gmfaAfter, 'k', 'LineWidth', 2.2);
	xline(axAfter, 0, 'k--');
	xlabel('Time (ms)');
	ylabel('Amplitude (\muV)');
	title(sprintf('After SSP-SIR (%d channels)', numCh));
	grid on;

	linkaxes([axBefore axAfter], 'xy');
	xlim(axBefore, [EEG.times(1) EEG.times(end)]);
	xlim(axAfter, [EEG.times(1) EEG.times(end)]);

	maxAbsY = max(abs([beforeMean(:); afterMean(:); gmfaBefore(:); gmfaAfter(:)]));
	if ~isfinite(maxAbsY) || maxAbsY <= 0
		maxAbsY = 1;
	end
	maxAbsY = ceil(maxAbsY / 2) * 2;
	yLimCommon = [-1 1] * maxAbsY;
	yTicksCommon = yLimCommon(1):2:yLimCommon(2);
	set(axBefore, 'YLim', yLimCommon, 'YTick', yTicksCommon);
	set(axAfter, 'YLim', yLimCommon, 'YTick', yTicksCommon);

	legend(axBefore, hGMFABefore, {'GMFA'}, 'Location', 'northwest');
	legend(axAfter, hGMFAAfter, {'GMFA'}, 'Location', 'northwest');

	drawnow;
	addSSPSIRInsetTopoplot(hf, axBefore, xCoords, yCoords, channelColors);
	addSSPSIRInsetTopoplot(hf, axAfter, xCoords, yCoords, channelColors);

	c_FigurePrinter.copyToFile(fullfile(s.outputDir, [s.outputFilePrefix '_QC_SSPSIR.png']),...
		'magnification', 2,...
		'resolution', 350,...
		'doTransparent', false);
	close(hf);

	md.didRunSSPSIR = true;
	md.didRunSecondICA = false;
	md.SSPSIR_timeRange = s.SSPSIRTimeRange;
	md.SSPSIR_PC = s.SSPSIRPC;

	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		intermediateLabels{end+1} = 'After SSP-SIR';
	end
else
	md.didRunSSPSIR = false;
end

%% ICA
if strcmp(s.postNotchCleanupMethod, 'ICA') || strcmp(s.postNotchCleanupMethod, 'BOTH')
	md.didRunSecondICA = true;

	c_say('Running ICA');
	EEG = c_EEG_ICA(EEG, 'method', s.ICAType);
	c_sayDone();

	c_say('Labeling ICs using ICLabel');
	[EEG, misc] = c_TMSEEG_runICLabel(EEG,...
		'brainComponentThreshold', s.brainComponentThreshold,...
		'doPlot', true,...
		'doRejection', false);
	c_sayDone();

	% EEG.reject.gcompreject should have been set to rejectComponents inside c_TMSEEG_runICLabel

	if s.doDebug
		tmpEEG = pop_subcomp(EEG,find(misc.rejectComponents));
		intermediateEEGs{end+1} = tmpEEG;
		intermediateLabels{end+1} = 'After ICA + ICLabel';
	end

	% inspired by tesa_compselect, classify as TMS-evoked muscle if very localized to short latency response
	tmsMuscleWinTimespan = [11 30]*1e-3; % in s
	icaact = eeg_getica(EEG);
	numComps = size(EEG.icaweights, 1);
	tmsMuscleRatio = nan(numComps, 1);
	for iC = 1:numComps
		indicesInWin = EEG.times >= tmsMuscleWinTimespan(1)*1e3 & EEG.times < tmsMuscleWinTimespan(2)*1e3;
		assert(any(indicesInWin));
		muscleScore = abs(mean(icaact(iC, :, :), 3));
		winScore = mean(muscleScore(:, indicesInWin), 2);
		tmsMuscleRatio(iC) = winScore / mean(muscleScore);
	end
	toRejectForMuscle = tmsMuscleRatio > s.stimMuscleComponentThreshold;

	numNewRejected = sum(~misc.rejectComponents & toRejectForMuscle);
	if numNewRejected > 0
		c_saySingle('Marked %d additional components for rejection as TMS-induced muscle artifacts');
		misc.rejectComponents = misc.rejectComponents | toRejectForMuscle;
		EEG.reject.gcompreject = misc.rejectComponents;
		didRejectExtra = true;
		
		% re-plot components after changing rejection labels
		if isfield(misc, 'hf') && ~isempty(misc.hf) && isgraphics(misc.hf, 'figure')
			close(misc.hf);
		end
		pop_viewprops(EEG, 0, 1:size(EEG.icawinv, 2), {}, {}, 0, 'ICLabel');
		misc.hf = gcf;
		if isgraphics(misc.hf, 'figure')
			adjustICAComponentFigureLayout(misc.hf);
		end
	else
		didRejectExtra = false;
	end

	if sum(~misc.rejectComponents) < 4
		% very small number of components remain after rejection
		msg = sprintf('Only %d/%d components remain after rejection', sum(~misc.rejectComponents), length(misc.rejectComponents));
		switch(s.onOverRejection)
			case 'error'
				error('%s', msg);
			case 'warn'
				warning('%s', msg);
			case 'pause'
				c_say('%s. Pausing', msg);
				keyboard
				c_sayDone();
			otherwise
				error('Not implemented');
		end
	end

	md.ICA_numComp = length(misc.rejectComponents);
	md.ICA_numRejComp = sum(misc.rejectComponents);

	if true
		% save EEG prior to component rejection to allow different component rejection choices later
		outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_preICARejection.mat']);
		c_say('Prior to ICA rejection, saving results to %s', outputPath);
		c_save(outputPath, 'EEG', 'md');
		c_sayDone();
	end

	if true
		if isfield(misc, 'hf') && ~isempty(misc.hf) && isgraphics(misc.hf, 'figure')
			hf = misc.hf;
		else
			hf = figure;
			misc.hf = hf;
		end
		figure(hf);
		adjustICAComponentFigureLayout(hf);
		c_FigurePrinter.copyToFile(fullfile(s.outputDir, [s.outputFilePrefix '_QC_ClassifiedComponents.png']), 'magnification', 2, 'resolution', 350, 'doTransparent', false);
		close(hf);
	end

	c_say('Rejecting %d/%d components', md.ICA_numRejComp, md.ICA_numComp);
	EEG = pop_subcomp(EEG,find(misc.rejectComponents));
	c_sayDone();

	if s.doDebug
		intermediateEEGs{end+1} = EEG;
		if didRejectExtra
			intermediateLabels{end+1} = 'After ICA + ICLabel + Extra';
		else
			intermediateLabels{end+1} = 'After ICA + ICLabel (+ Extra)';
		end
	end
else
	md.didRunSecondICA = false;
	md.ICA_numComp = nan;
	md.ICA_numRejComp = nan;
end

%% interpolate again
if s.doPostICAArtifactInterpolation && (strcmp(s.postNotchCleanupMethod, 'ICA') || strcmp(s.postNotchCleanupMethod, 'BOTH'))
	c_say('Interpolating artifact timespan');
	EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
		'timespanToReplace', s.artifactTimespan,...
		'method', 'ARExtrapolation',...
		'prePostFitDurations', [20 20]*1e-3);
	c_sayDone();
end

%% low-pass
c_say('Lowpass filtering');
EEG = c_EEG_filter_butterworth(EEG, [0 s.bandpassFreqSpan(2)]);
c_sayDone();

if s.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'Filtered';
end

%% rereference
c_say('Rereferencing');
EEG = pop_reref(EEG, []);
c_sayDone();

if s.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'Rereferenced';
end

%% save output
outputPath = fullfile(s.outputDir, [s.outputFilePrefix '.mat']);
c_say('Saving results to %s', outputPath);
c_save(outputPath, 'EEG', 'md');
c_sayDone();

%% plot
if s.doPlotFinalTimtopo
	c_say('Plotting final timtopo');
	outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_Plot_TEPTimtopo.png']);
	rejectedChannelsLegend = getRejectedChannelsLegendText(EEG, badChannels);
	hf = figure;
	tmp = c_EEG_plotTimtopo(EEG,...
		'doShowColorbars', false,... 
		'doPlotGMFA', true,...
		'TPOI', s.plotTPOIs,...
		'xlim', s.plotXLim*1e3);
	hf = gcf;
	hf.Position = [50 50 1800 1200];
	set(hf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'PaperPosition', [0 0 18 12]);
	allAxes = findall(hf, 'type', 'axes');
	if ~isempty(allAxes)
		set(allAxes, 'FontSize', 18, 'TitleFontSizeMultiplier', 1.3, 'LabelFontSizeMultiplier', 1.3);
		for iAx = 1:numel(allAxes)
			t = get(allAxes(iAx), 'Title');
			if isgraphics(t, 'text')
				set(t, 'FontSize', 24);
			end
		end
	end
	allText = findall(hf, 'type', 'text');
	if ~isempty(allText)
		allTitleText = gobjects(0);
		for iAx = 1:numel(allAxes)
			t = get(allAxes(iAx), 'Title');
			if isgraphics(t, 'text')
				allTitleText(end+1) = t;
			end
		end
		allText = setdiff(allText, allTitleText);
		if ~isempty(allText)
			set(allText, 'FontSize', 18);
		end
	end
	allLegends = findall(hf, 'type', 'legend');
	if ~isempty(allLegends)
		set(allLegends, 'FontSize', 18);
	end
	shiftFigureAxesLeftForTimtopo(hf, 0.03);
	addRejectedChannelsTextNote(hf, rejectedChannelsLegend, 'left');
	figure(hf);
	drawnow;
	c_FigurePrinter.copyToFile(outputPath,...
		'doCrop', false,...
		'doTransparent', false,...
		'magnification', 2,...
		'resolution', 450,...
		'hf', hf);
	close(hf);
	c_sayDone();
end

%%
if s.doDebug
    plotTEPComparison(s, intermediateEEGs, intermediateLabels);
    
    outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_Debug_IntermediateTEPs.png']);
    hf = gcf;
    c_fig_arrange('maximize', hf, 'monitor', s.maximizePlotsToMonitor);
    c_FigurePrinter.copyToFile(outputPath, 'magnification', 2, 'resolution', 450, 'doTransparent', false);
    close(hf);
	
	if ~isempty(s.plotChans)
		plotERPImComparison(s, intermediateEEGs, intermediateLabels);

		outputPath = fullfile(s.outputDir, [s.outputFilePrefix '_Debug_IntermediateERPIms.png']);
		hf = gcf;
		c_fig_arrange('maximize', hf, 'monitor',  s.maximizePlotsToMonitor);
		c_FigurePrinter.copyToFile(outputPath, 'magnification', 2, 'resolution', 450, 'doTransparent', false);
		close(hf);
	end
end

end

function addSSPSIRInsetTopoplot(hf, haMain, xCoords, yCoords, channelColors)
    if isempty(hf) || ~isgraphics(hf, 'figure') || isempty(haMain) || ~isgraphics(haMain, 'axes')
        return;
    end
    if isempty(xCoords) || isempty(channelColors)
        return;
    end

    % Normalise channel positions to a fixed head radius of 0.18
    rMax = max(sqrt(xCoords.^2 + yCoords.^2));
    if rMax > 0
        xN = xCoords / rMax * 0.18;
        yN = yCoords / rMax * 0.18;
    else
        xN = xCoords;
        yN = yCoords;
    end

    % Size the inset as a pixel-square so the head is circular
    figPos = hf.Position;  % [left bottom width_px height_px]
    insetPx = 120;
    insetW = insetPx / figPos(3);
    insetH = insetPx / figPos(4);

    mainPos = haMain.Position;
    insetLeft   = mainPos(1) + mainPos(3) * 0.02;
    insetBottom = mainPos(2) + mainPos(4) - insetH - mainPos(4) * 0.02;

    haInset = axes('Parent', hf, 'Position', [insetLeft insetBottom insetW insetH]);
    hold(haInset, 'on');

    % Head outline
    headR = 0.21;
    th = linspace(0, 2*pi, 200);
    plot(haInset, headR*cos(th), headR*sin(th), 'k-', 'LineWidth', 0.8);
    % Nose
    plot(haInset, [0 0.025 0 -0.025 0], [headR headR+0.04 headR+0.025 headR+0.04 headR], 'k-', 'LineWidth', 0.7);
    % Ears
    plot(haInset, [ headR  headR+0.025  headR+0.025], [0.04  0.025 -0.025], 'k-', 'LineWidth', 0.7);
    plot(haInset, [-headR -headR-0.025 -headR-0.025], [0.04  0.025 -0.025], 'k-', 'LineWidth', 0.7);

    % Channel dots — same index, same color as the line traces
    for iCh = 1:size(channelColors, 1)
        scatter(haInset, xN(iCh), yN(iCh), 18, channelColors(iCh,:), 'filled', 'MarkerEdgeColor', 'none');
    end

    axis(haInset, 'equal');
    axis(haInset, 'off');
    set(haInset, 'Color', 'none', 'XLim', [-0.27 0.27], 'YLim', [-0.27 0.27]);
	uistack(haInset, 'top');
end

function legendText = getRejectedChannelsLegendText(EEG, badChannels)
	if nargin < 2 || isempty(badChannels)
		legendText = 'Rejected channels (0): none';
		return;
	end

	badChannels = badChannels(:)';
	channelLabels = cell(1, numel(badChannels));
	for iCh = 1:numel(badChannels)
		idx = badChannels(iCh);
		if idx >= 1 && idx <= numel(EEG.chanlocs) && isfield(EEG.chanlocs, 'labels') && ~isempty(EEG.chanlocs(idx).labels)
			channelLabels{iCh} = formatRejectedChannelLabel(EEG.chanlocs(idx).labels);
		else
			channelLabels{iCh} = formatRejectedChannelLabel(sprintf('Ch%d', idx));
		end
	end

	maxLabelsPerLine = 5;
	lineTexts = {};
	for iStart = 1:maxLabelsPerLine:numel(channelLabels)
		iEnd = min(iStart + maxLabelsPerLine - 1, numel(channelLabels));
		lineTexts{end+1} = strjoin(channelLabels(iStart:iEnd), ', '); %#ok<AGROW>
	end

	if isempty(lineTexts)
		legendText = sprintf('Rejected channels (%d): none', numel(badChannels));
	elseif numel(lineTexts) == 1
		legendText = sprintf('Rejected channels (%d): %s', numel(badChannels), lineTexts{1});
	else
		legendText = sprintf('Rejected channels (%d): %s\n%s', numel(badChannels), lineTexts{1}, strjoin(lineTexts(2:end), '\n'));
	end
end

function shiftFigureAxesLeftForTimtopo(hf, shiftFraction)
	if nargin < 2 || isempty(shiftFraction)
		shiftFraction = 0.02;
	end
	if isempty(hf) || ~isgraphics(hf, 'figure')
		return;
	end

	allAxes = findall(hf, 'type', 'axes');
	for iAx = 1:numel(allAxes)
		pos = get(allAxes(iAx), 'Position');
		newPos = pos;
		newPos(1) = max(0, pos(1) - shiftFraction);
		newPos(3) = max(0, pos(3) - shiftFraction);
		set(allAxes(iAx), 'Position', newPos);
	end
end

function addRejectedChannelsTextNote(hObj, legendText, horizontalAlignment)
	hf = [];
	if isempty(hObj)
		return;
	end
	if nargin < 3 || isempty(horizontalAlignment)
		horizontalAlignment = 'left';
	end

	if isgraphics(hObj, 'figure')
		hf = hObj;
	elseif isgraphics(hObj)
		hf = ancestor(hObj, 'figure');
	end
	if isempty(hf) || ~isgraphics(hf, 'figure')
		return;
	end

	figure(hf);

	allAxes = findall(hf, 'type', 'axes');
	if isempty(allAxes)
		return;
	end

	oldNotes = findall(hf, 'type', 'text', 'Tag', 'RejectedChannelsTextNote');
	if ~isempty(oldNotes)
		delete(oldNotes);
	end

	% Use the upper-most axis so the note appears in a consistent location.
	axisTops = nan(numel(allAxes), 1);
	for iAx = 1:numel(allAxes)
		pos = get(allAxes(iAx), 'Position');
		axisTops(iAx) = pos(2) + pos(4);
	end
	[~, idxTop] = max(axisTops);
	targetAxis = allAxes(idxTop);

	textX = 0.5;
	if strcmpi(horizontalAlignment, 'left')
		textX = 0.02;
	end

	text(targetAxis, textX, 1.08, legendText, ...
		'Units', 'normalized', ...
		'Tag', 'RejectedChannelsTextNote', ...
		'Interpreter', 'tex', ...
		'HorizontalAlignment', horizontalAlignment, ...
		'VerticalAlignment', 'bottom', ...
		'FontSize', 10, ...
		'FontWeight', 'bold', ...
		'Color', [0 0 0], ...
		'Clipping', 'off');
end

function formattedLabel = formatRejectedChannelLabel(channelLabel)
	channelLabel = strtrim(channelLabel);
	escapedLabel = strrep(channelLabel, '_', '\_');
	if strcmpi(channelLabel, 'GND')
		formattedLabel = ['\color{red}' escapedLabel '\color{black}'];
	else
		formattedLabel = escapedLabel;
	end
end

function adjustICAComponentFigureLayout(hf)
    if isempty(hf) || ~isgraphics(hf, 'figure')
        return;
    end

    % Make the overall figure a bit larger so the dense ICA panels have room.
    figPos = get(hf, 'Position');
    if isempty(figPos) || numel(figPos) < 4
        figPos = [50 50 1400 900];
    end
    set(hf, 'Position', [figPos(1) figPos(2) max(figPos(3), 1400) max(figPos(4), 950)]);

    allAxes = findall(hf, 'type', 'axes');
    if ~isempty(allAxes)
        for iAx = 1:numel(allAxes)
            ax = allAxes(iAx);
            ax.FontSize = 9;
            ax.TitleFontSizeMultiplier = 1;
            ax.LabelFontSizeMultiplier = 1;

            % Give each subplot a bit more vertical space so labels do not collide.
            axPos = ax.Position;
            ax.Position = [axPos(1), axPos(2), axPos(3), axPos(4) * 0.88];
        end
    end

    % Fine-tune the component-number buttons so they sit close to the topography
    % plots but do not overlap the panel titles.
    allButtons = findall(hf, 'type', 'uicontrol', 'Style', 'pushbutton');
    if ~isempty(allButtons)
        for iBtn = 1:numel(allButtons)
            btn = allButtons(iBtn);
            if isempty(btn.Tag) || isempty(regexp(btn.Tag, '^comp\d+$', 'once'))
                continue;
            end
            btnPos = btn.Position;
            btn.Position = [btnPos(1), max(0, btnPos(2) - 0.01), btnPos(3), btnPos(4)];
        end
    end

    allText = findall(hf, 'type', 'text');
    if ~isempty(allText)
        set(allText, 'FontSize', 9);
    end

    allLegends = findall(hf, 'type', 'legend');
    if ~isempty(allLegends)
        set(allLegends, 'FontSize', 9);
    end
end

function plotTEPComparison(s, EEGs, datasetLabels)
c_say('Plotting TEP comparison');
hf = figure;
ht = c_GUI_Tiler();

ht.pauseAutoRetiling();
for iEEG = 1:length(EEGs)
    tmp = c_EEG_plotTimtopo(EEGs{iEEG},...
		'doShowColorbars', false,...
        'TPOI', s.plotTPOIs,...
		'doPlotGMFA', true,...
        'xlim', [-100 350],...
        'parent', ht.add(),...
		'doSymmetricYLim', true,...
        'title', datasetLabels{iEEG});
end

if length(EEGs) < 4
    hf.Position = [50 50 700*length(EEGs) 600];
    ht.numRows = 1;
else
    c_fig_arrange('maximize', hf, 'monitor',  s.maximizePlotsToMonitor);
end
ht.resumeAutoRetiling();

c_sayDone();
end

function plotERPImComparison(s, EEGs, datasetLabels)
c_say('Plotting ERPImage comparison');
hf = figure;
ht = c_GUI_Tiler();

ht.pauseAutoRetiling();
for iEEG = 1:length(EEGs)
	sht = c_GUI_Tiler('parent', ht.add(), 'title', datasetLabels{iEEG}, 'numRows', 1);
	has = gobjects(0);
	for iCh = 1:length(s.plotChans)
		if iscellstr(s.plotChans) && ~ismember(s.plotChans{iCh}, {EEGs{iEEG}.chanlocs.labels})
			ha = sht.addAxes(); % placeholder
		else
			ha = sht.addAxes();
			has(end+1) = ha;
			c_EEG_plotERPImage(EEGs{iEEG},...
				'COI', s.plotChans(iCh),...
				'axis', ha);
			xlim(ha, [-100 350]);
		end
		if iscellstr(s.plotChans)
			title(ha, sprintf('Chan %s', s.plotChans{iCh}));
		else
			title(ha, sprintf('Chan %d', s.plotChans(iCh)))
		end
	end
	ha = sht.addAxes('relWidth', 0.3, 'relHeight', 0.5);
	ignoreTimeIndices = EEGs{iEEG}.times > s.artifactTimespan(1)*2*1e3 & EEGs{iEEG}.times < s.artifactTimespan(2)*2*1e3; 
	c_plot_colorbar(...
		'axis', ha,...
		'clim', [-1 1]*max(abs(extrema(reshape(mean(EEGs{iEEG}.data(:, ~ignoreTimeIndices, :), 3), 1, []))))*2,...
		'linkedAxes', has,...
		'doForceSymmetric', true,...
		'clabel', 'Amplitude (\mu{V})');
	
end

if length(EEGs) < 4
    hf.Position = [50 50 700*length(EEGs) 600];
    ht.numRows = 1;
else
    c_fig_arrange('maximize', hf, 'monitor',  s.maximizePlotsToMonitor);
end
ht.resumeAutoRetiling();

c_sayDone();
end
