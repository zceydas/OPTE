function [EEG, misc] = c_TMSEEG_prepareForPreprocessing(varargin)
%% this function does several miscellaneous data input tasks to prepare for other preprocessing
% (implemented here to not duplicate this common code in multiple dataset-specific analysis scripts)
% Main steps include:
% - (If needed) Load data 
% - (If needed) Infer pulse event based on most frequent event
% - (If needed) Infer reasonable epoch timespan (not yet implemented)
% - (If needed) Concatenate multiple datasets into one EEG struct
%	- (If needed) label each dataset prior to concatenation to allow later identification of original dataset


c_EEG_openEEGLabIfNeeded();

%% parse inputs

p = inputParser();
p.addParameter('inputFilePaths', '', @(x) ischar(x) || iscellstr(x));
p.addParameter('inputEEGs',[],@(x) isstruct(x) || iscell(x));
p.addParameter('inputDatasetLabels', 'auto', @(x) ischar(x) || iscellstr(x));
p.addParameter('pulseEvent', 'auto', @ischar);
p.addParameter('pulseEvents', {}, @iscellstr);  % can specify multiple pulse event types to be treated as the same
p.addParameter('ifPulseTooCloseToBoundary', 'ignore', @ischar);
p.addParameter('epochTimespan', [], @c_isSpan);
p.parse(varargin{:});
s = p.Results;

assert(xor(isempty(s.inputFilePaths), isempty(s.inputEEGs)),'Must specify one of inputFilePaths or inputEEG');

misc = struct();

%%

if ~isempty(s.inputFilePaths)
	%% load EEG data
	assert(~isempty(s.inputFilePaths));
	if ischar(s.inputFilePaths)
		s.inputFilePaths = {s.inputFilePaths};
	end
	EEGs = {};
	prog = c_progress(length(s.inputFilePaths),'Loading EEG file %d/%d');
	prog.start('Loading EEG data');
	for iI = 1:length(s.inputFilePaths)
		prog.updateStart(iI);
		c_saySingle('Loading EEG from %s', s.inputFilePaths{iI});
		[~, ~, inputFileExt] = fileparts(s.inputFilePaths{iI});
		assert(c_exist(s.inputFilePaths{iI},'file')>0);
		switch(inputFileExt)
			case '.vhdr'
				EEGs{iI} = c_loadEEG_BrainProducts(s.inputFilePaths{iI});
			case '.mat'
				tmp = load(s.inputFilePaths{iI});
				if isfield(tmp, 'EEG')
					% mat file had an EEG struct, ignore anything else
					EEGs{iI} = tmp.EEG;
				elseif isfield(tmp, 'data')
					% if has a data field, assume mat file was a directly saved EEG struct
						EEGs{iI} = tmp;
				else
					error('EEG struct not found in loaded .mat file')
				end

			otherwise
				error('Unsupported input type: %s', ext);
		end
		prog.updateEnd(iI);
	end
	prog.stop()
	
	if isequal(s.inputDatasetLabels,'auto')
		[~, s.inputDatasetLabels] = c_str_findCommonPrefix(s.inputFilePaths);
		[~, s.inputDatasetLabels] = c_str_findCommonSuffix(s.inputFilePaths);
		c_saySingle('Using auto generated data labels: %s', c_toString(s.inputDatasetLabels));
	end
else
	if isstruct(s.inputEEGs)
		EEGs = {s.inputEEGs};
	else
		assert(all(cellfun(@isstruct, s.inputEEGs)));
		EEGs = s.inputEEGs;
	end
	if isequal(s.inputDatasetLabels,'auto')
		s.inputDatasetLabels = arrayfun(@(iD) sprintf('Dataset-%d', iD), 1:length(EEGs),'UniformOutput',false);
		c_saySingle('Using auto generated data labels: %s', c_toString(s.inputDatasetLabels));
	end
end

if ischar(s.inputDatasetLabels)
	assert(length(EEGs)==1);
	s.inputDatasetLabels = {s.inputDatasetLabels};
end
	
%% infer pulse event type if requested
if ~isempty(s.pulseEvents)
	assert(ismember('pulseEvent', p.UsingDefaults), 'Should not specify both pulseEvent and pulseEvents');
	c_saySingle('Specified pulse event types: %s', c_toString(s.pulseEvents));
	misc.pulseEvents = s.pulseEvents;
	pulseEvents = s.pulseEvents;
else
	if strcmpi(s.pulseEvent,'auto')
		args = cellfun(@(EEG) EEG.event, EEGs, 'UniformOutput',false);
		allEvents = cat(2,args{:});
		[counts, eventTypes] = c_countUnique({allEvents.type});
		[~,index] = max(counts);
		mostFrequentEvent = eventTypes{index};
		pulseEvent = mostFrequentEvent;
		assert(~ismember(pulseEvent,{'boundary'}));
		assert(c_str_matchRegex(pulseEvent, {'[RST][ 0-9]*','Pulse'}));
		s.pulseEvent = pulseEvent;
		c_saySingle('Inferred pulse event type: ''%s''', s.pulseEvent);
	else
		c_saySingle('Specified pulse event type: ''%s''', s.pulseEvent);
	end
	misc.pulseEvent = s.pulseEvent;
	pulseEvents = {s.pulseEvent};
end


%% infer reasonable epoch timespan if not specified
if isempty(s.epochTimespan)
	keyboard %TODO
	
end
misc.epochTimespan = s.epochTimespan; 

%% trim and concatenate

if ~isempty(s.inputDatasetLabels)
	% before concatenation, label individual epochs with their parent dataset so that 
	%  later analyses which epochs in concatenated data came from which original dataset
	assert(length(s.inputDatasetLabels)==length(EEGs));
	for iD = 1:length(EEGs)
		assert(~isfield(EEGs{iD}.event,'datasetLabel'));
		[EEGs{iD}.event.datasetLabel] = deal(s.inputDatasetLabels{iD});
	end
end

% trim any large excess of non-pulse data at beginning and end of each file
c_say('Trimming continuous data');
extraTime = s.epochTimespan*2; % time in seconds before first and after last pulse to keep
for iD = 1:length(EEGs)
	EEG = EEGs{iD};

	switch(s.ifPulseTooCloseToBoundary)
		case 'ignore'
			% do nothing
		case 'error'
			keyboard  % TODO
		case 'drop'
			for iD = 1:length(EEGs)
				EEG = EEGs{iD};
				pulseIndices = find(ismember({EEG.event.type}, pulseEvents));
				pulseTimes = [EEG.event(pulseIndices).latency] / EEG.srate;
				boundaryTimes = [EEG.event(ismember({EEG.event.type}, {'boundary'})).latency] / EEG.srate;
				boundaryTimes = [EEG.times(1)/1e3, boundaryTimes, EEG.times(end)/1e3];
				timesFromBoundaryToPulse = pulseTimes - boundaryTimes';
				timesFromPulseToBoundary = -timesFromBoundaryToPulse;
				timesFromBoundaryToPulse(timesFromBoundaryToPulse < 0) = NaN;
				timesFromPulseToBoundary(timesFromPulseToBoundary < 0) = NaN;
				toRemove = any(timesFromBoundaryToPulse < -s.epochTimespan(1), 1) | ...
					any(timesFromPulseToBoundary < s.epochTimespan(2), 1);
				if any(toRemove)
					c_saySingle('Removing %d pulse events too close to data boundary', sum(toRemove))
					EEG.event(pulseIndices(toRemove)) = [];
					EEGs{iD} = EEG;
				end
			end
		otherwise
			error('Not implemented')
	end


	firstEventIndex = find(ismember({EEG.event.type}, pulseEvents),1,'first');
	startTime = EEG.event(firstEventIndex).latency/EEG.srate + extraTime(1);
	lastEventIndex = find(ismember({EEG.event.type},pulseEvents),1,'last');
	endTime = EEG.event(lastEventIndex).latency/EEG.srate + extraTime(2);
	if startTime > 0 || (EEG.pnts-1)/EEG.srate - endTime > 0
		c_saySingle('Cutting %.2f s at beginning and %.2f s at end', startTime, (EEG.pnts-1)/EEG.srate - endTime);
		EEGs{iD} = pop_select(EEG,'time',[startTime, endTime]);
	else
		c_saySingle('Trim not needed');
	end
end
c_sayDone();

if length(EEGs) > 1
	c_say('Concatenating EEG data');
	EEG = pop_mergeset(cell2mat(EEGs), 1:length(EEGs), 0);
	c_sayDone();
else
	EEG = EEGs{1};
end
clearvars EEGs

%% [OLD] This particular dataset has Fpz instead of GND, the lines below correct this

%idx = strcmpi({EEG.chanlocs.labels}, 'Fpz');
%[EEG.chanlocs(idx).labels] = deal('GND');

%% [OLD] Two of the channels we have do not have coordinates, we will replace these

% Within EEG.chanlocs, copy the values from the 'Cz' row into the 'FCz' row
%src = find(strcmpi({EEG.chanlocs.labels}, 'Cz'), 1);
%dst = find(strcmpi({EEG.chanlocs.labels}, 'FCz'), 1);

%tmp = EEG.chanlocs(src);    % copy Cz struct
%tmp.labels = EEG.chanlocs(dst).labels;   % preserve 'FCz' label
%EEG.chanlocs(dst) = tmp;   % assign into FCz row

% Within EEG.chanlocs, copy the values from the 'Cz' row into the 'FCz' row
%src = find(strcmpi({EEG.chanlocs.labels}, 'Fz'), 1);
%dst = find(strcmpi({EEG.chanlocs.labels}, 'AFz'), 1);

%tmp = EEG.chanlocs(src);    % copy Cz struct
%tmp.labels = EEG.chanlocs(dst).labels;   % preserve 'FCz' label
%EEG.chanlocs(dst) = tmp;   % assign into FCz row

%% [NEW] Update the FCz and AFz chanlocs to match the coordinates found in the .ced file
% Find the index of the FCz channel
fczIdx = find(strcmpi({EEG.chanlocs.labels}, 'FCz'), 1);

% If FCz exists, replace the first 8 fields with the new values
if ~isempty(fczIdx)
    EEG.chanlocs(fczIdx).sph_radius = 1;  % Assuming the first value is for labels (as a string)
    EEG.chanlocs(fczIdx).sph_theta = 0;
    EEG.chanlocs(fczIdx).sph_phi = 67;
    EEG.chanlocs(fczIdx).theta = 0;
    EEG.chanlocs(fczIdx).radius = 0.12778;
    EEG.chanlocs(fczIdx).X = 0.39073;
    EEG.chanlocs(fczIdx).Y = 0;
    EEG.chanlocs(fczIdx).Z = 0.9205;
end

% Find the index of the AFz channel
afzIdx = find(strcmpi({EEG.chanlocs.labels}, 'AFz'), 1);

% If AFz exists, replace the first 8 fields with the new values
if ~isempty(afzIdx)
	EEG.chanlocs(afzIdx).sph_radius = 1;  % Assuming the first value is for labels (as a string)
	EEG.chanlocs(afzIdx).sph_theta = 0;
	EEG.chanlocs(afzIdx).sph_phi = 23;
	EEG.chanlocs(afzIdx).theta = 0;
	EEG.chanlocs(afzIdx).radius = 0.37222;
	EEG.chanlocs(afzIdx).X = 0.9205;
	EEG.chanlocs(afzIdx).Y = 0;
	EEG.chanlocs(afzIdx).Z = 0.39073;
end

%% [NEW] Some datasets won't have a GND channel in the data file
% This is usually because the GND electrode was used as a reference during recording, so it doesn't appear in the data.
% First we check if there is a GND channel in EEG object, and if not, we add a GND channel with all zero data
if ~any(strcmpi({EEG.chanlocs.labels}, 'GND'))
	c_saySingle('No GND channel found, adding GND channel with all zero data');
	EEG.nbchan = EEG.nbchan + 1;
	EEG.data(end+1,:) = 0; % add new row of zeros to data
	EEG.chanlocs(end+1) = EEG.chanlocs(1); % copy chanlocs from first channel (arbitrary since we will overwrite the label and coordinates)
	EEG.chanlocs(end).labels = 'GND'; % set correct label
end

% Now we update the chanlocs for the GND channel to match the coordinates to the closest channel found in the .ced file
gndIdx = find(strcmpi({EEG.chanlocs.labels}, 'GND'), 1);
if ~isempty(gndIdx)
	EEG.chanlocs(gndIdx).sph_radius = 1;
	EEG.chanlocs(gndIdx).sph_theta = -72;
	EEG.chanlocs(gndIdx).sph_phi = -23;
	EEG.chanlocs(gndIdx).theta = 72;
	EEG.chanlocs(gndIdx).radius = 0.62778;
	EEG.chanlocs(gndIdx).X = 0.28445;
	EEG.chanlocs(gndIdx).Y = -0.87545;
	EEG.chanlocs(gndIdx).Z = -0.39073;
end

%% make sure chanlocs are set (interpolation and plotting in ARTIST needs channel locations)
if length([EEG.chanlocs.X]) < EEG.nbchan ....
		&& length([EEG.chanlocs.sph_theta]) < EEG.nbchan ...
		&& length([EEG.chanlocs.theta]) < EEG.nbchan
	% missing chanlocs, load from default
	
	% require that labels are set
	assert(~any(cellfun(@isempty, {EEG.chanlocs.labels})));

	% drop EMG channels
	removeChanIndices = c_str_matchRegex({EEG.chanlocs.labels}, 'EMG.*');
	if any(removeChanIndices)
		c_say('Removing %d EMG channels', sum(removeChanIndices));
		EEG = pop_select(EEG, 'nochannel', find(removeChanIndices));
		c_sayDone();
	end
	
	switch(EEG.nbchan)
		case 95
			c_saySingle('No chanlocs set, loading ActiCAP-96 default locations');
			defaultChanlocsPath = fullfile(fileparts(which(mfilename)),'Resources','ActiCAP-96.ced');
		case 63
			c_saySingle('No chanlocs set, loading ActiCAP-64 default locations');
			defaultChanlocsPath = fullfile(fileparts(which(mfilename)),'Resources','ActiCAP-64.ced');
		otherwise
			error('No chanlocs template available for %d channel montage', EEG.nbchan);
	end
	
	chanlocs = readlocs(defaultChanlocsPath);
	% all our electrodes should be in the default chanlocs
	labelsAreNumeric = all(~isnan(cellfun(@str2double, {EEG.chanlocs.labels})));
	if labelsAreNumeric
		% assume channel order is the same as in template chanlocs
		indices = cellfun(@str2double, {EEG.chanlocs.labels});
	else
		assert(all(ismember({EEG.chanlocs.labels}, {chanlocs.labels})));
		indices = c_cell_findMatchingIndices({EEG.chanlocs.labels}, {chanlocs.labels});
	end
	 
end

if ~isfield(EEG.chanlocs, 'type') || all(cellfun(@isempty, {EEG.chanlocs.type}))
	% if no channel types are set, assume all channels are EEG channels
	c_saySingle('No channel types specified, assuming all are EEG.');
	[EEG.chanlocs.type] = deal('EEG');
end




end