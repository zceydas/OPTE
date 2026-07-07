function chan_labels = get_channel_labels_psd(EEG, nChannels)

chan_labels = strings(nChannels, 1);

if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= nChannels

    for ch = 1:nChannels

        if isfield(EEG.chanlocs(ch), 'labels') && ~isempty(EEG.chanlocs(ch).labels)
            chan_labels(ch) = string(EEG.chanlocs(ch).labels);
        else
            chan_labels(ch) = "Ch" + string(ch);
        end

    end

else

    for ch = 1:nChannels
        chan_labels(ch) = "Ch" + string(ch);
    end

end

end