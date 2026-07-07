function band_table = extract_band_power_table(psd_abs, psd_log, psd_rel, freqs, chan_labels, info, base_name, bands, band_names)

nChannels = size(psd_abs,1);
nBands = numel(band_names);
nRows = nChannels * nBands;

Participant = strings(nRows,1);
Session = strings(nRows,1);
Eyes = strings(nRows,1);
Epoch = strings(nRows,1);
File = strings(nRows,1);
Channel = strings(nRows,1);
Band = strings(nRows,1);
BandLowHz = nan(nRows,1);
BandHighHz = nan(nRows,1);
AbsolutePower = nan(nRows,1);
Log10Power = nan(nRows,1);
RelativePower = nan(nRows,1);

idx = 1;

for ch = 1:nChannels
    for b = 1:nBands

        band_name = band_names{b};
        fr = bands.(band_name);

        keep = freqs >= fr(1) & freqs < fr(2);

        Participant(idx) = string(info.participant);
        Session(idx) = string(info.session);
        Eyes(idx) = string(info.eyes);
        Epoch(idx) = string(info.epoch);
        File(idx) = string(base_name);
        Channel(idx) = string(chan_labels(ch));
        Band(idx) = string(band_name);
        BandLowHz(idx) = fr(1);
        BandHighHz(idx) = fr(2);

        AbsolutePower(idx) = mean(psd_abs(ch, keep), 'omitnan');
        Log10Power(idx) = mean(psd_log(ch, keep), 'omitnan');
        RelativePower(idx) = mean(psd_rel(ch, keep), 'omitnan');

        idx = idx + 1;

    end
end

band_table = table(Participant, Session, Eyes, Epoch, File, Channel, Band, ...
    BandLowHz, BandHighHz, AbsolutePower, Log10Power, RelativePower);

end