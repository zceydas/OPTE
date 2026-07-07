function psd_table = psd_to_long_table(psd_abs, psd_log, psd_rel, freqs, chan_labels, info, base_name)

[nChannels, nFreqs] = size(psd_abs);
nRows = nChannels * nFreqs;

Participant = strings(nRows,1);
Session = strings(nRows,1);
Eyes = strings(nRows,1);
Epoch = strings(nRows,1);
File = strings(nRows,1);
Channel = strings(nRows,1);
Frequency = nan(nRows,1);
AbsolutePower = nan(nRows,1);
Log10Power = nan(nRows,1);
RelativePower = nan(nRows,1);

idx = 1;

for ch = 1:nChannels
    for fi = 1:nFreqs

        Participant(idx) = string(info.participant);
        Session(idx) = string(info.session);
        Eyes(idx) = string(info.eyes);
        Epoch(idx) = string(info.epoch);
        File(idx) = string(base_name);
        Channel(idx) = string(chan_labels(ch));
        Frequency(idx) = freqs(fi);
        AbsolutePower(idx) = psd_abs(ch,fi);
        Log10Power(idx) = psd_log(ch,fi);
        RelativePower(idx) = psd_rel(ch,fi);

        idx = idx + 1;

    end
end

psd_table = table(Participant, Session, Eyes, Epoch, File, Channel, Frequency, ...
    AbsolutePower, Log10Power, RelativePower);

end