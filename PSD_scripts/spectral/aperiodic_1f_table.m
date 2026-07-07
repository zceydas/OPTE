function aperiodic_table = aperiodic_1f_table(psd_abs, freqs, chan_labels, info, base_name, fit_range, exclude_ranges)

nChannels = size(psd_abs,1);

Participant = strings(nChannels,1);
Session = strings(nChannels,1);
Eyes = strings(nChannels,1);
Epoch = strings(nChannels,1);
File = strings(nChannels,1);
Channel = strings(nChannels,1);
AperiodicSlope = nan(nChannels,1);
AperiodicIntercept = nan(nChannels,1);
AperiodicRSquared = nan(nChannels,1);

fit_keep = freqs >= fit_range(1) & freqs <= fit_range(2);

for r = 1:size(exclude_ranges,1)
    fit_keep = fit_keep & ~(freqs >= exclude_ranges(r,1) & freqs <= exclude_ranges(r,2));
end

x = log10(freqs(fit_keep));

for ch = 1:nChannels

    y = log10(psd_abs(ch, fit_keep) + eps);
    valid = isfinite(x) & isfinite(y);

    if sum(valid) >= 3

        p = polyfit(x(valid), y(valid), 1);
        yhat = polyval(p, x(valid));

        ss_res = sum((y(valid) - yhat).^2);
        ss_tot = sum((y(valid) - mean(y(valid))).^2);

        AperiodicSlope(ch) = p(1);
        AperiodicIntercept(ch) = p(2);

        if ss_tot ~= 0
            AperiodicRSquared(ch) = 1 - ss_res / ss_tot;
        end
    end

    Participant(ch) = string(info.participant);
    Session(ch) = string(info.session);
    Eyes(ch) = string(info.eyes);
    Epoch(ch) = string(info.epoch);
    File(ch) = string(base_name);
    Channel(ch) = string(chan_labels(ch));

end

aperiodic_table = table(Participant, Session, Eyes, Epoch, File, Channel, ...
    AperiodicSlope, AperiodicIntercept, AperiodicRSquared);

end