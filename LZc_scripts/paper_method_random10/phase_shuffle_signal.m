function x_surr = phase_shuffle_signal(x, seed)
% phase_shuffle_signal
% Creates a phase-randomized surrogate of a single-channel time series.
% The Fourier amplitude spectrum is preserved and the Fourier phases are randomized.
%
% Input:
%   x    : 1 x time or time x 1 signal
%   seed : optional RNG seed
%
% Output:
%   x_surr : phase-randomized surrogate signal with same length as x

if nargin >= 2 && ~isempty(seed)
    rng(seed);
end

x = double(x(:)');
n = numel(x);

x_mean = mean(x, 'omitnan');
x = x - x_mean;

Xf = fft(x);
Xf_surr = Xf;

if mod(n, 2) == 0
    % Even length: keep DC and Nyquist fixed.
    pos_idx = 2:(n/2);
    neg_idx = n:-1:(n/2+2);
else
    % Odd length: keep DC fixed.
    pos_idx = 2:((n+1)/2);
    neg_idx = n:-1:((n+3)/2);
end

random_phases = exp(1i * 2*pi * rand(1, numel(pos_idx)));

Xf_surr(pos_idx) = abs(Xf(pos_idx)) .* random_phases;
Xf_surr(neg_idx) = conj(Xf_surr(pos_idx));

x_surr = real(ifft(Xf_surr));
x_surr = x_surr + x_mean;

end
