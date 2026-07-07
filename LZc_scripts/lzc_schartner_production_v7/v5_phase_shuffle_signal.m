function x_surr = v5_phase_shuffle_signal(x, seed)

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
    pos_idx = 2:(n/2);
    neg_idx = n:-1:(n/2+2);
else
    pos_idx = 2:((n+1)/2);
    neg_idx = n:-1:((n+3)/2);
end

random_phases = exp(1i * 2*pi * rand(1, numel(pos_idx)));

Xf_surr(pos_idx) = abs(Xf(pos_idx)) .* random_phases;
Xf_surr(neg_idx) = conj(Xf_surr(pos_idx));

x_surr = double(real(ifft(Xf_surr)) + x_mean);

end
