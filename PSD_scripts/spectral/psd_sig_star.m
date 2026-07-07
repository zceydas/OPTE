function star = psd_sig_star(p)
% Convert corrected p-value to significance star label.
if ~isfinite(p)
    star = "";
elseif p < 0.001
    star = "***";
elseif p < 0.01
    star = "**";
elseif p < 0.05
    star = "*";
else
    star = "";
end
end
