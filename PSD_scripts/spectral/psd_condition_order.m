function condition_table = psd_condition_order()

Session = [
    "baseline"
    "dosing"
    "dosing"
    "dosing"
    "dosing"
    "1week"
    "2week"
    "1month"
];

Epoch = [
    "Epoch0"
    "Epoch1"
    "Epoch2"
    "Epoch3"
    "Epoch4"
    "Epoch0"
    "Epoch0"
    "Epoch0"
];

Condition = [
    "Baseline"
    "Dosing E1"
    "Dosing E2"
    "Dosing E3"
    "Dosing E4"
    "1 Week"
    "2 Week"
    "1 Month"
];

condition_table = table(Session, Epoch, Condition);

end
