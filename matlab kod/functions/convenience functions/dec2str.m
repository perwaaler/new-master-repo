function out = dec2str(decision)
% function that converts decision into a string
if    decision == 0
    out = 'default';
elseif decision == 1
    out = 'let other pass';
elseif decision == 2
    out = 'wait';
elseif decision == 3
    out = 'go first';
elseif decision == 4
    out = 'swerve and brake';
end

end