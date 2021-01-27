% Things left to do:

% I have noted that the extrapolations are sometimes not very reasonable,
% probably due to the randomization that occurs. I want to use averaging in
% order to set delta_theta and delta_speed, otherwise the extrapolations
% becomes to sensitive to the last step taken. Need smoothening over time.
% Also, I am considering having no randomization at all in the
% extrapolations, or at least I want to decrease the amount of random
% variation in the predictions/extrapolations.

% need to assign Tadv even to predicitions for which there is no collision
% of predicted trajectories.