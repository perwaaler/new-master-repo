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


% plan for tomorrow:
% I have noticed that there is a much easier way to check which side of
% eachother the road users are. I can simply draw a line through two
% neighbouring points and then compare whether the ball is above or below
% the line. Tomorrow my plan is to implement this way of cheking
% orientation.

% plan for tomorrow: The code for figuring out orientation mostly works
% well now, both for the case where there is path overlap, and the case
% when the paths diverge. The only remaining problem seems to be the case
% when RU A stops entirely, because then method for determining orientation
% does not work properly, since it needs movement to determine the tangent
% line. Will fix this tomorrow...