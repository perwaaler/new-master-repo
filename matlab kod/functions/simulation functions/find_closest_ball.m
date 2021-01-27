function candidates = find_closest_ball(candidates,...
                                            pathA,pathB,k,w)

% find index of closest neighbours in the other RUs trajectory.
% Arguments:
% candidates = set of neighbours that are closest in previous iteration
% path = the positions of the trajectory computed so far
% k    = iteration
% w    = determines number of candidates considered

% returns field that contains:
% index    = cell array with indeces for the closest candidates
% min_dist = vector with minimum distances for each set of candidates
% closest  = vecotr with indeces of closest candidates

candidatesA = candidates.index{1};
candidatesB = candidates.index{2};

%%%% path A
% remove elements with indeces that spill over 
candidatesB = candidatesB(candidatesB >= 1);
candidatesB = candidatesB(candidatesB <= k);

% find distance between last position of A and candidates of B
diffA = pathA(k) - pathB(candidatesB);
diffA = sqrt(diffA.*conj(diffA));

% update closest position in B
[min_distA,closest{2}] = min(diffA);
closest{2} = candidatesB(closest{2});

% recenter around new candidate
candidatesB = closest{2}-w:closest{2}+w;


%%%% path B
% remove elements with indeces that spill over 
candidatesA = candidatesA(candidatesA>=1);
candidatesA = candidatesA(candidatesA<=k);

% find distance between last position of A and candidates of B
diffB = pathB(k) - pathA(candidatesA);
diffB = sqrt(diffB.*conj(diffB));

% update closest position in B
[min_distB, closest{1}] = min(diffB);
closest{1} = candidatesA(closest{1});

% recenter around new candidate
candidatesA = closest{1}-w:closest{1}+w;

% save index of closest candidates, the minimum distances and index of
% closest cadidates:
candidates.index    = {candidatesA, candidatesB};
candidates.min_dist = [min_distA, min_distB];
candidates.closest  = closest;

end