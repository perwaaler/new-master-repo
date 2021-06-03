function [candidates,done_intersecting] = est_prox_points(candidates,...
                                                     pathA,pathB,k,w,d_thr)
% used to iteratively find the proximity points of RUA and RUB.
% Note
% Arguments:
% candidates = set of neighbours that are closest in previous iteration
%              candidates{1} is the index of the ball in pathA that lies
%              closest to the k:th ball of pathB.
% path = the positions of the trajectory computed so far
% k    = iteration
% w    = determines number of candidates considered

% returns field that contains:
% index    = cell array with indeces for the closest candidates
% min_dist = vector with minimum distances for each set of candidates
% closest  = vector with indeces of closest candidates. 
%            closest{1} = ball in pathA that lies closest to the proximity
%                         candidate of B.

done_intersecting.A = false;
done_intersecting.B = false; 

% candidates
canA = candidates.index{1};
canB = candidates.index{2};
prime_candA = candidates.prime_cand(1);
prime_candB = candidates.prime_cand(2); %#ok<NASGU>


%%%% path A
% remove elements with indeces that spill over 
canB = canB(canB >= 1);
canB = canB(canB <= k);

% find distance between the prime candidate of A and the candidates of B
d_pcA_2_cB = pathA(prime_candA) - pathB(canB);
d_pcA_2_cB = abs(d_pcA_2_cB);
tendrils(1,:) = [d_pcA_2_cB(1),d_pcA_2_cB(end)];

% update the prime candidate of B
[Da, I] = min(d_pcA_2_cB);
prime_candB = canB(I);

% recenter around new prime candidate
canB = prime_candB - w : prime_candB + w;


%%%% path B
% remove elements with indeces that spill over 
canA = canA(canA >= 1);
canA = canA(canA <= k);

% find distance between the prime candidate of B and the candidates of A
d_pcB_2_cA = pathB(prime_candB) - pathA(canA);
d_pcB_2_cA = abs(d_pcB_2_cA);
tendrils(2,:) = [d_pcB_2_cA(1),d_pcB_2_cA(end)];

% update the prime candidate of B
[Db, I] = min(d_pcB_2_cA);
prime_candA = canA(I);

% recenter around new prime candidate
canA = prime_candA - w : prime_candA + w;


% % neighbourhoodA format: who am I? which B-ball is my closest neighbour?
% neighbourhoodA = [k, prime_candA];
% % neighbourhoodB format: who am I? which A-ball is my closest neighbour?
% neighbourhoodB = [k, prime_candB];

% check if the paths are done intersecting
if sum(tendrils(1,:)>d_thr)==2 && Da<d_thr
    done_intersecting.A = true;
elseif prime_candB==1 && tendrils(1,2)>d_thr && Da<d_thr 
    done_intersecting.A = true;
end
if sum(tendrils(2,:)>d_thr)==2 && Db<d_thr
    done_intersecting.B = true;
elseif prime_candA==1 && tendrils(2,2)>d_thr && Db<d_thr 
    done_intersecting.B = true;
end

% save index of closest candidates, the minimum distances and index of
% closest cadidates:
candidates.index               = {canA, canB};
candidates.min_dist            = [Da, Db];
candidates.prime_cand          = [prime_candA, prime_candB];
candidates.intersection_status = done_intersecting;
% candidates.neighbours = [neighbourhoodA; neighbourhoodB];


end