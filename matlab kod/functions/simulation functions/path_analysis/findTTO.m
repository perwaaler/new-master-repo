function [t, d] = findTTO(d_prox_t,paths,r,Id,opt)
% computes the distance and time to reach the closest point in pathB

% extract variables from paths
ta = paths.t_prox(1);
tb = paths.t_prox(2);
% Ta_upper = paths.rel_time(1);
% Tb_upper = paths.rel_time(2);

g =@(X) d_prox_t(X);
t_low_bnd = find_t_first_contact(d_prox_t,[ta,tb],r);

if not(paths.contact)
    [t,d] = fminsearch(g,[ta,tb],opt);
else
%     % extract times of first path contact
%     ta_fc = paths.t_first_cont(1);
%     tb_fc = paths.t_first_cont(2);
    
    D  = 2*r;
    
    if Id==1
        %%%% Find proximity point A %%%%
        
        % search limits for ta:
%         al  = max(0,ta-5);
% %         al = t_low_bnd(1);
%         au  = ta;
% 
%         % search limits for ball in path B that A first touches:
%         bl = max(         tb - 5,0);
%         bu = min(Tb_upper,tb + 5);
%         
%         [t,d] = fminbnd(@(x) ...
%                  abs(g([x, fminbnd(@(y) g([x,y]),bl,bu,opt)]) - D),...
%                                                  al,au,opt);
                                             
                                             
        %%%% method 2 using minsearch %%%%                                          
        % search limits for tb:
        al = t_low_bnd(1);
        au = ta;

        [t,d] = fminbnd(@(x) ...
                 abs(g([x, fminsearch(@(y) g([x,y]),tb,opt)]) - D),...
                                                 al,au,opt);

    else

        % search limits for tb:
%         bl  = max(0,tb-5);
% %         bl = t_low_bnd(2);
%         bu = tb;
% 
% %         % search limits for ball in path B that A first touches:
%         al = max(         ta - 5,0);
%         au = min(Ta_upper,ta + 5);
% 
%         [t,d] = fminbnd(@(y) ...
%                         abs( g([fminbnd(@(x) g([x,y]),al,au,opt),y]) -D),...
%                                                       bl, bu, opt);
                                                  
                                                  
        %%%% method 2 using minsearch %%%%                                          
        % search limits for tb:
        bl = t_low_bnd(2);
        bu = tb;

        [t,d] = fminbnd(@(y) ...
                        abs( g([fminsearch(@(x) g([x,y]),ta,opt),y]) -D),...
                                                      bl, bu, opt);
                                                  
    end
end

end