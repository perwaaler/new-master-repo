function pc_string = pc_str(coll_ind)

if coll_type(coll_ind) == 1
    pc_string = 'P(C,NEA)';
elseif coll_type(coll_ind) == 2
    pc_string = 'P(C,EA)';
else
    pc_string = 'P(C)';
end

end