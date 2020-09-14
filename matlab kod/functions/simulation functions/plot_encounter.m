function plot_encounter(A_save, B_save, save_react_A, save_react_B, enc_id, pause_length, xinit, r)  

for i=1:length(enc_id)
    enc_length = sum(A_save(enc_id(i),:)<inf);

    for j=1:enc_length
        plot_pos(A_save(enc_id(i),j), B_save(enc_id(i),j), pause_length, xinit, r, 1, save_react_A(enc_id(i),j), save_react_B(enc_id(i),j))
    end
end

end