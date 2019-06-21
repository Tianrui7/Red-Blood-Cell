%this is to calculate total Energy 

H_curva
A_matrix

Ak = zeros(nv);
for i_ener = 1:nv
    con_ener = findface(v,nt,i_ener);
    if con_ener(6) == 0
        size =5;
    else
        size = 6;
    end
    for j_E = 1:size
        Ak(i_ener) = Ak(i_ener) + A_mat(con_ener(j_E))/3;
    end 
end


Energy_res = 0 ;
for i_ener = 1:nv
    Energy_res = Energy_res + ((norm(H_curvature(i_ener,:)))^2)/Ak(i_ener);
end
