H_curvature = zeros (nv,3);
%nv is the number of points
A_matrix

for i_cur = 1:nv
    connect_face = findface(v,nt,i_cur);
    if connect_face(6) == 0
        size = 5;
    else
        size = 6;
    end
    
    for j_cur = 1:size
        H_curvature (i_cur,:) = H_curvature (i_cur,:) - Derivative_Ai (x,v,connect_face(j_cur),i_cur,A_mat);
    end
end