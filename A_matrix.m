
%'A' matrix is zeros(nt), where nt refers to number of faces. 'A(i)' stores
%the value of ||Ai||

A_mat = zeros (nt);

for i_mat = 1:nt
    x1_mat = v(i_mat, 1);
    x2_mat = v(i_mat, 2);
    x3_mat = v(i_mat, 3);
    
    A_mat(i_mat) = norm(0.5*(cross(x(x1_mat,:),x(x2_mat,:)) + cross(x(x2_mat,:),x(x3_mat,:)) + cross(x(x3_mat,:),x(x1_mat,:))));
end



    
