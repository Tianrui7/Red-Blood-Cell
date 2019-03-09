%This program is to calculate the total area, for any given matrix x

function t_a = totalarea(v,nt,x)

temp_matrix_area = af(v,nt,x);
t_a = 0;
for i_t = 1:nt
    
    t_a = t_a + temp_matrix_area (i_t);
end
