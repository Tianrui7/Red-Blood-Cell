function connectf = findface(v,nt,n)
connectf = [ 0 0 0 0 0 0];
k_find=1;
%k is the index for connectf
for i_find= 1:nt
    %i is the index of the face that we are checking
    for j_find= 1:3
       if v(i_find,j_find)==n
           connectf(k_find) = i_find;
           k_find=k_find+1;
           break
       end
    end
    %There is probably a way to simplfy this program by changing 1:nt to a
    %smaller range
end
