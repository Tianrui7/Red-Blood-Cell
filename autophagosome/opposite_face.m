global opposite_faces
opposite_faces =zeros(nv,1); 
for i =1:nv
    count = 0;
    for j = 1:nt
        if (j == new(i,1)) || (j == new(i,2))||(j == new(i,3))||(j == new(i,4))||(j == new(i,5))||(j == new(i,6))
            
        else
            n_vector = inv([x(v(j,1),:).'-x(i,:).',x(v(j,2),:).'-x(i,:).',x(v(j,3),:).'-x(i,:).' ])*normal_matrix(i,:).';
            n1= n_vector(1);
            n2= n_vector(2);
            n3= n_vector(3);
            if (n1>=0) &(n2>=0)&(n3>=0)& (n1+n2+n3 ~=0)
                count = count+1;
                a = 1/(n1+n2+n3);
                if count >1
                    if a_old>a
                        opposite_faces(i,1) = j;     
                    end
                else
                    opposite_faces(i,1) = j;  
                    a_old = a;
                end
            end
        end
    end    
end
%the above is to find the opposite point