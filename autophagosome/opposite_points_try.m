%pass in adjoint_matrix, opposite_face_mat)
opposite_point=zeros(nv,4); 
for i =1:nv
    j_now = uint8(opposite_face_mat(i));
    n_vector = inv([x(v(j_now,1),:).'-x(i,:).',x(v(j_now,2),:).'-x(i,:).',x(v(j_now,3),:).'-x(i,:).' ])*normal_matrix(i,:).';
     
    n1= n_vector(1);
    n2= n_vector(2);
    n3= n_vector(3);
    if (n1>=0) &(n2>=0)&(n3>=0)& (n1+n2+n3 ~=0)
    else
        count = 0;
        for j =1:3
            j_now = adjoint_matrix (opposite_face_mat(i),j);
            n_vector = inv([x(v(j_now,1),:).'-x(i,:).',x(v(j_now,2),:).'-x(i,:).',x(v(j_now,3),:).'-x(i,:).' ])*normal_matrix(i,:).';
            n1= n_vector(1);
            n2= n_vector(2);
            n3= n_vector(3);
            if (n1>=0) &(n2>=0)&(n3>=0)& (n1+n2+n3 ~=0)
                count = count +1;
                global opposite_faces
                opposite_faces(i) = j_now;
                break
            end
            
        end    
        if count ==0
            new_option = zeros(6,1);
            pointer = 1;
            for j =1:3
                j_now = adjoint_matrix (opposite_face_mat(i),j);
                for j_second = 1:3
                    if adjoint_matrix (j_now,j_second) == opposite_face_mat(i)
                    else
                        new_option (pointer) = adjoint_matrix (j_now,j_second);
                        pointer = pointer+1;
                    end
                end
            end
            for j =1:6
                j_now = new_option(j);
                n_vector = inv([x(v(j_now,1),:).'-x(i,:).',x(v(j_now,2),:).'-x(i,:).',x(v(j_now,3),:).'-x(i,:).' ])*normal_matrix(i,:).';
                n1= n_vector(1);
                n2= n_vector(2);
                n3= n_vector(3);
                if (n1>=0) &(n2>=0)&(n3>=0)& (n1+n2+n3 ~=0)
                    count = count +1;
                    global opposite_faces
                    opposite_faces(i) = j_now;
                    break
                end
            end       
            if count==0
               % recheck all faces
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
                                    opposite_point(i,1:3) = x(i,:)+a*normal_matrix(i,:);
                                    opposite_point(i,4) =a;
                                    global opposite_faces
                                    opposite_faces(i) = j;
                                end
                            else
                                opposite_point(i,1:3) = x(i,:)+a*normal_matrix(i,:);
                                opposite_point(i,4) =a;
                                a_old = a;
                                global opposite_faces
                                opposite_faces(i) = j;
                            end
                        end
                    end
                end
                %re-check all faces 
            end
        end
    end
    a = 1/(n1+n2+n3);
    opposite_point(i,1:3) = x(i,:)+a*normal_matrix(i,:);
    opposite_point(i,4) =a;
end

