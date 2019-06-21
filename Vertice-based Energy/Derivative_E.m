Der_E = zeros(nv,3);
simple = [2 3 1];

H_curva
A_matrix

%this is to construct ak and d(ak)/d(xk) for each k
Ak = zeros(nv);
derAk = zeros (nv,3);
for i_E = 1:nv
    con_derE = findface(v,nt,i_E);
    if con_derE(6) == 0
        size =5;
    else
        size = 6;
    end
    for j_E = 1:size
        Ak(i_E) = Ak(i_E) + A_mat(con_derE(j_E))/3;
        derAk (i_E,:) = derAk (i_E,:) + Derivative_Ai (x,v,con_derE(j_E),i_E,A_mat)/3;
    end 
end
%this is to construct ak and d(ak)/d(xk) for each k
    
for i_E = 1:nv
    
    derHk = Derivative_Hk (x,v,nt,i_E,A_mat);
    Der_E (i_E,:) = (2*H_curvature(i_E,1)*derHk(1,:)+2*H_curvature(i_E,2)*derHk(2,:)+2*H_curvature(i_E,3)*derHk(3,:))/Ak(i_E) + (((norm(H_curvature(i_E,:)))^2)*derAk(i_E,:))/(Ak(i_E))^2;  
    
    con_derE = findface(v,nt,i_E);
    if con_derE(6) == 0
        size =5;
    else
        size = 6;
    end
    for j_E = 1:size
        for i_practice = 1:3
            if v(con_derE(j_E),i_practice) == i_E
                x_ = v(con_derE(j_E),simple(i_practice));
                break
            end
        end
        %x_is the index of point that is next to point x(i_E), on face(con_derE(j_E))
          
        %this is to calculate d(a(x_))/d(xk)
        con_derE_ = j_E -1;
        if j_E == 1
            con_derE_ = size;
        end
        %con_derE_ is the index of the face to the right of face con_derE(j_E)
        derAi = (Derivative_Ai (x,v,con_derE(j_E),i_E,A_mat) + Derivative_Ai (x,v,con_derE(con_derE_),i_E,A_mat))/3;
        %this is to calculate d(ai)/d(xk)
        
        derHi = Derivative_Hi (x,v,nt,i_E,x_,A_mat);
        Der_E (i_E,:) = Der_E (i_E,:) + (2*H_curvature(x_,1)*derHi(1,:)+2*H_curvature(x_,2)*derHi(2,:)+2*H_curvature(x_,3)*derHi(3,:))/Ak(x_)+(((norm(H_curvature(x_,:)))^2)*derAi)/(Ak(x_))^2;  
    end
end
