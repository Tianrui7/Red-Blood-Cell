%This is to calculate der|Ai|/derXk, given the ith triangle and kth vertex

function der_Ai = Derivative_Ai (x,v,i_face,k_derAi,A)

%nt is the number of faces 
simple = [2 3 1];
for j_der=1:3
    %to find the location of the point k in the ith face
    if k_derAi == v(i_face,j_der)
        x1_derAi = v(i_face, simple(j_der));
        x2_derAi = v(i_face, simple(simple(j_der)));
        %k,x1,x2 are indexes of vertices of the ith face
        break
    end
end

A_area = 0.5*(cross(x(k_derAi,:),x(x1_derAi,:))+cross(x(x1_derAi,:),x(x2_derAi,:))+cross(x(x2_derAi,:),x(k_derAi,:)));
der_Ai = (0.5/A(i_face))*(A_area(1)*[0,x(x1_derAi,3)-x(x2_derAi,3),-x(x1_derAi,2)+x(x2_derAi,2)] + A_area(2)*[-x(x1_derAi,3)+x(x2_derAi,3),0, x(x1_derAi,1)-x(x2_derAi,1)] + A_area(3)*[x(x1_derAi,2)-x(x2_derAi,2), -x(x1_derAi,1)+x(x2_derAi,1),0]);


