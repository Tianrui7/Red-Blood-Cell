function der_hk = Derivative_Hk (x,v,nt,k_Hk,A)
%this is to calculate derHk/derXk
%k is the index of given vertex 

der_hk = zeros(3,3);

note_face = findface(v,nt,k_Hk);
if note_face(6) == 0
    size = 5;
else
    size =6;
end
    
for i_hk = 1:size
    %note_face(i_hk) is the index of faces  
    
    simple = [2 3 1];
    for j_der=1:3
        %to find the location of the point k in the ith face
        if k_Hk == v(note_face(i_hk),j_der)
            x1_Hk = v(note_face(i_hk), simple(j_der));
            x2_Hk = v(note_face(i_hk), simple(simple(j_der)));
            %k,x1_Hk,x2_Hk are indexes of vertices of the ith face
            break
        end
    end
    
    A_area = 0.5*(cross(x(k_Hk,:),x(x1_Hk,:))+cross(x(x1_Hk,:),x(x2_Hk,:))+cross(x(x2_Hk,:),x(k_Hk,:)));
    D1 = A_area(1);
    D2 = A_area(2);
    D3 = A_area(3);
    
    
    der_hk(1,:) = der_hk(1,:)+ 0.5*(D2*(x(x2_Hk,3)-x(x1_Hk,3)) + D3*(x(x1_Hk,2) - x(x2_Hk,2)))*Derivative_Ai (x,v,note_face(i_hk),k_Hk,A)/(A(note_face(i_hk)))^2;
    der_hk(1,1) = der_hk(1,1) - ((x(x2_Hk,3)-x(x1_Hk,3))^2 + (x(x1_Hk,2) - x(x2_Hk,2))^2)/(4*A(note_face(i_hk)));
    der_hk(1,2) = der_hk(1,2) - ((x(x2_Hk,1)-x(x1_Hk,1))*(x(x1_Hk,2)- x(x2_Hk,2)))/(4*A(note_face(i_hk)));
    der_hk(1,3) = der_hk(1,3) - ((x(x1_Hk,1)-x(x2_Hk,1))*(x(x2_Hk,3)- x(x1_Hk,3)))/(4*A(note_face(i_hk)));
    
    der_hk(2,:) = der_hk(2,:)+ 0.5*(D1*(x(x1_Hk,3)-x(x2_Hk,3))+D3*(x(x2_Hk,1)-x(x1_Hk,1)))*Derivative_Ai (x,v,note_face(i_hk),k_Hk,A)/(A(note_face(i_hk)))^2;
    der_hk(2,1) = der_hk(2,1) - ((x(x2_Hk,1)-x(x1_Hk,1))*(x(x1_Hk,2)- x(x2_Hk,2)))/(4*A(note_face(i_hk)));
    der_hk(2,2) = der_hk(2,2) - ((x(x2_Hk,1)-x(x1_Hk,1))^2 + (x(x1_Hk,3) - x(x2_Hk,3))^2)/(4*A(note_face(i_hk)));
    der_hk(2,3) = der_hk(2,3) - ((x(x2_Hk,2)-x(x1_Hk,2))*(x(x1_Hk,3)- x(x2_Hk,3)))/(4*A(note_face(i_hk)));
    
    der_hk(3,:) = der_hk(3,:)+ 0.5*(D1*(x(x2_Hk,2)-x(x1_Hk,2))+D2*(x(x1_Hk,1)-x(x2_Hk,1)))*Derivative_Ai (x,v,note_face(i_hk),k_Hk,A)/(A(note_face(i_hk)))^2;
    der_hk(3,1) = der_hk(3,1) - ((x(x2_Hk,3)-x(x1_Hk,3))*(x(x1_Hk,1)- x(x2_Hk,1)))/(4*A(note_face(i_hk)));
    der_hk(3,2) = der_hk(3,2) - ((x(x1_Hk,3)-x(x2_Hk,3))*(x(x2_Hk,2)- x(x1_Hk,2)))/(4*A(note_face(i_hk)));
    der_hk(3,3) = der_hk(3,3) - ((x(x2_Hk,2)-x(x1_Hk,2))^2 + (x(x1_Hk,1) - x(x2_Hk,1))^2)/(4*A(note_face(i_hk)));
end
