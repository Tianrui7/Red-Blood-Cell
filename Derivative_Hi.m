function der_hi = Derivative_Hi (x,v,nt,k,k_,A)

%der hk_/der hk
der_hi = zeros(3,3);

%First, we want to find faces that are connected with both k and k_
connect = [0 0];
j_N=1;
%j is the index for connect
for i_N= 1:nt
    %i is the index of the face that we are checking
    if (v(i_N,1)==k||v(i_N,2)==k||v(i_N,3)==k)&&(v(i_N,1)==k_||v(i_N,2)==k_||v(i_N,3)==k_)
        connect(j_N) = i_N;
        j_N=j_N+1;
    end
end

x1 = 0;
x2 = 0;

for i_N=1:3
    if v(connect(1),i_N)== k
        simple = [2 3 1];
        if v(connect(1),simple(i_N)) == k_
            x1 = v(connect(1),simple(simple(i_N)));
            %x1 is the index of the point x1
        else
            temp = connect(1);
            connect(1) = connect(2);
            connect(2) = temp;
            x2 = v(connect(2),simple(i_N));
        end
    end  
end

if x1 == 0
    for i_N = 1:3
        if (v(connect(1),i_N)~=k)&&(v(connect(1),i_N)~=k_)
            x1 = v(connect(1),i_N);
        end
    end
end

if x2 == 0
    for i_N = 1:3
        if (v(connect(2),i_N)~=k)&&(v(connect(2),i_N)~=k_)
            x2 = v(connect(2),i_N);
        end
    end
end


    

%connect(1 or 2) is the index of faces  
A_area_Hi = 0.5*(cross(x(k,:),x(k_,:))+cross(x(k_,:),x(x1,:))+cross(x(x1,:),x(k,:)));    
B1 = A_area_Hi(1);
B2 = A_area_Hi(2);
B3 = A_area_Hi(3);
   
der_hi(1,:) = der_hi(1,:)+ 0.5*(B2*(x(k,3)-x(x1,3)) + B3*(x(x1,2) - x(k,2)))*Derivative_Ai (x,v,connect(1),k,A)/(A(connect(1)))^2;
der_hi(1,1) = der_hi(1,1) - ((x(x1,3)-x(k_,3))*(x(k,3)-x(x1,3)) + (x(k_,2)-x(x1,2))*(x(x1,2)-x(k,2)))/(4*A(connect(1)));
der_hi(1,2) = der_hi(1,2) - (- B3 + 0.5*(x(x1,1)-x(k_,1))*(x(x1,2)-x(k,2)) )/(2*A(connect(1)));
der_hi(1,3) = der_hi(1,3) - (B2 + 0.5*(x(k,3)-x(x1,3))*(x(k_,1)-x(x1,1)) )/(2*A(connect(1)));
    
der_hi(2,:) = der_hi(2,:)+ 0.5*(B1*(x(x1,3)-x(k,3))+B3*(x(k,1)-x(x1,1)))*Derivative_Ai (x,v,connect(1),k,A)/(A(connect(1)))^2;
der_hi(2,1) = der_hi(2,1) - (B3 + 0.5*(x(k,1)-x(x1,1))*(x(k_,2)-x(x1,2)) )/(2*A(connect(1)));
der_hi(2,2) = der_hi(2,2) - ((x(k_,3)-x(x1,3))*(x(x1,3)-x(k,3)) + (x(x1,1)-x(k_,1))*(x(k,1)-x(x1,1)) )/(4*A(connect(1)));
der_hi(2,3) = der_hi(2,3) - (-B1+ 0.5*(x(x1,3)-x(k,3))*(x(x1,2)-x(k_,2)))/(2*A(connect(1)));
    
der_hi(3,:) = der_hi(3,:)+ 0.5*(B1*(x(k,2)-x(x1,2))+B2*(x(x1,1)-x(k,1)))*Derivative_Ai (x,v,connect(1),k,A)/(A(connect(1)))^2;
der_hi(3,1) = der_hi(3,1) - (-B2+ 0.5*(x(x1,1)-x(k,1))*(x(x1,3)-x(k_,3)))/(2*A(connect(1)));
der_hi(3,2) = der_hi(3,2) - (B1 + 0.5*(x(k,2)-x(x1,2))*(x(k_,3)-x(x1,3)))/(2*A(connect(1)));
der_hi(3,3) = der_hi(3,3) - ((x(x1,2)-x(k_,2))*(x(k,2)-x(x1,2))+ (x(k_,1)-x(x1,1))*(x(x1,1)-x(k,1)))/(4*A(connect(1)));



A_area_Hi = 0.5*(cross(x(k,:),x(x2,:))+cross(x(x2,:),x(k_,:))+cross(x(k_,:),x(k,:))); 
C1 = A_area_Hi(1);
C2 = A_area_Hi(2);
C3 = A_area_Hi(3);

der_hi(1,:) = der_hi(1,:)+ 0.5*(C2*(x(x2,3)-x(k,3)) + C3*(x(k,2) - x(x2,2)))*Derivative_Ai (x,v,connect(2),k,A)/(A(connect(2)))^2;
der_hi(1,1) = der_hi(1,1) - ((x(k_,3)-x(x2,3))*(x(x2,3)-x(k,3)) + (x(x2,2)-x(k_,2))*(x(k,2)-x(x2,2)))/(4*A(connect(2)));
der_hi(1,2) = der_hi(1,2) - (C3 + 0.5*(x(k,2)-x(x2,2))*(x(k_,1)-x(x2,1)))/(2*A(connect(2)));
der_hi(1,3) = der_hi(1,3) - (-C2 + 0.5*(x(x2,1)-x(k_,1))*(x(x2,3)-x(k,3)))/(2*A(connect(2)));

der_hi(2,:) = der_hi(2,:)+ 0.5*(C1*(x(k,3)-x(x2,3)) + C3*(x(x2,1) - x(k,1)))*Derivative_Ai (x,v,connect(2),k,A)/(A(connect(2)))^2;
der_hi(2,1) = der_hi(2,1) - (-C3+ 0.5*(x(x2,1)-x(k,1))*(x(x2,2)-x(k_,2)))/(2*A(connect(2)));
der_hi(2,2) = der_hi(2,2) - ((x(x2,3)-x(k_,3))*(x(k,3)-x(x2,3))+(x(k_,1)-x(x2,1))*(x(x2,1)-x(k,1)))/(4*A(connect(2)));
der_hi(2,3) = der_hi(2,3) - (C1 + 0.5*(x(k,3)-x(x2,3))*(x(k_,2)-x(x2,2)))/(2*A(connect(2)));
    
der_hi(3,:) = der_hi(3,:)+ 0.5*(C1*(x(x2,2)-x(k,2)) + C2*(x(k,1) - x(x2,1)))*Derivative_Ai (x,v,connect(2),k,A)/(A(connect(2)))^2;
der_hi(3,1) = der_hi(3,1) - (C2+0.5*(x(k,1)-x(x2,1))*(x(k_,3)-x(x2,3)))/(2*A(connect(2)));
der_hi(3,2) = der_hi(3,2) - (-C1 + 0.5*(x(x2,2)-x(k,2))*(x(x2,3)-x(k_,3)))/(2*A(connect(2)));
der_hi(3,3) = der_hi(3,3) - ((x(k_,2)-x(x2,2))*(x(x2,2)-x(k,2)) + (x(x2,1)-x(k_,1))*(x(k,1)-x(x2,1)))/(4*A(connect(2)));
    
end
