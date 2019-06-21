function ode_calc = practice_ode_2(v,nt,nv,x,t,Vmin, Vmax, a_temp, new,adjoint_matrix, opposite_face_mat)

left  = zeros(2,2);
right = zeros(2,1);
simple = [2 3 1];

%This following is to find the corresponding point on the opposite side 
%first we want to calculate Nk = sum a*n/|sum a*n|
%new = [nv,6] is the findface matrix 
normal_matrix = zeros(nv,3);
%normal_matrix stores the value of Nk for every point
for i = 1:nv
    Note = new(i,:);
    if Note(6) == 0
        size =5;
    else
        size =6;
    end
    for j = 1:size
        %new(j) is the index of the face
        normal_matrix(i,:) = normal_matrix(i,:)+ (cross(x(v(Note(j),1),:),x(v(Note(j),2),:))+cross(x(v(Note(j),2),:),x(v(Note(j),3),:))+cross(x(v(Note(j),3),:),x(v(Note(j),1),:)));
    end
    normal_matrix(i,:) = -normal_matrix(i,:)/norm(normal_matrix(i,:));
    %there is a negative sign here since we want the pointing inward normal
    %inside of outward normal 
end

%opposite_point(i,1-3) stores coordinates of the opposite point of ith point
%and opposite_point(i,4) stores a
opposite_point=zeros(nv,4); 
for i =1:nv
    if opposite_face_mat(i) == 0
        count = 0;
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
      
    else
        j_now = opposite_face_mat(i);
        n_vector = inv([x(v(j_now,1),:).'-x(i,:).',x(v(j_now,2),:).'-x(i,:).',x(v(j_now,3),:).'-x(i,:).' ])*normal_matrix(i,:).';
        n1= n_vector(1);
        n2= n_vector(2);
        n3= n_vector(3);
        count = 0;
        if (n1>=0) &(n2>=0)&(n3>=0)& (n1+n2+n3 ~=0)
            count = count+1;
        else
            for j =1:3
                j_now = adjoint_matrix (opposite_face_mat(i),j);
                if (j_now == new(i,1)) || (j_now == new(i,2))||(j_now == new(i,3))||(j_now == new(i,4))||(j_now == new(i,5))||(j_now == new(i,6))
                
                else      
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
                    
                    if (j_now == new(i,1)) || (j_now == new(i,2))||(j_now == new(i,3))||(j_now == new(i,4))||(j_now == new(i,5))||(j_now == new(i,6))
                    else
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
                end
            end
        end
        
        if count ==0
            opposite_face_mat(i) = 0;
        else
            a = 1/(n1+n2+n3);
            opposite_point(i,1:3) = x(i,:)+a*normal_matrix(i,:);
            opposite_point(i,4) =a;
        end
    end
end
%the above is to find the opposite point

%This following is to calculate af, for every face f, for any given matrix x
area = zeros(1,nt);
for i_area=1:nt
    x1_hf = x(v(i_area,1),:);
    x2_hf = x(v(i_area,2),:);
    x3_hf = x(v(i_area,3),:);
    s_area = ((norm(x2_hf-x1_hf))+(norm(x3_hf-x2_hf))+(norm(x1_hf-x3_hf)))/2;
    %s_area = vpa (s_area,12);
    area(i_area) = sqrt(s_area*(s_area-(norm(x2_hf-x1_hf)))*(s_area-(norm(x3_hf-x2_hf)))*(s_area-(norm(x1_hf-x3_hf)))); 
   % area(i_area) = vpa (area(i_area),12);
end
%the above is to calculate triangle area

%this large for loop is to calculate derivative A, derivative V,
%derivative_E
DerA = zeros (nv,3);
DerV = zeros (nv,3);

Der_length = zeros(nv,3);

constant1 = 5; %A
constant2 = 6; %B
constant3 = 20; %k
bilayer_der = zeros(nv,3);
Area_around_point = zeros(nv,1);

for i=1:nv
    %i is the index of point/vertice
    New = new(i,:);
    if New(6) == 0
        size = 5;
    else
        size = 6;
    end
    for j=1:size
        %j is the counting of faces that are connected to the given point
        for k_ode=1:3
           %to find the location of the given point in the jth face
           if i == v(New(j),k_ode)
               x1 = v(New(j), simple(k_ode));
               x2 = v(New(j), simple(simple(k_ode)));
               %x1 and x2 are the index of two points that are on the jth face with the given point
               break
           end
        end
        a_ode=norm(x(x1,:)-x(i,:));
        b_ode=norm(x(x2,:)-x(x1,:));
        c_ode=norm(x(i,:)-x(x2,:));
        s_ode=(a_ode+b_ode+c_ode)/2;
        coeff_s = (s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode)+s_ode*(s_ode-b_ode)*(s_ode-c_ode)+s_ode*(s_ode-a_ode)*(s_ode-b_ode)+s_ode*(s_ode-a_ode)*(s_ode-c_ode);
        DerA(i,:)= DerA(i,:) + (1/(2*sqrt(s_ode*(s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode))))*(coeff_s*0.5*(-(x(x1,:)-x(i,:))/a_ode+(x(i,:)-x(x2,:))/c_ode)-(-s_ode*(s_ode-b_ode)*(s_ode-c_ode))*(x(x1,:)-x(i,:))/a_ode+(-s_ode*(s_ode-a_ode)*(s_ode-b_ode))*(x(i,:)-x(x2,:))/c_ode);
        DerV(i,:)= DerV(i,:) + (1/6)*cross(x(x1,:),x(x2,:));
        
        side_length = 1.5*norm(x(i,:)-x(x1,:));
        %The following is to calculate derivative of sum of squared length
        Der_length(i,:) = Der_length(i,:) + 2*side_length*(x(i,:)-x(x1,:));
        %The above is to calculate derivative of sum of squared length
        
        Area_around_point(i) = Area_around_point(i) + area(New(j));
        if (opposite_point(x1,4) == 0) & (opposite_point(x2,4) ~= 0)
            bilayer_der(i,:) = bilayer_der (i,:) + 0.5*(-constant1/(12*pi)*(1/(opposite_point(x2,4))^2) + constant2*constant3*(exp(-constant3*opposite_point(x2,4) )))*(1/(2*sqrt(s_ode*(s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode))))*(coeff_s*0.5*(-(x(x1,:)-x(i,:))/a_ode+(x(i,:)-x(x2,:))/c_ode)-(-s_ode*(s_ode-b_ode)*(s_ode-c_ode))*(x(x1,:)-x(i,:))/a_ode+(-s_ode*(s_ode-a_ode)*(s_ode-b_ode))*(x(i,:)-x(x2,:))/c_ode); 
        else
            if (opposite_point(x2,4) == 0) & (opposite_point(x1,4) ~= 0)
                bilayer_der(i,:) = bilayer_der (i,:) + 0.5*(-constant1/(12*pi)*(1/(opposite_point(x1,4))^2) + constant2*constant3*(exp(-constant3*opposite_point(x1,4) )))*(1/(2*sqrt(s_ode*(s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode))))*(coeff_s*0.5*(-(x(x1,:)-x(i,:))/a_ode+(x(i,:)-x(x2,:))/c_ode)-(-s_ode*(s_ode-b_ode)*(s_ode-c_ode))*(x(x1,:)-x(i,:))/a_ode+(-s_ode*(s_ode-a_ode)*(s_ode-b_ode))*(x(i,:)-x(x2,:))/c_ode); 
            else
                if (opposite_point(x2,4) == 0) & (opposite_point(x1,4) == 0)
         
                else
                    bilayer_der (i,:) = bilayer_der (i,:) + 0.5*(-constant1/(12*pi)*(1/(opposite_point(x1,4))^2+1/(opposite_point(x2,4))^2) + constant2*constant3*(exp(-constant3*opposite_point(x1,4))+exp(-constant3*opposite_point(x2,4) )))*(1/(2*sqrt(s_ode*(s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode))))*(coeff_s*0.5*(-(x(x1,:)-x(i,:))/a_ode+(x(i,:)-x(x2,:))/c_ode)-(-s_ode*(s_ode-b_ode)*(s_ode-c_ode))*(x(x1,:)-x(i,:))/a_ode+(-s_ode*(s_ode-a_ode)*(s_ode-b_ode))*(x(i,:)-x(x2,:))/c_ode);
                end
            end
        end
    end
    %the following is newly added for bilayer energy
    if opposite_point(i,4) == 0
    else
        bilayer_der (i,:) = bilayer_der (i,:)+ 0.5/6*Area_around_point(i)*(3*constant1*opposite_point(i,4)^(-4)*(x(i,:)-opposite_point(i,1:3))- constant2*(constant3^2)*exp(-constant3*opposite_point(i,4))*(x(i,:)-opposite_point(i,1:3))/opposite_point(i,4));
        bilayer_der (i,:) = bilayer_der (i,:) + 0.5/6*DerA(i,:)*(-constant1*opposite_point(i,4)^(-2)+ constant2*constant3*exp(-constant3*opposite_point(i,4)));
    end
end

Derivative_E
DerE = Der_E;


for i_calc = 1:nv
    left(1,1) = left(1,1) + dot(DerV(i_calc,:),DerV(i_calc,:));
end

for i_calc = 1:nv
    left(1,2) = left(1,2) + dot(DerV(i_calc,:),DerA(i_calc,:));
    left(2,1) = left(1,2);
end

for i_calc = 1:nv
    left(2,2) = left(2,2) + dot(DerA(i_calc,:),DerA(i_calc,:));
end


for i_calc = 1:nv
    right(1,1) = right(1,1) - 0.5*dot(DerV(i_calc,:),DerE(i_calc,:)+100*Der_length(i_calc,:)- bilayer_der(i_calc,:)/20 );
end

%right(1,1) = right(1,1) - Derivative_V0(Vmin, Vmax, a_temp, t);
right(1,1) = right(1,1) - Derivative_V1(a_temp);

for i_calc = 1:nv
    right(2,1) = right(2,1) - 0.5*dot(DerA(i_calc,:),DerE(i_calc,:)+100*Der_length(i_calc,:)- bilayer_der(i_calc,:)/20);   
end

mu_lambda = left\right;

lam = mu_lambda(1,1);
mu  = mu_lambda(2,1);

ode_calc = -0.5*DerE - 0.5*100*Der_length + 0.5*bilayer_der/20 - lam*DerV - mu*DerA;
%it is 100 before
%to t = 10, it is 500, now change to 100
%remember to change the above, there are 2 500 above
