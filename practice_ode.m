function ode_calc = practice_ode(v,nt,nv,x,t,Vmin, Vmax, a_temp, new)

left  = zeros(2,2);
right = zeros(2,1);

%The following is to calculate sum_an for every point
simple = [2 3 1];
Sum_an = zeros(nv,3);
for i = 1:nv
    Note = new(i,:);
    if Note(6) == 0
        size =5;
    else
        size =6;
    end
    an_xk = zeros(size,3);
    %This is a matrix that stores af*nf for every f that is connected to xk
    for i_SUM = 1:size
        for j_SUM=1:3
            %to find the location of the point i in the ith face
            if i == v(Note(i_SUM),j_SUM)
                x1_SUM = v(Note(i_SUM), simple(j_SUM));
                x2_SUM = v(Note(i_SUM), simple(simple(j_SUM)));
                %x1_SUM and x2_SUM are the index of two points that are on the jth face with the given point
                break
            end
        end
        an_xk(i_SUM,:) = 0.5*(cross(x(i,:),x(x1_SUM,:))+cross(x(x1_SUM,:),x(x2_SUM,:))+cross(x(x2_SUM,:),x(i,:)));
        Sum_an(i,:) = Sum_an(i,:) + an_xk(i_SUM,:);    
    end
end
%The above is to calculate sum_an for every point


%this large for loop is to calculate derivative A, derivative V,
%derivative_E
DerA = zeros (nv,3);
DerV = zeros (nv,3);

Der_length = zeros(nv,3);

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
        %{
        %Now we want to calculate derivative E
        %the index of the face is New(j)
        area_calc = sqrt(s_ode*(s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode)); 
        %area is the area for the New(j)th face
        coeff = (1/(2*sqrt(s_ode*(s_ode-a_ode)*(s_ode-b_ode)*(s_ode-c_ode))));
        coeff_a = -s_ode*(s_ode-b_ode)*(s_ode-c_ode);
        coeff_c = -s_ode*(s_ode-a_ode)*(s_ode-b_ode);
        Der_af = coeff*(coeff_s*0.5*(-(x(x1,:)-x(i,:))/a_ode+(x(i,:)-x(x2,:))/c_ode)-coeff_a*(x(x1,:)-x(i,:))/a_ode+coeff_c*(x(i,:)-x(x2,:))/c_ode);
   
        
        %THe following is to calculate derivative_ada
        Der_ada = zeros (1,3);
        Nfi = Sum_an(i,:);
        Nfi = Nfi/norm(Nfi);
        Nf1 = Sum_an(x1,:);
        Nf1 = Nf1/norm(Nf1);
        Nf2 = Sum_an(x2,:);
        Nf2 = Nf2/norm(Nf2);

        coeff1 = 0.25*(cross(Nfi,x(x1,:)-x(x2,:))+cross(Nf1,x(x2,:)-x(i,:))+cross(Nf2,x(i,:)-x(x1,:)));
        coeff2 = 0.25*(cross(x(i,:),x(x1,:))+cross(x(x1,:),x(x2,:))+cross(x(x2,:),x(i,:)));

        Der_nf1 = Derivative_N(v,nt,x,i,x1,Sum_an);
        Der_nf2 = Derivative_N(v,nt,x,i,x2,Sum_an);

        Der_ada(1) = dot(coeff1,cross([1 0 0],x(x1,:)-x(x2,:)))+dot(coeff2,cross(Nf2-Nf1,[1 0 0])+cross(Der_nf1(:,1).',x(x2,:)-x(i,:))+cross(Der_nf2(:,1).',x(i,:)-x(x1,:)));
        Der_ada(2) = dot(coeff1,cross([0 1 0],x(x1,:)-x(x2,:)))+dot(coeff2,cross(Nf2-Nf1,[0 1 0])+cross(Der_nf1(:,2).',x(x2,:)-x(i,:))+cross(Der_nf2(:,2).',x(i,:)-x(x1,:)));
        Der_ada(3) = dot(coeff1,cross([0 0 1],x(x1,:)-x(x2,:)))+dot(coeff2,cross(Nf2-Nf1,[0 0 1])+cross(Der_nf1(:,3).',x(x2,:)-x(i,:))+cross(Der_nf2(:,3).',x(i,:)-x(x1,:)));
        %The above is to calculate derivative_ada
        
        %The following is to calculate hf for the face New(j)
        xi_hf = x(i,:);
        x1_hf = x(x1,:);
        x2_hf = x(x2,:);
        h_f = (0.25*dot((cross(xi_hf,x1_hf)+cross(x1_hf,x2_hf)+cross(x2_hf,xi_hf)),(cross(Nfi,x1_hf-x2_hf)+cross(Nf1, x2_hf-xi_hf)+cross(Nf2, xi_hf-x1_hf))))/(area_calc)^2;
        %The above is to calculate hf for the face New(j)
        
        %}
        
        side_length = 1.5*norm(x(i,:)-x(x1,:));
        %The following is to calculate derivative of sum of squared length
        Der_length(i,:) = Der_length(i,:) + 2*side_length*(x(i,:)-x(x1,:));
        %The above is to calculate derivative of sum of squared length
        
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
    %right(1,1) = right(1,1) - 0.5*dot(DerV(i_calc,:),DerE(i_calc,:));
    right(1,1) = right(1,1) - 0.5*dot(DerV(i_calc,:),DerE(i_calc,:)+500*Der_length(i_calc,:));
end

right(1,1) = right(1,1) - Derivative_V0(Vmin, Vmax, a_temp, t);


for i_calc = 1:nv
    right(2,1) = right(2,1) - 0.5*dot(DerA(i_calc,:),DerE(i_calc,:)+500*Der_length(i_calc,:));
    %right(2,1) = right(2,1) - 0.5*dot(DerA(i_calc,:),DerE(i_calc,:));
    
end

mu_lambda = left\right;

lam = mu_lambda(1,1);
mu  = mu_lambda(2,1);

ode_calc = -0.5*DerE - 0.5*500*Der_length - lam*DerV - mu*DerA;
%ode_calc = -0.5*DerE - lam*DerV - mu*DerA;
%to t = 10, it is 500, now change to 100
%remember to change the above, there are 2 500 above
