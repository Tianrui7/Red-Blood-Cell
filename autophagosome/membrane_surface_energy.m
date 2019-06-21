total_surface_energy = 0;
surface_energy = zeros(1,nv);

opposite_point=zeros(nv,4); 
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
                        opposite_point(i,1:3) = x(i,:)+a*normal_matrix(i,:);
                        opposite_point(i,4) =a;
                    end
                else
                    opposite_point(i,1:3) = x(i,:)+a*normal_matrix(i,:);
                    opposite_point(i,4) =a;
                    a_old = a;
                end
            end
        end
    end    
end
constant1 = 5; %A
constant2 = 6; %B
constant3 = 20; %k

for i = 1: nv
    if opposite_point(i,4)==0
    else
        surface_energy (i) = - constant1/(opposite_point(i,4)^2) + constant2*constant3*exp(-constant3*opposite_point(i,4));
    end
    total_surface_energy = total_surface_energy + surface_energy (i);
end
total_surface_energy =  total_surface_energy/20;