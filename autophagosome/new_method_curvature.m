standard_deviation_error = zeros(1,6);
for i = 1:6
    nref = i;
   
    ellipsoid
    new = zeros(nv,6);
    for i_try1 = 1:nv
        new(i_try1,:) = findface(v,nt,i_try1);
    end
    simple = [2 3 1];

    area_vertex = zeros (nv,3);
    curvature_vertex = zeros (nv,3);
    normal_vertex = zeros (nv,3);
    
    for i_try1 = 1:nv
        Note = new(i_try1,:);
        if Note(6) == 0
            size =5;
        else
            size =6;
        end
        for j = 1:size

            for k_ode=1:3
               %to find the location of the given point in the jth face
               if i_try1 == v(Note(j),k_ode)
                   x1 = v(Note(j), simple(k_ode));
                   x2 = v(Note(j), simple(simple(k_ode)));
                   %x1 and x2 are the index of two points that are on the jth face with the given point
                   break
               end
            end
            area_vertex(i_try1,:) = area_vertex(i_try1,:)+ cross(x(x1,:),x(x2,:))/6;
            %area at k-th vertex
        end
        normal_vertex(i_try1,:) = area_vertex(i_try1,:);
        normal_vertex(i_try1,:) = normal_vertex(i_try1,:)/norm(normal_vertex(i_try1,:));
    end
    
    for i_try1 = 1:nv
        Note = new(i_try1,:);
        if Note(6) == 0
            size =5;
        else
            size =6;
        end
        for j = 1:size

            for k_ode=1:3
               %to find the location of the given point in the jth face
               if i_try1 == v(Note(j),k_ode)
                   x1 = v(Note(j), simple(k_ode));
                   x2 = v(Note(j), simple(simple(k_ode)));
                   %x1 and x2 are the index of two points that are on the jth face with the given point
                   break
               end
            end
            curvature_vertex(i_try1,:) = curvature_vertex(i_try1,:) + cross(normal_vertex(x2,:)+normal_vertex(x1,:),x(x2,:)-x(x1,:));
            
        end
        curvature_vertex(i_try1,:) = -curvature_vertex(i_try1,:)/(6*norm(area_vertex(i_try1,:)));
    end
    principal_cur = zeros(1,nv);
    for j_try = 1:nv
        principal_cur(j_try) = norm(curvature_vertex(j_try,:));
    end

    sphere
    continue_cur = zeros (1,nv);
    for i_ell = 1:nv
        x(i_ell,1) = x(i_ell,1)*2;
        x(i_ell,2) = x(i_ell,2)*1.5;
        x(i_ell,3) = x(i_ell,3)*3;
        continue_cur(i_ell) = abs(x(i_ell,1)^2+x(i_ell,2)^2+x(i_ell,3)^2-15.25)/(81*(x(i_ell,1)^2/16+x(i_ell,2)^2/5.0625+x(i_ell,3)^2/81)^1.5 );  
    end
    
    for j_try = 1:nv
        standard_deviation_error(i)=standard_deviation_error(i)+ (principal_cur(j_try)-continue_cur(j_try))^2;
    end
    
   instant_value = standard_deviation_error(i);
   
    standard_deviation_error(i) = (instant_value/nv)^0.5;
end
nref_matrix = [1, 2, 3, 4, 5,6];
plot(nref_matrix, standard_deviation_error,'b');

%histogram (principal_cur)