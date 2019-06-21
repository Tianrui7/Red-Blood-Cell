total_area_new_method = zeros(1,7);
facet_area = zeros (1,7);
for i = 7
    nref = i;
    %sphere
    ellipsoid
    
    new = zeros(nv,6);
    for i_try1 = 1:nv
        new(i_try1,:) = findface(v,nt,i_try1);
    end
    simple = [2 3 1];

    area_vertex = zeros (nv,3);
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
    end

    for i_try1 = 1:nv
        total_area_new_method (i) = total_area_new_method(i) + norm(area_vertex(i_try1,:));
    end
    
    %total_area_new_method (i) = total_area_new_method(i) - 4*pi;
    total_area_new_method (i) = total_area_new_method(i) - 4*pi*((3^1.6+6^1.6+4.5^1.6)/3)^(1/1.6);
    total_area_new_method(i) = log2(total_area_new_method(i));

    A_mat = zeros (1, nt);
    for i_mat = 1:nt
        x1_mat = v(i_mat, 1);
        x2_mat = v(i_mat, 2);
        x3_mat = v(i_mat, 3);
        A_mat(i_mat) = norm(0.5*(cross(x(x1_mat,:),x(x2_mat,:)) + cross(x(x2_mat,:),x(x3_mat,:)) + cross(x(x3_mat,:),x(x1_mat,:))));
    end
    for i_area = 1:nt
        facet_area(i) = facet_area(i) + A_mat(i_area);
    end
    %facet_area (i) = facet_area (i) - 4*pi;
    facet_area (i) = facet_area (i) - 4*pi*((3^1.6+6^1.6+4.5^1.6)/3)^(1/1.6);
    facet_area(i) = log2(facet_area(i));
end
hold on
nref_matrix = [1, 2, 3, 4, 5,6,7];
plot(nref_matrix, total_area_new_method,'b');
plot(nref_matrix, facet_area,'k');
%fplot(@(t) 4*pi,[1 6],'m' )
title('total surface area')
xlabel('number of refinements') 
ylabel('surface area')
%legend({'new method','facet based_area','continuous surface_area' },'Location','southwest')
legend({'new method','facet based_area'},'Location','southwest')
hold off