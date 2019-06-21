%find the three adjoint faces for each face
adjoint_face = zeros(nt,3);
for i_adjoint = 1:nt
    count_adjoint = 1;
    for j_adjoint = 1:nt
        if ismember([v(i_adjoint,1),v(i_adjoint,2),v(i_adjoint,3)],[v(j_adjoint,1),v(j_adjoint,2),v(j_adjoint,3)]) == [1,1,0]
            adjoint_face(i_adjoint,count_adjoint) = j_adjoint;
            count_adjoint = count_adjoint +1;
        else
            if ismember([v(i_adjoint,1),v(i_adjoint,2),v(i_adjoint,3)],[v(j_adjoint,1),v(j_adjoint,2),v(j_adjoint,3)]) == [1,0,1]
                adjoint_face(i_adjoint,count_adjoint) = j_adjoint;
                count_adjoint = count_adjoint +1;
            else
                if ismember([v(i_adjoint,1),v(i_adjoint,2),v(i_adjoint,3)],[v(j_adjoint,1),v(j_adjoint,2),v(j_adjoint,3)]) == [0,1,1]
                    adjoint_face(i_adjoint,count_adjoint) = j_adjoint;
                    count_adjoint = count_adjoint +1;
                end
            end
        end
    end    
end
