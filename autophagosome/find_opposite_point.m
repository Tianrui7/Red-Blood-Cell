%This program is to find the corresponding point on the opposite side 
%to find X' given X


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
    count = 0;
    for j = 1:nt
        if (j == new(i,1)) || (j == new(i,2))||(j == new(i,3))||(j == new(i,4))||(j == new(i,5))||(j == new(i,6))
            
        else
            n_vector = inv([x(v(j,1),:).'-x(i,:).',x(v(j,2),:).'-x(i,:).',x(v(j,3),:).'-x(i,:).' ])*normal_matrix(i,:).';
            n1= n_vector(1);
            n2= n_vector(2);
            n3= n_vector(3);
            if (n1>=0) &(n2>=0)&(n3>=0)
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

triplot2 
hold on
for i =1:nv
    if rem(i,15) == 0
        pts = [x(i,:);opposite_point(i,1:3)];
        plot3(pts(:,1), pts(:,2), pts(:,3),'b');
    end
end
hold off
    