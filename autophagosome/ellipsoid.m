%to get a ellipsoid from a sphere
sphere
%{
to get a ellipsoid from a sphere, by stretching it in the direction of <3,2,1> and <-1,1,1>
for i = 1:nv
    x(i,:) = x(i,:) + dot([3,2,1] , x(i,:))*[3,2,1]/7;
end

for i = 1:nv
    x(i,:) = x(i,:) + dot([-1,1,1] , x(i,:))*[-1,1,1]/2;
end
%}

for i_ell = 1:nv
    x(i_ell,1) = x(i_ell,1)*2;
    x(i_ell,2) = x(i_ell,2)*1.5;
    x(i_ell,3) = x(i_ell,3)*3;
end
alp = pi/6;
beta = pi/3;
sig = pi/2;
matrix_one = [1,0,0; 0, cos(alp), - sin(alp); 0, sin(alp),cos(alp)];
matrix_two = [cos(beta), 0, sin(beta); 0, 1,0;-sin(beta),0,cos(beta)];
matrix_three = [cos(sig), -sin(sig), 0; sin(sig), cos(sig), 0; 0,0,1];
for i_ell = 1:nv
    x(i_ell,:) = (matrix_three*matrix_two*matrix_one*(x(i_ell,:).')).';
end
triplot

