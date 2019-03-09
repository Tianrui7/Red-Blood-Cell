%test energy program
sphere

%z_test = [ rand rand rand];

ebs = 0.00001;
H_curva
%Hi is at 220th point
Fold = H_curvature(220,:);

derder = Derivative_Hi (x,v,nt,442,220,A_mat);

x(442,:) = x(442,:) + ebs*z_test;
H_curva
Fnew = H_curvature(220,:);

derderder=0;
for i_test = 1:3
    derderder =derderder + z_test(i_test)*(derder(:,i_test).');
end

result_test = Fnew - Fold - ebs*derderder

