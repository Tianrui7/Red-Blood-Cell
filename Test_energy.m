%test energy program
sphere
%{
z_test= zeros(nv,3);
for i_test = 1: nv
    z_test(i_test,:) = [ rand rand rand];
end
%}
ebs = 0.0000001;
Energy
Fold = Energy_res;

Derivative_E

x = x + ebs*z_test;
Energy
Fnew = Energy_res;

result_test=0;
for i_test = 1:nv
    result_test = result_test + dot (Der_E(i_test,:),z_test(i_test,:));
end

result_test = Fnew - Fold - ebs*result_test

