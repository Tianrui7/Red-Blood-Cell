%This program is to calculate af, for every face f, for any given matrix x
function area = af(v,nt,x)
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
