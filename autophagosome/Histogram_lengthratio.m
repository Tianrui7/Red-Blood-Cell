length_ratio = zeros(1,nt);
for i_ratio=1:nt
    x1_hf = x(v(i_ratio,1),:);
    x2_hf = x(v(i_ratio,2),:);
    x3_hf = x(v(i_ratio,3),:);
    length_1 = norm(x1_hf-x2_hf);
    length_2 = norm(x2_hf-x3_hf);
    length_3 = norm(x3_hf-x1_hf);
    
    if length_1<=length_2
        max_length = length_2;
        min_length = length_1;
    else
        max_length = length_1;
        min_length = length_2;
    end
    
    if max_length<=length_3
        max_length = length_3;
    end
    
    if min_length>=length_3
        min_length = length_3;
    end
    
    length_ratio(i_ratio) = max_length/min_length;
end

histogram(length_ratio,10)