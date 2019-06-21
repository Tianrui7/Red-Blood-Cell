%This program is to calculate total volume, for a given matrix x

function vol_curr = Volume(v,nt,x)
vol_curr = 0;
for i=1:nt
    x1_vol = x(v(i,1),:);
    x2_vol = x(v(i,2),:);
    x3_vol = x(v(i,3),:);  
    vol_curr = vol_curr +  (1/6)*dot(x1_vol, cross(x2_vol,x3_vol));
end
