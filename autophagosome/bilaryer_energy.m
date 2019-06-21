bilayer_energy = 0;
for i_try1= 1:301
    for j_try1=1:341
        if (x(points_above(i_try1),1) - x(points_below(j_try1),1)) <0.5
            distance = norm(x(points_above(i_try1),:) - x(points_below(j_try1),:));
            bilayer_energy = bilayer_energy - constant1/distance^2 + constant2*constant3*exp(-constant3*distance);
        end
    end
end