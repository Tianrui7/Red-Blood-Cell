sphere
Vmin =2*Volume(v,nt,x)/3;
Vmax =Volume(v,nt,x);
a_plot = 0.5;

new = zeros(nv,6);
for i_try1 = 1:nv
    new(i_try1,:) = findface(v,nt,i_try1);
end

%hold on
%fplot(@(t) Vmin + (Vmax-Vmin)*exp((-a_plot)*t),[0 22],'m' )
%ylim([0 5]);

%for i_plot= 1:1
    
 %   sphere
    
    colorplot = ['b','k','r'] ;
    % blue,black,red
    t_range = [0 22];
    %can change the number 5
    dt = (0.1)^(4);
    %can change the number 0.1
    
    %dimen = int16((t_range(2)/dt)+1);
    dimen = 220000;
    %to tract the change of volume and area 
    v_whole  = zeros(dimen ,1);
    a_whole  = zeros(dimen,1);
    %e_whole = zeros(dimen,1);
    t_matrix = zeros(dimen,1);

    t=0;
    f_RHS = @(x,t)practice_ode(v,nt,nv,x,t,Vmin, Vmax, a_plot, new);
    
    t_matrix(1,1) = t;
    v_whole(1,1) = Volume(v,nt,x);
    a_whole(1,1) = totalarea (v,nt,x);
   % e_whole (1,1) = energy(v,nt,nv,x,new);

    for i_try = 1:(dimen-1)

        x = x + dt*f_RHS(x,t);
        t = t+dt;
        t_matrix(i_try+1,1) = t;
        v_whole(i_try+1,1) = Volume(v,nt,x);
        a_whole(i_try+1,1) = totalarea(v,nt,x);
        %e_whole (i_try1+1,1) = energy(v,nt,nv,x,new);
    end
   
    %plot(t_matrix, v_whole, colorplot(i_plot))
    
   % plot(t_matrix, a_whole,colorplot(i_plot))
    
  %  plot(t_matrix, e_whole, colorplot(i_plot))   
%end

%hold off
  