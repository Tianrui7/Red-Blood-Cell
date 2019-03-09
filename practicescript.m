%profile on
t_range = [0 2];
%can change the number 5
dt = 0.0002;
%0.00001 is test to be a good time step given how it follows the volume and
%area constraint

%can change the number 0.1

sphere
%ellipse

Vmin =2*Volume(v,nt,x)/3 ;
Vmax =Volume(v,nt,x);
a_plot = 0.5;

new = zeros(nv,6);
for i_try1 = 1:nv
    new(i_try1,:) = findface(v,nt,i_try1);
end

t=0;
%}
f_RHS = @(x,t)practice_ode(v,nt,nv,x,t,Vmin, Vmax, a_plot, new);

%vidfile = VideoWriter('Sphere2');
%open(vidfile);
for i_try1 = 1:((t_range(2))/dt)
    x = x + dt*f_RHS(x,t);
    t = t+dt;
  %  triplot
 %   drawnow
  %  frame = getframe(gcf);
  %  writeVideo(vidfile,frame);
end
%close(vidfile)

triplot
%profile viewer