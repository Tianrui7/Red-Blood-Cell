%triplot.m   plots vertices and triangles

%figure(2)

%plot3(x(1:nv,1),x(1:nv,2),x(1:nv,3),'ko')
hold on
%axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
%axis equal manual

%plot3(x(1:nv,1),x(1:nv,2),x(1:nv,3),'ko')
%ko means black circle
%Vertices are all plotted. We are left to plot lines that connect points.

xplot=zeros(4,3);

for t_triplot=1:nt
  v1=v(t_triplot,1);
  v2=v(t_triplot,2);
  v3=v(t_triplot,3);
  %when t=1, v1,v2,v3 are original denotations (given by dodec) of the three vertices of face 1
  xplot(1,:)=x(v1,:);
  xplot(2,:)=x(v2,:);
  xplot(3,:)=x(v3,:);
  xplot(4,:)=x(v1,:);
  plot3(xplot(:,1),xplot(:,2),xplot(:,3))
end

view(0,0)
%view sets the viewing angle for a three-dimensional plot
