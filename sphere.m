%sphere.m   main program
nref=3;
nvmax=12*4^nref;
ntmax=20*4^nref;
x=zeros(nvmax,3);
v=zeros(ntmax,3);
dodec
for nr=1:nref
  refine
  %triplot
end

xnew= zeros(nv,3);
vnew=zeros(nt,3);
for i_la = 1:nv
    xnew(i_la,:) = x(i_la,:);
end
for i_la = 1:nt
     vnew(i_la,:) = v(i_la,:);
end

x=xnew;
v=vnew;
