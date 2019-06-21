
A= 5;
line(xlim,[0 0]);
hold on
ylim([-3*10,3*10 ]);
xlim([0,5]); 
B= 6;
k= 20;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'r' )
hold off
Energy 
Energy_res
D = 1.6;
-A/(12*pi*(D^2)) + B*k*exp(-k*D)

%114.5024

%w(D)*number of points is not a same unit with energy 