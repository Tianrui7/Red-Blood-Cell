A= 5*10^(-21);
B= 5*10^(-21);
k= 1 ;
hold on
%fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'m' )

line(xlim,[0 0]);
ylim([-3*10^-21,7*10^-21 ]);
xlim([0,5]); 

B= 2*10^(-21);
k= 5 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'r' )
%{
B= 1*10^(-21);
k= 5 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'b' )
%}
B= 5*10^(-21);
k= 5 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'r' )

%{
B= 5*10^(-22);
k= 10 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'g' )
%}
hold off