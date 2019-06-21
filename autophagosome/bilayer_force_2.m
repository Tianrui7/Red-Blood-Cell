A= 5;
B= 6;
k= 1 ;
hold on
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'m' )

line(xlim,[0 0]);
ylim([-6,5]);
xlim([0,1]); 

B= 6;
k= 20 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'r' )
%{
B= 1*10^(-21);
k= 5 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'b' )
%}
B= 5;
k= 50 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'b' )

%{
B= 5*10^(-22);
k= 10 ;
fplot(@(D) -A/(12*pi*(D^2)) + B*k*exp(-k*D),[0 5],'g' )
%}
hold off
legend({'A=5,B=6,k=1','A=5,B=5,k=50','A=5,B=6,k=20'})