mu  = inline('sin(2*pi*x)');

lam       = 1.233701;
barmu     = 0;
bari2mu   = 1/2; %integral(mu^2);
m        = 2*pi;% maximum of \mu'
for j=1:100
    delta = j/100;
     for i= 1:100
         t=0:.01:delta;
         a(i)=i/100;
         %Y1=gampdf(t, 1.088870, 0.810570).*(exp(barmu*a(i))+exp(-barmu*a(i))).*exp(-.5*bari2mu.*t-m*a(i)*(t+delta));
         Y1=lam/a(i)^2*exp(-lam/a(i)^2*t).*(exp(barmu*a(i))+exp(-barmu*a(i))).*exp(-.5*bari2mu.*t-m*a(i)*(t+delta));

         Y1=sum(t.*Y1*.01);
         Y2=delta*exp(-lam/a(i)^2*delta)*exp(barmu^2/2*delta)*exp(-.5*bari2mu.*delta-2*m*a(i)*(delta));

         U(i)=Y1+Y2;
     end
    [opt,Index]  = max(U);
    Delt(j)      = delta;
    opt_level(j) = a(Index);
    opt_val(j)   = opt;
end
plot(Delt,opt_level)

