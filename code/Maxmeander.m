function [M,step,ratio] = Maxmeander(t,X,Y )

% %Generate a sample of Brownian Menader such that 
% %over [0,t], B(s)>0
% %B(0)=x,B(t)=y
% The algorithm is Based on Luc Devroye (2009) 

x  =  X/sqrt(t);y=Y/sqrt(t);
n  =  -20:20;
 
if ( y>0.0 )
    r         = (y-x);
    Accepeted = 1;
    step      = 0;
    while Accepeted
        E      =  -log(rand);
        a      =  .5*(r+sqrt(r^2+2*E))+x;
        V      =   exp(-(x-y-2*n*a).^2/2).*(x-y-2.*n*a)+exp(-(x+y+2*n*a).^2/2).*(x+y+2*n*a);
        f      =  sum(n.*V)/(exp(-(x-y)^2/2)*(1-exp(-2*x*y)));
        g      =  exp(.5*(r^2-(2*(a-x)-r)^2)).*(2*a-(x+y));
        ratio  =  f/(2*g);
        U      =  rand;
        ccc    = max(1/abs(r),10);
        if f>ccc*g*U
           Accepeted = 0;
           M         = a;
           
        end

       step  =  step+1;
    end
elseif (y==0)
    
  if x>1.5
      zita   = 3.5789e-04;
      c      = 5*x/(10*x^2-8);
   
      
      Accepted = 1;
      step     = 0;
      while Accepted
       E    =  -log(rand);
       a    =   x+c*E;
       V    =   n.*(1-(2*n*a+x).^2).*exp(-(2*n*a+x).^2/2);
       f    =   2*sum(V)/(x*exp(-x^2/2));
       U    =   rand;
       g    =   10*x/(1-zita)*exp(-E);
       if f>g*U
           Accepted = 0;
           M        = a;
           
       end
       step   = step+1;
       ratio  = f/g;
       end

  end
  if x<1.5
      step  = 0;
      p     = 9.4505611743;tau=4*exp(-9);nu=4*tau; const=769.1762;
      q     =123*exp(3*x-4.5)/(2*(1-nu)*(4-2*x));
      Accepted=1;ratio=0;
      while Accepted
       step=step+1;
       U1=rand;U2=rand;
  
      if (U1<p/(p+q))
          
               E1=-log(rand);E2=-log(rand);Normalr=randn;
               a=pi/sqrt(Normalr^2+2*(E1+E2));
               %a=pi/sqrt(2*randg(2.5));
               V=n.*(1-(2*n*a+x).^2).*exp(-(2*n*a+x).^2/2);
               f=2*sum(V)/(x*exp(-x^2/2));
               g=const*exp(-pi^2/(2*a^2))/a^6;
               if (f>U2*g)&&(a<1.5)
                   Accepted=0;
                   M=a;
               end
    
                ratio=max(ratio,f/g);

      else
          step=step+1;

          E1  =  -log(rand);
          a   =   1.5+E1/(4-2*x);
          V   =   n.*(1-(2*n*a+x).^2).*exp(-(2*n*a+x).^2/2);
          f   =   2*sum(V)/(x*exp(-x^2/2));
          g   =   q*(4-2*x)*exp(-(4-2*x)*(a-1.5));
          if f>U2*g
                   Accepted = 0;
                   M        = a;
          end
                   ratio    = max(ratio,f/g);
          
      end
      end
      
end
end

M    =  M*sqrt(t);
end

