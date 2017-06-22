%Created on 03 Nov 2013
%@author: Mohammad Mousavi
function  [Max_path,T_x]= exactSampling_Max
%This fucntion is the main function 
% generate an excat sample path of a time-dependent B.M 
% output:
%    Max_path: the global maximum of this sample path. 
%    T_x:      the time at which maximum occurs
% The drift function mu (x) and its integral and derivative is defined in 
% the end of script  



global a 
a      = .4;
m      = 2*pi;         % maximum of \mu'
d      = 1/(2*pi);     % max of |\mu|

bargamma  = -.5;
cons      = 1.15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
T_x        = 0; 
current_t  = 0;
current_Y  = 0;
Max_path   = 0;
stop       = false;
while ~stop
     if a>0
         [MP,ltau,Wtau] = max_first_hit( current_t,.4 );
     else
         MP   = 0;
         ltau = 0;
         Wtau = 0;
    end
    Max_path = max( MP+current_Y,Max_path );
    if (current_Y+Wtau)> Max_path-cons *d       %It is not enough to jump
       
        current_Y = current_Y + Wtau;
        current_t = current_t + ltau;
        T_x       = T_x       + ltau;
   
    else 
        T_x       = T_x       + ltau;
        current_Y = current_Y + Wtau;
        current_t = current_t + ltau;
        level     = Max_path  - current_Y - d;
      
        [ infty,Idel_time] = firs_hiting_constant_drift( bargamma,level );

        if infty 
                stop        = true;
        else
                current_t   = current_t+Idel_time;
                delta_Gamma = level-Idel_time*bargamma+imu(current_t)-imu(current_t-Idel_time);
                current_Y   = current_Y  +delta_Gamma;
        end

    end
end


end


%-------------------------------------------------------------------
function [ infty,y] = firs_hiting_constant_drift(gamma,x)
%sample the first hiting time to the level x
% B_t+gamma t where gamma <0
%output 

mu    =  x/(-gamma); lambda=x^2;
U     =  rand;
infty =  0;
if U<exp(2*gamma*x)  %Generate a inverse gussian distribution
       infty = 0;
       v     = randn;   
       y     = v*v;
       x     = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda)) *sqrt(4*mu*lambda*y + mu*mu*y*y);
       test  =  rand; 
       if (test <= mu  /(mu + x) )
              y = x;
       else
              y = ( mu * mu )/x;
       end
else
       y     =  0;
       infty =  1;
end
end

%-------------------------------------------------------------------

function  [MP,ltau,Wtau]=max_first_hit(s,delta)

%Generate a sample path (B_{lata}-B_s)+\int _s^{ltau}\mu(u) du  starting at s 
%Input: 
  %delta: the  maximum length of horzin. 
  
%output:
%       ltau=the first hiting time to the a in abouslout vaule
%       Mp= Maximum of the sample path
%       Wtau= the value at ltau (+/- a or )
global a
    
m    =  2*pi;% maximum of \mu'
tm   =  1.5;%max of \mu
rate =  2*m*a;
Li   =  a;


accepted=false;

while   ( ~accepted )
    [ tau , BB ]  = BMatpossiontime( Li,rate,delta );
    if tau( length(tau)-1 )  ==  delta
        b     = length(tau)-2;
        ltau  = delta;
        Wtau  = BB(b+1);
        tau   = tau(1:(1+b));
    else
         b    = length(tau)-1;
         ltau = tau(b+1);
         Wtau = BB(b+1);
    end
    
    if (b<=1)
        U     =  rand(2,1);
        V     =  zeros(2,1);
        V(1)  =  exp(mu(ltau+s)*Wtau-tm*a);
        V(2)  =  exp(-1/2*(i2mu(ltau+s)-i2mu(s))+a*m*(ltau-delta));
     
        if (U(1)<V(1)) &&(U(2)<V(2))
             accepted = true;
        end
        
    
    
    else
        U            =  rand(b+1,1);
        V            =  zeros(b+1,1);
        V( 1:(b-1) ) = (dmu(tau(2:b)+s).*BB(2:(b))+m*a)/(2*m*a);
        V(b)         = exp(mu(ltau+s)*Wtau-tm*a);
        V(b+1)       = exp(-1/2*(i2mu(ltau+s)-i2mu(s))+a*m*(ltau-delta));
    

        if ( U(1:(b-1) )> V( 1:(b-1) )  )
           if (U(b)<V(b)) &&(U(b+1)<V(b+1))
             accepted  =  true;
           end
        end
    end
end

MP   =   MaxBRpath  ( tau,BB,Li );




end




function MP=MaxBRpath(tau,BB,level)
%sample of the max a path of zero drift BM condited to be less than a level
%in absoulute value
%Input :
% tau,BB= skleton of the path (the value of path in specfic point) 
% level

MP  = 0;
K   = length(tau);
if BB(K) == level
    MP = level;
else
    for i=1:(K-1)
        dt = tau(i+1)-tau(i);
        M  = level+1;
        while (M>level)
            M = Maxmeander(dt,BB(i)+ level,BB(i+1)+level )- level;
        end
        MP  = max(MP,M);
    end
end
end

function [tau,BB]= BMatpossiontime(Li,rate,delta)
%This function genrate a sample path of BM until it hits level Li 
%and deterimine the value of the path at random times coming from 
%possion distribution
%Input:
%        Li:hitting level |B_t|<Li
%        rate: the rate of possion distribution
%Output:
%        tau:sequance of random times incdluing the hiting time
%        BB: value of BM at times taue
%The output is truncted untill some time delta and the hitting time


varsigma=  GenerateNextExitTime(Li);
tau     =  Generatearrivaltime(rate,varsigma,delta);
H1      =  GenerateWsigma(Li); %generating W_varsigma candidate

if (length(tau)>2)
                    tau1 = SigmaNormalize(tau,Li);
                    I    = false;
                    while ~I
                        %OLD VERSION
                        BB=3;
                        while max(BB)>2
                             BB=GenerateBBTriplets(varsigma/Li^2,tau1);
                         %   BB=GenerateBB(varsigma/Li^2,tau1,B); %generating meander from triplets                    
                        end                        
                        
                        U  =  rand(1,length(tau)-1);
                        I  =  TestI(U,tau1,varsigma/Li^2,BB);
                        if I
                            if (H1<0)
                                
                                    BB=Li.*(-1+BB);
                                
                            else
                                
                                    BB=Li*(1-BB);
                            
                            end
                        end
                    end                  
   else
                    BB=[0 H1];
end   
      
    
end 
    

       
       
        
        

function I=TestI(U,tau1,varsigma,BB) %meander testing function
dif      = 0;
Accepted = true;
Resolved = false;
q        = 0;
n        = 0;
x        = BB(length(BB)-1);
s        = tau1(length(tau1))-tau1(length(tau1)-1);
while ~Resolved
    qprev = q;
    n     = n+1;
    if (mod(n,2))
        m  = ceil(n/2);
        q  = q-(1/x)*(4*m-x)*exp(-4*m*(2*m-x)/s);
    else
        m  = n/2;
        q  = q+(1/x)*(4*m+x)*exp(-4*m*(2*m+x)/s);
    end
    if n>=log(4)*(s)/8*1/x+2
        if ((qprev<=q)&&(q+1<U(length(U))))
            Resolved = true;
            Accepted = false;
        else
            if ((qprev>=q)&&(q>U(length(U))-1))
                Resolved = true;
                Accepted = Accepted&&true;
            end
        end
    end
end
for i   =   1:length(U)-1
    t   =   tau1(length(tau1))-tau1(length(tau1)-i-1);
    s   =   tau1(length(tau1))-tau1(length(tau1)-i);
    y   =   BB(length(BB)-i-1);
    x   =   BB(length(BB)-i);
    p1  =   1/(1-exp(-2*x*y/(t-s)));
	p   =   p1;
    n   =   0;
    Resolved  =  false;
    while ~Resolved
        pprev =  p;
        n     =  n+1;
        if ( mod(n,2) )
            m  =  ceil(n/2);
            p  =  p-p1*(exp(-2*(2*m-x)*(2*m-y)/(t-s))+exp(-2*(2*(m-1)+x)*(2*(m-1)+y)/(t-s)));
        else
            m  =  n/2;
            p  =  p+p1*(exp(-2*m*(4*m+2*(x-y))/(t-s))+exp(-2*m*(4*m+2*(x-y))/(t-s)));
        end
        if n>=log(3)*(t-s)/8*max(1/x,1/y)+1
            if ( (pprev<=p)&&(p<U(i)) )
                Resolved = true;
                Accepted = false;
            else
                if ((pprev>=p)&&(p>U(i)))
                    Resolved  =  true;
                    Accepted  =  Accepted&&true;
                end
            end
        end
    end
    if ~Accepted
        break;
    end
end
I=Accepted;
end
function b=GenerateBB(varsigma,tau,Triplets) %construct meander from triplets
for i=1:length(tau)
	b(i) = sqrt(((varsigma-tau(i))/varsigma+Triplets(i,1))^2+(Triplets(i,2))^2+(Triplets(i,3))^2);
end

b2 = sqrt(((varsigma-tau)/varsigma+Triplets(:,1)').^2+(Triplets(:,2)').^2+(Triplets(:,3)').^2);
end

function t=SigmaNormalize(tau,L) %scaling meander to standard BM with exit boundary L
t=zeros(1,length(tau));
for i=1:length(tau)
    t(i)=tau(i)/L^2;
end
end

function BM=GenerateBBTriplets(varsigma,tau) %generate a triplet of Brownian Bridges used in meander construction
len=length(tau);
b=zeros(len,3);
z=randn(3,length(tau));
for i=2:(len-1)
    for j=1:3
        b(i,j)=(varsigma-tau(i))*b(i-1,j)/(varsigma-tau(i-1))+z(j,i)*sqrt((varsigma-tau(i))*(tau(i)-tau(i-1))/(varsigma-tau(i-1)));
        
    end
end

Triplets=b;


BM=sqrt(((varsigma-tau)/varsigma+Triplets(:,1)').^2+(Triplets(:,2)').^2+(Triplets(:,3)').^2);

end

function tau=Generatearrivaltime(rate,varsigma,delta) %Poisson arrivals generation for acceptance over [0,min(delta,VS)] the output includes VS

tau=zeros(1,10);
i=1;
Lt=min(delta,varsigma);
while (tau(i)<Lt)
    dt=-1/rate*log(rand);
    if (tau(i)+dt<Lt)
        tau(i+1)=tau(i)+dt;
    else
        if (delta<varsigma)
        tau(i+1) = delta;
        tau(i+2 )= varsigma;
        tau      = tau(1:(i+2));
        
        else 
            
        tau(i+1)  = varsigma;
        tau       = tau(1:(i+1));
        
        end
        
    end
    i=i+1;
end

end
function H=GenerateWsigma(Li) %exit value generation
U=rand;
if U<=0.5
    H = -Li;
else
    H = Li;
end
end
function varsigma=GenerateNextExitTime(Li) %generate next exit time for standard BM
V        = xing();
varsigma = Li^2*V;
end

function X = xing()
%this function copied from Burq and Jones paper
accepted = 0;
while ~accepted
	X = gamrnd(1.088870, 0.810570);
    
    %  c1=.2737498;a=1.088870;lam=1.233701;
   % Y=rand*1.243707*exp((a-1)*log(X)-lam*X+c1);
	Y = rand*1.243707*gampdf(X, 1.088870, 0.810570);
    
    
	sqrt2piX3 = sqrt(2*pi*X^3);
    N         = max([ceil(0.275*X), 3]);
	K         = (1+2*(-N:N));
	fN0       = sum((-1).^(-N:N).*K.*exp(-K.^2./(2*X)))/sqrt2piX3;
    N         = N + 1;
    fN1       = fN0 + (-1)^N*((1-2*N)*exp(-(1-2*N)^2/(2*X))+(1+2*N)*exp(-(1+2*N)^2/(2*X)))/sqrt2piX3;
   while (((Y<fN0)&&(Y > fN1))||  ((Y>fN0)&&(Y < fN1)))
       fN0  = fN1;
       N    = N + 1;
       fN1  = fN0 + (-1)^N*((1-2*N)*exp(-(1-2*N)^2/(2*X)) + (1+2*N)*exp(-(1+2*N)^2/(2*X)))/sqrt2piX3;
    end
    if Y <= fN1
        accepted = 1;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = mu(x)

  y   =   cos(2*pi*x)-.5;
end
function y = imu(x)

  y   =   1/(2*pi)*sin(2*pi*x)-.5*x;
end

function y = dmu(x)
      y  =   -2*pi*sin(2*pi*x);
end
function y = i2mu(x)
      y  =    0.75* x- sin(2*pi *x)/(2*pi)+1/(8*pi) *sin(4*pi *x);
end

