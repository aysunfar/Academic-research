




%drwoing how the algorithm works
 function figure_generat
 
 
 delta=.001;
 cons=1.15;
 d=1/(2*pi);%max of |\mu|
 a=.4;
 bargamma=-.5;
 
 [sup_path,ltau,Wtau,path,time_seq]=Discritized_sample_path(delta,0,-.4,2)
 m=ceil(length(time_seq)/2);
 T1=time_seq(ceil(length(time_seq)/2));
 
 
 [sup_path,ltau,Wtau,path,time_seq]=Discritized_sample_path(delta,0,,1)
 plot(time_seq,path)
%     
%     
% current_t=0;
% current_Y=0;
% Max_path=0;
% stop=false;
% whole_path=[];
% while ~stop
%    % whole_path=[whole_path;current_t,current_Y,Max_path,a];
%   %  a=min(-Max_path+cons*d+current_Y,.4);
%     a=.4;
% 
%    [MP,ltau,Wtau,sup_path,time_seq]=Discritized_sample_path(.00001,current_t,a);
%    whole_path=[whole_path;time_seq',sup_path'+current_Y];
%    MP+current_Y
%     Max_path=max(MP+current_Y,Max_path);
%     if (current_Y+Wtau)> Max_path-cons *d       %It is not enough to jump
%         current_Y=current_Y+Wtau;
%         current_t=current_t+ltau;
%         
%     else 
%     whole_path=[whole_path;current_t,NaN];
%     current_Y=current_Y+Wtau;
%     current_t=current_t+ltau;
%     level=Max_path-current_Y-d;
%     [ infty,Idel_time] = firs_hiting_constant_drift(bargamma,level);
%     
%     if infty 
%             stop=true;
%         else
%             current_t=current_t+Idel_time;
%             delta_Gamma=level-Idel_time*bargamma+imu(current_t)-imu(current_t-Idel_time);
%             current_Y=current_Y+delta_Gamma;
%     end
% 
%     end
% end
% T_x=current_t;
% 
% 
% plot(whole_path(:,1),whole_path(:,2))
 end
function [sup_path,ltau,Wtau,path,time_seq]=Discritized_sample_path(delta,s,a,type)
%Generate a sample path of length [s,T+s] with step delta

T=10;
t=s+(0:delta:(T-s));
N=length(t)-1;
dW=randn(1,N)*sqrt(delta);
W=[0,cumsum(dW)];
I=imu(t)-imu(0);
y=I+W;
M=zeros(1,N);
E=-log(rand(1,N));
for i=1:N
   
   r=y(i+1)-y(i);
   M(i)=1/2 *r+1/2*sqrt(r^2+2*E(i)*delta);
   if (type==0)&&(abs(y(i+1))>a)
       hit=i+1;
       break;
       
       elseif (type==1)&&(y(i+1)>a)
           hit=i+1;
       break;
           
       elseif (type==2)&&(y(i+1)<a)
           hit=i+1;
       break;         
   end
end

[sup_path,ind]=max(M(1:(hit))+y(1:(hit)));
path=y(1:(hit));
Wtau=y(hit);
tim=t(ind);
ltau=t(hit);
time_seq=t(1:hit);
end
function [ infty,y] = firs_hiting_constant_drift(gamma,x)
%sample the first hiting time to the level x
% B_t+gamma t where gamma <0
%output 

mu=x/(abs(gamma)); lambda=x^2;

U=rand;
infty=0;
if U<exp(2*gamma*x)  %Generate a inverse gussian distribution
       infty=0;
       v = randn;   
       y = v*v;
       x = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda)) *sqrt(4*mu*lambda*y + mu*mu*y*y);
       test =rand; 
       if (test <= (mu)/(mu + x))
              y= x;
       else
              y=(mu*mu)/x;
       end
else
    y=0;
    infty=1;
end
end


function y=imu(x)

y=1/(2*pi)*sin(2*pi*x)-.5*x;
end