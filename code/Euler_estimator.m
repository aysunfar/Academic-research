function Est = Euler_estimator(delta,T)
     if nargin<2
         delta=.01;
         T=10;
     end

    [sup_path,W,tim]  = Discritized_sample_path(delta,T);
    Est.max           = (sup_path);
    Est.W             = W;
    Est.time          = tim;

end


function [sup_path,W_T,tim]=Discritized_sample_path(delta,T)
%Generate a sample path of length [0,T] with step delta
    t   =   0 : delta: T;
    N   =   length(t)-1;
    dW  =   randn(1,N)*sqrt(delta);
    W   =   [0,cumsum(dW)];
    I   =   imu(t)-imu(0);
    y   =   I+W;
    M   =   zeros(1,N);
    E   =   -log(rand(1,N));
    
    for i=1:N

        r        = y(i+1)-y(i);
        M(i)     = 1/2 *r+1/2*sqrt(r^2+2*E(i)*delta);
    end
    [sup_path,ind]  =  max(M+y(1:N));
    W_T             =  y(N);
    tim             =  t(ind);
end


function y=imu(x)

      y  =  1/(2*pi)*sin(2*pi*x)-.5*x;
end