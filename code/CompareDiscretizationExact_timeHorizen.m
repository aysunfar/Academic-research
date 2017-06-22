%Compute the estimator discretization for different time horizen (table 3) 
format long 
rand  ( 'twister' , sum(100*clock));
randn ( 'state'   , sum(100*clock));

%-------------------------------'discretizing'-------------------------
commenet ='obtimal budeget 1/5*se.^(-1/4) ';
Time_horizen        = [5,10,15,20,25,30,35,40];
N        = 1000;
delta    = 1/5*N.^(-1/4);

samples     = N;
M_discritiz = zeros(samples,7);
tim         = zeros(samples,7);


for j=1:length(Time_horizen)
    disp(['                  ', 'Time_horizen=',num2str(Time_horizen(j))])
    disp('     ')
    samples  =   N;
    T = Time_horizen(j)
    tic
    for i=1:samples

         est=Euler_estimator(delta(j),T);
         M_discritiz(i,j)=(est.max);
         tim(i,j)=est.time;
         if (~mod(i,round(samples/5)))
                disp(['Completed: ',num2str(20*round(i/round(samples/5))),'%']);
         end
    end

    time_D(j) = toc;
    mean_D(j) = mean(M_discritiz(1:samples,j));
    var_D(j)  = var(M_discritiz(1:samples,j));
end

 
save('longrun_discrite3');
 



