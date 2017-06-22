%Comparing the expectation estimator discretization and exact methods(table 2) 
format long 
rand  ( 'twister' , sum(100*clock));
randn ( 'state'   , sum(100*clock));


sampels         =  100000;    
Max_Sample_path =  zeros(sampels,1);
T_x             =  zeros(sampels,1);
tic
for i=1:sampels
    if (~mod(i,round(sampels/5)))
        disp( ['Completed: ',num2str(20*round(i/round(sampels/5))),'%'] );
        save( 'longrun_exact_N10_9' )
    end

[Max_Sample_path(i),T_x(i)]  =  exactSampling_Max;


end
time_I          =  toc;
running_lenngth =  mean(T_x);
mean_I          =  mean(Max_Sample_path);
var_I           =  var(Max_Sample_path);

save('longrun_exact_N09')
%-------------------------------'discretizing'-------------------------
commenet ='obtimal budeget 1/5*se.^(-1/4) ';
T        = 35;
N        = [10,20,40,50,100,150,200,250]*1000;
delta    = 1/5*N.^(-1/4);

samples     = 200000;
M_discritiz = zeros(samples,7);
tim         = zeros(samples,7);


for j=1:7
    disp(['                  ','delta=',num2str(delta(j))])
    disp('     ')
    samples  =   N(j);
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
 

for i=1:length(N)
    pd        =  fitdist(M_discritiz(1:N(i),i),'Kernel','Kernel','epanechnikov');
    x_values  =  0:.01:5;
    pdf1(i,:) =  pdf(pd,x_values);
end

pd                   =  fitdist(Max_Sample_path,'Kernel','Kernel','epanechnikov');
x_values             =  0:.01:5;
pdf1 (length(N)+1,:) =  pdf(pd,x_values);
plot(x_values,pdf1,'--')
%%%-------------------------


hold on;plot(x_values(1:500),pdf1(3,1:500),'black:','LineWidth',1.5)
hold on;plot(x_values(1:500),pdf1(5,1:500),'o','MarkerSize',1.5)
hold on;plot(x_values(1:500),pdf1(8,1:500),'r-.','LineWidth',1.1)
hold on;plot(x_values(1:500),pdf1(9,1:500),'g-','LineWidth',1.1)

legend('\delta=2^{-2}','\delta=2^{-6}','\delta=2^{-12}','Exact')

[h,p]=kstest2(M_discritiz(1:samples,8),Max_Sample_path)
end

