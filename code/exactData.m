%Comparing the exact estimator
format long 

sampels         =  200000;    
Max_Sample_path =  zeros(sampels,1);
T_x             =  zeros(sampels,1);
tic
for i=1:sampels
if (~mod(i,round(sampels/5)))
  disp( ['Completed: ',num2str(20*round(i/round(sampels/5))),'%'] );
  save( 'longrun_exact' )
end

[Max_Sample_path(i),T_x(i)]  =  exactSampling_Max;


end
time_I          =  toc;
running_lenngth =  mean(T_x);
mean_I          =  mean(Max_Sample_path);
var_I           =  var(Max_Sample_path);

save('longrun_exact_N09')


