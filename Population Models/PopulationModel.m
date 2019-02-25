%%ODE Solver for Population Model%%
%PopulationModel() - Takes inputs, calls appropriate functions, plots results
%exact() - Computes the exact value of the Population Function at each time step
%euler() - Calculates the approximate solution using Euler's Method
%heuns() - Calculates the approximate solution using Heun's Method
%rungekutta4() - Calculates the approximate solution using the 4th Order Runge Kutta Method
%error() - Computes the Global error as well as the error in each step

function PopulationModel() %Main function. Only for function calls and plots
    decision=1;
    while(decision==1)%This keeps the whole program running without an exit (unless selecting '0' at the promp).
        clc;
        clear all;
        disp('Enter 1 for Eulers Method, 2 for Heuns Method, 3 for Runge-Kutta Method ');
        choice=input('Enter Choice (1,2,3): ');
        step_size=[1 0.5 0.25 0.125];
        %Input variables
        y_initial=1;
        time=5;
          
        %y_initial=input('Enter initial population: ');
        %step=input('Enter Step-size: ');
        %time=input('Enter End Step (No. of years): ');
        switch choice
            case 1
                for i=1:4
                    steps=1:step_size(i):time;
                    y_exact=zeros(1,length(steps));
                    y_approximate=zeros(1,length(steps)); 
                    fprintf('\nEulers Method with Initial Population 1 and Time 5 years: ');
                    label='Eulers Method';
                    [t,y_exact]=exact(y_initial,step_size(i),time);%for exact values from the analytical solution
                    y_approximate=euler(y_initial,step_size(i),time);%for values obtained from Euler's Method
                    [error_method,step_error]=error(y_approximate,y_exact,step_size(i),time);%for calculation of error
                    fprintf('\nThe Cumulative error in Eulers Method for step size %f and time %f years is %f \n',step_size(i),time,error_method);
                    
                    %Graphs
                    
                    figure;
                    subplot(1,2,1);
                    plot(t,y_approximate,'-xr',t,y_exact,'-ob');
                    legend(label,'Exact Result');
                    xlim([0 1.2*time]);
                    ylim([0 12]);
                    hold on;
                    [T,P]=meshgrid(0:0.5:time,0:0.5*3:time*3);
                    dp=(1-P./10).*P;
                    dt=ones(length(dp));
                    L=sqrt(1+dp.^2);
                    quiver(T,P,dt./L,dp./L);
                    axis tight;
                    hold off;
                    subplot(1,2,2);
                    plot(t,step_error);
                    xlabel('Iteration No.');
                    ylabel('Step Error');                    
                end
            case 2
                for i=1:4
                    steps=1:step_size(i):time;
                    y_exact=zeros(1,length(steps));
                    y_approximate=zeros(1,length(steps));
                    fprintf('\nHeuns Method with Initial Population 1 and Time 5 years: ');
                    label='Heuns Method';
                    [t,y_exact]=exact(y_initial,step_size(i),time);%for exact values from the analytical solution
                    y_approximate=heuns(y_initial,step_size(i),time);%for values obtained from Heun's Method                
                    [error_method,step_error]=error(y_approximate,y_exact,step_size(i),time);%for calculation of error 
                    fprintf('\nThe Cumulative error in Heuns Method for step size %f and time %f years is %f \n',step_size(i),time,error_method);
                    
                    %Graphs
                    
                    figure;
                    subplot(1,2,1);
                    plot(t,y_approximate,'-xr',t,y_exact,'-ob');
                    legend(label,'Exact Result');
                    xlim([0 1.2*time]);
                    ylim([0 12]);
                    hold on;
                    [T,P]=meshgrid(0:0.5:time,0:0.5*3:time*3);
                    dp=(1-P./10).*P;
                    dt=ones(length(dp));
                    L=sqrt(1+dp.^2);
                    quiver(T,P,dt./L,dp./L);
                    axis tight;
                    hold off;
                    subplot(1,2,2);
                    plot(t,step_error);
                    xlabel('Iteration No.');
                    ylabel('Step Error');                    
                end
            case 3
                for i=1:4                  
                    steps=1:step_size(i):time;
                    y_exact=zeros(1,length(steps));
                    y_approximate=zeros(1,length(steps));
                    fprintf('\nRunge Kutta (4th Order) Method with Initial Population 1 and Time 5 years: ');
                    label='Runge Kutta (4th Order)';
                    [t,y_exact]=exact(y_initial,step_size(i),time);%for exact values from the analytical solution
                    y_approximate=rungekutta4(y_initial,step_size(i),time);%for values obtained from Runge Kutta's Method
                    [error_method,step_error]=error(y_approximate,y_exact,step_size(i),time);%for calculation of error
                    fprintf('\nThe Cumulative error in Runge Kutta (4th Order) Method for step size %f and time %f years is %f \n',step_size(i),time,error_method);
                    
                    %Graphs
                    
                    figure;
                    subplot(1,2,1);
                    plot(t,y_approximate,'-xr',t,y_exact,'-ob');
                    legend(label,'Exact Result');
                    xlim([0 1.2*time]);
                    ylim([0 12]);
                    hold on;
                    [T,P]=meshgrid(0:0.5:time,0:0.5*3:time*3);
                    dp=(1-P./10).*P;
                    dt=ones(length(dp));
                    L=sqrt(1+dp.^2);
                    quiver(T,P,dt./L,dp./L);
                    axis tight;
                    hold off;
                    subplot(1,2,2);
                    plot(t,step_error);
                    xlabel('Iteration No.');
                    ylabel('Step Error');                    
                end
            otherwise
                disp('Invalid Choice, Run again!!')
        end
            
        fprintf('\nWould you like to have another round? Press 1 for Yes and 0 for No.\n');
        decision=input('Enter Choice: ');
    end
end    

function [t,y_actual]=exact(y_initial,step,time)%for the calculation of exact solution at each point
    iter=(time/step);%number of iterations for the function
    t=0:step:time;
    y_actual=zeros(1,length(t));%pre-allocation of the vector
    y_actual(1)=y_initial;
    for i=1:iter
        y_actual(i+1)=10/(1+9*exp(-t(i+1)));
    end
end
    
function y_approx=euler(y_initial,step,time)%approximate solution by Euler's Method
    iter=(time/step);%number of iterations for the function
    t=0:step:time;
    y_approx=zeros(1,length(t));%pre-allocation of the vector
    y_approx(1)=y_initial;
    slope=@(y,t) (y-y^2/10);
    dydt=zeros(1,length(t));
    for i=1:iter
        %y_approx(i+1)=y_approx(i)+(step*((1-y_approx(i)/10)*y_approx(i)));%Euler's Method
        dydt(i)=slope(y_approx(i));
        y_approx(i+1)=y_approx(i)+(step*dydt(i));
    end
end

function y_approx=heuns(y_initial,step,time)%approximate solution by Heun's Method
    iter=(time/step);
    t=0:step:time;
    y_approx=zeros(1,length(t));
    y_approx(1)=y_initial;
    y_approx_predicted=zeros(1,length(t));%pre-allocation of the vector
    slope=@(y,t) (y-y^2/10);
    dydti=zeros(1,length(t));%pre-allocation of the vector
    dydtf=zeros(1,length(t));%pre-allocation of the vector
    for i=1:iter
        dydti(i)=slope(y_approx(i));%slope at the beginning of the interval
        y_approx_predicted(i+1)=y_approx(i)+(step*dydti(i));%Predictor Step (intermediate guess solution using Euler's Method)
        dydtf(i)=slope(y_approx_predicted(i+1));%slope at the end of the interval using the intermediate guess solution
        y_approx(i+1)=y_approx(i)+(0.5*step*(dydti(i)+dydtf(i)));%Corrector step (Heun's Method)
    end
end

function y_approx=rungekutta4(y_initial,step,time)%Approximate solution by Runge Kutta (4th order) Method
    iter=(time/step);
    t=0:step:time;
    y_approx=zeros(1,length(t));%pre-allocation of the vector
    y_approx(1)=y_initial;
    slope=@(y,t) (y-y^2/10);
    for i=1:iter
        k1=slope(y_approx(i));%slope at the beginning of the interval
        k2=slope((y_approx(i)+0.5*k1*step));%slope at the mid-point of the interval
        k3=slope((y_approx(i)+0.5*k2*step));%slope at the mid-point of the interval
        k4=slope((y_approx(i)+k3*step));%slope at the end of the interval
        phi=((k1+2*k2+2*k3+k4)/6);%weighted average of the slopes (increment function)
        y_approx(i+1)=y_approx(i)+phi*step;
    end
end
        
    
function [cumulative_error,step_error]=error(y_approx,y_actual,step,time)%calculation of error parameters
    step_error=(y_actual-y_approx);%error at each step
    squared_step_error=(step_error).^2;
    cumulative_error=sqrt((step/time)*sum(squared_step_error));%cumulative approximation error
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    