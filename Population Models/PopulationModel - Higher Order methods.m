%%ODE Solver for Population Model%%
%PopulationModel() - Takes inputs, calls appropriate functions, plots results
%exact() - Computes the exact value of the Population Function at each time step
%implicitEuler() - Calculates the approximate solution using the implicit Euler's Method
%explicitEuler() - Calculates the approximate solution using the explicit Euler's Method
%explicitHeuns() - Calculates the approximate solution using the explicit Heun's Method
%implicitAdamsmoulton() - Calculates the approximate solution using the Implicit Adams-Moulton Method
%implicitAdamsmoultonl1() - Calculates the approximate solution using the
%Implicit Adams-Moulton Method with Linearization Scheme 1
%implicitAdamsmoultonl2() - Calculates the approximate solution using the
%Implicit Adams-Moulton Method with Linearization Scheme 2
%solveNewton() - Finds the root of an equation using the Newton-Rhapson's
%Method
%error() - Computes the Global error as well as the error in each step
%stability() - Evaluates the stability of a method

function PopulationModel() %Main function. Only for function calls and plots
    decision=1;
    while(decision==1)%This keeps the whole program running without an exit (unless selecting '0' at the promp).
        clc;
        clear all;
        disp('Enter 1 for Eulers Method, 2 for Heuns & Adams-Moulton Method, 3 for Adams-Moulton (Linearized) Method ');
        choice=input('Enter Choice: ');
        %Input variables
        y_initial=20;
        time=5;       
        step_size=[1/2 1/4 1/8 1/16 1/32];
        %Some pre-allocations
        error_method_explicit=zeros(1,5);
        error_method_implicit=zeros(1,5);
        error_red_implicit=zeros(1,5);
        error_red_explicit=zeros(1,5);
        error_method_implicit_l1=zeros(1,5);
        error_method_implicit_l2=zeros(1,5);
        error_red_implicit_l1=zeros(1,5);
        error_red_implicit_l2=zeros(1,5);
        stability_explicit=["a","b","c","d","e"];
        stability_implicit=["a","b","c","d","e"];
        stability_implicit_l1=["a","b","c","d","e"];
        stability_implicit_l2=["a","b","c","d","e"];
          
        switch choice
            
            case 1
                fprintf('\nEulers Method (Implicit & Explicit 1st Order) with Initial Population 20 and Time 5 years:\n');
                for i=1:5
                    steps=1:step_size(i):time;
                    y_exact=zeros(length(steps));
                    y_approximate_implicit=zeros(length(steps)); 
                    [t,y_exact]=exact(y_initial,step_size(i),time);%for exact values from the analytical solution
                    y_approximate_implicit=implicitEuler(y_initial,step_size(i),time);%for values obtained from Euler's Method
                    y_approximate_explicit=explicitEuler(y_initial,step_size(i),time);
                    error_method_implicit(i)=error(y_approximate_implicit,y_exact,step_size(i),time);%for calculation of error
                    error_method_explicit(i)=error(y_approximate_explicit,y_exact,step_size(i),time);
                    if i==1
                        error_red_implicit(i)=0;
                        error_red_explicit(i)=0;
                    else
                        error_red_implicit(i)=error_method_implicit(i)/error_method_implicit(i-1);
                        error_red_explicit(i)=error_method_explicit(i)/error_method_explicit(i-1);
                    end
                    stability_implicit(i)=stability(y_approximate_implicit,step_size(i),time);
                    stability_explicit(i)=stability(y_approximate_explicit,step_size(i),time);
                    
                    
                    %Graphs
                    
                    figure;
                    plot(t,y_approximate_implicit,'-xr',t,y_exact,'-ob',t,y_approximate_explicit,'-*g');
                    title(['Step Size ',num2str(step_size(i))]);
                    legend('Implicit Euler','Exact Result','Explicit Euler');
                    xlim([0 5]);
                    ylim([0 20]);
                end
                
                    %Tables
                    
                dt=step_size';
                error_euler_implicit=error_method_implicit';
                error_red_euler_implicit=error_red_implicit';
                stability_euler_implicit=stability_implicit';
                T1=table(dt,error_euler_implicit,error_red_euler_implicit,stability_euler_implicit)
                Error_Euler_Explicit=error_method_explicit';
                Error_Red_Euler_Explicit=error_red_explicit';
                Stability_Euler_Explicit=stability_explicit';
                T2=table(dt,Error_Euler_Explicit,Error_Red_Euler_Explicit,Stability_Euler_Explicit)
                                
                
            case 2
                fprintf('\nHeuns Method & Adam-Moulton Method with Initial Population 20 and Time 5 years:\n');
                for i=1:5
                    steps=1:step_size(i):time;
                    y_exact=zeros(length(steps));
                    y_approximate_implicit=zeros(length(steps)); 
                    [t,y_exact]=exact(y_initial,step_size(i),time);%for exact values from the analytical solution
                    y_approximate_implicit=implicitAdamsmoulton(y_initial,step_size(i),time);%for values obtained from Euler's Method
                    y_approximate_explicit=explicitHeuns(y_initial,step_size(i),time);
                    error_method_implicit(i)=error(y_approximate_implicit,y_exact,step_size(i),time);%for calculation of error
                    error_method_explicit(i)=error(y_approximate_explicit,y_exact,step_size(i),time);
                    if i==1
                        error_red_implicit(i)=0;
                        error_red_explicit(i)=0;
                    else
                        error_red_implicit(i)=error_method_implicit(i)/error_method_implicit(i-1);
                        error_red_explicit(i)=error_method_explicit(i)/error_method_explicit(i-1);
                    end
                    stability_implicit(i)=stability(y_approximate_implicit,step_size(i),time);
                    stability_explicit(i)=stability(y_approximate_explicit,step_size(i),time);
                    
                    %Graphs
                    
                    figure;
                    plot(t,y_approximate_implicit,'-xr',t,y_exact,'-ob',t,y_approximate_explicit,'-*g');
                    title(['Step Size ',num2str(step_size(i))]);
                    legend('Implicit Adams-Moulton','Exact Result','Explicit Heuns');
                    xlim([0 5]);
                    ylim([0 20]);
                end
                
                    %Tables
                    
                dt=step_size';
                error_adamsmoulton_implicit=error_method_implicit';
                error_red_adamsmoulton_implicit=error_red_implicit';
                stability_adamsmoulton_implicit=stability_implicit';
                T1=table(dt,error_adamsmoulton_implicit,error_red_adamsmoulton_implicit,stability_adamsmoulton_implicit)
                Error_Heuns_Explicit=error_method_explicit';
                Error_Red_Heuns_explicit=error_red_explicit';
                Stability_Heuns_Explicit=stability_explicit';
                T2=table(dt,Error_Heuns_Explicit,Error_Red_Heuns_explicit,Stability_Heuns_Explicit)
                
                
            case 3
                fprintf('\nAdam-Moulton Method (Linearized) with Initial Population 20 and Time 5 years:\n');
                for i=1:5
                    steps=1:step_size(i):time;
                    y_exact=zeros(length(steps));
                    [t,y_exact]=exact(y_initial,step_size(i),time);%for exact values from the analytical solution
                    y_approximate_implicit_l1=implicitAdamsmoultonl1(y_initial,step_size(i),time);%for values obtained from Euler's Method
                    y_approximate_implicit_l2=implicitAdamsmoultonl2(y_initial,step_size(i),time);
                    error_method_implicit_l1(i)=error(y_approximate_implicit_l1,y_exact,step_size(i),time);%for calculation of error
                    error_method_implicit_l2(i)=error(y_approximate_implicit_l2,y_exact,step_size(i),time);
                    if i==1
                        error_red_implicit_l1(i)=0;
                        error_red_implicit_l2(i)=0;
                    else
                        error_red_implicit_l1(i)=error_method_implicit_l1(i)/error_method_implicit_l1(i-1);
                        error_red_implicit_l2(i)=error_method_implicit_l2(i)/error_method_implicit_l2(i-1);
                    end
                    stability_implicit_l1(i)=stability(y_approximate_implicit_l1,step_size(i),time);
                    stability_implicit_l2(i)=stability(y_approximate_implicit_l2,step_size(i),time);
                    
                    %Graphs
                    
                    figure;
                    plot(t,y_approximate_implicit_l1,'-xr',t,y_exact,'-ob',t,y_approximate_implicit_l2,'-*g');
                    title(['Step Size ',num2str(step_size(i))]);
                    legend('Implicit Adams-Moulton (Linearization 1)','Exact Result','Implicit Adams-Moulton (Linearization 2)');
                    xlim([0 5]);
                    ylim([0 20]);
                end
                
                    %Tables
                dt=step_size';
                Error_AdamsMoulton_Implicit_L1=error_method_implicit_l1';
                Error_Red_AdamsMoulton_Implicit_L1=error_red_implicit_l1';
                Stability_AdamsMoulton_Implicit_L1=stability_implicit_l1';
                T1=table(dt,Error_AdamsMoulton_Implicit_L1,Error_Red_AdamsMoulton_Implicit_L1,Stability_AdamsMoulton_Implicit_L1)
                Error_AdamsMoulton_Implicit_L2=error_method_implicit_l2';
                Error_Red_AdamsMoulton_Implicit_L2=error_red_implicit_l2';
                Stability_AdamsMoulton_Implicit_L2=stability_implicit_l2';
                T2=table(dt,Error_AdamsMoulton_Implicit_L2,Error_Red_AdamsMoulton_Implicit_L2,Stability_AdamsMoulton_Implicit_L2)
                
            otherwise
                fprintf('\nInvalid Choice, Try again!! ');
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
        y_actual(i+1)=200/(20-10*exp(-7*t(i+1)));
    end
end

function y_euler_implicit=implicitEuler(y_initial,step,time)%approximate solution by Implicit Euler's Method
    t=0:step:time;
    iter=time/step;
    y_euler_implicit=zeros(1,length(t));
    y_euler_implicit(1)=y_initial;
    for i=1:iter
        d=((10/(7*step))*(1-7*step))^2+((40/(7*step))*y_euler_implicit(i));%Checks the nature of roots using the Discriminant Value
        if d<0
            fprintf('Solution does not exist for Eulers Implicit Method for Step Size %f\n',step)
        else
            y_euler_implicit(i+1)=solveNewton(y_euler_implicit(i),step);%calcultes the next value using the Newton-Rhapson's Method
        end
    end
end

function y_approx=explicitEuler(y_initial,step,time)%approximate solution by Explicit Euler's Method
    iter=(time/step);%number of iterations for the function
    t=0:step:time;
    y_approx=zeros(1,length(t));%pre-allocation of the vector
    y_approx(1)=y_initial;
    slope=@(y) (7*y-0.7*(y^2));
    dydt=zeros(1,length(t));
    for i=1:iter
        dydt(i)=slope(y_approx(i));
        y_approx(i+1)=y_approx(i)+(step*dydt(i));%Euler's Method of finding the next value
    end
end

function y_approx=explicitHeuns(y_initial,step,time)%approximate solution by explicit Heun's Method
    iter=(time/step);
    t=0:step:time;
    y_approx=zeros(1,length(t));
    y_approx(1)=y_initial;
    y_approx_predicted=zeros(1,length(t));%pre-allocation of the vector
    slope=@(y,t) 7*(y-y^2/10); 
    for i=1:iter
        dydti=slope(y_approx(i));%slope at the beginning of the interval
        y_approx_predicted(i+1)=y_approx(i)+(step*dydti);%Predictor Step (intermediate guess solution using Euler's Method)
        dydtf=slope(y_approx_predicted(i+1));%slope at the end of the interval using the intermediate guess solution
        y_approx(i+1)=y_approx(i)+(0.5*step*(dydti+dydtf));%Corrector step (Heun's Method)
    end
end

function y_adamsmoulton_implicit=implicitAdamsmoulton(y_initial,step,time)%approximate solution by Implicit Adams-Moulton Method
    t=0:step:time;
    iter=time/step;
    y_adamsmoulton_implicit=zeros(1,length(t));
    y_adamsmoulton_implicit(1)=y_initial;
    for i=1:iter
        d=((1-3.5*step)^2)-1.4*step*(0.35*step*y_adamsmoulton_implicit(i)^2-3.5*step*y_adamsmoulton_implicit(i)-y_adamsmoulton_implicit(i));
        if d<0
            fprintf('Solution does not exist for Adams-Moulton Implicit Method for Step Size %f\n',step);
            break;
        else
            y_adamsmoulton_implicit(i+1)=solveNewton(y_adamsmoulton_implicit(i),step);%calcultes the next value using the Newton-Rhapson's Method
        end
    end
end

function y_adamsmoulton_implicit_l1=implicitAdamsmoultonl1(y_initial,step,time)%approximate solution by Implicit Adams-Moulton Method (Linearization 1)
    t=0:step:time;
    iter=time/step;
    y_adamsmoulton_implicit_l1=zeros(1,length(t));
    y_adamsmoulton_implicit_l1(1)=y_initial;
    for i=1:iter
        y_adamsmoulton_implicit_l1(i+1)=y_adamsmoulton_implicit_l1(i)*(1+7*step-((7/20)*step*y_adamsmoulton_implicit_l1(i)))/(1+((7/20)*step*y_adamsmoulton_implicit_l1(i)));
    end
end

function y_adamsmoulton_implicit_l2=implicitAdamsmoultonl2(y_initial,step,time)%approximate solution by Implicit Adams-Moulton Method (Linearization 2)
    t=0:step:time;
    iter=time/step;
    y_adamsmoulton_implicit_l2=zeros(1,length(t));
    y_adamsmoulton_implicit_l2(1)=y_initial;
    for i=1:iter
        y_adamsmoulton_implicit_l2(i+1)=(y_adamsmoulton_implicit_l2(i)*(1+(7/2)*step-((7/20)*step*y_adamsmoulton_implicit_l2(i))))/(1+((7/20)*step*y_adamsmoulton_implicit_l2(i))-(7/2)*step);
    end
end

function y_approx=solveNewton(y_euler_implicit,step)
    guess_solution=zeros(1,2);
    guess_solution(1)=y_euler_implicit;
    h=step;
    y=@(y) y^2+((10/(7*h))*(1-7*h))*y-(10/(7*h))*y_euler_implicit;
    y_prime=@(y) 2*y+(10/(7*h));
    while (abs(guess_solution(1)-guess_solution(2))>=1e-4)%Tolerance criterion
        guess_solution(2)=guess_solution(1)-y(guess_solution(1))/y_prime(guess_solution(1));%Newton-Rhapson's Method for finding roots
        guess_solution(1)=guess_solution(2);
    end
    y_approx=guess_solution(1);
end

function cumulative_error=error(y_approx,y_actual,step,time)%calculation of error parameters
    step_error=(y_actual-y_approx);%error at each step
    squared_step_error=(step_error).^2;
    cumulative_error=sqrt((step/time)*sum(squared_step_error));%cumulative approximation error
end

function stb=stability(y_approximate,step,time)
    iter=time/step;
    counter=0;
    for i=1:iter
        if y_approximate(i)<=0%Negative population is the defining criteria for unstable soulution
            counter=counter+1;
        end
    end
    if counter>1
        stb="unstable";
    else
        stb="stable";
    end
end    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        