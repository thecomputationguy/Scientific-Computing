function HeatEquation()
    clc;
    clear all;
    
    % Initializing values
    
    N=[3 7 15 31];
    dt=[1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096];
    limit=[1/8 2/8 3/8 4/8];
    
    label=string(zeros(4,7,4));
        
    % Explicit Method
    % Iteration over the pre-defined values
    for i=1:4
        fig=figure('Position', get(0, 'Screensize')); %for full-screen figure
        n=1;
        for j=1:4 % For the Number of points
            x=0:1/(N(j)+1):1;
            y=0:1/(N(j)+1):1;
            [X,Y]=meshgrid(x,y);
            for k=1:7 % For the different time-steps
                [T,l]=ExplicitEuler(N(j),N(j),dt(k),limit(i));
                label(j,k,i)=l;
                subplot(4,7,n);
                surf(X,Y,T);
                set(gca,'fontsize',8);
                title(['dt=',num2str(dt(k)),' & N=',num2str(N(j))]);
                xlabel('X','fontsize',8);
                ylabel('Y','fontsize',8);
                n=n+1;               
            end
        end
        saveas(fig,sprintf('Explicit Time t=%d.png',limit(i)));
    end
    
    % Tables
    label=label(:,:,4);
    fprintf("Table for solution via Explicit Method \n");
    tab=["dt=1/64" "dt=1/128" "dt=1/256" "dt=1/512" "dt=1/1024" "dt=1/2048" "dt=1/4096";label];
    tab=[["Nx=Ny";"3";"7";"15";"31"] tab]
    
    % Implicit Method
    % Iteration over the pre-defined values and fixed time-step of 1/64
    
    label=string(zeros(4,4));
    
    for i=1:4
        fig2=figure('Position', get(0, 'Screensize')); % For full-screen figure
        n=1;
        for j=1:4
            x=0:1/(N(j)+1):1;
            y=0:1/(N(j)+1):1;
            [X,Y]=meshgrid(x,y);
            [T,l]=ImplicitEuler(N(j),N(j),1/64,limit(i));
            label(j,i)=l;
            subplot(2,2,n);
            surf(X,Y,T);
            set(gca,'fontsize',10);
            title(['dt=',num2str(dt(k)),' & N=',num2str(N(j))]);
            xlabel('X','fontsize',10);
            ylabel('Y','fontsize',10);
            n=n+1;
        end
        saveas(fig2,sprintf('Implicit Time t=%d.png',limit(i)));
    end
    
    % Tables
    fprintf("Table for solution via Implicit Method \n");
    tab=["t=1/8" "t=2/8" "t=3/8" "t=4/8";label];
    tab=[["Nx=Ny";"3";"7";"15";"31"] tab]
end
    
    
    