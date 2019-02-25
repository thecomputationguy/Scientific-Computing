% Main function implementing the solution of 2-D Poisson's Equation 


function StationaryPoisson()
    clc;
    clear all;
    N=[7 15 31 63 127];
    fprintf("Enter 1 for Solution using full system matrix\n");
    fprintf("Enter 2 for Solution using sparse matrix\n");
    fprintf("Enter 3 for Solution using Gauss Seidel Iteration\n");
    choice=input('Enter Choice: ');
    switch(choice)
        case 1
            runtime_in_s=zeros(5,1);
            storage_in_bytes=zeros(5,1);
            for i=1:5
                Nx=N(i);
                Ny=N(i);
                A=full_matrix_formulation(Nx,Ny);
                size_A=length(A);
                x=0:1/(Nx+1):1;
                y=0:1/(Ny+1):1;
                [X,Y]=meshgrid(x,y);
                b=(-2*(pi^2))*sin(pi*X).*sin(pi*Y);
                b=b(2:end-1,2:end-1);
                size_b=length(b);
                c=reshape(b,[],1);
                tic
                T=A\c; %solving directly
                runtime_in_s(i)=toc;
                storage_in_bytes(i)=(size_A^2+size_b^2)*8; %total storage in bytes assuming each element has a size of 8 bytes
                x=x(2:end-1);
                y=y(2:end-1);
                z=analytical_solution(Nx-2,Ny-2);
                [X,Y]=meshgrid(x,y);
                T=reshape(T,Nx,Ny);
                figure
                subplot(2,2,1);
                surf(X,Y,T);
                zlim([0 1]);
                title(["Direct (Full Matrix) solution with no of interior points = ",num2str(N(i))]);
                subplot(2,2,2);
                contour(T);
                zlim([0 1]);
                title(["Direct (Full Matrix) solution with no of interior points = ",num2str(N(i))]);
                subplot(2,2,3);
                surf(X,Y,z);
                zlim([0 1]);
                title(["Analytical Solution with no of interior points = ",num2str(N(i))]);
            end
            N=N';
            t=table(N,runtime_in_s,storage_in_bytes)
        case 2
            runtime_in_s=zeros(5,1);
            storage_in_bytes=zeros(5,1);
            for i=1:5
                Nx=N(i);
                Ny=N(i);
                A=sparse_matrix_formulation(Nx,Ny);
                size_A=nnz(A); %calculates the number of non-zero elements in the sparse matrix
                x=0:1/(Nx+1):1;
                y=0:1/(Ny+1):1;
                [X,Y]=meshgrid(x,y);
                b=(-2*(pi^2))*sin(pi*X).*sin(pi*Y);
                b=b(2:end-1,2:end-1);
                size_b=length(b);
                c=reshape(b,[],1);
                tic
                T=A\c; % solving directly
                runtime_in_s(i)=toc;
                storage_in_bytes(i)=(size_A+size_b^2)*8; %total storage in bytes assuming each element has a size of 8 bytes
                x=x(2:end-1);
                y=y(2:end-1);
                z=analytical_solution(Nx-2,Ny-2);
                [X,Y]=meshgrid(x,y);
                T=reshape(T,Nx,Ny);
                figure
                subplot(2,2,1);
                surf(X,Y,T);
                zlim([0 1]);
                title(["Direct (Full Matrix) solution with no of interior points = ",num2str(N(i))]);
                subplot(2,2,2);
                contour(T);
                zlim([0 1]);
                title(["Direct (Full Matrix) solution with no of interior points = ",num2str(N(i))]);
                subplot(2,2,3);
                surf(X,Y,z);
                zlim([0 1]);
                title(["Analytical Solution with no of interior points = ",num2str(N(i))]);
            end
            N=N';
            t=table(N,runtime_in_s,storage_in_bytes)
        case 3
            runtime_in_seconds=zeros(5,1);
            storage_in_bytes=zeros(5,1);
            error=zeros(5,1);
            error_red=zeros(5,1);
            for i=1:5
                Nx=N(i);
                Ny=N(i);
                hx=1/(Nx+1);
                hy=1/(Ny+1);
                x=0:hx:1;
                y=0:hy:1;
                [X,Y]=meshgrid(x,y);
                b=(-2*(pi^2))*sin(pi*X).*sin(pi*Y);
                z=analytical_solution(Nx,Ny);
                [Temperature,time]=gauss_seidel(b,Nx,Ny);
                runtime_in_seconds(i)=time;
                length_Temperature=length(Temperature);
                length_b=length(b);
                storage_in_bytes(i)=(length_Temperature^2+length_b^2)*8;
                error(i)=error_estimation_gauss_seidel(Temperature,z,Nx,Ny);
                if i~=1
                    error_red(i)=error(i-1)/error(i);
                end
                figure
                subplot(2,2,1);
                surf(X,Y,Temperature);
                zlim([0 1]);
                title(["Gauss-Seidel solution with no of interior points = ",num2str(N(i))]);
                subplot(2,2,2);
                contour(Temperature);
                zlim([0 1]);
                title(["Gauss-Seidel solution with no of interior points = ",num2str(N(i))]);
                subplot(2,2,3);
                surf(X,Y,z);
                zlim([0 1]);
                title(["Analytical Solution with no of interior points = ",num2str(N(i))]);
            end
            N=N';
            t=table(N,runtime_in_seconds,storage_in_bytes)
            e=table(N,error,error_red)
        otherwise
            fprintf("Invalid Choice!! Run Again");
    end
end
    
    
    