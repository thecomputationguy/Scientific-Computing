% Function implementing the Gauss-Seidel Iterative solver

function [T,time_gaussseidel]=gauss_seidel(b,Nx,Ny)
    T=zeros(Nx+2,Ny+2);
    hx=1/(Nx+1);
    hy=1/(Ny+1);
    c1=(Nx+1)^2;
    c2=(Ny+1)^2;
    c3=-2*(c1+c2);
    %Boundary Conditions
    %%%%%%%%%%%%%%%%%%%
    T(1,:)=0;
    T(end,:)=0;
    T(:,1)=0;
    T(:,end)=0;
    %%%%%%%%%%%%%%%%%%%
    residual=1;
    tic
    while (residual>1e-4)
        for i=2:Nx+1
            for j=2:Ny+1
                T(i,j)=(b(i,j)-c1*(T(i-1,j)+T(i+1,j))-c2*(T(i,j+1)+T(i,j-1)))/c3;
            end
        end
        sum=0;
        for i=2:Nx+1
            for j=2:Ny+1
                sum=sum+((b(i,j)-c1*(T(i-1,j)+T(i+1,j))-c2*(T(i,j+1)+T(i,j-1)))-c3*T(i,j))^2;
            end
        end
        residual=sqrt(sum/(Nx*Ny));
    end
        time_gaussseidel=toc;
end