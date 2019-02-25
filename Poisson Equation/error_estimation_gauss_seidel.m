% Function for the evaluation of the error in the Gauss-Seidel method

function error=error_estimation_gauss_seidel(T,Z,Nx,Ny)
    sum=0;
    for i=2:Nx+1
        for j=2:Ny+1
            squared_error=(T(i,j)-Z(i,j))^2;
            sum=sum+squared_error;
        end
    end
    error=sqrt(sum/(Nx*Ny));
end