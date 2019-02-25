% Function for creating a full matrix of the coefficients

function A=full_matrix_formulation(Nx,Ny)
    hx=1/(Nx+1);
    hy=1/(Ny+1);
    hx=1/(Nx+1);
    hy=1/(Ny+1);
    A=zeros(Nx*Ny);
    for i=1:Ny*Nx %diagonal elements
        A(i,i)=-2*(1/hx^2 + 1/hy^2);
    end
    for i=1:Nx-1
        for j=1:Ny
            A(i+(j-1)*Nx, i+(j-1)*Nx+1)=1/hx^2; %upper main diagonal elements
            A(i+(j-1)*Nx+1, i+(j-1)*Nx)=1/hy^2; %lower main diagonal elements
        end
    end
    for i=1:Nx
        for j=1:Ny-1
            A(i+(j-1)*Nx, i+j*Nx)=1/hy^2; %upper off-diagonal elements
            A(i+j*Nx, i+(j-1)*Nx)=1/hy^2; %lower off-diagonal elements
        end
    end
    
end