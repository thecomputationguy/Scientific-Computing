% Function for the evaluation of the Analytical Solution

function exact=analytical_solution(Nx,Ny)
    hx=1/(Nx+1);
    hy=1/(Ny+1);
    x=0:hx:1;
    y=0:hy:1;
    [X,Y]=meshgrid(x,y);
    exact=sin(pi*X).*sin(pi*Y);
end