function [T, label] = ExplicitEuler(Nx, Ny, dt, limit)
    T = ones(Nx+2,Ny+2); % Pre-allocating the geometry
    
    % Defining the co-efficients of terms
    
    c1 = dt*((Nx+1)^2);
    c2 = dt*((Ny+1)^2);
    
    % Boundary conditions
    
    T(1,:)=0;
    T(end,:)=0;
    T(:,1)=0;
    T(:,end)=0;
     
    t=0;
    counter=0;
    
    % Numerical solution via Explicit Euler method
    
    while (t <= limit)
        T_previous = T;
        for i = 2:Nx+1
            for j = 2:Ny+1
                T(i,j) = T_previous(i,j) + c1*(T_previous(i-1,j) - 2*T_previous(i,j) + T_previous(i+1,j)) + c2*(T_previous(i,j-1) - 2*T_previous(i,j) + T_previous(i,j+1)); % update rule of PDE
                if T(i,j)<0
                    counter = counter+1;
                end
            end
        end
        t = t+dt;
    end
    
    % Criteria for stability
    
    if counter >= 1
        label = "Unstable";
    else
        label = "Stable";
    end
 
end