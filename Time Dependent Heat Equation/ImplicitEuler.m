function [T, label] = ImplicitEuler(Nx, Ny, dt, limit)
    T = ones(Nx+2,Ny+2); % Pre-allocating the geometry
    
    % Defining the co-efficients of terms
    
    c1 = dt*((Nx+1)^2);
    c2 = dt*((Ny+1)^2);
    c3 = (1+2*c1+2*c2);
    
    % Boundary conditions
    
    T(1,:)=0;
    T(end,:)=0;
    T(:,1)=0;
    T(:,end)=0;
    
    counter=0;
    t=0;
    
    % Numerical solution via Explicit Euler method
    
    while (t <= limit)
        residual = 1;
        T_previous = T;
        while (residual >= 1e-4)
            for i = 2:Nx+1
                for j = 2:Ny+1
                    T(i,j) = (T_previous(i,j) + c1*(T(i-1,j) + T(i+1,j)) + c2*(T(i,j-1) + T(i,j+1))) / c3; % Update rule of the PDE
                    if T(i,j) < 0
                        counter = counter+1;
                    end
                end
            end
            sum=0;
        
        % Calculation of Residual (for convergence and termination)
        
            for i = 2:Nx+1
                for j = 2:Ny+1
                    sum = sum + ((c3*T(i,j) - c1*(T(i-1,j) + T(i+1,j)) - c2*(T(i,j-1) + T(i,j+1))) - T_previous(i,j))^2;
                end
            end
            residual = sqrt(sum/(Nx*Ny));
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