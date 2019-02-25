function [T,t]=tempExplicitEuler(T,Nx,Ny,dt)
    c1=dt*((Nx+1)^2);
    c2=dt*((Ny+1)^2);
    T(1,:)=0;
    T(end,:)=0;
    T(:,1)=0;
    T(:,end)=0;
    tic
    n=1;
    frame=0;
    x=0:1/(Nx+1):1;
    y=0:1/(Ny+1):1;
    [X,Y]=meshgrid(x,y);
    residual=1;
    while (residual>=0.0001)
        T_previous=T;
        for i=2:Nx+1
            for j=2:Ny+1
                T(i,j)=T_previous(i,j)+c1*(T_previous(i-1,j)-2*T_previous(i,j)+T_previous(i+1,j))+c2*(T_previous(i,j-1)-2*T_previous(i,j)+T_previous(i,j+1));
            end
        end
        residual=norm(T-T_previous,2);
        n=n+1;
        if(mod(n,16)==0)
            frame=frame+1;
            surf(X,Y,T);
            axis([0 1 0 1 0 1]);
            h=gca;
            get(h,'FontSize');
            set(h,'FontSize',12);
            colorbar('location','eastoutside','fontsize',12);
            caxis([0.0001 1]);
            xlabel('X','fontSize',12);
            ylabel('Y','fontSize',12);
            title(['Temperature Distribution (Tt=Txx+Tyy), Iteration No. ',num2str(n)],'fontsize',12);
            fh = figure(1);
            set(fh, 'color', 'white'); 
            F=getframe;
        end     
    end
    t=toc;
    movie(F);
end