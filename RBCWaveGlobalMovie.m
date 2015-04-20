
    % Shape parameter
    e = 4;
    
    %nrm2 = @(TI, TJ, PI, PJ) abs(2*(1 - cos(PI).*cos(PJ).*cos(TI - TJ) - sin(PI).*sin(PJ)));

    % Set time-step
    dt =  1.e-2;
    
    % Set number of time-steps.
    T = 5000;
    
%     % Nodes from sphere_pts.
%     [t, p] = sphere_pts(.1);
%     
%     
%     
%     % Set node spacing (interpolation)
%     N = length(t);
%     %n = sqrt(N);
    
    x=load('md044.02025'); N = length(x(:,1));

x = x(:,1:3);

[t, p] = cart2sph(x(:,1),x(:,2),x(:,3));
    
    %Constants for Red Blood Cell
    
    r0=3.39; c0=0.81/r0; c2=7.83/r0; c4=-4.39/r0; a=3.91/r0;
    
     % Parameterization of Red Blood Cell.
    x = a*cos(t).*cos(p);
    y = a*sin(t).*cos(p);
    z = 0.5*sin(p).*(c0+c2*(cos(p)).^2+c4*(cos(p)).^4);
   
   rbc = @(la,th) [a*cos(la).*cos(th) a*sin(la).*cos(th) ...
   0.5*sin(th).*(c0+c2*(cos(th)).^2+c4*(cos(th)).^4)];
        
 sz=[101 201]; M = prod(sz);  % surface grid parameters
[ll,tt]=meshgrid(linspace(-pi,pi,sz(2)),linspace(-pi/2,pi/2,sz(1)));
xx=rbc(ll(:),tt(:));
% Interpolate to the grid
re2=(repmat(xx(:,1),[1 N])-repmat(x.',[M 1])).^2;
re2=re2+(repmat(xx(:,2),[1 N])-repmat(y.',[M 1])).^2;
re2=re2+(repmat(xx(:,3),[1 N])-repmat(z.',[M 1])).^2;

% % Computing normals symbolically
% syms t p; % We will find the surface normal vectors symbolically
% nr = inline(cross(diff(rbc(t,p),t),diff(rbc(t,p),p)));
% [t, p] = sphere_pts(.1);
% nr= nr(t,p); 
% nr=nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);

 syms th ph; % We will find the surface normal vectors symbolically
%nr = inline(cross(diff(rbc(th,ph),th),diff(rbc(th,ph),ph)));
nr = @(th,ph) [-(391.*cos(ph).*cos(th).*((sin(ph).*((522.*cos(ph).*sin(ph))./113 - (1756.*cos(ph).^3.*sin(ph))./339))./2 - (cos(ph).*((261.*cos(ph).^2)./113 - (439.*cos(ph).^4)./339 + 27./113))./2))./339, -(391.*cos(ph).*sin(th).*((sin(ph).*((522.*cos(ph).*sin(ph))./113 - (1756.*cos(ph).^3.*sin(ph))./339))./2 - (cos(ph).*((261.*cos(ph).^2)./113 - (439.*cos(ph).^4)./339 + 27./113))./2))./339, (152881.*cos(ph).*cos(th).^2.*sin(ph))./114921 + (152881.*cos(ph).*sin(ph).*sin(th).^2)./114921];
nr= nr(t,p); 
nr=nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);
% Change to surf size to Plot the results
yy=reshape(xx(:,2),sz); zz=reshape(xx(:,3),sz); xx=reshape(xx(:,1),sz);

   
    
    %Set initial conditions
    
    u1 = .5*exp(-10*((x+1).^2+y.^2+z.^2))+0*exp(-10*((x).^2+y.^2+(z-.2).^2));
    u2 = zeros(N,1);
    u = [u1;u2];
    
    % Implement Grady's Method
    
    xmat = repmat(x, [1,N]);
    xmat = xmat-xmat.';
    nx = repmat(nr(:,1),[1,N]);
    
    ymat = repmat(y, [1,N]);
    ymat = ymat-ymat.';
    ny = repmat(nr(:,2),[1,N]);
    
    zmat = repmat(z, [1,N]);
    zmat = zmat-zmat.';
    nz = repmat(nr(:,3),[1,N]);
    
    phi=@(r2) 1./sqrt(1+e^2*r2);         % IMQ
    dphi=@(r2) -e^2./sqrt(1+e^2*r2).^3; % Derivative of IMQ over r
    
%     phi = @(r2) exp(-e^2*r2);            %Gaussian
%     dphi = @(r2) -e^2*exp(-e*r2);      % Derivative of Gaussian
    
    r2 = xmat.^2 + ymat.^2 +zmat.^2;
    A2 = dphi(r2);
    
    Dx = ((1-nx.^2).*xmat - nx.*ny.*ymat - nx.*nz.*zmat).*A2;
    Dy = (-nx.*ny.*xmat + (1-ny.^2).*ymat - ny.*nz.*zmat).*A2;
    Dz = (-nx.*nz.*xmat - ny.*nz.*ymat + (1-nz.^2).*zmat).*A2;
   
    
    A2 = chol(phi(r2)); Dx = (Dx/A2)/A2.'; Dy = (Dy/A2)/A2.'; Dz = (Dz/A2)/A2.';
    L2 = Dx*Dx + Dy*Dy  + Dz*Dz; %This is the surface Laplacian
    
    AD = phi(r2);
    
%      for jj =0:5
         gam1 = 1/(N^(1.3));
    gam2 = gam1;
    
    
    
    D1 = [zeros(N,N),eye(N,N);L2, -(gam1*AD^-1)];
%      D2 =  [zeros(N,N),eye(N,N);L2, zeros(N,N)];
%     D3 = [zeros(N,N)-(gam1*AD^-1),eye(N,N);L2, -(gam2*AD^-1)];
%     
%     
% %     e1 = eig(L2);
%     e2 = eig(D1);
%     maxe2 = [maxe2 max(real(e2))];
%      e3 = eig(D2);
% %     e4 = eig(D3);
% %     figure(1)
% %     plot(real(e1),imag(e1),'b*')
% %     title('Evals of Laplace-Beltrami')
%     figure
%     plot(real(e2),imag(e2),'r*')
%     title('Evals of D with gam2/gam1')
%     end
%     figure
%     plot(real(e3),imag(e3),'g*')
%     title('Evals of D without added hyperviscosity')
% %     figure(4)
% %     plot(real(e4),imag(e4),'m*')
% %     title('Evals of D with gam1/gam2')
% %     hold on
% %     
%     
    
    
    
%     % fine grid for plots
%     xplot = a*cos(tmatplot).*cos(pmatplot);
%     yplot =  a*sin(tmatplot).*cos(pmatplot);
%     zplot = 0.5*sin(pmatplot).*(c0+c2*(cos(pmatplot)).^2+c4*(cos(pmatplot)).^4);
%     
%     % Debugging plots below (why initial condition is so badly
%     % approxiamted?)
%     
%     % plot initial condition

%     figure(100)
%     u1plot = .5*exp(-10*((xx).^2+yy.^2+(zz-.2).^2));
%     surf(xx,yy,zz,u1plot);
%     shading interp; daspect([1 1 1]); axis tight;
%     colorbar
%     drawnow
% 
%     figure(200)
%     lambda = A2\(A2.'\u1);
%     uplot = reshape(phi(re2)*lambda, sz);
%     surf(xx, yy, zz, abs(uplot-u1plot));
%     shading interp; daspect([1 1 1]); axis tight;
% %     hold on
% %     plot3(x,y,z,'.k','markersize',20)
%     colorbar
%     drawnow
    
    
% %     

% Movie Test.


 
%% Set up the movie.
% writerObj = VideoWriter('rbcwave2.avi'); % Name it.
% writerObj.FrameRate = 30; % How many frames per second.
% open(writerObj); 
 

clear mov;
    Aplot= phi(re2);
        for i = 1:T
               
        % Plot solution
        if mod(i, 5) == 0
   
 
            
            % Compute solution on plotting mesh
            lambda = A2\(A2.'\u1);
            uplot = reshape(Aplot*lambda, sz);
            
            % Plot interpolated solution
            surf(xx, yy, zz, uplot);
            shading interp; daspect([1 1 1]); axis tight;
            colorbar
            drawnow
            mov(i/5)=getframe(gcf);
%             frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
%             
            
        end

        % Fourth-Order Runge-Kutta
        d1 = dt*D1*u;         
        d2 = dt*D1*(u + (1/2)*d1);    
        d3 = dt*D1*(u + (1/2)*d2);     
        d4 = dt*D1*(u + d3);            
        u = u + (1/6)*(d1 + 2*d2 + 2*d3 + d4);
        u1 = u(1:N);
        i
        
        end
%         hold off
% close(writerObj); % Saves the movie.
%      
        
