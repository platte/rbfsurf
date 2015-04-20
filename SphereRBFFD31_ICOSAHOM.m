  % Shape parameter
    %e = 3;
    KT = 10^12;
    
    %Size of local RBF-FD stencil
    size = 31;
    
       % Set time-step
    dt =  5e-3;
    
    % Set number of time-steps.
    T = 200;
    
    time = (0:.005:1);

    x=load('md077.06084'); N = length(x(:,1));
    x = x(:,1:3);
    nr= x./repmat(sqrt(sum(x.^2,2)),[1 3]);
    
    
    %Obtain the size nearest neighbors to each node
    
    [IDX,Dist] = knnsearch(x,x,'K',size);
    
    LX = zeros(N,N);
    F1 = zeros(N,N);
    F2 = zeros(N,N);
    
    
    Y11 = -(1/2)*sqrt(3/(2*pi))*(x(:,1)+1i*x(:,2));
   
    
    %Set initial conditions
    
    u1 = Y11;
    u2 = -sqrt(-2)*Y11;
    u = [u1;u2];
        
 gam1 = 1/N^3;
    gam2 = gam1;
    
    for i = 1:N
        
        
        indxP = IDX(i,:);
        Px = x(indxP,1);
        Py = x(indxP,2);
        Pz = x(indxP,3);
      
        
        
    % Implement Grady's Method locally
    
    xmat2 = repmat(Px, [1 size]);
    ymat2 = repmat(Py, [1 size]);
    zmat2 = repmat(Pz, [1 size]);
    
    
    xmat = xmat2-xmat2.';
    
    
    nx = repmat(nr(indxP,1),[1 size]);
    
    
    ymat = ymat2-ymat2.';
    ny = repmat(nr(indxP,2),[1 size]);
    
   
    zmat = zmat2-zmat2.';
    nz = repmat(nr(indxP,3),[1 size]);
%     

r2 = xmat.^2 + ymat.^2 +zmat.^2; %r^2 for kernel

  %Choose kernel here-Gaussian or IMQ     
%     phi=@(r2,e) 1./sqrt(1+e^2*r2);         % IMQ
%     
%     fobj = @(e) log(cond(phi(r2,e))/KT);
%     
%     eloc = fzero(fobj,2);
   eloc = 2.5;
    
%     phi=@(r22) 1./sqrt(1+eloc^2*r22);         % IMQ local
%     dphi=@(r22) -eloc^2./sqrt(1+eloc^2*r22).^3; % Derivative of IMQ local over r
    
     phi = @(r2) exp(-eloc^2*r2);            %Gaussian
     dphi = @(r2) -2*eloc^2*exp(-eloc^2*r2);      % Derivative of Gaussian
      
    phif= @(r2) exp(-eloc^2*r2);
    
    A2 = dphi(r2);
    
    Dx = ((1-nx.^2).*xmat - nx.*ny.*ymat - nx.*nz.*zmat).*A2;
    Dy = (-nx.*ny.*xmat + (1-ny.^2).*ymat - ny.*nz.*zmat).*A2;
    Dz = (-nx.*nz.*xmat - ny.*nz.*ymat + (1-nz.^2).*zmat).*A2;
    
    phimod = [phi(r2), ones(size,1); ones(1,size), 0];
    Dxmod = [Dx zeros(size,1)];
    Dymod = [Dy zeros(size,1)];
    Dzmod = [Dz zeros(size,1)];
    
%     Amod = phimod^-1;
%     Amod = Amod(1:size,1:size);
% %     
%     A2 = chol(phi(r2)); Dx = (Dx/A2)/A2.'; Dy = (Dy/A2)/A2.'; Dz = (Dz/A2)/A2.';
%     Lloc = Dx*Dx + Dy*Dy  + Dz*Dz; %This is the surface Laplacian
    
     [L,U] = lu(phimod); Dxmod = (Dxmod/U)/L; Dymod = (Dymod/U)/L; Dzmod = (Dzmod/U)/L; %Solve modified system
    Dxmod(:, end) = []; Dymod(:, end) = []; Dzmod(:, end) = []; %Delete last row to get correct sizes
    
    Llu= Dxmod*Dxmod + Dymod*Dymod  + Dzmod*Dzmod; %This is the local surface Laplacian
    
    
    
    phifmod = [phif(r2), ones(size,1); ones(1,size), 0];
   
    [Lfilt,Ufilt] = lu(phifmod);

%     filt = phi(r2); %Hyper-viscosity with local Lap^p
     filt = eloc^8*((256*eloc^8*r2.^4) - 4096*(eloc^6*r2.^3) +18432*(eloc^4*r2.^2)-24576*(eloc^2 * r2)+6144).*phif(r2);
     filt = [filt zeros(size,1)];
     filt = (filt/Ufilt)/Lfilt;
     filt(:,end) =[];
      F1(i,indxP) = -gam1*filt(1,:); %Filter 1 in hyp system, see derivation.
    F2(i,indxP) = -gam2*filt(1,:); %Filter 2 in hyp system, see derivation.
    Xjwght = Llu(1,:);
    LX(i,indxP) = Xjwght;
    
   
    
     end
      D2 =  [zeros(N,N),eye(N,N);LX, F2];
%        e3 = eig(D2);
% % %     
%      maxe6 = [maxe6 max(real(e3))];
% % % % 
%     figure
%     plot(real(e3),imag(e3),'r*')
%     title('Evals of D2')
%       end
% %       end
%     r = symrcm(LX);
%     LXrcm = LX(r,r);
%     e1 = eig(LXrcm);
%     figure
%     plot(real(e1),imag(e1),'b*')
%     title('Evals of Laplace-Beltrami')
%         display('entering time step')
%   
% % 
         for i = 1:T
% % % %                
% % % % %         % Plot solution
% % % % %         if mod(i, 5) == 0
% % % % %             
% % % % %             % Compute solution on plotting mesh
% % % % %             lambda = A2\(A2.'\u1);
% % % % %             uplot = reshape(Aplot*lambda, sz);
% % % % %             
% % % % %             % Plot interpolated solution
% % % % %             surf(xx, yy, zz, uplot);
% % % % %             shading interp; daspect([1 1 1]); axis tight;
% % % % %             colorbar
% % % % %             drawnow
% % % % %             
% % % % %             
% % % % %         end
% % 
        % Fourth-Order Runge-Kutta
        d1 = dt*D2*u;         
        d2 = dt*D2*(u + (1/2)*d1);    
        d3 = dt*D2*(u + (1/2)*d2);     
        d4 = dt*D2*(u + d3);            
        u = u + (1/6)*(d1 + 2*d2 + 2*d3 + d4);
        u1 = u(1:N);
        i
        
        
        end
        
        uexact = exp(-sqrt(-2)*time(end)).*Y11;
        uerror = uexact-u1;
        erropt = [erropt norm(uerror,2)];
%         
%        