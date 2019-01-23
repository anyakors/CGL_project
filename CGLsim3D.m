% CGLsim2D.m
% Copyright David M. Winterbottom 2005

% ************************************************************
% !Edited code! original copyright saved
% This piece is built on the basis of the 2D simulation
% but simulates 3D GL eqn with the same methods.
% Edited by Anna Korsakova, 2018
%
% -----------------------------------------------------------
%
% Simulating the complex Ginzburg-Landau equation in 2D using
% pseudo-spectral code and ETD2 exponential time-stepping (See
% Matthews and Cox, J. Comp. Phys., 176 430-455, 2002).
% Periodic boundary conditions are used.
% ************************************************************

disp('*** 2D CGL SIMULATION ***');

% Set equation coefficients
a = 1.5;
b = 1;

% Set system parameters
L    = 50;        % Domain size (assume square container)
Tmax = 50;        % End time
N    = 128;       % Number of grid points
dT   = 0.05;      % Timestep
dps  = 100;       % Number of stored times
co   = 0;         % Whether to continue simulation from final values of previous

% Calculate some further parameters
nmax  = round(Tmax/dT);
XX    = (L/N)*(-N/2:N/2-1)'; 
[X,Y,Z] = meshgrid(XX);
nplt  = floor(nmax/dps);

% Define initial conditions
if co == 0 
	Tdata = zeros(1,dps+1);
	Tdata(1) = 0;
	A = zeros(size(X)) + 10^(-4)*randn(size(X));
else
	A         = Adata(:,:,:,end);
	starttime = Tdata(end);
	Tdata     = zeros(1,dps+1);
	Tdata(1)  = starttime;
	disp('    CARRYING OVER...')
end	

% Set wavenumbers.
k  				= [0:N/2-1 0 -N/2+1:-1]'*(2*pi/L);
k2 				= k.*k;  
k2(N/2+1) = ((N/2)*(2*pi/L))^2;
[k2x,k2y,k2z] = meshgrid(k2); 
del2 			= k2x+k2y+k2z;
Adata     = zeros(N,N,N,dps+1);
A_hatdata = zeros(N,N,N,dps+1);
A_hat            = fftn(A);
Adata(:,:,:,1)     = A;
A_hatdata(:,:,:,1) = A_hat;

% Compute exponentials and nonlinear factors for ETD2 method
cA 	    	= 1 - del2*(1+i*a);
expA 	  	= exp(dT*cA);
nlfacA  	= (exp(dT*cA).*(1+1./cA/dT)-1./cA/dT-2)./cA;
nlfacAp 	= (exp(dT*cA).*(-1./cA/dT)+1./cA/dT+1)./cA;

% Solve PDE
dataindex = 2;
for n = 1:nmax
	T = Tdata(1) + n*dT;
	A = ifftn(A_hat);
	
	% Find nonlinear component in Fourier space
	nlA	= -(1+i*b)*fftn(A.*abs(A).^2);
	
	% Setting the first values of the previous nonlinear coefficients
	if n == 1
		nlAp = nlA;
	end
	
	% Time-stepping
	A_hat = A_hat.*expA + nlfacA.*nlA + nlfacAp.*nlAp;
	nlAp  = nlA;
	
	% Saving data
	if mod(n,nplt) == 0 
		A = ifftn(A_hat);
		Adata(:,:,:,dataindex)     = A; 
		A_hatdata(:,:,:,dataindex) = A_hat; 
		Tdata(dataindex)         = T;
		dataindex                = dataindex + 1;
	end
	
	% Commenting on time elapsed
	if mod(n,floor(nmax/10)) == 0
		outp = strcat('  n= ', num2str(n), ' completed'); disp(outp);
	end
end

% Plot evolution
%figure('position', [200 200 300 300])
%A = reshape(Adata(:,:,:,end),[1024,1024,3]);
%plot3(A(:,1),A(:,2),A(:,3),'.');
A1 = abs(Adata(:,:,:,end));
%A2 = real(Adata(:,:,:,end))
figure,
isosurf(X,Y,Z,A1);
%surf(X,Y, abs(Adata(:,:,:,end)))
%view(0,90), shading interp, axis tight

