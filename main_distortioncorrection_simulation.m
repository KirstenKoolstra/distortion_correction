clear all;
close all;

global m n mu lambda RcEc ER;

% Scan parameters
TE1     = 1e-6;       % shift of readout gradient in s, for first acquisition
TE2     = 50e-6;      % shift of readout gradient in s, for second acquisition
BW      = 20000;      % BW of readout gradient in Hz
fovx    = 0.192;      % FOV in readout direction
fovy    = 0.192;      % FOV in phase enc direction

% Recon parameters
m       = 128;        % Not tested for non-square
n       = 128;
gamma   = 2e-10;      % Regularization for B0 mapping
fac0    = 1;          
mu      = 1/fac0;     % Regularization for image recon: data fidelity
lambda  = 5e-9/fac0;  % Regularization for image recon: total variation
nriters = 2;          % Nr of outer iters: number of times the B0 map and images are updated

addpath('utils');

% Prepare phantom
im = phantom('Shepp-Logan',m); 
im(im>=0.1)=0.04;
im = imresize(im,[m n]);
imagemask=(im~=0);
imr = im(:);

%% Read in B0 
load('data/B0map.mat'); %in Hz
b = Post_shimmed(:,:,round(end/4));           % Take a slice
b(23,23) = b(23,24);
gyro   = 42.57608*10^6;                       % In Hz/T
B0     = fliplr(b);                           % In mT
B0     = B0*gyro*10^(-3);                     % in Hz

modulfreq = B0(round(end/2),round(end/2));
% Replace NAN by center frequency
B0(isnan(B0)) = modulfreq;

%% Define FOV 
fovx = 5*10^-3*size(B0,1);                    % in m 
fovy = 5*10^-3*size(B0,2);                    % in m 

B0   = imresize(B0,[m n]);
% Demodulate B0 field map
B0   = B0-modulfreq;                          % Demodulation with center frequency

SH = prepareSH(m,n);
B0 = decompose_field(SH,B0.*imagemask);

figure; imagesc(B0,[-1500 1500]);

% Convert to rad/s
B0    = B0*2*pi;
B0_GT = B0; % save true B0 for checking accuracy

%% Define gradients 
dx   = fovx/m;                          % in m
dy   = fovy/n;                          % in m
dxdy = dx*dy;                           % in m^2
gx   = BW/fovx;                         % In Hz/m
gy   = BW/fovx;                         % In Hz/m
rvecx          = ((1:m)-round(m/2))*dx; % In m
rvecy          = ((1:n)-round(n/2))*dy; % In m
Gradvec_read   = rvecx*gx; 
Gradvec_phase  = rvecy*gy; 
Grad_read      = -rot90(repmat(Gradvec_read,[size(B0,2) 1]),1);   % Readout gradient in Hz/m
Grad_phase     = -rot90(repmat(Gradvec_phase,[size(B0,1) 1]),2);  % Phase enc gradient in Hz/m

figure;
subplot(1,3,1); imagesc(B0/(2*pi));     colorbar; title('B0 field');
subplot(1,3,2); imagesc(Grad_read);     colorbar; title('Readout gradient');
subplot(1,3,3); imagesc(Grad_phase);    colorbar; title('Phase enc gradient');

% Convert to rad/s
Grad_phase = Grad_phase*2*pi;
Grad_read  = Grad_read*2*pi;

%% Experimental setup
nrpoints = size(B0,1);
dt       = 1/BW;
tvec     = ((-m/2*dt):dt:((m/2-1)*dt))';
tvec_ph  = ((-n/2*dt):dt:((n/2-1)*dt))';

shiftvec1 = TE1*ones(length(tvec),1);
shiftvec2 = TE2*ones(length(tvec),1);

%% Prepare the encoding matrix E (without phase steps)
disp('Rotating the B0 map + CSMs, and building E...');
E1 = dxdy*exp(-1i*tvec*reshape(B0+Grad_read,1,m*n)-1i*shiftvec1*reshape(B0,1,m*n));
E2 = dxdy*exp(-1i*tvec*reshape(B0+Grad_read,1,m*n)-1i*shiftvec2*reshape(B0,1,m*n));

%% Prepare the encoding matrix for phase encoding steps
R   = zeros(length(tvec),m*n,n);
for index_ph = 1:n
     phvec = tvec_ph(index_ph)*ones(size(tvec));
     R(:,:,index_ph)  = exp(-1i*phvec*reshape(Grad_phase,1,m*n));
end

%% Simulate the data
disp('Simulate the data...');
for index_ph = 1:n
    y1(:,index_ph) = (E1.*R(:,:,index_ph))*imr;
    y2(:,index_ph) = (E2.*R(:,:,index_ph))*imr;
end
 
mm = 0;
varr = 1e-13;
realy1 = real(y1)+sqrt(varr)*randn(size(y1)) + mm;
imagy1 = imag(y1)+sqrt(varr)*randn(size(y1)) + mm;
realy2 = real(y2)+sqrt(varr)*randn(size(y2)) + mm;
imagy2 = imag(y2)+sqrt(varr)*randn(size(y2)) + mm;

y1 = realy1+1i*imagy1;
y2 = realy2+1i*imagy2;

%% Create system matrix for B0 mapping
A1 = spdiags(TE1*ones(m*n,1),0,m*n,m*n);
A2 = spdiags(TE2*ones(m*n,1),0,m*n,m*n);
A3 = spdiags(ones(m*n,1),0,m*n,m*n);

e           = ones(m,1);
dx          = spdiags([-e e], -1:1:0, m,m);
dx(1,end)   = -1;
e           = ones(n,1);
dy          = spdiags([-e e], -1:1:0, n,n);
dy(1,end)   = -1;
Gx          = sparse(kron(eye(n),dx));
Gy          = sparse(kron(dy,eye(m)));
tv          = (Gx'*Gx+Gy'*Gy);

A  = [A1 A3 ; A2 A3];
TV = [tv zeros(size(tv)) ; zeros(size(tv)) tv];

%% Start the reconstruction
SH     = prepareSH(m,n);
images = zeros(m,n,nriters);
B0maps = zeros(m,n,nriters);
B0     = zeros(size(B0));

for iter = 1:nriters
    
    % Update encoding matrix
    disp(['Update encoding matrix, ITERATION ',num2str(iter)])
    E  = dxdy*exp(-1i*tvec*reshape(B0+Grad_read,1,m*n));

    % Reconstruct 2 images
    [ER, RcEc, y1t, y2t] = prepare_SB(E,R,y1,y2);
    disp('Reconstruct image 1...')
    u1 = SB(ER,RcEc,y1t,3,1,zeros(m*n,1),1e-2,500);
    disp('Reconstruct image 2...')
    u2 = SB(ER,RcEc,y2t,3,1,zeros(m*n,1),1e-2,500);
     
    % B0 mapping
    rhs = -[angle(u1); angle(u2)];
    disp('Reconstruct B0 map...')
    sol = pcg((A'*A+gamma*TV),A'*rhs);
    B0  = reshape(sol(1:m*n),m,n);
    L   = reshape(sol((m*n+1):end),m,n);
        
    % Crop mask for spherical harmonics
    filtersize = 3;
    mask_image = abs(reshape(u1,m,n))>0.18*max(max(abs(reshape(u1,m,n))));
    mask_B0    = imfilter(mask_image,1/(filtersize^2)*ones(filtersize,filtersize));
    mask_B0    = (mask_B0>=1);
    
    % Decompose field
    B0 = decompose_field(SH,B0/(2*pi).*mask_B0)*2*pi;
    
    % Save the data
    images(:,:,iter) = reshape(u1,m,n);    
    B0maps(:,:,iter) = B0;

    figure; 
    subplot(1,3,1); imagesc(mask_image);                                    axis image;  colorbar; title('Mask');
    subplot(1,3,2); imagesc(mask_B0);                                       axis image;  colorbar; title('Cropped mask');
    subplot(1,3,3); imagesc(B0maps(:,:,iter)/(2*pi),[-1800 1800]);          axis image;  colorbar; title('New B0 map');
    
end

figure;
subplot(2,2,1);imagesc(abs(images(:,:,1)));                          axis image; colorbar; title('Begin image')
subplot(2,2,2);imagesc(abs(images(:,:,nriters)));                    axis image; colorbar; title('Final image')
subplot(2,2,3);imagesc(zeros(size(B0)),[-1500 1500]);                axis image; colorbar; title('Begin B0'); 
subplot(2,2,4);imagesc(B0maps(:,:,nriters)/(2*pi),[-1500 1500]);     axis image; colorbar; title('Final B0'); 
colormap gray

figure;
subplot(1,3,1); imagesc(B0_GT/(2*pi),[-1500 1500]);                     axis image; colorbar; title('Ground truth B0')
subplot(1,3,2); imagesc(B0maps(:,:,nriters)/(2*pi),[-1500 1500]);       axis image; colorbar; title('Estimated B0')
subplot(1,3,3); imagesc((B0_GT-B0maps(:,:,nriters))/(2*pi),[-500 500]); axis image; colorbar; title('Error')
