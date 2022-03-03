clear all
close all
clc

% Scan parameters
TE1     = 1e-6;       % shift of readout gradient in s, for first acquisition
TE2     = 50e-6;      % shift of readout gradient in s, for second acquisition
BW      = 20000;      % BW of readout gradient in Hz

% Recon parameters
slicenr   = 6;           % slice nr to perform correction on
filt      = 0.5;         % strength of k-space filter (0-1)
gamma     = 2e-10;       % Regularization for B0 mapping
nriters   = 2;           % Nr of outer iters: number of times the B0 map and images are updated
freqrange = [-2000 2000];% Min and max frequency range used for CPR (Hz)
stepsize  = 3;           % stepsize used for CPR (Hz)

% Data location
Datafolder1 = 'data/1/'; % .mat data stored as y1: [readout,phase1,phase2]
Datafolder2 = 'data/2/'; % .mat data stored as y2: [readout,phase1,phase2]

addpath('utils');

%% Read the data  
load([Datafolder1,'data.mat']);
load([Datafolder2,'data.mat']);

% Parameters needed for recon
m         = size(y1,1);
n         = size(y1,2);
dt        = 1/BW;
tvec      = ((-m/2*dt):dt:((m/2-1)*dt))';

%% Filter the k-space data
y1 = filter_kspace(y1,filt);
y2 = filter_kspace(y2,filt);
    
%% Select slice from 3D acquisition
IM1    = fftshift(fftshift(ifft3cc(y1),3),2);
slice1 = squeeze(IM1(:,:,slicenr));
S      = angle(slice1(round(end/2),round(end/2)));
slice1 = slice1.*exp(-1i*S); % Phase correction to avoid phase wraps
y1     = fftshift(fft2(fftshift(slice1)));

IM2    = fftshift(fftshift(ifft3cc(y2),3),2);
slice2 = squeeze(IM2(:,:,slicenr));
slice2 = slice2.*exp(-1i*S); % Phase correction to avoid phase wraps
y2     = fftshift(fft2(fftshift(slice2)));

DFTx = 1/sqrt(m)*dftmtx(m);
DFTy = 1/sqrt(n)*dftmtx(n);
DFTx = circshift(DFTx,ceil(size(DFTx)/2));
DFTy = circshift(DFTy,ceil(size(DFTy)/2));

kspace1 = squeeze(y1);
image1  = rot90(sqrt(m*n)*DFTx*kspace1*DFTy',3);
kspace2 = squeeze(y2);
image2  = rot90(sqrt(m*n)*DFTx*kspace2*DFTy',3);

imagemask = abs(image1)>0.18*max(max(abs(image1)));

figure;
subplot(2,2,1); imagesc(abs(image1).*imagemask);            colormap gray; axis image; title('Image 1 magnitude')
subplot(2,2,2); imagesc(angle(image2).*imagemask,[-pi pi]); colormap gray; axis image; title('Image 1 phase')
subplot(2,2,3); imagesc(abs(image2).*imagemask);            colormap gray; axis image; title('Image 2 magnitude')
subplot(2,2,4); imagesc(angle(image2).*imagemask,[-pi pi]); colormap gray; axis image; title('Image 2 phase')

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
B0     = zeros(m,n);

data_in(:,:,1) = image1;
data_in(:,:,2) = image2;

for iter = 1:nriters
    % Correcting images
    disp(['Performing CPR, ITERATION ',num2str(iter)])
    data_out = CPR_correct(data_in,B0,tvec,freqrange,stepsize);
    u1 = reshape(data_out(:,:,1),[],1);
    u2 = reshape(data_out(:,:,2),[],1);
    
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
subplot(2,2,3);imagesc(B0maps(:,:,1)/(2*pi),[-3000 3000]);           axis image; colorbar; title('Begin B0'); 
subplot(2,2,4);imagesc(B0maps(:,:,nriters)/(2*pi),[-3000 3000]);     axis image; colorbar; title('Final B0'); 
colormap gray


