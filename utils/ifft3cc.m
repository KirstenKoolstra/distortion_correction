function res = ifft3cc(x)

% original:
res = sqrt(length(x(:)))*fftshift(ifftn(ifftshift(x)),1);   

% gives the same:
% res = fftshift(ifft(ifftshift(x,1),[],1),1);
% res = rot90_3D(res,2,1);
% res = ifftshift(ifftshift(ifftshift(ifft2(fftshift(res))),2),1);  
% res = rot90_3D(res,2,3);
% res = sqrt(length(x(:)))*res;

% But rotation makes it difficult to write as one function, so can also do
% this, but is a little different in phase part.
% res = fftshift(ifft(ifftshift(x,1),[],1),1); % along 1 dimension
% res = sqrt(length(x(:)))*ifftshift(ifftshift(ifftshift(ifft(ifft(fftshift(res),[],3),[],2),2),3)); % along dimension 2 and 3

end

 