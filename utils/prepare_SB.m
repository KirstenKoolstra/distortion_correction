function [ ER, RcEc, y1t, y2t ] = prepare_SB(E,R,y1,y2)

global m n 

% Prepare matrices and rhs
ER  = zeros(size(E,1)*m,size(E,2));
y1t = zeros(m*n,1);
y2t = zeros(m*n,1);
for ind = 1:n
    x = (ind-1)*size(E,1)+1;
    ER(x:(x+size(E,1)-1),:)  = E.*R(:,:,ind);
    y1t(x:(x+size(E,1)-1),:) = y1(:,ind);
    y2t(x:(x+size(E,1)-1),:) = y2(:,ind);
end
RcEc   = ER';

end

