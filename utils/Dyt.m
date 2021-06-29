function d = Dyt(u)

    global m n 
    ur         = reshape(u,m,n);
    d          = zeros(m,n);
    d(1:m-1,:) = ur(1:m-1,:)-ur(2:m,:);
    d(m,:)     = ur(m,:)-ur(1,:);
    d          = reshape(d,[],1);

end

