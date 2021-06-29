function d = Dy(u)

    global m n 
    ur       = reshape(u,m,n);
    d        = zeros(m,n);
    d(2:m,:) = ur(2:m,:)-ur(1:m-1,:);
    d(1,:)   = ur(1,:)-ur(m,:);
    d        = reshape(d,[],1);

end

