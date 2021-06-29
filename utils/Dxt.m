function d = Dxt(u)
    global m n
    ur         = reshape(u,m,n);
    d          = zeros(m,n);
    d(:,1:n-1) = ur(:,1:n-1)-ur(:,2:n);
    d(:,n)     = ur(:,n)-ur(:,1);
    d          = reshape(d,[],1);
end
