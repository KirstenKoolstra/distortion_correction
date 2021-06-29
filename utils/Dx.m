function d = Dx(u)

    global m n
    ur          = reshape(u,m,n);
    d           = zeros(m,n);
    d(:,2:n)    = ur(:,2:n)-ur(:,1:n-1);
    d(:,1)      = ur(:,1)-ur(:,n);
    d           = reshape(d,[],1);
end

