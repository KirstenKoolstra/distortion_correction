function u = SB(ER,RcEc,y1t,nBreg,nInner,u0,tol,maxit)

global mu lambda m n

teller = 1;
u  = u0;
xx = zeros(m*n,1);
yy = zeros(m*n,1);
bx = zeros(m*n,1);
by = zeros(m*n,1);

figure();       
colormap(gray(256));
h=imagesc(abs( reshape(u,m,n)));     axis equal tight;      title('Reconstructed Image');
drawnow

f1 = y1t;
for outer = 1:nBreg
    % put noise back in rhs
    f1=f1+y1t-ER*u;  
    for inner = 1:nInner

        rhs1 = RcEc*f1;
        rhs2 = Dxt(xx-bx)+Dyt(yy-by);
        rhs  = mu*rhs1 + lambda*rhs2;
        
        [u,flag,~,iter(teller)]   = pcg(@funcSB,rhs,tol,maxit,[],[],u); 
        disp(['Outer:',num2str(outer),', PCG iterations:' num2str(iter(teller))]);
        if flag~= 0
            warning(['PCG flag: ' num2str(flag)]);
        end
        iter(teller);
        
        T = abs(reshape(u,m,n));
        set(h ,'CData',T);   
        drawnow
        
        % update x and y
        dxu     = Dx(u);
        dyu     = Dy(u);
        xx      = shrink( dxu+bx, 1/lambda ); 
        yy      = shrink( dyu+by, 1/lambda ); 

        % update Bregman params
        bx = bx+(dxu-xx);
        by = by+(dyu-yy);
        
        teller = teller + 1;
    end
end

end

