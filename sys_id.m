function [A, B, C, Dd, K] =  sys_id(u,y,n, k)  
    % u = input measurements with most recent measure at right
    % y = output measurements with most recent measure at right
    % n = dimension of deterministic subspace
    % k = 
   %% Deterministic Subspace
    % ORT Method by Picci and Katayama
    [p,Ndat] = size(y); [m,Ndat] = size(u);
    N = Ndat-2*(k+1);
    %form hankel matrices from u,y
            
    U0_k = hankel(u(:,1:k),u(:,k:k+N+2));
    Uk_2k = hankel(u(:,k+1:2*k),u(:,2*k:end));
    U = [Uk_2k;U0_k;];
   
    Y0_k = hankel(y(:,1:k),y(:,k:k+N+2));
    Yk_2k = hankel(y(:,k+1:2*k),y(:,2*k:end));
    Y = [Y0_k;Yk_2k];
    km = size(U,1)/2; kp = size(Y,1)/2;
    %LQ Decomp
    L = triu(qr([U;Y]'))';
    L11 = L(1:km,1:km);

    L21 = L(km+1:2*km,1:km);
    L22 = L(km+1:2*km,km+1:2*km);

    L31 = L(2*km+1:2*km+kp,1:km);
    L32 = L(2*km+1:2*km+kp,km+1:2*km);
    L33 = L(2*km+1:2*km+kp,2*km+1:2*km+kp);

    L41 = L(2*km+kp+1:2*km+2*kp,1:km);
    L42 = L(2*km+kp+1:2*km+2*kp,km+1:2*km);
    L43 = L(2*km+kp+1:2*km+2*kp,2*km+1:2*km+kp);
    L44 = L(2*km+kp+1:2*km+2*kp,2*km+kp+1:2*km+2*kp);
    % SVD
    [UU,SS,VV] = svd([L42 L43]);
    % [U,S,~] = svd(L42);
    U1 = UU(:,1:n);
    Ok = U1*sqrtm(SS(1:n,1:n));
    
    %A and C
    Cd = Ok(1:p,1:n);
    Ad = pinv(Ok(1:p*(k-1),1:n))*Ok(p+1:k*p,1:n);
    
    %B and D (from PO-MOESP algorithm, gave better estimates of B and D in
    %testing when D = 0)
    U2 = UU(:,n+1:size(UU',1));
    Z = (U2'*[L31 L32 L41])/[L21 L22 L11];
    XX = [];
    RR = [];
    for j = 1:k
        XX = [XX;Z(:,m*(j-1)+1:m*j)];
        Okj = Ok(1:p*(k-j),:);
        Rj = [zeros(p*(j-1),p),zeros(p*(j-1),n);
             eye(p), zeros(p,n);
             zeros(p*(k-j),p),Okj];
        RR = [RR;U2'*Rj];
    end
    DB = pinv(RR)*XX;
    Dd = DB(1:p,:);
    Bd = DB(p+1:end,:);
   
   %% Stochastic Realization (from ORT Algorithm)
     Sfp = (L43*L33')/N;
     Sff = (L43*L43'+L44*L44')/N;
     Spp = (L33*L33')/N;

     L = sqrtm(Sff);
     M = sqrtm(Spp);
     [U,S,V] = svd(L\Sfp/M');
     U1 = U(:,1:n);
     Ok = L*U1*sqrtm(S(1:n,1:n));
     Ck = sqrtm(S(1:n,1:n))*V(1:n,:)*M';
    
     As = pinv(Ok(1:p*(k-1),1:n))*Ok(p+1:kp,:);
     Cs = Ok(1:p,:);
     Cs_bar = Ck(:,1:p);
     lbd_0 = Sff(1:p,1:p);
     Ks = (Cs_bar-As*S(1:n,1:n)*Cs')/(lbd_0-Cs*S(1:n,1:n)*Cs');
     
     A = blkdiag(Ad, As);
     B = [Bd; zeros(length(As),1)];
     K = [zeros(length(Ad),1);Ks];
     C = [Cd Cs]; 
end
