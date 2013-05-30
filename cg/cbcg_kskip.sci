function [x, hist] = cbcg_kskip(A, b, Min, Max, k, Size)
    rho0 = Max-Min
    rho1 = Max+Min
    rho2 = rho0^2
    rho3 = rho0*rho1
    rho4 = rho1^2
    x = zeros(b)
    r = b
    p = r
    hist(1,:) = [0,(r' * r)]
    for i=1:(floor(2*Size/(k+1))+10)
        for j=1:k+2
            if j == 1 then
                S(:,1) = r
                Q(:,1) = p
            elseif j == 2 then
                S(:,2) = (2/rho0)* A * r - (rho1/rho0) * r
                Q(:,2) = (2/rho0)* A * p - (rho1/rho0) * p
            elseif j == k+2 then
                Q(:,j) = (4/rho0) * A * Q(:,j-1) - (2*rho1/rho0) * Q(:,j-1) - Q(:,j-2)
            else
                S(:,j) = (4/rho0) * A * S(:,j-1) - (2*rho1/rho0) * S(:,j-1) - S(:,j-2)
                Q(:,j) = (4/rho0) * A * Q(:,j-1) - (2*rho1/rho0) * Q(:,j-1) - Q(:,j-2)
            end
        end
        for j=1:k
            Delta(1, j) = (S(:,j))'*S(:,j)
            Zeta(1, j)  = (S(:,j))'*S(:,j+1)
            Eta(1, j)   = (S(:,j))'*Q(:,j)
            Theta(1, j) = (S(:,j))'*Q(:,j+1)
            Iota(1, j)  = (Q(:,j))'*Q(:,j)
            Kappa(1, j) = (Q(:,j))'*Q(:,j+1)
        end
        Delta(1, k+1) = (S(:,k+1))'*S(:,k+1)
        Eta(1, k+1)   = (S(:,k+1))'*Q(:,k+1)
        Theta(1, k+1) = (S(:,k+1))'*Q(:,k+2)
        Iota(1, k+1)  = (Q(:,k+1))'*Q(:,k+1)
        Iota(1, k+2)  = (Q(:,k+2))'*Q(:,k+2)
        Kappa(1, k+1) = (Q(:,k+1))'*Q(:,k+2)

        for j=1:k+1
          Alpha = 2*Delta(j, 1) / (rho0*Kappa(j,1)+rho1*Iota(j,1))
          Delta(j+1, 1) = Delta(j, 1)...
                          - Alpha*(rho0*Theta(j, 1)+rho1*Eta(j, 1))...
                          + (Alpha^2)*((rho2/4)*Iota(j,2)+(rho3/2)*Kappa(j,1)+(rho4/4)*Iota(j,1))
          Beta = Delta(j+1, 1) / Delta(j, 1)
          Eta(j+1, 1) = Delta(j+1, 1)
          Iota(j+1, 1) = Delta(j+1, 1) + (Beta^2)*Iota(j,1)
          for l=2:k+2-j
            Delta(j+1, l) = Delta(j, l)...
                            - Alpha*((rho0/2)*(Theta(j, l-1)+Theta(j, l))+rho1*Eta(j, l))...
                            + (Alpha^2)*((rho2/16)*(Iota(j,l-1)+2*Iota(j,l)+Iota(j,l+1)+2*Iota(j,2)-2*Iota(j,1))...
                                         +(rho3/4)*(Kappa(j,l-1)+Kappa(j,l))...
                            +(rho4/4)*Iota(j,l))
            Eta(j+1, l) = Delta(j+1, l)...
                          + Beta*Eta(j,l)...
                          -Alpha*Beta*((rho0/4)*(Kappa(j,l-1)+Kappa(j,l))+(rho1/2)*Iota(j,l))
//            Iota(j+1, l) = Eta(j+1, l)...
//                          + Beta*Eta(j,l) + (Beta^2)*Iota(j,l) ...
//                          -Alpha*Beta*((rho0/4)*(Kappa(j,l-1)+Kappa(j,l))+(rho1/2)*Iota(j,l))
            Iota(j+1, l) = 2*Eta(j+1, l) - Delta(j+1, l) + (Beta^2)*Iota(j,l)
          end
          for l=1:k+1-j
            if l==1 then
              Zeta(j+1, 1) = Zeta(j, 1)...
                             - Alpha*(rho0*Eta(j,2)+rho1*Theta(j,1))...
                             +(Alpha^2)*((rho2/8)*(Kappa(j,1)+Kappa(j,2))...
                                         +(rho3/2)*Iota(j,2)+(rho4/4)*Kappa(j,1))
              Theta(j+1, 1) = Zeta(j+1, 1) + Beta*Theta(j, 1)...
                              - Alpha*Beta * ((rho0/2)*Iota(j,2)...
                                              +(rho1/2)*Kappa(j,1))
//              Kappa(j+1, 1) = Theta(j+1, 1) -(Beta^2)*(rho1/rho0)*Iota(j,1)
//              Kappa(j+1, 1) = Theta(j+1, 1) + Beta*Theta(j,1) + (Beta^2)*Kappa(j,1)...
//                             - Alpha * Beta * ((rho0/2) * Iota(j, 2)+ (rho1/2)*Kappa(j,1))
              Kappa(j+1, 1) = 2*Theta(j+1, 1) -Zeta(j+1,1) + (Beta^2)*Kappa(j,1)
            else
              Zeta(j+1, l) = Zeta(j,l)...
                             - Alpha * ((rho0/2)*(Eta(j,l)+Eta(j,l+1)+Eta(j,2)-Eta(j,1))+rho1*Theta(j,l))...
                             + (Alpha^2) * ((rho2/16)*(Kappa(j,l-1)+2*Kappa(j,l)+Kappa(j,l+1)+Kappa(j,2)-Kappa(j,1))+(rho3/4)*(Iota(j,l)+Iota(j,l+1)+Iota(j,2)-Iota(j,1))+(rho4/4)*Kappa(j,l))
              Theta(j+1, l) = Zeta(j+1, l) + Beta*Theta(j,l) - Alpha*Beta*((rho0/4)*(Iota(j,l)+Iota(j,l+1)+Iota(j,2)-Iota(j,1))+(rho1/2)*Kappa(j,l))
              Kappa(j+1, l) = 2*Theta(j+1,l) - Zeta(j+1,l) + (Beta^2)*Kappa(j,l)
            end
          end
          x = x + Alpha * p
          r = r - Alpha * A * p
//          r = b - A*x
          p = r + Beta * p
        end
        hist((i+1),:) = [i*(k+1),(r' * r)]
//        disp([i*(k+1),r' * r])
        if (norm(r) / norm(b)< 1e-10) then
            break
        end
    end
endfunction

function cbcg_kskip_main(filename, A, Min, Max, Size, b, n, k)
    [x,hist]=cbcg_kskip(A,A*b(1:Size,1),Min,Max,k, Size);
    for i = 2:n
        [x,tmp]=cbcg_kskip(A,A*b(1:Size,i),Min,Max,k, Size);
	if hist($,1) < tmp($,1) then
            hist = tmp
	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
            hist = tmp
        end
    end
    str=sprintf("cbcg_kskip_%s_k=%d_sample=%d.txt",filename,k,n)
    fprintfMat(str,hist,"%10.7e")
    disp([k,hist($,1)])
endfunction
