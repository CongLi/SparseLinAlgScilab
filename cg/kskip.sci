function [x, hist] = kskip(A, b, k, Size)
    x = zeros(b)
    r = b
    p = r
    hist(1,:) = [0,(r' * r)]
    for i=1:(floor(2*Size/(k+1))+10)
        for j=1:k+1
            if j == 1 then
                R(:,1) = A*r
		P(:,1) = A*p
            else
                R(:,j) = A * R(:,j-1)
                P(:,j) = A * P(:,j-1)
            end
        end
	Gamma = r' * r
	Delta(1, 1:k) = r' * R(:,1:k)
	Delta(1, k+1:2*k) = (R(:,k))' * R(:,1:k)
	Eta(1, 1:k) = r' * P(:,1:k)
	Eta(1, k+1:2*k+1) = (R(:,k))' * P(:,1:k+1)
	Zeta(1, 1:k+1) = p' * P(:,1:k+1)
	Zeta(1, k+2:2*k+2) = (P(:,k+1))' * P(:,1:k+1)
	for j=1:k+1
	    Alpha = Gamma / Zeta(1,1)
	    Beta = (Gamma-Alpha*Eta(1,1))-Alpha*(Eta(1,1)-Alpha*Zeta(1,2))/Gamma
	    Gamma = (Gamma-Alpha*Eta(1,1))-Alpha*(Eta(1,1)-Alpha*Zeta(1,2))
	    for l=1:2*(k+1-j)
	        Delta(1,l) = Delta(1,l)-2*Alpha*Eta(1,l+1)+Alpha*Alpha*Zeta(1,l+2)
		tmp = Eta(1,l)
		Eta(1,l) = Delta(1,l)+Beta*Eta(1,l)-Alpha*Beta*Zeta(1,l+1)
		Zeta(1,l) = Eta(1,l)+Beta*tmp+Beta*Beta*Zeta(1,l)-Alpha*Beta*Zeta(1,l+1)
	    end
	    x = x + Alpha * p
	    r = r - Alpha * A * p
	    p = r + Beta * p
        end
        hist((i+1),:) = [i*(k+1),(r' * r)]
        if (norm(r) / norm(b) < 1e-10) then
            break
        end
    end
endfunction

function kskip_main(filename, A, Size, b, n, k)
    [x,hist]=kskip(A,A*b(1:Size,1),k, Size);
    for i = 2:n
        [x,tmp]=kskip(A,A*b(1:Size,i),k, Size);
	if hist($,1) < tmp($,1) then
            hist = tmp
	elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
            hist = tmp
        end
    end
    str=sprintf("kskip_%s_k=%d_sample=%d.txt",filename,k,n)
    fprintfMat(str,hist,"%10.7e")
    disp([k,hist($,1)])
endfunction
