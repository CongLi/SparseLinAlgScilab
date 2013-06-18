// function 11 : mminfo.sci in other file
exec("mminfo.sci")
// funcition 12 : mmread in other file
exec("mmread.sci")
// function 13 : mmwrite in other file
exec("mmwrite.sci")

function Matrix = MakeMat(Min, Max, rho, n)
    Matrix = eye(n,n)
    for i=1:n
    	Matrix(i,i) = Min + (rho^(n-i))*(i-1)*(Max-Min)/(n-1)
    end
endfunction

// function 9 : MakeTri
function A = MakeTri(Diag, n)
    A = Diag * eye(n,n)
    for i=1:n-1
    	A(i,i+1) = -1
	A(i+1,i) = -1
    end
endfunction

// function 10 : MakeP
function P = MakeP(Diag, n)
    P = (1/Diag) * eye(n,n)
endfunction

// function 14 : Matrix
//function [A,rows,cols,entries,rep,field,symm] = mmread(filename)
function [A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename)
    if filename=="TRI"
       A = MakeTri(2.0, 100)
       Size = 100
    else
       [a0,num_rows,num_cols,entries,a4,a5,a6] = mmread(filename);
       A = full(a0);
    end
    [m,n] = size(A)
    d = diag(A)
    for i=1:m
        A(i,:) = A(i,:) ./ sqrt(d(i))
        A(:,i) = A(:,i) ./ sqrt(d(i))
    end
endfunction

function [A,num_rows,num_cols,entries] = Matrix_nonprecondtioned(filename)
    if filename=="TRI"
       A = MakeTri(2.0, 100)
       Size = 100
    else
       [a0,num_rows,num_cols,entries,a4,a5,a6] = mmread(filename);
       A = full(a0);
    end
endfunction
