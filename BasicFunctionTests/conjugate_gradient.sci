function hist = cg_1(A, x, b, toleration)
    r = b - A * x;
    p = r;
    iter=1;
    hist(iter,:) = [iter,(r' * r),(x'*x)];
    while hist(iter,2) > toleration,
        iter=iter+1;
        rrip = (r' * r);
        Ap=A * p;
        alp = rrip / (p' * Ap);
        r = r - alp * Ap;
        x = x + alp * p;
        hist(iter,:) = [iter,(r' * r),(x'*x)];
        r_temp=(r' * r);
        //printf("r_norm=%f\n",r_temp);
        //printf("x_nrom=%f\n",x'*x);
        printf("iter=%f\n",iter);
        bet = r_temp / rrip
        p = r + bet * p 
        if  iter > 100000
            break;
        end
    end
endfunction

function hist = cg_2(A, x, b, toleration)
    r = b - A * x;
    p = r;
    iter=1;
    hist(iter,:) = [iter,(r' * r)];
    bip = b' * b;  // wrong when b=0 vector
    while hist(iter,2) > toleration,
        iter=iter+1;
        rrip = (r' * r);
        Ap=A * p;
        alp = rrip / (p' * Ap);
        r = r - alp * Ap;
        x = x + alp * p;
        hist(iter,:) = [iter,(r' * r)];
        r_temp=(r' * r);
        printf("r_norm=%f\n",r_temp);
        printf("x_nrom=%f\n",x'*x);
        bet = r_temp / rrip
        p = r + bet * p 
        if  r_temp/bip < 1.0e-5 then
            break;
        end
        printf("r_temp/bip=%f\n",r_temp/bip);
    end
endfunction

stacksize('max')
//cd D:\WorkSpace\SparseLinAlgScilab\BasicFunctionTests
// function 11 : mminfo.sci in other file
cd C:\Users\sc2012\Documents\GitHub\SparseLinAlgScilab\BasicFunctionTests
exec("mminfo.sci")
// funcition 12 : mmread in other file
exec("mmread.sci")
// function 13 : mmwrite in other file
exec("mmwrite.sci")

//filename="s1rmq4m1.mtx";
//filename_b="s1rmq4m1_coord.mtx";
filename="bcsstk26.mtx";

[A,rows,cols,entries,rep,field,symm] = mmread(filename);
//[b,rows_b,cols_b,entries_b,rep_b,field_b,symm_b] = mmread(filename_b);
A=full(A);
//x=ones(cols,2);
//x(:,2)=x(:,2)+1;
x=ones(cols,1);
//x=x*1e-4;
b=zeros(cols,1);
tol = 1.0e-10;
//r_history = cg_2(A, x, b(:,2), tol)
r_history = cg_1(A, x, b, tol)

