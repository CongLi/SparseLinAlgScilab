// block conjugate gradient method implementations
//
function [converge_label, max_difference] =  Converge_Checking (R , B , epsilon)
    converge_label = %t;
    max_difference = -%inf;
    
    assert_checktrue (size(R) == size(B));
    
    [num_rows, num_cols] = size(R);
    assert_checktrue(num_cols>0);
//    disp (num_cols);
    if num_cols>1 then
        for i = 1: num_cols
            norm_list_R(i) = norm ( R(:,i) );
//          disp (norm_list_R(i));
            norm_list_B(i) = norm ( B(:,i) );
//          disp (norm_list_B(i));
            check_sign = norm_list_R(i) < (norm_list_B(i) * epsilon);
            if ~check_sign then
                converge_label = %f;
                if ( norm_list_R(i) / (norm_list_B(i)) - epsilon ) > max_difference then
                    max_difference = norm_list_R(i) / (norm_list_B(i)) - epsilon;
                end
            end
        end
    else
        r_norm = norm(R);
        check_sign = r_norm < epsilon;
        if ~check_sign then
            converge_label = %f;
            max_difference = r_norm;
        end
    end

endfunction
//
function [X, hist] = bcg(A, B, max_iters,epsilon)
//    disp ("size of b_rhs_m is :", size(b_rhs_m));
   
    X = ones(B);
    R = B - A * X;
//    disp(R);
    P = R;

    A_P = A * P;
    
    Ptranspose_A_P = P' * A_P;
//    disp(Ptranspose_A_P);
    Ptranspose_A_P_inverse = inv(Ptranspose_A_P);
//    disp(Ptranspose_A_P_inverse);
    Rtranspose_R = R' * R;
//    disp (Rtranspose_R);
    
    for i=1:max_iters
        LAMBDA = Ptranspose_A_P_inverse * Rtranspose_R;
//        disp (LAMBDA)
        X = X + P * LAMBDA;
        R = R - A_P * LAMBDA;
        
        [converge_label , max_difference] = Converge_Checking (R,B,epsilon);
        disp (max_difference);
        hist(i,:) = [i,max_difference];
        if converge_label then
            disp ("converged ...........");
            break;
        end
        
        PSI = inv (Rtranspose_R);
//        disp (PSI);
        Rtranspose_R = R' * R;
//        disp (sqrt(Rtranspose_R));
         
        PSI = PSI * Rtranspose_R;
//        disp (PSI);
        P = R + P * PSI;
        
        A_P = A * P;
        Ptranspose_A_P = P' * A_P;
        Ptranspose_A_P_inverse = inv(Ptranspose_A_P);
        
//        disp (P);
    end
endfunction

function bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon)
    [num_rows, num_cols] = size(A);
//    disp (size(b(1:num_rows,1:rhs_m)));
    [x,hist]=bcg( A , b(1:num_rows,1:rhs_m) , max_iters,epsilon);
    for i = 2:num_samples
        [x,tmp]=bcg( A , b(1:num_rows, ((i-1)*rhs_m + 1):(i* rhs_m)), max_iters,epsilon); 
//      if hist($,1) < tmp($,1) then
//            hist = tmp
//	     elseif hist($,1) == tmp($,1) & hist($,2) < tmp($,2) then
//              hist = tmp
//      end
    end
//    str=sprintf("cbcg_%s_k=%d_sample=%d.txt",filename,k,n)
//    fprintfMat(str,hist,"%10.7e")
//    disp([k,hist($,1)])
endfunction

stacksize('max')
//cd D:\WorkSpace\SparseLinAlgScilab\cg
//cd /home/sc2012/SparseLinAlgScilab-gh/cg
cd /home/scl/SparseLinAlgScilab-gh/cg
exec('Matrix.sci');

epsilon = 1e-20;
max_iters = 100;
rhs_m = 4;
num_samples = 1;
//b=fscanfMat("/home/skkk/ExperimentsRandom/Random");
b=rand(5000,rhs_m * num_samples);
//b=zeros(5000,rhs_m * num_samples);

//filename="/home/sc2012/MStore/SPD/bcsstk26.mtx";
//filename="/home/sc2012/MStore/SPD/shallow_water2.mtx";
//filename="/home/sc2012/MStore/SPD/nasasrb.mtx";

//filename="/home/scl/MStore/SPD/bcsstk26.mtx";
//filename="/home/scl/MStore/SPD/sts4098.mtx";
filename="/home/scl/MStore/SPD/crystm01.mtx";

[A,num_rows,num_cols,entries] = Matrix_precondtioned_1(filename); // the returned matrix is preconditioed
//[A,num_rows,num_cols,entries] = Matrix_nonprecondtioned(filename); // the returned matrix is nonpreconditioed

bcg_main(filename, A, b, rhs_m, num_samples, max_iters,epsilon);
