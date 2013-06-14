//spmv test
//#1
stacksize('max')
//cd D:\WorkSpace\SparseLinAlgScilab\BasicFunctionTests
// function 11 : mminfo.sci in other file
//cd C:\Users\sc2012\Documents\GitHub\SparseLinAlgScilab\BasicFunctionTests
cd /home/scl/SparseLinAlgScilab-gh/BasicFunctionTests
exec("mminfo.sci")
// funcition 12 : mmread in other file
exec("mmread.sci")
// function 13 : mmwrite in other file
exec("mmwrite.sci")

//filename="/home/scl/MStore/SPD/shallow_water2.mtx"
filename="/home/scl/MStore/SPD/bcsstk26.mtx"

[A,rows,cols,entries,rep,field,symm] = mmread(filename);
A=full(A);
//x=ones(cols,2);
//x(:,2)=x(:,2)+1;
x=ones(rows,1);
//b=x;
Ax=A*x;
//r0=b-Ax;
//ip= r0'* r0
