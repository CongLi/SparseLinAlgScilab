//spmv test
//#1
stacksize('max')
cd D:\WorkSpace\SparseLinAlgScilab\BasicFunctionTests
// function 11 : mminfo.sci in other file
exec("mminfo.sci")
// funcition 12 : mmread in other file
exec("mmread.sci")
// function 13 : mmwrite in other file
exec("mmwrite.sci")

filename="bcsstk26.mtx"

[A,rows,cols,entries,rep,field,symm] = mmread(filename);
A=full(A);
x=ones(cols,1)
Ax=A*x
ip= Ax'* Ax
