stacksize('max')
//cd D:\WorkSpace\SparseLinAlgScilab\BasicFunctionTests
// function 11 : mminfo.sci in other file
//Windows 1
//cd C:\Users\sc2012\Documents\GitHub\SparseLinAlgScilab\BasicFunctionTests
//ubuntu 1
cd /home/scl/SparseLinAlgScilab-gh/BasicFunctionTests
exec("mminfo.sci");
// funcition 12 : mmread in other file
exec("mmread.sci");
// function 13 : mmwrite in other file
exec("mmwrite.sci");
exec("gerschgorin.sci");

//filename="s1rmq4m1.mtx";
//filename_b="s1rmq4m1_coord.mtx";
filename="/home/scl/MStore/SPD/bcsstk26.mtx";

[A,rows,cols,entries,rep,field,symm] = mmread(filename);
//[b,rows_b,cols_b,entries_b,rep_b,field_b,symm_b] = mmread(filename_b);
A=full(A);
[Min, Max] = gerschgorin(A, rows);
