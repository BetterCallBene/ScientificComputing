rvecA = 	[1, 1,  1,  2, 5, 3, 5, 5,  7,  7];
cvecA = 	[2, 1,  3,  1, 2, 4, 1, 6,  7,  6];
vvecA = [4, -1, -4, 1, 2, 6, 4, 10, 11, 20];

rvecB = [1, 1, 2, 2,  3, 5, 7];
cvecB = [3, 6, 3, 2, 4, 6, 7];
vvecB = [-1, -4, 1, 5, 2, 4, -5];

A = sparse(rvecA, cvecA, vvecA, 8, 8);
B = sparse(rvecB, cvecB, vvecB, 8, 8);