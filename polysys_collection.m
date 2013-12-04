%%katsura3 8 solutions
% truncated co-rank stabilizes to 8 from degree d = 6
% all 8 solutions are found from d >= 6
% [x,y,z,t]
% 2*t^2 + 2*z^2+2*y^2+2*x^2-x
polysys{1,1} = [2 2 2 2 -1];
polysys{1,2} = [0 0 0 2;0 0 2 0; 0 2 0 0;2 0 0 0;1 0 0 0];
% 2*z*t+2*y*z+2*x*y-y
polysys{2,1} = [2 2 2 -1];
polysys{2,2} = [0 0 1 1;0 1 1 0;1 1 0 0;0 1 0 0];
% 2*y*t+2*x*z+y^2-z
polysys{3,1} = [2 2 1 -1];
polysys{3,2} = [0 1 0 1;1 0 1 0;0 2 0 0;0 0 1 0];
% 2*t+2*z+2*y+x-1
polysys{4,1} = [2 2 2 1 -1];
polysys{4,2} = [0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0; 0 0 0 0];

% polynomial optimization problems are eigenvalue problems
polysys{1,1} = [3 2 -1 -1 1];
polysys{1,2} = [0 0; 0 1; 1 1; 0 2; 2 0];
polysys{2,1} = [5 4 3 1 -2 1];
polysys{2,2} = [0 0; 0 1; 1 0; 1 1; 0 2; 2 0];

%%y=x^2;
%%y=-2;
polysys{1,1} = [1 -1];
polysys{1,2} = [0 1;2 0];
polysys{2,1} = [1 2];
polysys{2,2} = [0 1;0 0];

% numerical grobner basis, nagasaka k.
polysys{1,1}=[2,3];polysys{1,2}=[1 0;0 1];
polysys{2,1}=[1,-2];polysys{2,2}=[1 1;0 0];
% added perturbations of 1e-6, set tolerance to 1e-5
polysys{1,1}=[2.000005,3.000001];polysys{1,2}=[1 0;0 1];
polysys{2,1}=[.999999,-2.000003];polysys{2,2}=[1 1;0 0];

% original polynomial system Sturmfels IWT
polysys{1,1} = [-0.2825625e-1 0.004575 0.00027156 -0.179706e-2 -0.3825e-4 0.4252e-3 -0.15195e-2 0.0038808 -0.1338e-4 0.01359375];
polysys{1,2} = [1 0;0 1;1 2;2 1;0 2;1 1;2 0;3 0;0 3; 0 0];
polysys{2,1} = [.004575 -.000775 -.00004052 .00027332 .0000065 -.0000655 .0002115 -.000601 .00000196 -.0021875];
polysys{2,2} = [1 0;0 1;1 2;2 1;0 2; 1 1;2 0;3 0;0 3;0 0];

% % explicit polynomial system CpG (simplify in maple)
% polysys{1,1} = [5.8458450 +3.1088125 +3.1386525 +2.8657500 +2.9781828 2.5099650 2.5767972 .6985440 1.0547406 1.0250000];
% polysys{1,2} = [1 1 1;0 1 2;0 2 1;1 0 2;1 2 0;2 0 1;2 1 0;3 0 0;0 3 0;0 0 3];
% polysys{2,1} = [11.2895325 5.6906875 5.6973150 5.5863750 5.6966598 5.2204275 5.3163882 1.5400800 1.8987696 1.8921875];
% polysys{2,2} = [1 1 1;0 1 2;0 2 1;1 0 2;1 2 0;2 0 1;2 1 0;3 0 0;0 3 0;0 0 3];
% polysys{3,1} = [4.538175 2.305625 2.311350  2.242500  2.293212  2.076975  2.117628 .608040 .771324 .765625];
% polysys{3,2} = [1 1 1;0 1 2;0 2 1;1 0 2;1 2 0;2 0 1;2 1 0;3 0 0;0 3 0;0 0 3];
% polysys{4,1} = [1 1 1 -1];
% polysys{4,2} = [1 0 0;0 1 0;0 0 1;0 0 0];
% 
% % explicit polynomial system CpG (factor in maple)
% polysys{1,1} = [1.169169000 .6217625000 0.6277305001 .573150000 0.5956365600 0.5019930000 .5153594400 0.1397088000 0.2109481201 .205000000];
% polysys{1,2} = [1 1 1;0 1 2;0 2 1;1 0 2;1 2 0;2 0 1;2 1 0;3 0 0;0 3 0;0 0 3];
% polysys{2,1} = [1.128953250 0.56906875 0.5697315 0.5586375 0.5696659801 0.52204275 0.53163882 0.154008 .1898769 0.18921875];
% polysys{2,2} = [1 1 1;0 1 2;0 2 1;1 0 2;1 2 0;2 0 1;2 1 0;3 0 0;0 3 0;0 0 3];
% polysys{3,1} = [1.134543750 .5764062499 .5778374999 .5606249999 0.5733029999 0.5192437499 0.5294069999 0.15201 0.192831 .1914062 ];
% polysys{3,2} = [1 1 1;0 1 2;0 2 1;1 0 2;1 2 0;2 0 1;2 1 0;3 0 0;0 3 0;0 0 3];
% polysys{4,1} = [1 1 1 -1];
% polysys{4,2} = [1 0 0;0 1 0;0 0 1;0 0 0];

% explicit polynomial system CpG (factor in maple)
polysys{1,1} = [1.5625 1.366056 4.755075 1.5625 1.55625 1.59309 2.229375 1.3347 2.300562 2.340625 -118.480275 4.531875001 -14.256 -20.734375 -21.310092 -62.860625 -63.43829999 4.724375001 .3906249999 .387504 .28512 -58.11187499 -60.32910599 -51.05879999 -52.388784];
polysys{1,2} = [0 1 3 1;3 1 0 1;1 2 1 1;1 0 3 1;0 3 1 1;1 3 0 1;2 0 2 1;3 0 1 1;2 2 0 1;0 2 2 1;1 1 1 0;2 1 1 1;3 0 0 0;0 0 3 0;0 3 0 0;0 1 2 0;0 2 1 0;1 1 2 1;0 0 4 1;0 4 0 1;4 0 0 1;1 0 2 0;1 2 0 0;2 0 1 0;2 1 0 0];
polysys{2,1} = [1.5625 1.366056 4.755075 1.5625 1.55625 1.59309 2.229375 1.3347 2.300562 2.340625 -115.7094 4.531875001 -15.914016 -19.34375 -19.3752 -58.13875 -58.17059999 4.724375001 .3906249999 .387504 .28512 -57.298125 -58.344408 -53.75137499 -54.69899401];
polysys{2,2} = [0 1 3 1;3 1 0 1;1 2 1 1;1 0 3 1;0 3 1 1;1 3 0 1;2 0 2 1;3 0 1 1;2 2 0 1;0 2 2 1;1 1 1 0;2 1 1 1;3 0 0 0;0 0 3 0;0 3 0 0;0 1 2 0;0 2 1 0;1 1 2 1;0 0 4 1;0 4 0 1;4 0 0 1;1 0 2 0;1 2 0 0;2 0 1 0;2 1 0 0];
polysys{3,1} = [1.5625 1.366056 4.755075 1.5625 1.55625 1.59309 2.229375 1.3347 2.300562 2.340625 -116.06 4.531875001 -15.6762 -19.53125 -19.6419 -58.78124999 -58.89249999 4.724375001 .3906249999 .387504 .28512 -57.390625 -58.60604999 -53.356875 -54.36210001];
polysys{3,2} = [0 1 3 1;3 1 0 1;1 2 1 1;1 0 3 1;0 3 1 1;1 3 0 1;2 0 2 1;3 0 1 1;2 2 0 1;0 2 2 1;1 1 1 0;2 1 1 1;3 0 0 0;0 0 3 0;0 3 0 0;0 1 2 0;0 2 1 0;1 1 2 1;0 0 4 1;0 4 0 1;4 0 0 1;1 0 2 0;1 2 0 0;2 0 1 0;2 1 0 0];
polysys{4,1} = [1 1 1 -1];
polysys{4,2} = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 0];
    
% small instance of Jukes-Cantor model phylogenetics
polysys{1,1} = [-24 9 9 9 -3 -3 -3 1];polysys{1,2}= [1 1 1;1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
polysys{2,1} = [-48 6 6 6];polysys{2,2} = [1 1 1;1 1 0;1 0 1;0 1 1];
polysys{3,1} = [24 3 -9 -9 3];polysys{3,2} = [1 1 1;1 1 0;1 0 1;0 1 1;0 0 1];
polysys{4,1} = [24 -9 3 -9 3];polysys{4,2} = [1 1 1;1 1 0;1 0 1;0 1 1;0 1 0];
polysys{5,1} = [24 -9 -9 3 3];polysys{5,2} = [1 1 1;1 1 0;1 0 1;0 1 1;1 0 0];

% has 13 solutions: 4 at inf, (0,0)^8 and (1,2)
% multiplicity structure = [d00_(0,0) d10_(0,0) d01_(0,0) d11_(0,0) d02_(0,0) d12_(0,0) d03_(0,0) d13_(0,0) d00_(1,2)]
% y^4*x + 3*x^3 - y^4 - 3*x^2
polysys{1,1} = [1 3 -1 -3];
polysys{1,2} = [1 4;3 0;0 4;2 0];
% x^2*y - 2*x^2
polysys{2,1} = [1 -2];
polysys{2,2} = [2 1;2 0];
% 2*y^4*x - x^3 - 2*y^4 + x^2
polysys{3,1} = [2 -1 -2 1];
polysys{3,2} = [1 4;3 0;0 4;2 0];

% radical of upper example by adding two extra polynomials
% has 2 solutions: (0,0) and (1,2)
% y^4*x + 3*x^3 - y^4 - 3*x^2
polysys{1,1} = [1 3 -1 -3];
polysys{1,2} = [1 4;3 0;0 4;2 0];
% x^2y - 2x^2
polysys{2,1} = [1 -2];
polysys{2,2} = [2 1;2 0];
% 2y^4x - x^3 - 2y^4 + x^2
polysys{3,1} = [2 -1 -2 1];
polysys{3,2} = [1 4;3 0;0 4;2 0];
% x^2-x
polysys{4,1} = [1 -1];
polysys{4,2} = [2 0;1 0];
% y^2-2*y
polysys{5,1} = [1 -2];
polysys{5,2} = [0 2;0 1];

% dayton & zeng: (0,0) multiplicity 12
polysys{1,1}=1;polysys{1,2}=[3 0];
polysys{2,1}=[1,1];polysys{2,2}=[2 1;0 4];

% cox: Ideals, varieties and algorithms, pg 115
polysys{1,1} = [1 1 1 -1];
polysys{1,2} = [2 0 0;0 1 0;0 0 1;0 0 0];
polysys{2,1} = [1 1 1 -1];
polysys{2,2} = [1 0 0;0 2 0;0 0 1;0 0 0];
polysys{3,1} = [1 1 1 -1];
polysys{3,2} = [1 0 0;0 1 0;0 0 2;0 0 0];

% PNLA test system
polysys{1,1} = [2.13, 6, 1, -1];
polysys{1,2} = [2 1 0;0 1 0;0 0 1;0 0 0];
polysys{2,1} = [1, 8, 1.7, 3];
polysys{2,2} = [1 0 0;0 2 0;0 0 1;0 0 0];
polysys{3,1} = [7.9, 32.6, 1, 9];
polysys{3,2} = [1 0 0;0 1 0;2 0 2;0 0 0];

% cox: using algebraic geometry, pg27
% x^2 + y^2 + z^2 - 4
polysys{1,1} = [1 1 1 -4];polysys{1,2} = [2 0 0;0 2 0;0 0 2;0 0 0];
% x^2 + 2* y^2 - 5
polysys{2,1} = [1 2 -5];polysys{2,2} = [2 0 0;0 2 0;0 0 0];
% x*z - 1
polysys{3,1} = [1 -1];polysys{3,2} = [1 0 1;0 0 0]; 

% groebner basis for pg27 from above
gsys{1,1} = [-3 1 2];gsys{1,2} = [0 0 0;2 0 0;0 0 2];
gsys{2,1} = [1 2 -5];gsys{2,2} = [2 0 0;0 2 0;0 0 0];
gsys{3,1} = [1 -1];gsys{3,2} = [1 0 1;0 0 0]; 
gsys{4,1} = [2 -3 1];gsys{4,2} = [0 0 1;1 0 0;3 0 0];

% (x-2) * (x-3)
polysys{1,1} = [1 -5 6];
polysys{1,2} = [2 0;1 0;0 0];
% y + 5x - 7
polysys{2,1} = [1 5 -7];
polysys{2,2} = [0 1;1 0;0 0];

% multiplicity root at infinity of 2
polysys{1,1} = [1 -7 10];polysys{1,2} = [2 1;1 1;0 1];
polysys{2,1} = [1 -2];polysys{2,2} = [1 1;1 0];

% multiple roots x = 2 en x = 3, y = 5 en y = 4
% (x-2) * (x-3) * y
polysys{1,1} = [1 -5 6];
polysys{1,2} = [2 1;1 1;0 1];
% (y-4) * (y-5)
polysys{2,1} = [1 -9 20];
polysys{2,2} = [0 2;0 1;0 0];

% (x-2) * y
polysys{1,1} = [1 -2];
polysys{1,2} = [1 1;0 1];
% (y-3)
polysys{2,1} = [1 -3];
polysys{2,2} = [0 1;0 0];

% (x-2) * y
polysys{1,1} = [1 -2];
polysys{1,2} = [1 1;0 1];
% (y-3)*x
polysys{2,1} = [1 -3];
polysys{2,2} = [1 1;1 0];

% x = 2, y=3, root at infinity (0,1,0) with multiplicity of 2 
% (x-2) * y^2
polysys{1,1} = [1 -2];
polysys{1,2} = [1 2;0 2];
% (y-3) *
polysys{2,1} = [1 -3];
polysys{2,2} = [0 1;0 0];

% (x-2) * y^2
polysys{1,1} = [1 -2];
polysys{1,2} = [1 2;0 2];
% (x-y)
polysys{2,1} = [1 -1];
polysys{2,2} = [1 0;0 1];
% ksb = getKSB(5,1,2,[0 0]);
% root = [ksb(:,1) ksb(:,2)+ksb(:,3) makeRoot(5,[2 2])];
% M = getM(polysys,5);

% (x-2)^2 y^2
polysys{1,1} = [1 -4 4];
polysys{1,2} = [2 2;1 2;0 2];
% y-3 = 0
polysys{2,1} = [1 -3];
polysys{2,2} = [0 1;0 0];

% 1 affine root (2,3) with multiplicity 2, no roots @ inf
% (y-3)^2 = 0
polysys{1,1} = [1 -6 9];polysys{1,2} = [0 2;0 1;0 0];
% x+1 - y = 0
polysys{2,1} = [1 1 -1];polysys{2,2} = [1 0;0 0;0 1];

% 1 affine root (2,3) with multiplicity 4, no roots @ inf
% (y-3)^2 = 0
polysys{1,1} = [1 -6 9];polysys{1,2} = [0 2;0 1;0 0];
% (x+1-y)^2 = 0
polysys{2,1} = [1 1 1 2 -2 -2];polysys{2,2} = [2 0;0 0;0 2;1 0;1 1;0 1];

% (X-2)^2 * y = 0
polysys{1,1} = [1 -4 4];polysys{1,2} = [2 1;1 1; 0 1];
% y^2-9 = 0
polysys{2,1} = [1 -9];polysys{2,2} = [0 2; 0 0];

% http://www.math.uic.edu/~jan/Demo/noon3.html
% V. W. Noonburg:
% A neural network modeled by an adaptive Lotka-Volterra system",
% SIAM J. Appl. Math., Vol. 49, No. 6, 1779-1792, 1989.
% PHCpack finds 8 solutions
% x1*x2^2 + x1*x3^2 - 1.1*x1 + 1;
polysys{1,1} = [1 1 -1.1 1];
polysys{1,2} = [1 2 0;1 0 2;1 0 0;0 0 0];
% x2*x1^2 + x2*x3^2 - 1.1*x2 + 1;
polysys{2,1} = [1 1 -1.1 1];
polysys{2,2} = [2 1 0;0 1 2;0 1 0;0 0 0];
% x3*x1^2 + x3*x2^2 - 1.1*x3 + 1;
polysys{3,1} = [1 1 -1.1 1];
polysys{3,2} = [2 0 1;0 2 1;0 0 1;0 0 0];

% zeng, numerical elimination method for polynomial computations
% very well-conditioned polynomial system for elimination
polysys{1,1} = [-310 959 774 1389 1313];polysys{1,2} = [0 0 0;0 2 0;0 0 2;0 1 1;0 2 2];
polysys{2,1} = [-365 755 917 1451 1269];polysys{2,2} = [0 0 0;0 0 2;2 0 0;1 0 1;2 0 2];
polysys{3,1} = [-413 837 838 1655 1352];polysys{3,2} = [0 0 0;2 0 0;0 2 0;1 1 0;2 2 0];

% http://www.math.uic.edu/~jan/Demo/lorentz.html
% equilibrium points of a 4-dimensional Lorentz attractor
% Tien-Yien Li : "Solving polynomial systems",
% The Mathematical Intelligencer 9(3):33-39, 1987.
% finds all 11 solutions, tcr stabilizes
% x1*x2-x1*x3-x4+ 1;
polysys{1,1} = [1 -1 -1 1];
polysys{1,2} = [1 1 0 0;1 0 1 0; 0 0 0 1;0 0 0 0];
% x2*x3-x2*x4-x1+ 1;
polysys{2,1} = [1 -1 -1 1];
polysys{2,2} = [0 1 1 0;0 1 0 1;1 0 0 0;0 0 0 0];
% -x1*x3+x3*x4-x2+ 1;
polysys{3,1} = [-1 1 -1 1];
polysys{3,2} = [1 0 1 0;0 0 1 1;0 1 0 0;0 0 0 0];
% x1*x4-x2*x4-x3+ 1;
polysys{4,1} = [1 -1 -1 1];
polysys{4,2} = [1 0 0 1;0 1 0 1;0 0 1 0;0 0 0 0];

% http://www.math.uic.edu/~jan/Deflate/mth191.html
% problem with 15 regular roots and 3 roots of multiplicity 4
% There are 27 solutions, counted with multiplicities.
% Besides 9 real regular roots and 6 conjugated complex ones,
% the three solutions (1,0,0), (0,1,0) and (0,0,1) each have
% multiplicity four, which all adds up to 27.
% You can count the multiplicities simply by counting the paths
% that end up close to these roots.
% This problem comes from mth191, a course taught by Bernd Sturmfels.
%  x^3+y^2+z^2-1;
polysys{1,1} = [1 1 1 -1];
polysys{1,2} = [3 0 0;0 2 0;0 0 2;0 0 0];
%  x^2+y^3+z^2-1;
polysys{2,1} = [1 1 1 -1];
polysys{2,2} = [2 0 0;0 3 0;0 0 2;0 0 0];
%  x^2+y^2+z^3-1;
polysys{3,1} = [1 1 1 -1];
polysys{3,2} = [2 0 0;0 2 0;0 0 3;0 0 0];

% sparf does not find the correct solutions because the pure components are
% wrong! Second QR delivers wrong permutation
% http://www.math.uic.edu/~jan/Demo/eco5.html
% 5-dimensional economics problem
%  (x1 + x1*x2 + x2*x3 + x3*x4)*x5 - 1;
polysys{1,1} = [1 1 1 1 -1];
polysys{1,2} = [1 0 0 0 1;1 1 0 0 1;0 1 1 0 1;0 0 1 1 1;0 0 0 0 0];
%          (x2 + x1*x3 + x2*x4)*x5 - 2;
polysys{2,1} = [1 1 1 -2];
polysys{2,2} = [0 1 0 0 1;1 0 1 0 1;0 1 0 1 1;0 0 0 0 0];
%                  (x3 + x1*x4)*x5 - 3;
polysys{3,1} = [1 1 -3];
polysys{3,2} = [0 0 1 0 1;1 0 0 1 1;0 0 0 0 0];
%                            x4*x5 - 4;
polysys{4,1} = [1 -4];
polysys{4,2} = [0 0 0 1 1;0 0 0 0 0];
%                x1 + x2 + x3 + x4 + 1;
polysys{5,1} = [1 1 1 1 1];
polysys{5,2} = [1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 0];

% http://www.math.uic.edu/~jan/Demo/quadfor2
% Jan Verschelde and Karin Gatermann:
%`Symmetric Newton Polytopes for Solving Sparse Polynomial Systems',
%Adv. Appl. Math., 16(1): 95-127, 1995.
% [x1 x2 w1 w2]
%  w1 + w2 - 1;
polysys{1,1} = [1 1 -1];
polysys{1,2} = [0 0 1 0;0 0 0 1; 0 0 0 0];
%  w1*x1 + w2*x2;
polysys{2,1} = [1 1];
polysys{2,2} = [1 0 1 0;0 1 0 1];
%  w1*x1**2 + w2*x2**2 - 2/3;
polysys{3,1} = [1 1 -2/3];
polysys{3,2} = [2 0 1 0;0 2 0 1;0 0 0 0];
%  w1*x1**3 + w2*x2**3;
polysys{4,1} = [1 1];
polysys{4,2} = [3 0 1 0;0 3 0 1];

 strings{1} = 'a[1]^2+a[1]*a[2]+a[2]^2-2*a[1]*a[3]-4*a[2]*a[3]+3*a[3]^2-3*a[1]*a[4]+2*a[2]*a[4]+a[4]^2-3*a[1]-2*a[2]+3*a[3]-2*a[4]-2';
 strings{2} = '2*a[1]^2-a[1]*a[2]+a[2]^2-a[1]*a[3]-a[2]*a[3]-6*a[3]^2-a[1]*a[4]+a[2]*a[4]-5*a[3]*a[4]-3*a[4]^2-5*a[1]+a[2]+5*a[3]+2*a[4]+5';
 strings{3} ='-3-3*a[1]*a[2]+2*a[1]*a[3]+a[1]*a[4]^2-5*a[1]*a[3]^2-5*a[3]^2*a[4]-3*a[1]*a[4]-2*a[3]*a[4]+a[1]*a[2]*a[3]+a[1]*a[2]*a[4]-a[1]^2*a[3]+a[1]^2 -a[2]^2+2*a[3]^2+11*a[3]-2*a[4]-a[1]+a[2]+a[1]^3+a[2]^3-3*a[3]^3+2*a[4]^3-3*a[4]^2-5*a[2]^2*a[3]+7*a[2]*a[3]^2';
 strings{4} = '-15+2*a[1]*a[2]+11*a[1]*a[4]^2+5*a[1]*a[3]^2-a[3]*a[4]-4*a[1]*a[2]*a[3]+6*a[1]*a[2]*a[4]-a[1]^2*a[3]+3*a[1]^2+2*a[2]^2-a[3]^2+4*a[3] -10*a[4]-35*a[1]-14*a[2]-a[1]^3+6*a[2]^3+15*a[3]^3+4*a[4]^3+5*a[4]^2+6*a[2]^2*a[3]+4*a[2]*a[3]^2-a[1]*a[3]*a[4]+6*a[1]^2*a[2] -12*a[1]*a[2]^2-7*a[2]^2*a[4]+2*a[2]*a[4]';
polysys = lti2polysys(strings,[],[],'na',4);

% Patrick Vandewalle's MultiChannel (38), p.10
strings{1}='-11+a[1]+a[2]+a[3]+a[4]+a[5]';
strings{2}='-2.4142135623730950488*b[1]^2+a[1]+a[2]*b[1]+a[3]*b[1]^2+a[4]*b[1]^3+a[5]*b[1]^4';
strings{3}='7-1*a[1]-i*a[2]+a[3]+i*a[4]-1*a[5]';
strings{4}='3.2426406871192851464*b[1]^2-1*a[1]-i*a[2]*b[1]+a[3]*b[1]^2+i*a[4]*b[1]^3-1*a[5]*b[1]^4';
strings{5}='-3+a[1]-1*a[2]+a[3]-1*a[4]+a[5]';
strings{6}='.4142135623730950488*b[1]^2+a[1]-1*a[2]*b[1]+a[3]*b[1]^2-1*a[4]*b[1]^3+a[5]*b[1]^4';
strings{7}='3-1*a[1]+i*a[2]+a[3]-i*a[4]-1*a[5]';
strings{8}='-5.2426406871192851464*b[1]^2-1*a[1]+i*a[2]*b[1]+a[3]*b[1]^2-i*a[4]*b[1]^3-1*a[5]*b[1]^4';
polysys=lti2polysys(strings,[],[],'na',5,'nb',1);

% Patrick Vandewalle's MultiChannel noisy example (40), p.11
strings{1}='-3.48449+a[1]+a[2]+a[3]+a[4]+a[5]';
strings{2}='-2.09168*b[1]^2+a[1]+a[2]*b[1]+a[3]*b[1]^2+a[4]*b[1]^3+a[5]*b[1]^4';
strings{3}='21.2468-1*a[1]-i*a[2]+a[3]+i*a[4]-1*a[5]';
strings{4}='7.44795*b[1]^2-1*a[1]-i*a[2]*b[1]+a[3]*b[1]^2+i*a[4]*b[1]^3-1*a[5]*b[1]^4';
strings{5}='-3.66716+a[1]-1*a[2]+a[3]-1*a[4]+a[5]';
strings{6}='-.729979*b[1]^2+a[1]-1*a[2]*b[1]+a[3]*b[1]^2-1*a[4]*b[1]^3+a[5]*b[1]^4';
strings{7}='9.53102-1*a[1]+i*a[2]+a[3]-i*a[4]-1*a[5]';
strings{8}='19.3078*b[1]^2-1*a[1]+i*a[2]*b[1]+a[3]*b[1]^2-i*a[4]*b[1]^3-1*a[5]*b[1]^4';
polysys=lti2polysys(strings,[],[],'na',5,'nb',1);


% http://www.math.uic.edu/~jan/Demo/wright.html
% system of A.H. Wright, 32 solutions
% x1^2 - x1 + x2 + x3 + x4 + x5 - 10;
polysys{1,1} = [1 -1 1 1 1 1 -10];
polysys{1,2} = [2 0 0 0 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;zeros(1,4) 1;zeros(1,5)];
% x2^2 + x1 - x2 + x3 + x4 + x5 - 10;
polysys{2,1} = [1 1 -1 1 1 1 -10];
polysys{2,2} = [0 2 0 0 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;zeros(1,4) 1;zeros(1,5)];
% x3^2 + x1 + x2 - x3 + x4 + x5 - 10;
polysys{3,1} = [1 1 1 -1 1 1 -10];
polysys{3,2} = [0 0 2 0 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;zeros(1,4) 1;zeros(1,5)];
% x4^2 + x1 + x2 + x3 - x4 + x5 - 10;
polysys{4,1} = [1 1 1 1 -1 1 -10];
polysys{4,2} = [0 0 0 2 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;zeros(1,4) 1;zeros(1,5)];
% x5^2 + x1 + x2 + x3 + x4 - x5 - 10;
polysys{5,1} = [1 1 1 1 1 -1 -10];
polysys{5,2} = [0 0 0 0 2;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;zeros(1,4) 1;zeros(1,5)];

% cannot find correct canonical eigenvectors
polysys{1,1} = [2 1 -6];polysys{1,2} = [2 0;0 2;0 0];
polysys{2,1} = [1 2 -9];polysys{2,2} = [2 0;0 2;0 0];

polysys{1,1} = [4 1 36 4 -24 -12];polysys{1,2} = [4 0;0 4;0 0;2 2;2 0;0 2];
polysys{2,1} = [1 2 -9];polysys{2,2} = [2 0;0 2;0 0];

polysys{1,1} = [2 1 -6];polysys{1,2} = [2 0;0 2;0 0];
polysys{2,1} = [1 4 81 4 -18 -36];polysys{2,2} = [4 0;0 4;0 0;2 2;2 0;0 2];


% cannot find correct canonical eigenvectors ^2
polysys{1,1} = [4 1 36 4 -24 -12];polysys{1,2} = [4 0;0 4;0 0;2 2;2 0;0 2];
polysys{2,1} = [1 4 81 4 -18 -36];polysys{2,2} = [4 0;0 4;0 0;2 2;2 0;0 2];

% test
% y = 2
polysys{1,1} = [1 -2];polysys{1,2} = [0 1 0;0 0 0];
% x+z = 0;
polysys{2,1} = [1 1];polysys{2,2} = [1 0 0;0 0 1];
% x^2 -1 = 0
polysys{3,1} = [1 -1];polysys{3,2} = [2 0 0;0 0 0];

% eigenvector is [1 x x^2] so need to use linear function of x and y as shift to construct the eigenvalue problem
polysys{1,1} = [1 -6 11 -6];polysys{1,2} = [3 0;2 0;1 0;0 0];
polysys{2,1} = [1 -1 -1];polysys{2,2} = [0 1;1 0;0 0];


polysys{1,1} = [8 2 4 1 2 2 2 2];polysys{1,2} = [4 1 0 0 3];
polysys{2,1} = [4 2 1 1 1 2 2 1 1];polysys{2,2} = [0 0 2 0 3;2 1 0 0 1;0 0 1 0 1;0 0 0 0 1;0 1 0 0 0;8 0 4 0 0;0 0 2 0 0;1 0 0 0 0;0 0 1 0 0];
polysys{3,1} = [1 2]; polysys{3,2} = [2 0 1 0 0;1 0 3 0 3];
polysys{4,1} = [8 1 1 4 1 2];polysys{4,2} = [2 2 2 1 0;1 0 0 1 0;0 1 1 0 0;2 1 1 3 0;0 0 0 3 0;0 0 5 0 0;0 0 0 1 0];
polysys{5,1} = [16 2 4 1 2];polysys{5,2} = [3 0 3 6 2;4 0 0 0 2;4 0 5 0 3;1 0 1 0 0;0 0 1 0 1];

% gcd(f1,f2) = gcd(f1,f3) = gcd(f2,f3) = 1
% x^2 +y^2 +z^2 -9 = 0
polysys{1,1} = [1 1 1 -9];polysys{1,2} = [2 0 0;0 2 0;0 0 2;0 0 0];
% x + y + z = 0
polysys{2,1} = [1 1 1];polysys{2,2} = [1 0 0;0 1 0;0 0 1];
% x^2z + yz -5
polysys{3,1} = [1 1 -5];polysys{3,2} = [2 0 1;0 1 1;0 0 0];

% cox, ideals, varieties and algo's, pg 459
polysys{1,1} = [1 3 1 1];polysys{1,2} = [3 2;2 2;0 3;0 0];


% cox, ideals, varieties and algo's pg 105
polysys{1,1} = 1;polysys{1,2} = [1 0 0];
polysys{2,1} = [1 1];polysys{2,2} = [2 0 0;0 0 1];
polysys{3,1} = [1 1];polysys{3,2} = [0 1 0;0 0 1]; 

polysys{1,1} = [1 2 1];polysys{1,2} = [0;1;2];
polysys{2,1} = [1 1 1];polysys{2,2} = [0;1;2];

% Cox, ideals, varieties, pg 107
% gcd(f2,f3) = y, deg(f1f2) = 7, deg(f1f3)=7, deg(f2f3) = 6, maar
% deg(f2f3/gcd) = 5, niet-triviale syzygies
polysys{1,1} = [1 1];polysys{1,2} = [2 2 0;0 0 1];
polysys{2,1} = [1 -1];polysys{2,2} = [1 2 0;0 1 0];
polysys{3,1} = [1 1];polysys{3,2} = [2 1 0;0 1 1];

polysys{1,1} = [1 -1];polysys{1,2} = [1 2 0;0 1 0];
polysys{2,1} = [1 1];polysys{2,2} = [2 2 0;0 0 1];
polysys{3,1} = [1 1];polysys{3,2} = [2 1 0;0 1 1];

% divided y from f2 and f3 from the system above, non-trivial syzygies
polysys{1,1} = [1 1];polysys{1,2} = [2 2 0;0 0 1];
polysys{2,1} = [1 -1];polysys{2,2} = [1 1 0;0 0 0];
polysys{3,1} = [1 1];polysys{3,2} = [2 0 0;0 0 1];

polysys{1,1} = [1 -2];polysys{1,2} = [0 1 0;0 0 0];
polysys{2,1} = [1 1];polysys{2,2} = [1 0 0;0 0 1];
polysys{3,1} = [1 -1];polysys{3,2} = [2 0 0;0 0 0];


% (1+x) (x+y-z)
polysys{1,1} = [1 1 -1 1 1 -1];polysys{1,2} = [1 0 0;0 1 0;0 0 1;2 0 0;1 1 0;1 0 1];
% (1+x) (9+5x-y^2)
polysys{2,1} = [9 14 -1 5 -1];polysys{2,2} = [0 0 0;1 0 0;0 2 0;2 0 0;1 1 0];
% (1+x) (z-2y+7xyz)
polysys{3,1} = [1 -2 7 1 -2 7];polysys{3,2} = [0 0 1;0 1 0;1 1 1;1 0 1;1 1 0;2 1 1];

% non-trivial linear dependencies in rows, {x,y,x+y}
polysys{1,1} = 1;polysys{1,2} = [1 0];
polysys{2,1} = 1;polysys{2,2} = [0 1];
polysys{3,1} = [1 1];polysys{3,2} = [1 0;0 1];

% intersection of sphere with xz plane
polysys{1,1} = [-9 1 1 1];polysys{1,2} = [0 0 0;2 0 0;0 2 0;0 0 2];
polysys{2,1} = [1];polysys{2,2} = [0 1 0];

% x^2 + y^2 + z^2 - 9
polysys{1,1} = [1 1 1 -9];polysys{1,2} = [2 0 0;0 2 0;0 0 2;0 0 0];
% x^2 + y^2 + z^2 - 9
polysys{2,1} = [1 1 1 ];polysys{2,2} = [1 0 0;0 1 0;0 0 1];
% x^2*z + y*z - 5
polysys{3,1} = [1 1 -5];polysys{3,2} = [2 0 1;0 1 1;0 0 0];

% rock-sciccors-paper game
polysys{1,1} = [1 -1];polysys{1,2} = [1 1 0;1 0 1];
polysys{2,1} = [1 -1];polysys{2,2} = [0 1 1;1 1 0];
polysys{3,1} = [1 -1];polysys{3,2} = [1 0 1;0 1 1];
polysys{4,1} = [1 1 1 -1];polysys{4,2} = [1 0 0;0 1 0;0 0 1;0 0 0];

% I^h does not equal ideal with basis f_i^h, Cox p 387
polysys{1,1}=[1 -1];polysys{1,2}=[0 1 0;2 0 0];
polysys{2,1}=[1 -1];polysys{2,2}=[0 0 1;3 0 0];
psys{1,1}=[1 -1];psys{1,2}=[2 0 0 1;1 1 1 0];


% http://www.math.uic.edu/~jan/Demo/conform1.html
% conformal analysis of cyclic molecules
% 16 solutions
% need to go to getDreg(polysys)+1 to find a Grobner basis
% [t1 t2 t3]
%  -9 - t2**2 - t3**2 - 3*t2**2*t3**2 + 8*t2*t3;
polysys{1,1} = [-9 -1 -1 -3 8];
polysys{1,2} = [0 0 0;0 2 0;0 0 2;0 2 2;0 1 1];
%  -9 - t3**2 - t1**2 - 3*t3**2*t1**2 + 8*t3*t1;
polysys{2,1} = [-9 -1 -1 -3 8];
polysys{2,2} = [0 0 0;0 0 2;2 0 0;2 0 2;1 0 1];
%  -9 - t1**2 - t2**2 - 3*t1**2*t2**2 + 8*t1*t2;
polysys{3,1} = [-9 -1 -1 -3 8];
polysys{3,2} = [0 0 0;2 0 0;0 2 0;2 2 0;1 1 0];

% http://www.math.uic.edu/~jan/Demo/boon.html
% neurophysiology, sjirk boon
% total degree = 1024
% 3-homogeneous bezout number = 344
% generalized Bezout number = 216
% here are only 8 finite solutions for general values of
% the constant terms.
% It can be proved that it is equivalent to a quadrature formula
% problem, so that there is only one solution upon symmetry.
% [ s1 s2 g1 g2 C1 C2]
%  s1**2+g1**2 - 1;
polysys{1,1} = [1 1 -1];
polysys{1,2} = [2 0 0 0 0 0;0 0 2 0 0 0;0 0 0 0 0 0];
%  s2**2+g2**2 - 1;
polysys{2,1} = [1 1 -1];
polysys{2,2} = [0 2 0 0 0 0;0 0 0 2 0 0;0 0 0 0 0 0];
%  C1*g1**3+C2*g2**3 - 1.2;
polysys{3,1} = [1 1 -1.2];
polysys{3,2} = [0 0 3 0 1 0; 0 0 0 3 0 1;0 0 0 0 0 0];
%  C1*s1**3+C2*s2**3 - 1.2;
polysys{4,1} = [1 1 -1.2];
polysys{4,2} = [3 0 0 0 1 0; 0 3 0 0 0 1;0 0 0 0 0 0];
%  C1*g1**2*s1+C2*g2**2*s2 - 0.7;
polysys{5,1} = [1 1 -.7];
polysys{5,2} = [1 0 2 0 1 0; 0 1 0 2 0 1;0 0 0 0 0 0];
%  C1*g1*s1**2+C2*g2*s2**2 - 0.7;
polysys{6,1} = [1 1 -.7];
polysys{6,2} = [2 0 1 0 1 0; 0 2 0 1 0 1;0 0 0 0 0 0];

% no zero-dimensional variety!, inf oplossingen
% (x-2) * (y-5) * (y-9) 
polysys{1,1} = [1 -14 45 -2 28 -90];
polysys{1,2} = [1 2; 1 1; 1 0; 0 2; 0 1; 0 0];
% (x-2) * (y-5) * (x+3)
polysys{2,1} = [1 1 -5 -5 -6 30];
polysys{2,2} = [2 1;1 1;2 0; 1 0;0 1; 0 0];

% Non zero-dimensional variety, however, for low degrees valid solutions are
% found
polysys{1,1} = [1 2 -2];
polysys{1,2} = [1 0;1 1;0 1];
polysys{2,1} = [1 -2 1];
polysys{2,2} = [2 0;1 0;0 0];

% M.C. Steenkamp :
% Die numeriese oplos van stelsels polinoomvergelykings.
% Technical report, Nasionale Navorsingsinstituut vir Wiskundige Wetenskappe,
% Pretoria, 1982.
% 10-dimensional system of Ku, 2 solutions
%  5*x1*x2+ 5*x1+ 3*x2+ 55;
polysys{1,1} = [5 5 3 55];
polysys{1,2} = [1 1 zeros(1,8);1 zeros(1,9);0 1 zeros(1,8);zeros(1,10)];
%  7*x2*x3+ 9*x2+ 9*x3+ 19;
polysys{2,1} = [7 9 9 19];
polysys{2,2} = [0 1 1 zeros(1,7);0 1 zeros(1,8);0 0 1 zeros(1,7); zeros(1,10)];
%  3*x3*x4+ 6*x3+ 5*x4-4;
polysys{3,1} = [3 6 5 -4];
polysys{3,2} = [0 0 1 1 zeros(1,6);0 0 1 zeros(1,7);zeros(1,3) 1 zeros(1,6);zeros(1,10)];
%  6*x4*x5+ 6*x4+ 7*x5+ 118;
polysys{4,1} = [6 6 7 118];
polysys{4,2} = [zeros(1,3) 1 1 zeros(1,5);zeros(1,3) 1 zeros(1,6);zeros(1,4) 1 zeros(1,5);zeros(1,10)];
% x5*x6+ 3*x5+ 9*x6+ 27;
polysys{5,1} = [1 3 9 27];
polysys{5,2} = [zeros(1,4) 1 1 zeros(1,4);zeros(1,4) 1 zeros(1,5); zeros(1,5) 1 zeros(1,4);zeros(1,10)];
%  6*x6*x7+ 7*x6+x7+ 72;
polysys{6,1} = [6 7 1 72];
polysys{6,2} = [zeros(1,5) 1 1 zeros(1,3);zeros(1,5) 1 zeros(1,4);zeros(1,6) 1 zeros(1,3);zeros(1,10)];
%  9*x7*x8+ 7*x7+x8+ 35;
polysys{7,1} = [9 7 1 35];
polysys{7,2} = [zeros(1,6) 1 1 0 0;zeros(1,6) 1 0 0 0;zeros(1,7) 1 0 0;zeros(1,10)];
%  4*x8*x9+ 4*x8+ 6*x9+ 16;
polysys{8,1} = [4 4 6 16];
polysys{8,2} = [zeros(1,7) 1 1 0;zeros(1,7) 1 0 0;zeros(1,8) 1 0;zeros(1,10)];
%  8*x9*x10+ 4*x9+ 3*x10-51;
polysys{9,1} = [8 4 3 -51];
polysys{9,2} = [zeros(1,8) 1 1;zeros(1,8) 1 0;zeros(1,9) 1;zeros(1,10)];
%  3*x1*x10-6*x1+x10+ 5;
polysys{10,1} = [3 -6 1 5];
polysys{10,2} = [1 zeros(1,8) 1;1 zeros(1,9);zeros(1,9) 1;zeros(1,10)];

% twisted cubic
polysys{1,1}=[1 -1];polysys{1,2}=[0 1 0;2 0 0];
polysys{2,1}=[1 -1];polysys{2,2}=[0 0 1;3 0 0];

% computing floating point grobner basis stably - sasaki, kako 1)
polysys{1,1}=[1/10 3 1];polysys{1,2}=[3 0;2 1;0 2];
polysys{2,1}=[1 -3 -1];polysys{2,2}=[2 2;1 2;1 1];
polysys{3,1}=[1/10 2];polysys{3,2}=[0 3;2 0];

% computing floating point grobner basis stably - sasaki, kako 2)
polysys{1,1}=[2 1];polysys{1,2}=[2 1 1;2 0 0];
polysys{2,1}=[3 -2];polysys{2,2}=[1 1 2;0 1 1];
polysys{3,1}=[3 -2 -2 -1];polysys{3,2}=[2 0 1;1 0 0;0 1 1;0 0 0];

polysys{1,1}=[1 4 -6 6 -18 13];polysys{1,2}=[2 0;1 1;1 0;0 2;0 1;0 0];
polysys{2,1}=[1 16 -7 118 -286 147 -1 6 1 5];polysys{2,2}=[3 0;2 1;2 0;1 2;1 1;1 0;0 3;0 2;0 1;0 0];
polysys{3,1}=[1 10 -5 72 -176 91 -1 4 1 3];polysys{3,2}=[3 0;2 1;2 0;1 2;1 1;1 0;0 3;0 2;0 1;0 0];


% Sparse 5
% http://homepages.math.uic.edu/~jan/Demo/sparse5.html
strings{1}='a[1]^2*a[2]^2*a[3]^2*a[4]^2*a[5]^2+3*a[1]^2+a[2]^2+a[3]^2+a[4]^2+a[5]^2+a[1]*a[2]*a[3]*a[4]*a[5]+5';
strings{2}='a[1]^2*a[2]^2*a[3]^2*a[4]^2*a[5]^2+a[1]^2+3*a[2]^2+a[3]^2+a[4]^2+a[5]^2+a[1]*a[2]*a[3]*a[4]*a[5]+5';
strings{3}='a[1]^2*a[2]^2*a[3]^2*a[4]^2*a[5]^2+a[1]^2+a[2]^2+3*a[3]^2+a[4]^2+a[5]^2+a[1]*a[2]*a[3]*a[4]*a[5]+5';
strings{4}='a[1]^2*a[2]^2*a[3]^2*a[4]^2*a[5]^2+a[1]^2+a[2]^2+a[3]^2+3*a[4]^2+a[5]^2+a[1]*a[2]*a[3]*a[4]*a[5]+5';
strings{5}='a[1]^2*a[2]^2*a[3]^2*a[4]^2*a[5]^2+a[1]^2+a[2]^2+a[3]^2+a[4]^2+3*a[5]^2+a[1]*a[2]*a[3]*a[4]*a[5]+5';
polysys=lti2polysys(strings,[],[],'na',5);

% Cyc6
% http://homepages.math.uic.edu/~jan/Demo/redcyc6.html
strings{1}='1+a[1]+a[2]+a[3]+a[4]+a[5]';
strings{2}='a[1]+a[1]*a[2]+a[2]*a[3]+a[3]*a[4]+a[4]*a[5]+a[5]';
strings{3}='a[1]*a[2]+a[1]*a[2]*a[3]+a[2]*a[3]*a[4]+a[3]*a[4]*a[5]+a[4]*a[5]+a[5]*a[1]';
strings{4}='a[1]*a[2]*a[3]+a[1]*a[2]*a[3]*a[4]+a[2]*a[3]*a[4]*a[5]+a[3]*a[4]*a[5]+a[4]*a[5]*a[1]+a[5]*a[1]*a[2]';
strings{5}='a[1]*a[2]*a[3]*a[4]+a[1]*a[2]*a[3]*a[4]*a[5]+a[2]*a[3]*a[4]*a[5]+a[3]*a[4]*a[5]*a[1]+a[4]*a[5]*a[1]*a[2]+a[5]*a[1]*a[2]*a[3]';
strings{6}='a[6]^6*a[1]*a[2]*a[3]*a[4]*a[5]-1';
polysys=lti2polysys(strings,[],[],'na',6);

% s9_1
% http://homepages.math.uic.edu/~jan/Demo/s9_1.html
strings{1}='-a[5]*a[7]-2*a[4]*a[8]';
strings{2}='9*a[5]+4*a[2]';
strings{3}='-4*a[3]*a[8]-2*a[5]*a[6]-3*a[4]*a[7]';
strings{4}='-7*a[3]+9*a[1]-8*a[6]';
strings{5}='-4*a[4]*a[6]-5*a[3]*a[7]-6*a[8]-3*a[5]';
strings{6}='-5*a[4]-6*a[3]*a[6]-7*a[7]+9*a[2]';
strings{7}='9*a[4]+6*a[1]-5*a[2]';
strings{8}='9*a[3]-7*a[1]+8';
polysys=lti2polysys(strings,[],[],'na',8);

% geneig
% http://homepages.math.uic.edu/~jan/Demo/geneig.html
strings{1}='-10*a[1]*a[6]^2+2*a[2]*a[6]^2-a[3]*a[6]^2+a[4]*a[6]^2+3*a[5]*a[6]^2+a[1]*a[6]+2*a[2]*a[6]+a[3]*a[6]+2*a[4]*a[6]+a[5]*a[6]+10*a[1]+2*a[2]-a[3]+2*a[4]-2*a[5]';
strings{2}='2*a[1]*a[6]^2-11*a[2]*a[6]^2+2*a[3]*a[6]^2-2*a[4]*a[6]^2+a[5]*a[6]^2+2*a[1]*a[6]+a[2]*a[6]+2*a[3]*a[6]+a[4]*a[6]+3*a[5]*a[6]+2*a[1]+9*a[2]+3*a[3]-a[4]-2*a[5]';
strings{3}='-a[1]*a[6]^2+2*a[2]*a[6]^2-12*a[3]*a[6]^2-a[4]*a[6]^2+a[5]*a[6]^2+a[1]*a[6]+2*a[2]*a[6]-2*a[4]*a[6]-2*a[5]*a[6]-a[1]+3*a[2]+10*a[3]+2*a[4]-a[5]';
strings{4}='a[1]*a[6]^2-2*a[2]*a[6]^2-a[3]*a[6]^2-10*a[4]*a[6]^2+2*a[5]*a[6]^2+2*a[1]*a[6]+a[2]*a[6]-2*a[3]*a[6]+2*a[4]*a[6]+3*a[5]*a[6]+2*a[1]-a[2]+2*a[3]+12*a[4]+a[5]';
strings{5}='3*a[1]*a[6]^2+a[2]*a[6]^2+a[3]*a[6]^2+2*a[4]*a[6]^2-11*a[5]*a[6]^2+a[1]*a[6]+3*a[2]*a[6]-2*a[3]*a[6]+3*a[4]*a[6]+3*a[5]*a[6]-2*a[1]-2*a[2]-a[3]+a[4]+10*a[5]';
strings{6}='a[1]+a[2]+a[3]+a[4]+a[5]-1';
polysys=lti2polysys(strings,[],[],'na',6);

% cyclic5
% http://homepages.math.uic.edu/~jan/Demo/cyclic5.html
strings{1}='a[1]+a[2]+a[3]+a[4]+a[5]';
strings{2}='a[1]*a[2]+a[2]*a[3]+a[3]*a[4]+a[4]*a[5]+a[5]*a[1]';
strings{3}='a[1]*a[2]*a[3]+a[2]*a[3]*a[4]+a[3]*a[4]*a[5]+a[4]*a[5]*a[1]+a[5]*a[1]*a[2]';
strings{4}='a[1]*a[2]*a[3]*a[4]+a[2]*a[3]*a[4]*a[5]+a[3]*a[4]*a[5]*a[1]+a[4]*a[5]*a[1]*a[2]+a[5]*a[1]*a[2]*a[3]';
strings{5}='a[1]*a[2]*a[3]*a[4]*a[5]-1';
polysys=lti2polysys(strings,[],[],'na',5);

univarstring={'a[1]^15+122*a[1]^10-122*a[1]^5-1'};
p1sys=lti2polysys(univarstring,[],[],'na',1);

% noon4
% http://homepages.math.uic.edu/~jan/Demo/noon4.html
strings{1}='a[1]*a[2]^2+a[1]*a[3]^2+a[1]*a[4]^2-1.1*a[1]+1';
strings{2}='a[2]*a[1]^2+a[2]*a[3]^2+a[2]*a[4]^2-1.1*a[2]+1';
strings{3}='a[3]*a[1]^2+a[3]*a[2]^2+a[3]*a[4]^2-1.1*a[3]+1';
strings{4}='a[4]*a[1]^2+a[4]*a[2]^2+a[4]*a[3]^2-1.1*a[4]+1';
polysys=lti2polysys(strings,[],[],'na',4);

% noon3
% http://homepages.math.uic.edu/~jan/Demo/noon3.html
strings{1}='a[1]*a[2]^2+a[1]*a[3]^2-1.1*a[1]+1';
strings{2}='a[2]*a[1]^2+a[2]*a[3]^2-1.1*a[2]+1';
strings{3}='a[3]*a[1]^2+a[3]*a[2]^2-1.1*a[3]+1';
polysys=lti2polysys(strings,[],[],'na',3);


% eco6
% http://homepages.math.uic.edu/~jan/Demo/eco6.html
strings{1}='a[1]*a[6]+a[2]*a[1]*a[6]+a[6]*a[2]*a[3]+a[3]*a[6]*a[4]+a[6]*a[5]*a[4]-1';
strings{2}='a[2]*a[6]+a[3]*a[1]*a[6]+a[2]*a[6]*a[4]+a[3]*a[6]*a[5]-2';
strings{3}='a[3]*a[6]+a[1]*a[6]*a[4]+a[2]*a[6]*a[5]-3';
strings{4}='a[4]*a[6]+a[1]*a[6]*a[5]-4';
strings{5}='a[5]*a[6]-5';
strings{6}='a[1]+a[2]+a[3]+a[4]+a[5]+1';
polysys=lti2polysys(strings,[],[],'na',6);

% puma
strings{1}='a[1]^2+a[2]^2-1';
strings{2}='a[3]^2+a[4]^2-1';
strings{3}='a[5]^2+a[6]^2-1';
strings{4}='a[7]^2+a[8]^2-1';
strings{5}='0.004731*a[1]*a[3]-0.3578*a[2]*a[3]-0.1238*a[1]-0.001637*a[2]-0.9338*a[4]+a[7]-0.3571';
strings{6}='0.2238*a[1]*a[3]+0.7623*a[2]*a[3]+0.2638*a[1]-0.07745*a[2]-0.6734*a[4]-0.6022';
strings{7}='a[6]*a[8]+0.3578*a[1]+0.004731*a[2]';
strings{8}='-0.7623*a[1]+0.2238*a[2]+0.3461';
polysys=lti2polysys(strings,[],[],'na',8);

% cyclic 10
strings{1}='a[1]+a[2]+a[3]+a[4]+a[5]+a[6]+a[7]+a[8]+a[9]+a[10]';
strings{2}='a[1]*a[2]+a[2]*a[3]+a[3]*a[4]+a[4]*a[5]+a[5]*a[6]+a[6]*a[7]+a[7]*a[8]+a[8]*a[9]+a[9]*a[10]+a[10]*a[1]';
strings{3}='a[1]*a[2]*a[3]+a[2]*a[3]*a[4]+a[3]*a[4]*a[5]+a[4]*a[5]*a[6]+a[5]*a[6]*a[7]+a[6]*a[7]*a[8]+a[7]*a[8]*a[9]+a[8]*a[9]*a[10]+a[9]*a[10]*a[1]+a[10]*a[1]*a[2]';
strings{4}='a[1]*a[2]*a[3]*a[4]+a[2]*a[3]*a[4]*a[5]+a[3]*a[4]*a[5]*a[6]+a[4]*a[5]*a[6]*a[7]+a[5]*a[6]*a[7]*a[8]+a[6]*a[7]*a[8]*a[9]+a[7]*a[8]*a[9]*a[10]+a[8]*a[9]*a[10]*a[1]+a[9]*a[10]*a[1]*a[2]+a[10]*a[1]*a[2]*a[3]';
strings{5}='a[1]*a[2]*a[3]*a[4]*a[5]+a[2]*a[3]*a[4]*a[5]*a[6]+a[3]*a[4]*a[5]*a[6]*a[7]+a[4]*a[5]*a[6]*a[7]*a[8]+a[5]*a[6]*a[7]*a[8]*a[9]+a[6]*a[7]*a[8]*a[9]*a[10]+a[7]*a[8]*a[9]*a[10]*a[1]+a[8]*a[9]*a[10]*a[1]*a[2]+a[9]*a[10]*a[1]*a[2]*a[3]+a[10]*a[1]*a[2]*a[3]*a[4]';
strings{6}='a[1]*a[2]*a[3]*a[4]*a[5]*a[6]+a[2]*a[3]*a[4]*a[5]*a[6]*a[7]+a[3]*a[4]*a[5]*a[6]*a[7]*a[8]+a[4]*a[5]*a[6]*a[7]*a[8]*a[9]+a[5]*a[6]*a[7]*a[8]*a[9]*a[10]+a[6]*a[7]*a[8]*a[9]*a[10]*a[1]+a[7]*a[8]*a[9]*a[10]*a[1]*a[2]+a[8]*a[9]*a[10]*a[1]*a[2]*a[3]+a[9]*a[10]*a[1]*a[2]*a[3]*a[4]+a[10]*a[1]*a[2]*a[3]*a[4]*a[5]';
strings{7}='a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]+a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]+a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]+a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]+a[5]*a[6]*a[7]*a[8]*a[9]*a[10]*a[1]+a[6]*a[7]*a[8]*a[9]*a[10]*a[1]*a[2]+a[7]*a[8]*a[9]*a[10]*a[1]*a[2]*a[3]+a[8]*a[9]*a[10]*a[1]*a[2]*a[3]*a[4]+a[9]*a[10]*a[1]*a[2]*a[3]*a[4]*a[5]+a[10]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]';
strings{8}='a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]+a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]+a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]+a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]*a[1]+a[5]*a[6]*a[7]*a[8]*a[9]*a[10]*a[1]*a[2]+a[6]*a[7]*a[8]*a[9]*a[10]*a[1]*a[2]*a[3]+a[7]*a[8]*a[9]*a[10]*a[1]*a[2]*a[3]*a[4]+a[8]*a[9]*a[10]*a[1]*a[2]*a[3]*a[4]*a[5]+a[9]*a[10]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]+a[10]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]';
strings{9}='a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]+a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]+a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]*a[1]+a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]*a[1]*a[2]+a[5]*a[6]*a[7]*a[8]*a[9]*a[10]*a[1]*a[2]*a[3]+a[6]*a[7]*a[8]*a[9]*a[10]*a[1]*a[2]*a[3]*a[4]+a[7]*a[8]*a[9]*a[10]*a[1]*a[2]*a[3]*a[4]*a[5]+a[8]*a[9]*a[10]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]+a[9]*a[10]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]+a[10]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]';
strings{10}='a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7]*a[8]*a[9]*a[10]-1';
polysys=lti2polysys(strings,[],[],'na',10);


% moment matrices, trace matrices and the radical of ideals: (-1,3)^3 and
% (2,2)^2
% K=[d{00}(-1,3) d{10}(-1,3) (2*d{20}+d{01})(-1,3) d{00}(2,2) d{10}(2,2)]
polysys{1,1}=[3,18,-48,21,-114,156];polysys{1,2}=[2 0;1 1;1 0;0 2;0 1;0 0];
polysys{2,1}=[1,-259/4,493/4,-611/4,2423/4,-1175/2,-5,6,1,5];polysys{2,2}=[3 0;2 1;2 0;1 2;1 1;1 0;0 3;0 2;0 1;0 0];
polysys{3,1}=[1,81/4,-163/4,21/4,87/4,-151/2,-1,4,2,3];polysys{3,2}=[3 0;2 1;2 0;1 2;1 1;1 0;0 3;0 2;0 1;0 0];

% approximate radical for clusters: a global approach using Gaussian
% Elimination or SVD; (1,1)^3 and (-1,2)^2: cannot be solved using matrix
% methods!
polysys{1,1}=[1 4 -6 6 -18 13];polysys{1,2}=[2 0;1 1;1 0;0 2;0 1;0 0];
polysys{2,1}=[1 16 -7 118 -286 147 -1 6 1 5];polysys{2,2}=[3 0;2 1;2 0;1 2;1 1;1 0;0 3;0 2;0 1;0 0];
polysys{3,1}=[1 10 -5 72 -176 91 -1 4 1 3];polysys{3,2}=[3 0;2 1;2 0;1 2;1 1;1 0;0 3;0 2;0 1;0 0];
% matrix of traces
R=[5 1 7 -1 5;1 5 -1 7 1;7 -1 11 -5 7;-1 7 -5 11 -1;5 1 7 -1 5];

% groebner basis of above polynomial system: can be solved perfectly
strings{1}='a[1]^2+4*a[1]*a[2]-6*a[1]+6*a[2]^2-18*a[2]+13';
strings{2}='3*a[1]^3-a[1]^2-4*a[1]*a[2]-3*a[1]-12*a[2]+17';
strings{3}='3*a[1]^2*a[2]-5-2*a[1]*a[2]-5*a[1]^2+3*a[2]+6*a[1]';
polysys=lti2polysys(strings,[],[],'na',2);

polysys{1,1}=[8 1 3];polysys{1,2}=[7 0 0 0;0 1 2 2;0 1 0 0];
polysys{2,1}=[8 1 3 1];polysys{2,2}=[0 7 0 0;1 0 2 2;1 0 0 0;0 0 1 0];
polysys{3,1}=[8 2 1 1 2];polysys{3,2}=[0 0 7 0;1 1 1 2;0 0 0 2;0 1 0 0;0 0 0 1];
polysys{4,1}=[8 2 2 2 1];polysys{4,2}=[0 0 0 7;1 1 2 1;0 0 1 1;0 0 1 0;0 0 0 0];

% toy problem
strings{1}='a[1]*a[7]';
strings{2}='a[1]*a[8]+a[2]*a[7]-4*a[8]+a[9]';
strings{3}='a[1]*a[9]+a[3]*a[7]-4*a[9]-a[8]';
strings{4}='a[1]*a[10]+a[2]*a[8]+a[4]*a[7]-8*a[10]+a[11]';
strings{5}='a[1]*a[11]+a[2]*a[9]+a[3]*a[8]+a[5]*a[7]-8*a[11]-2*a[10]+2*a[12]';
strings{6}='a[1]*a[12]+a[3]*a[9]+a[6]*a[7]-8*a[12]-a[11]';
strings{7}='a[2]*a[10]+a[4]*a[8]+a[8]';
strings{8}='a[2]*a[11]+a[3]*a[10]+a[4]*a[9]+a[5]*a[8]+a[9]';
strings{9}='a[2]*a[12]+a[3]*a[11]+a[5]*a[9]+a[6]*a[8]+a[8]';
strings{10}='a[3]*a[12]+a[6]*a[9]+a[9]';
strings{11}='a[4]*a[10]+2*a[10]';
strings{12}='a[4]*a[11]+a[5]*a[10]+2*a[11]';
strings{13}='a[4]*a[12]+a[5]*a[11]+a[6]*a[10]+2*a[10]+2*a[12]';
strings{14}='a[5]*a[12]+a[6]*a[11]+2*a[11]';
strings{15}='a[6]*a[12]+2*a[12]';
polysys=lti2polysys(strings,[],[],'na',12);

% toy problem 2
strings{1}='a[6]-4*a[1]+a[2]';
strings{2}='a[7]-a[1]-4*a[2]';
strings{3}='a[4]-8*a[3]+a[6]*a[1]-2';
strings{4}='2*a[5]-8*a[4]-2*a[3]+a[6]*a[2]+a[7]*a[1]';
strings{5}='a[7]*a[2]-8*a[5]-a[4]-2';
strings{6}='a[6]*a[3]-a[1]';
strings{7}='a[6]*a[4]-a[2]+a[7]*a[3]';
strings{8}='a[6]*a[5]-a[1]+a[7]*a[4]';
strings{9}='a[7]*a[5]-a[2]';
polysys=lti2polysys(strings,[],[],'na',7);

% toy problem deg(p)=1
strings{1}='a[4]';
strings{2}='a[4]+a[5]-3';
strings{3}='a[4]+a[6]-5';
strings{4}='a[5]+a[7]';
strings{5}='a[5]+a[6]+a[8]';
strings{6}='a[6]+a[9]';
strings{7}='a[7]+1';
strings{8}='a[7]+a[8]+1';
strings{9}='a[8]+a[9]+1';
strings{10}='a[9]+1';
polysys=lti2polysys(strings,[],[],'na',9);

% simple dynamical system
dynsys{1,1}=[1 -1 -1 4];dynsys{1,2}=[0 1;3 0;1 2;1 0];
dynsys{2,1}=[-1 -1 -1 4];dynsys{2,2}=[1 0;2 1;0 3;0 1];

% simple dynamical system, different degrees
dynsys{1,1}=[1 -1 -1 4];dynsys{1,2}=[0 1;3 0;1 2;1 0];
dynsys{2,1}=[-1 -1 -1 4];dynsys{2,2}=[1 0;2 1;0 4;0 1];

% Wien Bridge oscillator
strings{1}='3.334*a[4]^3*a[2]*a[1]+10.002*a[4]*a[2]^3*a[1]+3.334*a[4]*a[2]*a[3]^2*a[1]-6.585*a[4]*a[2]*a[1]-1.0*a[1]^2+1.0';
strings{2}='13.336*a[4]^4+6.668*a[4]^2*a[2]^2+3.334*a[4]^2*a[2]*a[3]+26.672*a[4]^2*a[3]^2-13.17*a[4]^2+6.668*a[2]^4+16.67*a[2]^3*a[3]+6.668*a[2]^2*a[3]^2-6.585*a[2]^2+3.334*a[2]*a[3]^3-6.585*a[2]*a[3]+13.336*a[3]^4-13.17*a[3]^2+0.234';
strings{3}='10.002*a[4]^5*a[1]+40.008*a[4]^3*a[2]^2*a[1]+20.004*a[4]^3*a[3]^2*a[1]-19.755*a[4]^3*a[1]+43.342*a[4]*a[2]^4*a[1]+6.668*a[4]*a[2]^3*a[3]*a[1]+40.008*a[4]*a[2]^2*a[3]^2*a[1]-39.51*a[4]*a[2]^2*a[1]+10.002*a[4]*a[3]^4*a[1]-19.755*a[4]*a[3]^2*a[1]+0.702*a[4]*a[1]-9.0*a[3]*a[1]^2+a[3]';
strings{4}='10.002*a[4]^4*a[3]*a[1]+3.334*a[4]^2*a[2]^3*a[1]+40.008*a[4]^2*a[2]^2*a[3]*a[1]+20.004*a[4]^2*a[3]^3*a[1]-19.755*a[4]^2*a[3]*a[1]+9.0*a[4]*a[1]^2-1.0*a[4]+10.002*a[2]^5*a[1]+43.342*a[2]^4*a[3]*a[1]+10.002*a[2]^3*a[3]^2*a[1]-6.585*a[2]^3*a[1]+40.008*a[2]^2*a[3]^3*a[1]-39.51*a[2]^2*a[3]*a[1]+10.002*a[3]^5*a[1]-19.755*a[3]^3*a[1]+0.702*a[3]*a[1]';
polysys=lti2polysys(strings,[],[],'na',4);

% Van der Pol oscillator
% [omega c1 c3i cr3]
strings{1}='-a[1]^2-3*a[2]*a[3]*a[1]+1';
strings{2}='6*a[4]^2+3*a[2]*a[4]+6*a[3]^2+3*a[2]^2-1';
strings{3}='-9*a[4]*a[1]^2-9*a[3]*a[4]^2*a[1]-9*a[3]^3*a[1]-18*a[2]^2*a[3]*a[1]+3*a[3]*a[1]+a[4]';
strings{4}='9*a[4]^3*a[1]-9*a[3]*a[1]^2+9*a[3]^2*a[4]*a[1]+18*a[2]^2*a[4]*a[1]-3*a[4]*a[1]+3*a[2]^3*a[1]+a[3]';
polysys=lti2polysys(strings,[],[],'na',4);

% Van der Pol oscillator
% Chebyshev basis, N=3
strings{1}='';
strings{2}='';
strings{3}='';
strings{4}='';
polysys=lti2polysys(strings,[],[],'na',4);

% Van der Pol oscillator
% Chebyshev basis, N=5
strings{1}='576*a[1]*a[3]*a[5]+288*a[2]*a[3]*a[4]-1152*a[3]^2*a[5]+10368*a[3]*a[4]*a[6]-2304*a[1]*a[5]^2-4608*a[1]*a[6]*a[4]-17280*a[3]*a[6]^2+288*a[1]*a[4]^2-2016*a[2]*a[3]*a[6]+7776*a[2]*a[5]*a[6]+576*a[1]*a[6]*a[2]+5184*a[3]*a[5]^2+16*a[6]+144*a[2]^2*a[5]-1152*a[3]*a[4]^2-31968*a[4]*a[5]*a[6]+48*a[3]^3-5376*a[5]^3+10080*a[1]*a[6]^2+5040*a[4]^2*a[5]-2016*a[2]*a[4]*a[5]+42480*a[5]*a[6]^2';
strings{2}='540*a[4]^3+60*a[2]*a[3]^2+1200*a[2]*a[5]^2-420*a[3]^2*a[4]-4560*a[4]*a[5]^2+60*a[2]^2*a[4]-360*a[2]*a[4]^2+240*a[1]*a[2]*a[5]+240*a[1]*a[3]*a[4]-1680*a[1]*a[4]*a[5]-720*a[2]*a[3]*a[5]+3360*a[3]*a[4]*a[5]+8*a[5]+240*a[1]^2*a[6]+1740*a[3]^2*a[6]+11040*a[5]^2*a[6]-3000*a[6]^2*a[2]+10500*a[6]^2*a[4]-300*a[2]^2*a[6]-4500*a[4]^2*a[6]+2400*a[2]*a[6]*a[4]+6480*a[1]*a[5]*a[6]-1680*a[1]*a[3]*a[6]-10080*a[3]*a[5]*a[6]-7500*a[6]^3-80*a[6]';
strings{3}='12*a[3]^3-5856*a[5]^3+96*a[1]^2*a[5]+48*a[1]*a[3]^2-1920*a[1]*a[5]^2-912*a[3]^2*a[5]+5136*a[3]*a[5]^2+144*a[1]*a[3]*a[5]+96*a[1]*a[2]*a[4]+120*a[2]*a[3]*a[4]-1848*a[2]*a[4]*a[5]+4*a[4]-32*a[5]+24*a[2]^2*a[3]+84*a[2]^2*a[5]+72*a[1]*a[4]^2-936*a[3]*a[4]^2+5148*a[4]^2*a[5]-18600*a[3]*a[6]^2+48300*a[5]*a[6]^2+10200*a[1]*a[6]^2+240*a[1]*a[6]*a[2]-3840*a[1]*a[6]*a[4]-35160*a[4]*a[5]*a[6]+8280*a[2]*a[5]*a[6]+10320*a[3]*a[4]*a[6]-1800*a[2]*a[3]*a[6]+320*a[6]';
strings{4}='3*a[2]^3+459*a[4]^3+36*a[1]^2*a[4]+24*a[2]*a[3]^2+1056*a[2]*a[5]^2-276*a[3]^2*a[4]-4092*a[4]*a[5]^2+33*a[2]^2*a[4]-279*a[2]*a[4]^2+36*a[1]*a[2]*a[3]+96*a[1]*a[2]*a[5]+60*a[1]*a[3]*a[4]-1176*a[1]*a[4]*a[5]-540*a[2]*a[3]*a[5]+2748*a[3]*a[4]*a[5]-12*a[4]+96*a[5]+2*a[3]+60*a[1]^2*a[6]+1380*a[3]^2*a[6]+10140*a[5]^2*a[6]-2775*a[6]^2*a[2]+9825*a[6]^2*a[4]-255*a[2]^2*a[6]-4095*a[4]^2*a[6]+2130*a[2]*a[6]*a[4]+5400*a[1]*a[5]*a[6]-1140*a[1]*a[3]*a[6]-8820*a[3]*a[5]*a[6]-7125*a[6]^3-20*a[6]';
strings{5}='6*a[3]^3-2760*a[5]^3+12*a[1]^2*a[3]+24*a[1]^2*a[5]+12*a[1]*a[3]^2-816*a[1]*a[5]^2-396*a[3]^2*a[5]+2340*a[3]*a[5]^2+48*a[1]*a[3]*a[5]+36*a[1]*a[2]*a[4]+36*a[2]*a[3]*a[4]-792*a[2]*a[4]*a[5]+24*a[4]-8*a[5]+a[2]-4*a[3]+6*a[1]*a[2]^2+12*a[2]^2*a[3]+24*a[2]^2*a[5]+18*a[1]*a[4]^2-396*a[3]*a[4]^2+2340*a[4]^2*a[5]-8700*a[3]*a[6]^2+23100*a[5]*a[6]^2+4650*a[1]*a[6]^2+60*a[1]*a[6]*a[2]-1620*a[1]*a[6]*a[4]-16560*a[4]*a[5]*a[6]+3840*a[2]*a[5]*a[6]+4680*a[3]*a[4]*a[6]-780*a[2]*a[3]*a[6]+120*a[6]';
strings{6}='3*a[2]^3+324*a[4]^3+6*a[1]^2*a[2]+18*a[1]^2*a[4]+15*a[2]*a[3]^2+762*a[2]*a[5]^2-189*a[3]^2*a[4]-2970*a[4]*a[5]^2+18*a[2]^2*a[4]-189*a[2]*a[4]^2+24*a[1]*a[2]*a[3]+48*a[1]*a[2]*a[5]+36*a[1]*a[3]*a[4]-792*a[1]*a[4]*a[5]-372*a[2]*a[3]*a[5]+1944*a[3]*a[4]*a[5]-6*a[4]+64*a[5]+2*a[1]-2*a[2]+8*a[3]+30*a[1]^2*a[6]+975*a[3]^2*a[6]+7410*a[5]^2*a[6]-2025*a[6]^2*a[2]+7200*a[6]^2*a[4]-180*a[2]^2*a[6]-2970*a[4]^2*a[6]+1530*a[2]*a[6]*a[4]+3840*a[1]*a[5]*a[6]-780*a[1]*a[3]*a[6]-6360*a[3]*a[5]*a[6]-5250*a[6]^3-10*a[6]';
polysys=lti2polysys(strings,[],[],'na',6);

% Van der Pol oscillator
% Haotian's equations, N=3, has 2D variety
strings{1}='-9*a[4]^3*a[1]-18*a[4]*a[2]^2*a[1]-9*a[4]*a[3]^2*a[1]+3*a[4]*a[1]-9*a[3]*a[1]^2+a[3]';
strings{2}='-3*a[4]*a[2]^2*a[1]-a[2]*a[1]^2+a[2]';
strings{3}='-9*a[4]^2*a[3]*a[1]+9*a[4]*a[1]^2-a[4]-3*a[2]^3*a[1]-18*a[2]^2*a[3]*a[1]-9*a[3]^3*a[1]+3*a[3]*a[1]';
strings{4}='-6*a[1]*a[4]^2*a[2]-3*a[1]*a[2]^3-3*a[1]*a[2]^2*a[3]-6*a[1]*a[2]*a[3]^2+a[1]*a[2]';
polysys=lti2polysys(strings,[],[],'na',4);

% Van der Pol oscillator
% Haotian's equations, N=3, simplified
strings{1}='-9*a[4]^3*a[1]-18*a[4]*a[2]^2*a[1]-9*a[4]*a[3]^2*a[1]+3*a[4]*a[1]-9*a[3]*a[1]^2+a[3]';
strings{2}='3*a[1]*a[2]*a[4]+a[1]^2-1';
strings{3}='-9*a[4]^2*a[3]*a[1]+9*a[4]*a[1]^2-a[4]-3*a[2]^3*a[1]-18*a[2]^2*a[3]*a[1]-9*a[3]^3*a[1]+3*a[3]*a[1]';
strings{4}='6*a[4]^2+3*a[2]^2+3*a[2]*a[3]+6*a[3]^2-1';
polysys=lti2polysys(strings,[],[],'na',4);

% Van der Pol oscillator
% Haotian's equations, N=5, simplified
strings{1}='-30*a[1]*a[4]^2*a[6]-15*a[1]*a[2]^2*a[4]-30*a[1]*a[2]*a[3]*a[4]-15*a[1]*a[6]^3-30*a[1]*a[2]^2*a[6]-30*a[1]*a[3]^2*a[6]-15*a[1]*a[5]^2*a[6]+5*a[1]*a[6]-25*a[1]^2*a[5]+a[5]';
strings{2}='-9*a[1]*a[4]^3-18*a[1]*a[4]*a[6]^2-18*a[1]*a[2]^2*a[4]+18*a[1]*a[2]*a[4]*a[5]-9*a[1]*a[3]^2*a[4]-18*a[1]*a[4]*a[5]^2+3*a[1]*a[4]-9*a[1]*a[2]^2*a[6]-18*a[1]*a[2]*a[3]*a[6]-9*a[1]^2*a[3]+a[3]';
strings{3}='-3*a[1]*a[4]^2*a[6]-3*a[1]*a[2]^2*a[4]+6*a[1]*a[2]*a[4]*a[5]-6*a[1]*a[3]*a[4]*a[5]-6*a[1]*a[2]*a[3]*a[6]-a[1]^2*a[2]+a[2]+3*a[1]*a[3]^2*a[6]';
strings{4}='15*a[1]*a[2]*a[4]^2-30*a[1]*a[4]^2*a[5]-15*a[1]*a[5]*a[6]^2+25*a[1]^2*a[6]-a[6]-15*a[1]*a[2]^2*a[3]-30*a[1]*a[2]^2*a[5]-15*a[1]*a[2]*a[3]^2-30*a[1]*a[3]^2*a[5]-15*a[1]*a[5]^3+5*a[1]*a[5]';
strings{5}='-9*a[1]*a[3]*a[4]^2-18*a[1]*a[2]*a[4]*a[6]+9*a[1]^2*a[4]-a[4]-18*a[1]*a[3]*a[6]^2-3*a[1]*a[2]^3-18*a[1]*a[2]^2*a[3]-9*a[1]*a[2]^2*a[5]-18*a[1]*a[2]*a[3]*a[5]-9*a[1]*a[3]^3-18*a[1]*a[3]*a[5]^2+3*a[1]*a[3]';
strings{6}='6*a[2]*a[4]^2-3*a[4]^2*a[5]+6*a[2]*a[4]*a[6]+6*a[3]*a[4]*a[6]+6*a[2]*a[6]^2+3*a[2]^3+3*a[2]^2*a[3]+6*a[2]*a[3]^2+6*a[2]*a[3]*a[5]+6*a[2]*a[5]^2-a[2]+3*a[3]^2*a[5]';
polysys=lti2polysys(strings,[],[],'na',6);

% Van der Pol oscillator
% Haotian's equations, N=5, simplified, 0-dim affine solution set
% [omega cr1 cr3 ci3 cr5 ci5 z]
strings{1}='-30*a[1]*a[4]^2*a[6]-15*a[1]*a[2]^2*a[4]-30*a[1]*a[2]*a[3]*a[4]-15*a[1]*a[6]^3-30*a[1]*a[2]^2*a[6]-30*a[1]*a[3]^2*a[6]-15*a[1]*a[5]^2*a[6]+5*a[1]*a[6]-25*a[1]^2*a[5]+a[5]';
strings{2}='-9*a[1]*a[4]^3-18*a[1]*a[4]*a[6]^2-18*a[1]*a[2]^2*a[4]+18*a[1]*a[2]*a[4]*a[5]-9*a[1]*a[3]^2*a[4]-18*a[1]*a[4]*a[5]^2+3*a[1]*a[4]-9*a[1]*a[2]^2*a[6]-18*a[1]*a[2]*a[3]*a[6]-9*a[1]^2*a[3]+a[3]';
strings{3}='-3*a[1]*a[4]^2*a[6]-3*a[1]*a[2]^2*a[4]+6*a[1]*a[2]*a[4]*a[5]-6*a[1]*a[3]*a[4]*a[5]-6*a[1]*a[2]*a[3]*a[6]-a[1]^2*a[2]+a[2]+3*a[1]*a[3]^2*a[6]';
strings{4}='15*a[1]*a[2]*a[4]^2-30*a[1]*a[4]^2*a[5]-15*a[1]*a[5]*a[6]^2+25*a[1]^2*a[6]-a[6]-15*a[1]*a[2]^2*a[3]-30*a[1]*a[2]^2*a[5]-15*a[1]*a[2]*a[3]^2-30*a[1]*a[3]^2*a[5]-15*a[1]*a[5]^3+5*a[1]*a[5]';
strings{5}='-9*a[1]*a[3]*a[4]^2-18*a[1]*a[2]*a[4]*a[6]+9*a[1]^2*a[4]-a[4]-18*a[1]*a[3]*a[6]^2-3*a[1]*a[2]^3-18*a[1]*a[2]^2*a[3]-9*a[1]*a[2]^2*a[5]-18*a[1]*a[2]*a[3]*a[5]-9*a[1]*a[3]^3-18*a[1]*a[3]*a[5]^2+3*a[1]*a[3]';
strings{6}='6*a[2]*a[4]^2-3*a[4]^2*a[5]+6*a[2]*a[4]*a[6]+6*a[3]*a[4]*a[6]+6*a[2]*a[6]^2+3*a[2]^3+3*a[2]^2*a[3]+6*a[2]*a[3]^2+6*a[2]*a[3]*a[5]+6*a[2]*a[5]^2-a[2]+3*a[3]^2*a[5]';
strings{7}='a[2]*a[7]-1';
polysys=lti2polysys(strings,[],[],'na',7);

% Van der Pol, N=5
% additional pure powers added
% [omega cr1 cr3 ci3 cr5 ci5 z]
strings{1}='-30*a[1]*a[4]^2*a[6]-15*a[1]*a[2]^2*a[4]-30*a[1]*a[2]*a[3]*a[4]-15*a[1]*a[6]^3-30*a[1]*a[2]^2*a[6]-30*a[1]*a[3]^2*a[6]-15*a[1]*a[5]^2*a[6]+5*a[1]*a[6]-25*a[1]^2*a[5]+a[5]';
strings{2}='-9*a[1]*a[4]^3-18*a[1]*a[4]*a[6]^2-18*a[1]*a[2]^2*a[4]+18*a[1]*a[2]*a[4]*a[5]-9*a[1]*a[3]^2*a[4]-18*a[1]*a[4]*a[5]^2+3*a[1]*a[4]-9*a[1]*a[2]^2*a[6]-18*a[1]*a[2]*a[3]*a[6]-9*a[1]^2*a[3]+a[3]';
strings{3}='-3*a[1]*a[4]^2*a[6]-3*a[1]*a[2]^2*a[4]+6*a[1]*a[2]*a[4]*a[5]-6*a[1]*a[3]*a[4]*a[5]-6*a[1]*a[2]*a[3]*a[6]-a[1]^2*a[2]+a[2]+3*a[1]*a[3]^2*a[6]';
strings{4}='15*a[1]*a[2]*a[4]^2-30*a[1]*a[4]^2*a[5]-15*a[1]*a[5]*a[6]^2+25*a[1]^2*a[6]-a[6]-15*a[1]*a[2]^2*a[3]-30*a[1]*a[2]^2*a[5]-15*a[1]*a[2]*a[3]^2-30*a[1]*a[3]^2*a[5]-15*a[1]*a[5]^3+5*a[1]*a[5]';
strings{5}='-9*a[1]*a[3]*a[4]^2-18*a[1]*a[2]*a[4]*a[6]+9*a[1]^2*a[4]-a[4]-18*a[1]*a[3]*a[6]^2-3*a[1]*a[2]^3-18*a[1]*a[2]^2*a[3]-9*a[1]*a[2]^2*a[5]-18*a[1]*a[2]*a[3]*a[5]-9*a[1]*a[3]^3-18*a[1]*a[3]*a[5]^2+3*a[1]*a[3]';
strings{6}='6*a[2]*a[4]^2-3*a[4]^2*a[5]+6*a[2]*a[4]*a[6]+6*a[3]*a[4]*a[6]+6*a[2]*a[6]^2+3*a[2]^3+3*a[2]^2*a[3]+6*a[2]*a[3]^2+6*a[2]*a[3]*a[5]+6*a[2]*a[5]^2-a[2]+3*a[3]^2*a[5]';
strings{7}='a[2]*a[7]-1';
strings{8}='365709686019733887557659344611377654194124181137740190953517251343960003454823499304254133993222058389966004757408632223639171770288*a[7]-5012764663565731524354295846079126387093690662717827467714627057849326286960844967319278088890023944514160651033438659749332031653624*a[1]^2*a[7]+7689668318518853666056310840671125268972715095310981001493159783219424105888941860245180250579522301975880299006708524723553233775440*a[4]*a[6]*a[3]+21719671128167137155495616589349841096777848307999027998495554368741355264719282557243448122663233160281279233324805376631855236277056*a[2]*a[3]*a[5]+19939789058713865162598471374266016066921620260966576538096877734369483009483018852819083012617120465313483697106827908984594712130766*a[4]*a[6]*a[2]+36769957880759680904944303906219980442790476030914797433284345306498536420859844512536677900454034125804998903014754157166073149760*a[4]*a[6]*a[5]-9811503350460026996525352273323230940749546024572335357793375254310670790458305576567777663683383731783130178268844508225526119537522*a[2]-3273199154323968375133644504715064332115321050146391514844230554199159658075806284195346852482401274930832634799658610978627535339440*a[3]+1871567118857543529560932155080392759279762598141989528285321960887392366533409730622827885672239922643730318323949941157453840189253*a[5]+12715082202042404333730087410843775811339895401662032846354562837960655921728314906887911666075923967359115114323410301624252588566590*a[2]*a[1]^2+6657469457836479091911762718026985844167615127241888768340692087986297564809882806040485265869518816207186766649380201919187778366321*a[6]*a[1]-9480620695348464714577861049479752824395437238287167604758130907278709921011723890860777041974893448805677995746523301523835665988605*a[5]*a[1]^2+7855764963060713297954420840765405455138671431703822278415167599677383572445229312519905906302137303092255922812026815594559463108742*a[4]*a[1]+15376635072246500116812681167961758482840792356710005003177980347462524122225125803713274813980610038338030891731923296298575477504520*a[3]*a[1]^2-19449697643446907752811658196394528744357579231609270537644320923090510401368627549002039169748182823245704852682171773676996496450188*a[2]^2*a[5]+1520312773346429994441697399579839424367856445247677234801664989625396564777718704294632094186329747279960904204926296128772531431712*a[4]^2*a[7]+476482254167411683715397039312431597059958318078327280870366361193554970683511154671458586438078778944909735711628651324272921112816*a[4]*a[6]*a[7]+43747405878566340367270156945649140991348561547564509400083683335772267038727490631363690592708284694922369204759860704809750133598440*a[4]^2*a[2]-7642575730758553823450112831303763981929473616911105091323805857933604626310029660058533306659432196022534277317566964994643686872384*a[4]^2*a[5]+66511843173185657046460477960474178140129994234981350667715334214221043055811840156539902249130040195692411218895591608447075994133320*a[6]^2*a[2]+6937008729061313741217184159292840816917857472633406543783028154083892926025288896441696379265733084846923706068692844549738724795144*a[2]^2*a[3]+45559820524205570808928073628881587376810590069571742995886980150733122066983565763330749048948903417017801815002242509577609614636900*a[2]*a[3]^2+64505196806612993053604904553879328700966995231289483471762344140470785117953841259031035479600631627905819237813409861848863229430984*a[2]*a[5]^2+29385201087283681656138973366909432857901571216691801230142049949052229320174711532868748873477007744772626086820237420012996732946496*a[2]^3+2327712302185496520872594124730132541201899855653121595030509789745619704081234727379806012908843141073045023022805251110676583716180*a[3]^3+2914655281004700350312187938477619798988119642787598331345319657869252995210573321315252527592258371655742929503779852414976592879772*a[4]^2*a[3]';
strings{9}='-233751906421527571130085891238538894190594192512772915380165127280974600558387787194001893604036460742224910952894453852487290543427599084520878063795099128*a[1]*a[7]+10220610001733227931539882716887742077564517602105893881499153788035972691049985503039174902477694648844265370566553335512805021392803039503508182776455415181*a[4]^2*a[6]-39816031177152590027751459710400763002100437106921011840191563360797494637479401816424286162641488177104786644625827557037907755724266453311960864892141937774*a[4]*a[6]^2-17950995002897925581192423424617083769886129181889382692837389597901331251748854386059569301666131097546383183725387802301064074825564618022326397799428624638*a[4]*a[3]*a[5]-27349383632756199586020926535706421930531388106920815255612509645993525066204844333752884364778929739331250101997844686462210069817944979212288728193094543778*a[4]*a[5]^2+6480726125870785198080925061764433203224692150224002016360585331801506567718321763832452361258620782804209977795060078808734301969880281062750678581123553624*a[4]-3641234110759816933122133759097503287817827720742649819035159782512036673479806510224911289573881494367664308468709159406683611680808268400093178389328625540*a[6]+630914962118206033994983656170761300662654106638289836835284764602772420044370037484285354711906828382034994628957257177918952653182575538233347649651999304*a[2]*a[1]-2021701162065951229397770764795907034019557800466532832347830367856776427710244281232497323490282845232053027197938935410513587106589570117365460962433936160*a[6]*a[1]^2+24259406521931067797366232765375402559327399906350198365713701035644275847954390127123874237055899483050173758040794456694527022397417714871832589743173875424*a[5]*a[1]-243472385669905269255172825522316305255227791924261015570594808321450192700551580643045925071964647432230813517084287678117641767577385550644361082266389680*a[4]*a[1]^2-7226147073679488465817094878642357841371104687634497372066031329570643792147447019891188046875975111448431985973452201121314232811912943109112456071294041664*a[3]*a[1]-9892235919735649466071509607318545088285897940667889375994231933143911033624016796859821942109618294395721681756390032361776388368032911252247551198765729783*a[4]*a[2]^2+25397604571499245121528705792600014951697138850375803976462937853451999879426896241736802965383608555649012025436781233707465122810300740102228110125668352960*a[4]*a[2]*a[3]+22063757696942183479777851994631975241013770096914535791870180458469568479722439317069621636343449509670036755274296904143985815737819548523638661585659522218*a[6]*a[2]^2-8218831281451964074777230769378863493850937681794393415108241105994055617785574098072599253423657106500800426378140359312494920081128642996863793486067240427*a[6]*a[2]*a[3]+2226193531588354921928548261857937641606495803210625814170252916601526153482192765953042318291171808572549352709806979149300627636100495151725290654166139827*a[2]*a[4]*a[5]+1850541455196761684776793048780866816896526857654121012618042086004091405098422247998533670663843773344855819746797801991072345083454340015755918064037021040*a[6]*a[2]*a[5]+626752454483932912297189136481361868757438573724181872605768616760667196375416901849615955687647529672806250507897968780954564415773574419843093829937094944*a[4]*a[5]*a[7]-83522945803218927143272697172304715702051516370421287077764861105797208257751172710043033982749732502673864122104745285694793848397146863804094508017786880*a[4]*a[7]^2+54148365713476352702026175006026820955664692693816358939561719366927259063497175211305130766993642399789393490162814579799098690810422835364603797996702160*a[4]*a[3]*a[7]-12416572431355783184146115277275541538215592296672337078836445756090524221703009481595411561277545637633257613562067323010331486676369113414366273036078461580*a[6]*a[3]*a[5]+25342033359690345055178729954858472888193856941199156563042027336413840936492137839840593101659135973490964450595217949119561806138168029781783656764939617739*a[6]*a[3]^2-23913072926086702376732546554384253629975805304704202482207160564669138660369215350484364998877253462539563653989021628304381081925205997834340414900754118322*a[4]*a[3]^2+17277017186149648209717672869536428267283400772879522367378957148030439492021840093175514917060902506322266073203907685241290745616567323488189138064774400*a[6]^3-23383277126811421453133277857335344234937843331037379520756912391102600952926599024670131546092037120371992884332927423154084266490274033213464651535955365282*a[4]^3';
strings{10}='-379645740934800599093813900258495007790808659117633078639390033045038592960269816603200067770356415337613595926273037574043393598977251429307741449324086053717844232*a[7]+7116649425285170650241541402395956517230130547759777762572854642768676754358996192896728539353813815955332451394433617898293002072442124676792253525604518484579744896*a[1]^2*a[7]-12174680873348267228330766657237287495901382162086569870153306351720165676884292341274551147182294663806118536687294931875405869062798923742526559235563616412094232320*a[4]*a[6]*a[3]-30462031168016927554908562163107575819987601674808765734321732328571358628282509795629348377680341588009365887042824620384100794356905601663841899021074339686195143264*a[2]*a[3]*a[5]-28808671084118009688539887610952731325711854738763711091638734783170513946340040912851429122915931008816770436180754815125497468094571696813460376120712359839916338834*a[4]*a[6]*a[2]-41969782539143012916250765091715122602536945750902542328128129021855849275635741596435132319529751514375023360930601293479478470091537542460674236357270927220312480*a[4]*a[6]*a[5]+11931257169140904451148111386890372282973005384590944762083839037697135972149736280098690962158536765139878180313628272936996936090028718584578520824714936460118117778*a[2]+2158010734342862861993294432177218836590956578220757235666874156349216757574398314168424741112738625763457030341376495262939109205504138270213173006275413419326693100*a[3]-2231124210602888526942414159679834783386875563293143754975822365585937286088147441554342507629892069793587162622674074915111672478053624953432031049159956227254015437*a[5]-18621431744501800587260656798894454931556745693302123353233920317185230808117346293366141548237520973748324274868022036420959632330824489576623915650924866893806147910*a[2]*a[1]^2-6660309163042718381625360074451489456517556734308848847317122084712155931786623039990758591118141421186429285696229338931870813291381887277283907981775445525738973049*a[6]*a[1]+12328557977472576333313585678674718548607827433721603540136062874780152185588993834853733504908748202415092545002352503441358623271900107244231821759547957085595019645*a[5]*a[1]^2-10673750902808724830521828687555990335778762359089617548646624953751123401216645535023151792177359009625004697193665138484937163005426912198027753055879541885008443478*a[4]*a[1]-17502317118702154922182820727680996796692112141857359565601070357006039759118113609072779878613923275975763445126164031345070391638697557600105215165817209012507047280*a[3]*a[1]^2+22544725605446658428627402935076330271429151054461236300674528576880150471690175357686880190039816695955323842817536888817162291998208379916524400338167366725557103712*a[2]^2*a[5]-2032194523416482408527035545170065490082502372787887693405332533620419624276283663525862379022607785359853708117930171638124627184938024110093440537187731900789336368*a[4]^2*a[7]-472734484128133744441634459348596212839343114412534749171193166284028990916037521716523039814094440129493947773628102686314662744903570710079619939820014665645428864*a[4]*a[6]*a[7]-54551071398905658840023797448624937823094666764143611085742874116438456983523799854471409647481768197292677846104339606278541588986197497871657737713955061205809834560*a[4]^2*a[2]+11269764093990807444925365632027806874125135712006041164056623562128331171413129267530709655611349405869280849151720484746283001136513667030145779921449820886506359456*a[4]^2*a[5]-82388864052586966948361435191220491085625310078751813483247654220606967689975565193614649391385823688267289243076075014410899006010122500741846859461299665879140252880*a[6]^2*a[2]-6334494101862561858706533612165189586948841625897213206103602996913487024702743738283067757775923470387636068387199213184253618639598621686329040842154906279154255016*a[2]^2*a[3]-58151407782285814490148082649503709256553622593180482645220503294183348395318392912431691800792236886792650227942325536739192829524565351078873883665370749388556978220*a[2]*a[3]^2-79730735386766725616935148673840551367723652897173142444903921460481510401274247784651444587190417433034210289507056583626816844361823118400631761131000973248727465816*a[2]*a[5]^2+5079696573628221019337611570395476760067829492649520642766896451343442987968490126608019015329256296113908535731508887760781392045020335568840463733579096282717120*a[5]^3-35478442336894132248765824211713786859684547153486035647595775868860834350037661488574218695635001260252149168033938240636291240459537690732875842680871741667684992404*a[2]^3-719138390108916529317347196731153654566357290625816108545620095431620299536752701778729307290764499345250370950022294056664083195696531844817913817572644227056109728*a[4]^2*a[3]';
strings{11}='155230499639888045922124736056355864425328007716427184749047123591250744212694750190370979296369443609338323794686380262016813943526717744*a[7]-1291509357655320342125948317568985010804696492655180484653760984914945236056951616340954762197258090039612596823091784231676473234565139816*a[1]^2*a[7]+5042480054441067307675173382141317047477134881378104475729500199937236971993695804094588954541965621633554599207210000569001078182440064*a[7]^3+7735049756426146460044334043524978514925418041210430470446729863877758048003759506999477481488592347207177506480525733343553568776301130448*a[4]*a[6]*a[3]+18765729201967968935250827469040551376189074731043035202220604477457286441974831330475701541483047108616855926639486988124087728099085695248*a[2]*a[3]*a[5]+19248832836694663542515280948650794480355092787293625407759719419218509777620426702134491864823871236502571306140515390285335508309364878350*a[4]*a[6]*a[2]+6388421815534506406457964245928076397193094582614308805638571548566425010281434638280977080942741015759147487662940753896943417163485440*a[4]*a[6]*a[5]-5980674903791713926742141402727507188931943235442454276327479662114497234609912644547611442263246904015073166659993199706917712691123483162*a[2]-178518231672769329958964061485742628596933029807643781871140115696473480686142387859239438241131753065658714983896334506275414894941428652*a[3]+465864743504307326633590600833376356337982141961654661161854003674143698865135358806275234360291518228473418185872912139481377841523565885*a[5]+3162569291604994354524234438484165685222157650935979168469261633695991183430168276803386940864031941429027726714897697414531051972765638270*a[2]*a[1]^2+2967153306115603635010562623942906104183990258261255956514725350946443198410274635196640083367225550679875176723780139415324595126951248929*a[6]*a[1]-687499399600445372145270561106401990704329243702807817251152168406131757775530194012155004836105400223393095384500418377383651450609990365*a[5]*a[1]^2+2039775880974975406745693899552668874166102290931694877333486450072118732714560024767368000747079460499085620119956407028650849034868096662*a[4]*a[1]+3891909280278342376697468247893089553918615583451717772749067859872290718958037817664132706679899578066343883003497851868605650083726118400*a[3]*a[1]^2-7806931494077588724338489819854461405073579795720354230752799929315634072215045216489197968222677688023568665828298187255308021094993517808*a[2]^2*a[5]+722776804417887026370528815889541146964522117927516067033978070099474575097031646914025103047903392094948391684700116612979187800909892328*a[4]^2*a[7]+145447246926645564969156703725893251855485518542123260982000788276549014767087401392244281126281960701948388439789940074078518972809042096*a[4]*a[6]*a[7]+26874439344985833805861319412505197816882619519707785576970141422096361715074419703820285391200744551467749903569997462139328936304808308452*a[4]^2*a[2]-7696619667670175584527009582425197036114658814199613955950197512230025143274660330317206201130709590132757427013373605273778641943494195808*a[4]^2*a[5]+39691672850132744640498563302434657550647522147356335535535396863485047091477128025201347663518519368995611976539084619384921703713400932192*a[6]^2*a[2]+8028766526953250701201453034420163064963569044305643124430021690704226552356151157603151163446096338445042392572997929344669541522712697864*a[2]^2*a[3]+25363932594450676654628123767050678607415975184277601903104809288544359743833410432444976949907948262958023453514339954370482542631229065964*a[2]*a[3]^2+39032877309335499190127868254462481888083101372718462754901231900627079359528728449415638745186771110363692829434650729649382517039411757320*a[2]*a[5]^2+17608224754338456492275199805773765568032105534879936299886382615840591748539857840256193130273427189722090461675931778391193947278615229844*a[2]^3-713776546077684505237128621518217653656049030696156937294644558960730091008310313543801620865809302305366087503165463603197399532468770344*a[4]^2*a[3]';
strings{12}='-619893619465247420436135430928272738901189478251345388610852489417690339731131642919688885386066177657173580619304685*a[1]^2+657622643031288158713049261852525515209978507125161585541359341952552575235130269533682081377197441650003335058235325*a[2]^2*a[1]^2+669597946540750422934587100577052790390585914142220130003067989272929869890414998645764140926028754028336830678843000*a[1]^4+1484473583697209951794366175973316050203147345494323628518757949148513014652349001128468332821278575539463215135631465*a[2]*a[3]*a[1]^2+39245689071899880251005824656110392803484307138502037263881771361635751388348761703770122748516956549270247354094000*a[2]^4+14735432968025339034482407392908449856308841679582040864675536946580089756032363317600211037061319079404156139412976*a[5]*a[7]+103256913765975509146072443415126650541662530858973590676496609716670315145241656754939312740556399270060481977715360*a[3]*a[7]-141844494512406829451353839680175786405370573057391260316280452549191050840179973735543296458790647371557910064576440*a[1]*a[4]*a[3]+302857012422988088285398353180008483426763538133966604834309553322240903101641096530423022586770652483750968508366320*a[1]*a[4]*a[7]-99619503770278319710659981930841133826422679331520210225888862926900470680579049436360299385871942047076738750543240*a[6]*a[1]*a[7]+384871333567926505170003180711311249849269569135773359216509984003554855289602279964297285160657743158314070599084800*a[4]*a[6]-7490491472158042351613698561979478264809936145915983028440088921953228403318734633848646661607573422671361161935680*a[7]^2-1964220902310056162432899131731163483239047357209251455183201628637725423319588944160292611747189646207189711672195425*a[1]*a[4]*a[5]-405637060215586488169773100343733022687389743055958051739424277982169055672840405279978803467675864689210022069751945*a[2]*a[3]-308625633767261614294543066305149727953878570085986588531464080622265592163362578810086833975618996818104982784510112*a[2]*a[5]+712650398804411376362521788503791411908047222509344679115627821791985174768781725285287708086989959408196287690670155*a[1]*a[2]*a[4]+564576184479218131174947875039304894890280062154300927757496142165538345852672086617143211656170684915306705872799960*a[6]*a[2]*a[1]+2343544375892343144748833684640179472878366881711825783536889748602779196204996476999146658023818214470851551625485105*a[6]*a[3]*a[1]+28967317434599229738780711266751327077374654601389775118760802251756010946736042563633248261827265369761559405340960*a[6]*a[5]*a[1]+277615485111778005610545596398320432845232630333573761570671238335263992613937568056659787238319107304368941157380032*a[3]*a[5]+1013897471536721994924282435882601386656165436315523622428367152212911503024214190285210087425914489151971673316637440*a[4]^2+690905468611932639739754561581933277367264103009796066731521773265616144315341222453663530223023393917014034976777355*a[2]^2+622685602754300753036145442148769874496471500002173042771206710383561342787011124710049950191796851019742908247282240*a[3]^2+1541850784868895249250876542347419333387318793994296946440026644362217036242794199687456993428652936996970107348240112*a[5]^2+1537250210974527644492537062783985008716694398755615277764473821957377065389743208196905663475012963710138290880200880*a[6]^2-438299661292026045773437838972869511374624353988681708793116592302218259111611436562632942397976315826352468628611115';
strings{13}='-631185729965067022105435762414669162520720405453326851923387327123597300795647717841624543730299527714270050694357359181411736513734670672295*a[1]^2+556484306586584872544663919955983091893936362212401226171861875378951597607411284483843043930731494981227710121035862218165389131520298341025*a[2]^2*a[1]^2+261233048198790492098186276095276839471293237041897300823187973154619298913071710283046679884498407862664173742492038456579433091773060261000*a[1]^4+1920362594122707598795127638012390320373601481925609577839991288176649561578788222920801533487666997258500525978245304886509086867800897511005*a[2]*a[3]*a[1]^2-7690086806145547764182915371791011590555693388345129975949510463216428536388832294245130015898602785256071872635058497351290471631126679448*a[5]*a[7]+59595843772198612109451602623217250632348106826200387219397397868590623595723639753351833868038722532202273538073683696298877340646670839720*a[3]*a[7]-187983333441690135591593761218847496591933176055426972624151927779320982993446481494939596284565487136045428346013200340968815396362316985580*a[1]*a[4]*a[3]-121220896980430129423230460200956848437045137565147185637275401250266895904930914796617763617336388229934982121533721649463826464882634360560*a[1]*a[4]*a[7]-84228881005813340086626521446530934891242838725972742017008794560560142853251099268169722626317319236588147737449750870387275880057366593880*a[6]*a[1]*a[7]+88534473613567332977231036880288063489386164627090927706811220387508639879448832951728571757789066261511750737926257668386488545830973684100*a[4]*a[6]+557426251061641449268412919769362920309193464753692588780860843858317447648198994098542281468634472616469625304632175269468479283025128640*a[7]^2+1410554045830375576148458791296397931321091913812879119054162313866883725469631261832832109926727546092107215298634165224043831950478081162775*a[1]*a[4]*a[5]-288289495141703451904718204956115338286278893649080551046500667679667935047599207961946000960589104839838245933696009891553808356542530440365*a[2]*a[3]-82164809148276914805692006482746842244083223724439994124505867373953724787841733702196303792221448050178109216724227339115516438086656329384*a[2]*a[5]-367849975064581533202325973717732398314429221867047144525941214008551207467146220862490621407776110035417105329508802066119665018391943799265*a[1]*a[2]*a[4]+394210462183522448283354456529128022356306601462740269340664337368206648658083918678591280403925757579614429926678553277886391849372523085370*a[6]*a[2]*a[1]-1125689595966180149562727020329737572443525881771794512224086533233669373680543883849073531178006524562417432655995465455180097880441461782615*a[6]*a[3]*a[1]+21147610375950360257884883496924489164431732075758224340500035523692139136652892678511015141603619533590647000333457602894585934747868151920*a[6]*a[5]*a[1]-2892406678491635790511664760279664936174984357782184206016873299744468504615078182048117519793351977715709963147587420004542938827760628576*a[3]*a[5]+93268371813469745650217686407936784273661772356714001584038978671845430789408604088658539725312282177692344380636824171910127692363330308495+539033595046288411846837582246789328313735197011857701600240757205186620951622707654870619074765937885631824365716953806048023452807463482080*a[4]^2+411765686028175204747837870596297674301610237458151355417777372611098692785774467624097546055868624565984927172497853332584173969686400214735*a[2]^2+132891513814603914011432977924235375508109567262132465522191332921840564141014002689178875527050082838538406704779221541463817350575838822680*a[3]^2+631824705353334694064885071003109798723913720911837140946204077117476513746904959510461162800357746181221226161629577970431389076459330082984*a[5]^2+613443690259932987638435753116659666372385654264641274047816718984749146219352228040599194338925856004721695140237340843140648470471063594160*a[6]^2+7740926594389689924282044012320850097225128555965981722744494691171590520004365409564174497217378005891513967078702884773178786611796552000*a[4]^4';
strings{14}='-120259471260867014005758668628783144121186051241992529830196524079799527336800468170794263277689697697429576644501697495*a[1]+43287966136759031691157641744890892479993645746568522252677960282082237962111864032947844345714103035771971906959129950*a[6]*a[3]-15221721996858738213426929505727238302499825357010752996059898468815782319958119823182108233065122760248731409917577275*a[6]*a[2]*a[1]^2-1390415643953631756603404735588948448771941902073523708359389776250736028004492597154057297243003524650696337544480240*a[1]*a[7]^2+272821649945872064282178321709952995690260056578130233104989271934972227793149054230335870393639301167340596334578878360*a[1]^2*a[2]*a[4]-275473797620681921700098123201498987024512900632885928897148197840022497715057113618020113465202379316169668669715912*a[6]*a[7]-95791169778902231522117081546288465582699971582976944773728195119542977671624291911683470531971843419735571978213048*a[6]*a[5]-65406073905321203773603056675641701511978569888562268538765565524746580761574108087491522544897070477010546346239336*a[5]*a[1]*a[7]+414407201183734254484648709218581973395884684999117239866259723457066558020146522566867295371286478439850110754574055820*a[1]*a[4]^2+412122540795748350357215394637409078746748029696673830917268394362560457914141179313452959164537729620876378162153313660*a[3]^2*a[1]+600834878171028417399331953959290999897823057256691405809381633807726732961798506436216104364779020666038801399021165800*a[6]^2*a[1]+11961890406600332213381138648937477951060399050828685286163009124830784525252722324146971903347068990420204180231622100*a[3]*a[1]*a[7]-320108690195131738036902199051570406471682852949961466588500702488881378431323417278700066510291717582829769455110458345*a[1]^3-197697329969579962587398940451576154528469365604207962091674071706509936171897401071045548725470944323388613878376202570*a[3]*a[4]*a[1]^2+73788293989245700942588084114339134303987209245738765257472234121683421998284918589276654614789261326222619036015484200*a[2]*a[5]*a[1]+643388739159785214245012564912112390563861174137504004543780990597276342708281546726051795438698999845606481831804170760*a[3]*a[5]*a[1]+9479103516059937036048373593755231546499273307118771086373876044345307885032501917395228463241682798265230844733875750*a[2]*a[4]+21269901737792631442026897536552656900148732671815129329124702947480840525176283338062029717022797709379557827131548630*a[6]*a[2]+247401778375910846010282344153945333819447170165632299424641857426788506459615393305110555611517469439972763200738852340*a[2]^2*a[1]+370846681243928358598268442997638185568998410457388226415984072594697790312376808174391897057082407385344140649624427160*a[2]*a[3]*a[1]+651413769639635074189046430323046473941754264616221004451913093458755058562223414804033378056074733785455901596013068600*a[1]*a[4]*a[6]-38723759236176514372250878249956615390010470538893558086926458084175254217165891325621097353989554716638248816813614220*a[4]*a[3]+21893774285447506553351342336708302185919829915070261317420476666981119439134391252414947813978500854727317918497994660*a[4]*a[7]+607080002281621505744282450597983182519295155587532828324682247132259650465307196813918204002683688519535978548271263840*a[5]^2*a[1]-22836304684864365661608991805990475406728174813214548943728021676115340991394201796675606041449429634408243687143736622*a[4]*a[5]+400752226681299895659356941054550738456734511698961697781036574547435744035611304871650660160751573043992501209072499500*a[1]^5';
polysys=lti2polysys(strings,[],[],'na',7);

% Van der Pol oscillator
% Haotian's equations, N=7, simplified
strings{1}='-42*a[4]^2*a[8]*a[1]-21*a[4]*a[6]^2*a[1]-42*a[4]*a[2]*a[3]*a[1]-42*a[4]*a[2]*a[5]*a[1]+21*a[4]*a[5]^2*a[1]-42*a[6]^2*a[8]*a[1]-21*a[6]*a[2]^2*a[1]-42*a[6]*a[2]*a[3]*a[1]-42*a[6]*a[3]*a[5]*a[1]-21*a[8]^3*a[1]-42*a[8]*a[2]^2*a[1]-42*a[8]*a[3]^2*a[1]-42*a[8]*a[5]^2*a[1]-21*a[8]*a[7]^2*a[1]+7*a[8]*a[1]-49*a[7]*a[1]^2+a[7]';
strings{2}='-30*a[4]^2*a[6]*a[1]-30*a[4]*a[6]*a[8]*a[1]-15*a[4]*a[2]^2*a[1]-30*a[4]*a[2]*a[3]*a[1]+30*a[4]*a[2]*a[7]*a[1]-30*a[4]*a[5]*a[7]*a[1]-15*a[6]^3*a[1]-30*a[6]*a[8]^2*a[1]-30*a[6]*a[2]^2*a[1]-30*a[6]*a[3]^2*a[1]+30*a[6]*a[3]*a[7]*a[1]-15*a[6]*a[5]^2*a[1]-30*a[6]*a[7]^2*a[1]+5*a[6]*a[1]-15*a[8]*a[2]^2*a[1]-30*a[8]*a[2]*a[3]*a[1]-30*a[8]*a[3]*a[5]*a[1]-25*a[5]*a[1]^2+a[5]';
strings{3}='-9*a[4]^3*a[1]-18*a[4]*a[6]^2*a[1]-18*a[4]*a[8]^2*a[1]-18*a[4]*a[2]^2*a[1]+18*a[4]*a[2]*a[5]*a[1]+18*a[4]*a[2]*a[7]*a[1]-9*a[4]*a[3]^2*a[1]-18*a[4]*a[5]^2*a[1]-18*a[4]*a[7]^2*a[1]+3*a[4]*a[1]-9*a[6]^2*a[8]*a[1]-9*a[6]*a[2]^2*a[1]-18*a[6]*a[2]*a[3]*a[1]+18*a[6]*a[2]*a[7]*a[1]-18*a[6]*a[5]*a[7]*a[1]-18*a[8]*a[2]*a[3]*a[1]-18*a[8]*a[2]*a[5]*a[1]+9*a[8]*a[5]^2*a[1]-9*a[3]*a[1]^2+a[3]';
strings{4}='a[2]-a[2]*a[1]^2-3*a[4]^2*a[6]*a[1]+3*a[4]^2*a[8]*a[1]-3*a[4]*a[2]^2*a[1]+3*a[6]*a[3]^2*a[1]-3*a[8]*a[3]^2*a[1]-6*a[4]*a[6]*a[8]*a[1]+6*a[4]*a[2]*a[5]*a[1]-6*a[6]*a[2]*a[3]*a[1]-6*a[4]*a[3]*a[5]*a[1]+6*a[4]*a[3]*a[7]*a[1]+6*a[6]*a[2]*a[7]*a[1]-6*a[8]*a[2]*a[5]*a[1]-6*a[4]*a[5]*a[7]*a[1]-6*a[6]*a[3]*a[7]*a[1]+6*a[8]*a[3]*a[5]*a[1]';
strings{5}='21*a[4]^2*a[2]*a[1]-42*a[4]^2*a[7]*a[1]+42*a[4]*a[6]*a[2]*a[1]-42*a[4]*a[6]*a[5]*a[1]+21*a[6]^2*a[3]*a[1]-42*a[6]^2*a[7]*a[1]-21*a[8]^2*a[7]*a[1]+49*a[8]*a[1]^2-a[8]-21*a[2]^2*a[5]*a[1]-42*a[2]^2*a[7]*a[1]-21*a[2]*a[3]^2*a[1]-42*a[2]*a[3]*a[5]*a[1]-42*a[3]^2*a[7]*a[1]-21*a[3]*a[5]^2*a[1]-42*a[5]^2*a[7]*a[1]-21*a[7]^3*a[1]+7*a[7]*a[1]';
strings{6}='15*a[4]^2*a[2]*a[1]-30*a[4]^2*a[5]*a[1]-30*a[4]*a[6]*a[7]*a[1]-30*a[4]*a[8]*a[2]*a[1]+30*a[4]*a[8]*a[5]*a[1]-15*a[6]^2*a[5]*a[1]-30*a[6]*a[8]*a[3]*a[1]+25*a[6]*a[1]^2-a[6]-30*a[8]^2*a[5]*a[1]-15*a[2]^2*a[3]*a[1]-30*a[2]^2*a[5]*a[1]-15*a[2]^2*a[7]*a[1]-15*a[2]*a[3]^2*a[1]-30*a[2]*a[3]*a[7]*a[1]-30*a[3]^2*a[5]*a[1]-30*a[3]*a[5]*a[7]*a[1]-15*a[5]^3*a[1]-30*a[5]*a[7]^2*a[1]+5*a[5]*a[1]';
strings{7}='-9*a[4]^2*a[3]*a[1]-18*a[4]*a[6]*a[2]*a[1]-18*a[4]*a[8]*a[2]*a[1]+9*a[4]*a[1]^2-a[4]-18*a[6]^2*a[3]*a[1]+9*a[6]^2*a[7]*a[1]-18*a[6]*a[8]*a[2]*a[1]-18*a[6]*a[8]*a[5]*a[1]-18*a[8]^2*a[3]*a[1]-3*a[2]^3*a[1]-18*a[2]^2*a[3]*a[1]-9*a[2]^2*a[5]*a[1]-18*a[2]*a[3]*a[5]*a[1]-18*a[2]*a[3]*a[7]*a[1]-18*a[2]*a[5]*a[7]*a[1]-9*a[3]^3*a[1]-18*a[3]*a[5]^2*a[1]-18*a[3]*a[7]^2*a[1]+3*a[3]*a[1]-9*a[5]^2*a[7]*a[1]';
strings{8}='-6*a[4]^2*a[2]+3*a[4]^2*a[5]+3*a[4]^2*a[7]-6*a[4]*a[6]*a[2]-6*a[4]*a[6]*a[3]+6*a[4]*a[6]*a[7]-6*a[4]*a[8]*a[3]-6*a[4]*a[8]*a[5]-6*a[6]^2*a[2]-6*a[6]*a[8]*a[2]-6*a[6]*a[8]*a[3]-6*a[8]^2*a[2]-3*a[2]^3-3*a[2]^2*a[3]-6*a[2]*a[3]^2-6*a[2]*a[3]*a[5]-6*a[2]*a[5]^2-6*a[2]*a[5]*a[7]-6*a[2]*a[7]^2+a[2]-3*a[3]^2*a[5]-3*a[3]^2*a[7]-6*a[3]*a[5]*a[7]';
polysys=lti2polysys(strings,[],[],'na',8);

% Biochemical example, N=1,
strings{1}='a[4]-6.2831853072*a[7]+a[6]+0.5';
strings{2}='a[5]-1.0*a[4]-6.2831853072*a[8]';
strings{3}='a[6]-1.0*a[4]+10.0*a[1]*a[6]+10.0*a[4]*a[3]';
strings{4}='a[1]+a[3]';
strings{5}='a[2]-1.0*a[1]';
strings{6}='a[3]-1.0*a[1]+20.0*a[7]*a[9]+10.0*a[1]*a[3]+20.0*a[4]*a[6]';
strings{7}='-1.0*a[7]-1.0*a[9]-6.2831853072*a[4]';
strings{8}='a[7]-1.0*a[8]-6.2831853072*a[5]';
strings{9}='a[7]-1.0*a[9]-10.0*a[7]*a[3]-10.0*a[9]*a[1]';
polysys=lti2polysys(strings,[],[],'na',9);

% Biochemical example,n=3, N=2
strings{1}='a[7]-12.566370614*a[13]+a[9]';
strings{2}='a[8]-1.0*a[7]-12.566370614*a[14]';
strings{3}='a[9]-1.0*a[7]-10.0*a[10]*a[12]+10.0*a[1]*a[9]+10.0*a[4]*a[6]+10.0*a[7]*a[3]';
strings{4}='a[4]-6.2831853072*a[10]+a[6]+0.5';
strings{5}='a[5]-1.0*a[4]-6.2831853072*a[11]';
strings{6}='a[6]-1.0*a[4]+10.0*a[10]*a[15]+10.0*a[13]*a[12]+10.0*a[1]*a[6]+10.0*a[4]*a[3]+10.0*a[4]*a[9]+10.0*a[7]*a[6]';
strings{7}='a[1]+a[3]';
strings{8}='a[2]-1.0*a[1]';
strings{9}='a[3]-1.0*a[1]+20.0*a[10]*a[12]+20.0*a[13]*a[15]+10.0*a[1]*a[3]+20.0*a[4]*a[6]+20.0*a[7]*a[9]';
strings{10}='-1.0*a[13]-1.0*a[15]-12.566370614*a[7]';
strings{11}='a[13]-1.0*a[14]-12.566370614*a[8]';
strings{12}='a[13]-1.0*a[15]-10.0*a[10]*a[6]-10.0*a[13]*a[3]-10.0*a[12]*a[4]-10.0*a[15]*a[1]';
strings{13}='-1.0*a[10]-1.0*a[12]-6.2831853072*a[4]';
strings{14}='a[10]-1.0*a[11]-6.2831853072*a[5]';
strings{15}='a[10]-1.0*a[12]-10.0*a[10]*a[3]-10.0*a[12]*a[1]+10.0*a[10]*a[9]-10.0*a[13]*a[6]+10.0*a[12]*a[7]-10.0*a[15]*a[4]';
polysys=lti2polysys(strings,[],[],'na',15);

% Biochemical example,n=2, N=2
strings{1}='a[5]-12.566370614*a[9]+a[6]';
strings{2}='a[6]-1.0*a[5]-10.0*a[7]*a[8]+10.0*a[1]*a[6]+10.0*a[3]*a[4]+10.0*a[5]*a[2]';
strings{3}='a[3]-6.2831853072*a[7]+a[4]+0.5';
strings{4}='a[4]-1.0*a[3]+10.0*a[7]*a[10]+10.0*a[9]*a[8]+10.0*a[1]*a[4]+10.0*a[3]*a[2]+10.0*a[3]*a[6]+10.0*a[5]*a[4]';
strings{5}='a[1]+a[2]';
strings{6}='a[2]-1.0*a[1]+20.0*a[7]*a[8]+20.0*a[9]*a[10]+10.0*a[1]*a[2]+20.0*a[3]*a[4]+20.0*a[5]*a[6]';
strings{7}='-1.0*a[9]-1.0*a[10]-12.566370614*a[5]';
strings{8}='a[9]-1.0*a[10]-10.0*a[7]*a[4]-10.0*a[9]*a[2]-10.0*a[8]*a[3]-10.0*a[10]*a[1]';
strings{9}='-1.0*a[7]-1.0*a[8]-6.2831853072*a[3]';
strings{10}='a[7]-1.0*a[8]-10.0*a[7]*a[2]-10.0*a[8]*a[1]+10.0*a[7]*a[6]-10.0*a[9]*a[4]+10.0*a[8]*a[5]-10.0*a[10]*a[3]';
polysys=lti2polysys(strings,[],[],'na',10);
% Robiano's example from presentation
polysys{1,1}=[1/4 1 -1];polysys{1,2}=[2 0;0 2;0 0];
polysys{2,1}=[1 1/4 -1];polysys{2,2}=[2 0;0 2;0 0];
% noisy conics, 1e-5
polysys2{1,1}=[1/4 1 -1 1e-5];polysys2{1,2}=[2 0;0 2;0 0;1 1];
polysys2{2,1}=[1 1/4 -1 1e-5];polysys2{2,2}=[2 0;0 2;0 0;1 1];
% noisy conics, 1e-15
polysys3{1,1}=[1/4 1 -1 1e-15];polysys3{1,2}=[2 0;0 2;0 0;1 1];
polysys3{2,1}=[1 1/4 -1 1e-15];polysys3{2,2}=[2 0;0 2;0 0;1 1];

% example discontinuous dependence groebner basis, stetter pg 323
% epsilon = 0
polysys{1,1}=[1 1 -1];polysys{1,2}=[2 0;0 2;0 0];
polysys{2,1}=[1 -3];polysys{2,2}=[0 3;2 1];
% epsilon = 1e-5
polysys2{1,1}=[1 1e-5 1 -1];polysys2{1,2}=[2 0;1 1;0 2;0 0];
polysys2{2,1}=[1 -3];polysys2{2,2}=[0 3;2 1];


polysys{1,1}=[1 3 -7];polysys{1,2}=[0 0 2;0 1 0;0 0 1];
polysys{2,1}=[1 -4];polysys{2,2}=[0 1 1;0 1 0];
polysys{3,1}=[1 -4];polysys{3,2}=[1 0 1;0 1 0];
polysys{4,1}=[1 -4];polysys{4,2}=[0 2 0;0 1 0];
polysys{5,1}=[1 -4];polysys{5,2}=[1 1 0;0 1 0];
polysys{6,1}=[1 -8 14 8 -15 15];polysys{6,2}=[5 0 0;4 0 0;3 0 0;2 0 0;1 0 0;0 1 0];

% example 19 from "Computing Border Bases" by Achim Kehrein & Martin Kreuzer
polysys{1,1}=[1 -1];polysys{1,2}=[3 0;1 0];
polysys{2,1}=[1 -1];polysys{2,2}=[0 3;0 1];
polysys{3,1}=[1 -.5 -.5];polysys{3,2}=[2 1;0 1;0 2];
polysys{4,1}=[1 -1 -.5 1 -.5];polysys{4,2}=[1 1;1 0;0 1;2 0;0 2];
polysys{5,1}=[1 -1 -.5 1 -.5];polysys{5,2}=[1 2;1 0;0 1;2 0;0 2];

% stetter - numerical polynomial algebra - pg 434
polysys{1,1}=[4.831 4.597 .417 1.688 .351 -1.428 3.344 .64 3.308 2.728];
polysys{1,2}=[2 0 0;1 1 0;0 2 0;1 0 1;0 1 1;0 0 2;1 0 0;0 1 0;0 0 1;0 0 0];
polysys{2,1}=[4.036 3.655 2.988 2.19 1.473 1.96 4.27 3.572 .853 .239];
polysys{2,2}=[2 0 0;1 1 0;0 2 0;1 0 1;0 1 1;0 0 2;1 0 0;0 1 0;0 0 1;0 0 0];
polysys{3,1}=[4.229 1.95 2.988 1.298 4.86 1.249 3.056 1.267 2.887 3.853];
polysys{3,2}=[2 0 0;1 1 0;0 2 0;1 0 1;0 1 1;0 0 2;1 0 0;0 1 0;0 0 1;0 0 0];

%% Difficult problems

% Cassou
strings{1}='15*a[1]^4*a[2]*a[3]^2+6*a[1]^4*a[2]^3+21*a[1]^4*a[2]^2*a[3]-144*a[1]^2*a[2]-8*a[1]^2*a[2]^2*a[4]-28*a[1]^2*a[2]*a[3]*a[4]-648*a[1]^2*a[3]+36*a[1]^2*a[3]^2*a[4]+9*a[1]^4*a[3]^3-120';
strings{2}='30*a[2]^3*a[1]^4*a[3]-32*a[3]*a[4]^2*a[2]-720*a[3]*a[1]^2*a[2]-24*a[2]^3*a[1]^2*a[4]-432*a[2]^2*a[1]^2+576*a[4]*a[2]-576*a[3]*a[4]+16*a[2]*a[1]^2*a[3]^2*a[4]+16*a[3]^2*a[4]^2+16*a[4]^2*a[2]^2+9*a[2]^4*a[1]^4+5184+39*a[3]^2*a[1]^4*a[2]^2+18*a[3]^3*a[1]^4*a[2]-432*a[3]^2*a[1]^2+24*a[3]^3*a[1]^2*a[4]-16*a[2]^2*a[1]^2*a[3]*a[4]-240*a[2]';
strings{3}='216*a[3]*a[1]^2*a[2]-162*a[3]^2*a[1]^2-81*a[2]^2*a[1]^2+5184+1008*a[4]*a[2]-1008*a[3]*a[4]+15*a[2]^2*a[1]^2*a[3]*a[4]-15*a[2]^3*a[1]^2*a[4]-80*a[3]*a[4]^2*a[2]+40*a[3]^2*a[4]^2+40*a[4]^2*a[2]^2';
strings{4}='261+4*a[3]*a[1]^2*a[2]-3*a[3]^2*a[1]^2-4*a[2]^2*a[1]^2+22*a[4]*a[2]-22*a[3]*a[4]';
polysys=lti2polysys(strings,[],[],'na',4);

% STLS Philippe Dreesen, dreg = 21, 18 affine solutions for d=23
d1=5;
d2=3;
d3=5;
d4=7;
d5=2;
d6=2;
polysys{1,1} = [-4*d2 4 -8*d3 8 -12*d4 12 -16*d5 16 -10*d6 10];
polysys{1,2} = [ 0 1; 1 2; 1 1; 3 2; 2 1; 5 2; 3 1; 7 2; 4 1; 9 2];
polysys{2,1} = [-2*d1 +2 -4*d2 +4 -4*d3 +4 -4*d4 +4 -4*d5 +4 -2*d6 +2];
polysys{2,2} = [ 0 0; 0 1; 1 0;2 1; 2 0; 4 1; 3 0; 6 1; 4 0; 8 1; 5 0; 10 1];

% Kanno
polysys{1,1}=[2,-2,-1,10,1,8,5,-10];polysys{1,2}=[1 0 0 1;0 1 1 0;0 1 0 2;0 0 0 0;0 0 0 1;0 1 0 0;1 0 0 0;0 1 0 1];
polysys{2,1}=[-2,-1,2,17,15,1,-86,33,-66];polysys{2,2}=[1 0 1 0;1 0 0 2;0 1 1 1;0 0 0 0;0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0;0 1 0 1];
polysys{3,1}=[2,-1,50,15,-345,-79,79,-158];polysys{3,2}=[1 0 1 1;0 1 2 0;0 0 0 1;0 0 1 0;0 1 0 0;0 0 0 0;1 0 0 0;0 1 0 1];
polysys{4,1}=[50,-1,-50,50,-100,-250];polysys{4,2}=[0 0 1 0;1 0 2 0;0 0 0 0;1 0 0 0;0 1 0 1;0 1 0 0];

% http://www.math.uic.edu/~jan/Demo/trinks.html
% TITLE : system of Trinks from the PoSSo test suite
% [x y z u v t], d=7, truncated co-rank stabilizes at 10, finds 10
% solutions
% 45*y + 35*u - 165*v - 36
polysys{1,1} = [45 35 -165 -36];
polysys{1,2} = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 0];
% 35*y + 25*z + 40*t - 27*u
polysys{2,1} = [35 25 40 -27];
polysys{2,2} = [0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 1;0 0 0 1 0 0];
%  25*y*u - 165*v**2 + 15*x - 18*z + 30*t
polysys{3,1} = [25 -165 15 -18 30];
polysys{3,2} = [0 1 0 1 0 0; 0 0 0 0 2 0;1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 1];
% 15*y*z + 20*t*u - 9*x
polysys{4,1} = [15 20 -9];
polysys{4,2} = [0 1 1 0 0 0; 0 0 0 1 0 1;1 0 0 0 0 0];
% -11*v**3 + x*y + 2*z*t
polysys{5,1} = [-11 1 2];
polysys{5,2} = [0 0 0 0 3 0; 1 1 0 0 0 0;0 0 1 0 0 1];
% -11*u*v + 3*v**2 + 99*x
polysys{6,1} = [-11 3 99];
polysys{6,2} = [0 0 0 1 1 0; 0 0 0 0 2 0;1 0 0 0 0 0];

% sturmfels IWT, Grobner basis: 
polysys{1,1} = [13003050 2744 -2116125 -6290625 ];
polysys{1,2} = [1 0;0 2;0 1; 0 0];
polysys{2,1} = [134456 -10852275 -4304728125 935718750];
polysys{2,2} = [0 3;0 2;0 1;0 0]; 

% reduced Asys contains all pure components at getDreg(polysys)+2!!
% basis syzygies are all found at sum(getDorig(polysys))=getDreg(polysys)+2
% gcd's = 1, deg(f1) = 2, deg(f2) = 4, deg(f3) = 4, deg(f1f2) = 6,
% deg(f1f3) = 6, deg(f2f3) = 8
% x^2 + x*z - 2*y + 5
polysys{1,1} = [1 1 -2 5]; polysys{1,2} = [2 0 0;1 0 1;0 1 0;0 0 0];
% 2*x^3*y + 7*y*z^2 - 4*x*y*z + 3*x - 2
polysys{2,1} = [2 7 -4 3 -2];polysys{2,2} = [3 1 0;0 1 2;1 1 1;1 0 0;0 0 0];
% y^4 + 2*y*z + 5*x^2 -5
polysys{3,1} = [1 2 5 -5];polysys{3,2} = [0 4 0;0 1 1;2 0 0;0 0 0];

% x^2*z^2 + x*z - 2*y + 5
polysys{1,1} = [1 1 -2 5]; polysys{1,2} = [2 0 2;1 0 1;0 1 0;0 0 0];
% 2*x^3*y + 7*y*z^2 - 4*x*y*z + 3*x - 2
polysys{2,1} = [2 7 -4 3 -2];polysys{2,2} = [3 1 0;0 1 2;1 1 1;1 0 0;0 0 0];
% y^4 + 2*y*z + 5*x^2 -5
polysys{3,1} = [1 2 5 -5];polysys{3,2} = [0 4 0;0 1 1;2 0 0;0 0 0];

% (1-x+y) (x^2+5-6y)
polysys{1,1} = [1 5 -1 -1 -5 6 1 -6];polysys{1,2} = [2 0 0;0 0 0;0 1 0;3 0 0;1 0 0;1 1 0;2 1 0;0 2 0];
% (1-x+y) (2x+5)
polysys{2,1} = [-3 5 -2 2 5];polysys{2,2} = [1 0 0;0 0 0;2 0 0;1 1 0;0 1 0];
% x+y+z
polysys{3,1} = [1 1 1];polysys{3,2} = [1 0 0;0 1 0;0 0 1];

% flip f2 with f3
% (1-x+y) (x^2+5-6y)
polysys{1,1} = [1 5 -1 -1 -5 6 1 -6];polysys{1,2} = [2 0 0;0 0 0;0 1 0;3 0 0;1 0 0;1 1 0;2 1 0;0 2 0];
% x+y+z
polysys{2,1} = [1 1 1];polysys{2,2} = [1 0 0;0 1 0;0 0 1];
% (1-x+y) (2x+5)
polysys{3,1} = [-3 5 -2 2 5];polysys{3,2} = [1 0 0;0 0 0;2 0 0;1 1 0;0 1 0];

% flip f1 with f3
% x+y+z
polysys{1,1} = [1 1 1];polysys{1,2} = [1 0 0;0 1 0;0 0 1];
% (1-x+y) (2x+5)
polysys{2,1} = [-3 5 -2 2 5];polysys{2,2} = [1 0 0;0 0 0;2 0 0;1 1 0;0 1 0];
% (1-x+y) (x^2+5-6y)
polysys{3,1} = [1 5 -1 -1 -5 6 1 -6];polysys{3,2} = [2 0 0;0 0 0;0 1 0;3 0 0;1 0 0;1 1 0;2 1 0;0 2 0];

% Sturmfels IWT, gxel ordering
% 291950 x - 42225 y + 2744 x^1  - 143125,
polysys{1,1} = [291950 42225 2744  -143125]; polysys{1,2} = [1 0;0 1;2 0;0 0];
% 1894650 x - 291125 y - 920625 + 2744 x y,
polysys{2,1} = [1894650 -291125  -920625 2744]; polysys{2,2} = [1 0;0 1;0 0;1 1];
% 13003050 x + 2744 y^2  - 2116125 y - 6290625
polysys{3,1} = [13003050 2744 -2116125 -6290625]; polysys{3,2} = [1 0;0 2;0 1;0 0];

% http://www.math.uic.edu/~jan/Demo/caprasse.html, needs 
% d = sum(getDorig(polysys))+1 in order to see al basis-syzygies!
% the system caprasse of the PoSSo test suite
% There are 54 isolated solutions, so 6 ones with zero
% components which are not counted by the mixed volume(=48)
% [x y z t]
% y**2*z + 2*x*y*t - 2*x-z;
polysys{1,1} = [1 2 -2 -1];
polysys{1,2} = [0 2 1 0;1 1 0 1;1 0 0 0;0 0 1 0];
% -x**3*z + 4*x*y**2*z + 4*x**2*y*t + 2*y**3*t + 4*x**2 - 10*y**2 + 4*x*z - 10*y*t + 2;
polysys{2,1} = [-1 4 4 2 4 -10 4 -10 2];
polysys{2,2} = [3 0 1 0;1 2 1 0;2 1 0 1;0 3 0 1;2 0 0 0;0 2 0 0;1 0 1 0;0 1 0 1;0 0 0 0]; 
% 2*y*z*t + x*t**2 - x - 2*z;
polysys{3,1} = [2 1 -1 -2];
polysys{3,2} = [0 1 1 1;1 0 0 2;1 0 0 0;0 0 1 0];
% -x*z**3 + 4*y*z**2*t + 4*x*z*t**2 + 2*y*t**3 + 4*x*z + 4*z**2 - 10*y*t - 10*t**2 + 2;
polysys{4,1} = [-1 4 4 2 4 4 -10 -10 2];
polysys{4,2} = [1 0 3 0;0 1 2 1;1 0 1 2;0 1 0 3;1 0 1 0;0 0 2 0;0 1 0 1;0 0 0 2;0 0 0 0];

% http://www.math.uic.edu/~jan/Demo/reimer5.html
% The 5-dimensional system of Reimer.
% The system has 12 roots with all components positive,
% but any fine mixed subdivision has only one mixed cell.
% total degree : 720
% mixed volume : 720
% [x y z t u]
% -1 + 2*x**2 - 2*y**2 + 2*z**2 - 2*t**2 + 2*u**2;
polysys{1,1} = [-1 2 -2 2 -2 2];
polysys{1,2} = [0 0 0 0 0;2 0 0 0 0;0 2 0 0 0;0 0 2 0 0;0 0 0 2 0;0 0 0 0 2];
% -1 + 2*x**3 - 2*y**3 + 2*z**3 - 2*t**3 + 2*u**3;
polysys{2,1} = [-1 2 -2 2 -2 2];
polysys{2,2} = [0 0 0 0 0;3 0 0 0 0;0 3 0 0 0;0 0 3 0 0;0 0 0 3 0;0 0 0 0 3];
% -1 + 2*x**4 - 2*y**4 + 2*z**4 - 2*t**4 + 2*u**4;
polysys{3,1} = [-1 2 -2 2 -2 2];
polysys{3,2} = [0 0 0 0 0;4 0 0 0 0;0 4 0 0 0;0 0 4 0 0;0 0 0 4 0;0 0 0 0 4];
% -1 + 2*x**5 - 2*y**5 + 2*z**5 - 2*t**5 + 2*u**5;
polysys{4,1} = [-1 2 -2 2 -2 2];
polysys{4,2} = [0 0 0 0 0;5 0 0 0 0;0 5 0 0 0;0 0 5 0 0;0 0 0 5 0;0 0 0 0 5];
% -1 + 2*x**6 - 2*y**6 + 2*z**6 - 2*t**6 + 2*u**6;
polysys{5,1} = [-1 2 -2 2 -2 2];
polysys{5,2} = [0 0 0 0 0;6 0 0 0 0;0 6 0 0 0;0 0 6 0 0;0 0 0 6 0;0 0 0 0 6];

% http://www.math.uic.edu/~jan/Demo/virasoro.html
% the construction of Virasoro algebras
% Schrans, S. and Troost, W.:
% `Generalized Virasoro Constructions for SU(3)',
%  Nuclear Phys. B, Vol. 345, No. 2--3, 1990, pp. 584--606.
% 76 solutions
% 8*x1^2 + 8*x1*x2 + 8*x1*x3 + 2*x1*x4 + 2*x1*x5 + 2*x1*x6 + 2*x1*x7 - 8*x2*x3 - 2*x4*x7 - 2*x5*x6 - x1;
polysys{1,1} = [8 8 8 2 2 2 2 -8 -2 -2 -1];
polysys{1,2} = [2 zeros(1,7); 1 1 zeros(1,6); 1 0 1 zeros(1,5); 1 zeros(1,2) 1 zeros(1,4) ; 1 zeros(1,3) 1 zeros(1,3);1 zeros(1,4) 1 zeros(1,2); 1 zeros(1,5) 1 0; 0 1 1 zeros(1,5); zeros(1,3) 1 0 0 1 0;zeros(1,4) 1 1 0 0;1 zeros(1,7)]; 
%  8*x1*x2 - 8*x1*x3 + 8*x2^2 + 8*x2*x3 + 2*x2*x4 + 2*x2*x5 + 2*x2*x6 + 2*x2*x7 -2*x4*x6 - 2*x5*x7 - x2;
polysys{2,1} = [8 -8 8 8 2 2 2 2 -2 -2 -1];
polysys{2,2} = [1 1 zeros(1,6); 1 0 1 zeros(1,5);0 2 zeros(1,6); 0 1 1 zeros(1,5); 0 1 0 1 zeros(1,4); 0 1 0 0 1 zeros(1,3);0 1 zeros(1,3) 1 0 0;0 1 zeros(1,4) 1 0;zeros(1,3) 1 0 0 1 0;zeros(1,4) 1 0 1 0; 0 1 zeros(1,6)]; 
 % -8*x1*x2 + 8*x1*x3 + 8*x2*x3 + 8*x3^2 + 2*x3*x4 + 2*x3*x5 + 2*x3*x6 + 2*x3*x7 - 2*x4*x5 - 2*x6*x7 - x3;
polysys{3,1} = [-8 8 8 8 2 2 2 2 -2 -2 -1];
polysys{3,2} = [1 1 zeros(1,6); 1 0 1 zeros(1,5);0 1 1 zeros(1,5); 0 0 2 zeros(1,5); 0 0 1 1 zeros(1,4); 0 0 1 0 1 zeros(1,3);0 0 1 zeros(1,2) 1 0 0;0 0 1 zeros(1,3) 1 0;zeros(1,3) 1 1 0 0 0;zeros(1,5) 1 1 0; 0  0 1 zeros(1,5)];
%  2*x1*x4 - 2*x1*x7 + 2*x2*x4 - 2*x2*x6 + 2*x3*x4 - 2*x3*x5 + 8*x4^2 + 8*x4*x5 + 2*x4*x6 + 2*x4*x7 + 6*x4*x8 - 6*x5*x8 - x4;
polysys{4,1} = [2 -2 2 -2 2 -2 8 8 2 2 6 -6 -1];
polysys{4,2} = [1 0 0 1 zeros(1,4); 1 zeros(1,5) 1 0;0 1 0 1 zeros(1,4); 0 1 zeros(1,3) 1 0 0; 0 0 1 1 zeros(1,4); 0 0 1 0 1 zeros(1,3);0 0 0 2 zeros(1,4);0 0 0 1 1 zeros(1,3);zeros(1,3) 1 0 1 0 0;zeros(1,3) 1 0 0 1 0; zeros(1,3) 1 0 0 0 1;zeros(1,4) 1 0 0 1; 0 0 0 1 zeros(1,4)];
%  2*x1*x5 - 2*x1*x6 + 2*x2*x5 - 2*x2*x7 - 2*x3*x4 + 2*x3*x5 + 8*x4*x5 - 6*x4*x8 + 8*x5^2 + 2*x5*x6 + 2*x5*x7 + 6*x5*x8 -x5;
polysys{5,1} = [2 -2 2 -2 -2 2 8 -6 8 2 2 6 -1];
polysys{5,2} = [1 0 0 0 1 zeros(1,3); 1 zeros(1,4) 1 0 0;0 1 0 0 1 zeros(1,3); 0 1 zeros(1,4) 1 0; 0 0 1 1 zeros(1,4); 0 0 1 0 1 zeros(1,3);0 0 0 1 1 zeros(1,3);0 0 0 1 0 0 0 1;zeros(1,4) 2 0 0 0;zeros(1,4) 1 1 0 0; zeros(1,4) 1 0 1 0;zeros(1,4) 1 0 0 1; 0 0 0 0 1 zeros(1,3)];
% -2*x1*x5 + 2*x1*x6 - 2*x2*x4 + 2*x2*x6 + 2*x3*x6 - 2*x3*x7 + 2*x4*x6 + 2*x5*x6 + 8*x6^2 + 8*x6*x7 + 6*x6*x8 - 6*x7*x8 - x6;
polysys{6,1} = [-2 2 -2 2 2 -2 2 2 8 8 6 -6 -1];
polysys{6,2} = [1 0 0 0 1 zeros(1,3); 1 zeros(1,4) 1 0 0;0 1 0 1 zeros(1,4); 0 1 zeros(1,4) 1 0; 0 0 1 0 0 1 0 0; 0 0 1 0 0 0 1 0;0 0 0 1 0 1 zeros(1,2);zeros(1,4) 1 1 0 0;zeros(1,5) 2 0 0;zeros(1,5) 1 1 0; zeros(1,5) 1 0 1;zeros(1,6) 1 1; zeros(1,5) 1 0 0];
% -2*x1*x4 + 2*x1*x7 - 2*x2*x5 + 2*x2*x7 - 2*x3*x6 + 2*x3*x7 + 2*x4*x7 + 2*x5*x7 + 8*x6*x7 - 6*x6*x8 + 8*x7^2 + 6*x7*x8 -x7;
polysys{7,1} = [-2 2 -2 2 -2 2 2 2 8 -6 8 6 -1];
polysys{7,2} = [1 0 0 1 zeros(1,4); 1 zeros(1,5) 1 0;0 1 0 0 1 zeros(1,3); 0 1 zeros(1,4) 1 0; 0 0 1 0 0 1 0 0; 0 0 1 0 0 0 1 0;0 0 0 1 0 0 1 0;zeros(1,4) 1 0 1 0;zeros(1,5) 1 1 0;zeros(1,5) 1 0 1; zeros(1,6) 2 0;zeros(1,6) 1 1; zeros(1,6) 1 0];
% -6*x4*x5 + 6*x4*x8 + 6*x5*x8 - 6*x6*x7 + 6*x6*x8 + 6*x7*x8 + 8*x8^2 - x8;
polysys{8,1} = [-6 6 6 -6 6 6 8 -8];
polysys{8,2} = [0 0 0 1 1 zeros(1,3); zeros(1,3) 1 0 0 0 1;zeros(1,4) 1 0 0 1; zeros(1,5) 1 1 0; zeros(1,5) 1 0 1; zeros(1,6) 1 1;zeros(1,7) 2;zeros(1,7) 1];

% http://www.math.uic.edu/~jan/Demo/cohn3.html = nonzero-dimensional
% variety!! Groebner basis contains only 1 pure component
%
% The reduced polynomial system has norms which are beyond the range of
% double precision!!
%
% mixed volume : 213, 110 solutions
% [x y z t]
% -x**3*y**2 + 2*x**2*y**2*z - x**2*y*z**2 - 144*x**2*y**2 - 207*x**2*y*z + 288*x*y**2*z +
%    78*x*y*z**2 + x*z**3 - 3456*x**2*y - 5184*x*y**2 - 9504*x*y*z - 432*x*z**2
%    - 248832*x*y + 62208*x*z - 2985984*x;
polysys{1,1} = [-1 2 -1 -144 -207 288 78 1 -3456 -5184 -9504 -432 -248832 62208 -2985984];
polysys{1,2} = [3 2 0 0;2 2 1 0;2 1 2 0;2 2 0 0;2 1 1 0;1 2 1 0;1 1 2 0;1 0 3 0;2 1 0 0;1 2 0 0;1 1 1 0;1 0 2 0;1 1 0 0;1 0 1 0;1 0 0 0];
%  -x**3*z*t**2 + x**2*z**2*t**2 - 6*x**3*z*t + 4*x**2*z**2*t + 32*x**3*t**2 -
%    72*x**2*z*t**2 - 87*x*z**2*t**2 - z**3*t**2 - 8*x**3*z - 432*x**2*z*t - 414*x*z**2*t +
%    2592*x*z*t**2 + 864*z**2*t**2 - 1728*x**2*z - 20736*x*z*t + 3456*z**2*t
%    - 186624*z*t**2 - 124416*x*z - 1492992*z*t - 2985984*z;
polysys{2,1} = [-1 1 -6 4 32 -72 -87 -1 -8 -432 -414 2592 864 -1728 -20736 3456 -186624 -124416 -1492992 -2985984];
polysys{2,2} = [3 0 1 2;2 0 2 2;3 0 1 1;2 0 2 1;3 0 0 2;2 0 1 2;1 0 2 2;0 0 3 2;3 0 1 0;2 0 1 1;1 0 2 1; 1 0 1 2;0 0 2 2;2 0 1 0;1 0 1 1;0 0 2 1;0 0 1 2;1 0 1 0;0 0 1 1;0 0 1 0]; 
%  x**2*y*t**3 - 2*x*y**2*t**3 + y**3*t**3 + 8*x**2*y*t**2 - 12*x*y**2*t**2 + 4*y**3*t**2 -
%    24*x*y*t**3 + 24*y**2*t**3 + 20*x**2*y*t - 20*x*y**2*t - 160*x*y*t**2 + 96*y**2*t**2 +
%    128*x*t**3 + 16*x**2*y + 96*x*y*t + 2304*x*t**2 + 1152*x*y + 13824*x*t + 27648*x;
polysys{3,1} = [1 -2 1 8 -12 4 -24 24 20 -20 -160 96 128 16 96 2304 1152 13824 27648];
polysys{3,2} = [2 1 0 3;1 2 0 3;0 3 0 3;2 1 0 2;1 2 0 2;0 3 0 2;1 1 0 3;0 2 0 3;2 1 0 1;1 2 0 1;1 1 0 2;0 2 0 2;1 0 0 3;2 1 0 0;1 1 0 1;1 0 0 2;1 1 0 0;1 0 0 1;1 0 0 0];
%  y**3*t**3 - y**2*z*t**3 + 4*y**3*t**2 - 2*y**2*z*t**2 + 72*y**2*t**3 + 71*y*z*t**3 +
%    288*y**2*t**2 + 360*y*z*t**2 + 6*z**2*t**2 + 1728*y*t**3 - 464*z*t**3 + 432*y*z*t
%    + 8*z**2*t + 6912*y*t**2 - 4320*z*t**2 + 13824*t**3 + z**2 - 13824*z*t + 55296*t**2
%    - 13824*z;
polysys{4,1} = [1 -1 4 -2 72 71 288 360 6 1728 -464 432 8 6912 -4320 13824 1 -13824 55296 -13824];
polysys{4,2} = [0 3 0 3;0 2 1 3;0 3 0 2;0 2 1 2;0 2 0 3;0 1 1 3;0 2 0 2;0 1 1 2;0 0 2 2;0 1 0 3;0 0 1 3; 0 1 1 1;0 0 2 1;0 1 0 2;0 0 1 2;0 0 0 3;0 0 2 0;0 0 1 1;0 0 0 2;0 0 1 0];

strings{1}='-a[1]^3*a[2]^2+2*a[1]^2*a[2]^2*a[3]-a[1]^2*a[2]*a[3]^2-144*a[1]^2*a[2]^2-207*a[1]^2*a[2]*a[3]+288*a[1]*a[2]^2*a[3]+78*a[1]*a[2]*a[3]^2+a[1]*a[3]^3-3456*a[1]^2*a[2]-5184*a[1]*a[2]^2-9504*a[1]*a[2]*a[3]-432*a[1]*a[3]^2-248832*a[1]*a[2]+62208*a[1]*a[3]-2985984*a[1]';
strings{2}='-a[1]^3*a[3]*a[4]^2+a[1]^2*a[3]^2*a[4]^2-6*a[1]^3*a[3]*a[4]+4*a[1]^2*a[3]^2*a[4]+32*a[1]^3*a[4]^2-72*a[1]^2*a[3]*a[4]^2-87*a[1]*a[3]^2*a[4]^2-a[3]^3*a[4]^2-8*a[1]^3*a[3]-432*a[1]^2*a[3]*a[4]-414*a[1]*a[3]^2*a[4]+2592*a[1]*a[3]*a[4]^2+864*a[3]^2*a[4]^2-1728*a[1]^2*a[3]-20736*a[1]*a[3]*a[4]+3456*a[3]^2*a[4]-186624*a[3]*a[4]^2-124416*a[1]*a[3]-1492992*a[3]*a[4]-2985984*a[3]';
strings{3}='a[1]^2*a[2]*a[4]^3-2*a[1]*a[2]^2*a[4]^3+a[2]^3*a[4]^3+8*a[1]^2*a[2]*a[4]^2-12*a[1]*a[2]^2*a[4]^2+4*a[2]^3*a[4]^2-24*a[1]*a[2]*a[4]^3+24*a[2]^2*a[4]^3+20*a[1]^2*a[2]*a[4]-20*a[1]*a[2]^2*a[4]-160*a[1]*a[2]*a[4]^2+96*a[2]^2*a[4]^2+128*a[1]*a[4]^3+16*a[1]^2*a[2]+96*a[1]*a[2]*a[4]+2304*a[1]*a[4]^2+1152*a[1]*a[2]+13824*a[1]*a[4]+27648*a[1]';
strings{4}='a[2]^3*a[4]^3-a[2]^2*a[3]*a[4]^3+4*a[2]^3*a[4]^2-2*a[2]^2*a[3]*a[4]^2+72*a[2]^2*a[4]^3+71*a[2]*a[3]*a[4]^3+288*a[2]^2*a[4]^2+360*a[2]*a[3]*a[4]^2+6*a[3]^2*a[4]^2+1728*a[2]*a[4]^3-464*a[3]*a[4]^3+432*a[2]*a[3]*a[4]+8*a[3]^2*a[4]+6912*a[2]*a[4]^2-4320*a[3]*a[4]^2+13824*a[4]^3+a[3]^2-13824*a[3]*a[4]+55296*a[4]^2-13824*a[3]';
polysys = lti2polysys(strings,[],[],'na',4);
