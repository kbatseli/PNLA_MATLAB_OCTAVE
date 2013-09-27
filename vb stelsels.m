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

% stelsel polynomial optimization problems are eigenvalue problems
polysys{1,1} = [3 2 -1 -1 1];
polysys{1,2} = [0 0; 0 1; 1 1; 0 2; 2 0];
polysys{2,1} = [5 4 3 1 -2 1];
polysys{2,2} = [0 0; 0 1; 1 0; 1 1; 0 2; 2 0];

% problematisch stelsel, werkt
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
%
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
% has 2 solutions: (0,0) and (1,2), taking radical therefore also
% eliminates solutions at inf for this example. Is this a general thing?
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

% eenvoudig stelsel, werkt
% (x-2) * (x-3)
polysys{1,1} = [1 -5 6];
polysys{1,2} = [2 0;1 0;0 0];
% y + 5x - 7
polysys{2,1} = [1 5 -7];
polysys{2,2} = [0 1;1 0;0 0];

% multipliciteit root op oneindig van 2
polysys{1,1} = [1 -7 10];polysys{1,2} = [2 1;1 1;0 1];
polysys{2,1} = [1 -2];polysys{2,2} = [1 1;1 0];

% eenvoudig stelsel, meervoudige wortels x = 2 en x = 3, y = 5 en y = 4
% (x-2) * (x-3) * y
polysys{1,1} = [1 -5 6];
polysys{1,2} = [2 1;1 1;0 1];
% (y-4) * (y-5)
polysys{2,1} = [1 -9 20];
polysys{2,2} = [0 2;0 1;0 0];

% eenvoudig stelsel, x = 2 
% (x-2) * y
polysys{1,1} = [1 -2];
polysys{1,2} = [1 1;0 1];
% (y-3) *
polysys{2,1} = [1 -3];
polysys{2,2} = [0 1;0 0];

% eenvoudig stelsel, x = 2, y=3, root at infinity (0,1,0) with multiplicity of 2 
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

% eigenvector van dit probleem is [1 x x^2], moet dus slimmere manier
% vinden om lin. independent rows te kiezen. suggestie: check eerst alle
% zuivere componenten
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

% y weggedeeld uit f2 en f3 van hierboven, niet-triviale syzygies
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

% I^h is niet gelijk aan ideaal met basis f_i^h, Cox p 387
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

% geen zero-dimensional variety!, inf oplossingen
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

% Ngai's toy problem
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

% Ngai's toy problem 2
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

% Ngai's toy problem deg(p)=1
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

% Ngai's simple dynamical system
dynsys{1,1}=[1 -1 -1 4];dynsys{1,2}=[0 1;3 0;1 2;1 0];
dynsys{2,1}=[-1 -1 -1 4];dynsys{2,2}=[1 0;2 1;0 3;0 1];

% Ngai's simple dynamical system, different degrees
dynsys{1,1}=[1 -1 -1 4];dynsys{1,2}=[0 1;3 0;1 2;1 0];
dynsys{2,1}=[-1 -1 -1 4];dynsys{2,2}=[1 0;2 1;0 4;0 1];


% Ngai's Van der Pol oscillator
% [omega c1 c3i cr3]
strings{1}='-a[1]^2-3*a[2]*a[3]*a[1]+1';
strings{2}='6*a[4]^2+3*a[2]*a[4]+6*a[3]^2+3*a[2]^2-1';
strings{3}='-9*a[4]*a[1]^2-9*a[3]*a[4]^2*a[1]-9*a[3]^3*a[1]-18*a[2]^2*a[3]*a[1]+3*a[3]*a[1]+a[4]';
strings{4}='9*a[4]^3*a[1]-9*a[3]*a[1]^2+9*a[3]^2*a[4]*a[1]+18*a[2]^2*a[4]*a[1]-3*a[4]*a[1]+3*a[2]^3*a[1]+a[3]';
polysys=lti2polysys(strings,[],[],'na',4);

% Ngai's Van der Pol oscillator
% Haotian's equations, N=3, has 2D variety
strings{1}='-9*a[4]^3*a[1]-18*a[4]*a[2]^2*a[1]-9*a[4]*a[3]^2*a[1]+3*a[4]*a[1]-9*a[3]*a[1]^2+a[3]';
strings{2}='-3*a[4]*a[2]^2*a[1]-a[2]*a[1]^2+a[2]';
strings{3}='-9*a[4]^2*a[3]*a[1]+9*a[4]*a[1]^2-a[4]-3*a[2]^3*a[1]-18*a[2]^2*a[3]*a[1]-9*a[3]^3*a[1]+3*a[3]*a[1]';
strings{4}='-6*a[1]*a[4]^2*a[2]-3*a[1]*a[2]^3-3*a[1]*a[2]^2*a[3]-6*a[1]*a[2]*a[3]^2+a[1]*a[2]';
polysys=lti2polysys(strings,[],[],'na',4);

% Ngai's Van der Pol oscillator
% Haotian's equations, N=3, simplified
strings{1}='-9*a[4]^3*a[1]-18*a[4]*a[2]^2*a[1]-9*a[4]*a[3]^2*a[1]+3*a[4]*a[1]-9*a[3]*a[1]^2+a[3]';
strings{2}='3*a[1]*a[2]*a[4]+a[1]^2-1';
strings{3}='-9*a[4]^2*a[3]*a[1]+9*a[4]*a[1]^2-a[4]-3*a[2]^3*a[1]-18*a[2]^2*a[3]*a[1]-9*a[3]^3*a[1]+3*a[3]*a[1]';
strings{4}='6*a[4]^2+3*a[2]^2+3*a[2]*a[3]+6*a[3]^2-1';
polysys=lti2polysys(strings,[],[],'na',4);

% Ngai's Van der Pol oscillator
% Haotian's equations, N=5, simplified
strings{1}='-30*a[1]*a[4]^2*a[6]-15*a[1]*a[2]^2*a[4]-30*a[1]*a[2]*a[3]*a[4]-15*a[1]*a[6]^3-30*a[1]*a[2]^2*a[6]-30*a[1]*a[3]^2*a[6]-15*a[1]*a[5]^2*a[6]+5*a[1]*a[6]-25*a[1]^2*a[5]+a[5]';
strings{2}='-9*a[1]*a[4]^3-18*a[1]*a[4]*a[6]^2-18*a[1]*a[2]^2*a[4]+18*a[1]*a[2]*a[4]*a[5]-9*a[1]*a[3]^2*a[4]-18*a[1]*a[4]*a[5]^2+3*a[1]*a[4]-9*a[1]*a[2]^2*a[6]-18*a[1]*a[2]*a[3]*a[6]-9*a[1]^2*a[3]+a[3]';
strings{3}='-3*a[1]*a[4]^2*a[6]-3*a[1]*a[2]^2*a[4]+6*a[1]*a[2]*a[4]*a[5]-6*a[1]*a[3]*a[4]*a[5]-6*a[1]*a[2]*a[3]*a[6]-a[1]^2*a[2]+a[2]+3*a[1]*a[3]^2*a[6]';
strings{4}='15*a[1]*a[2]*a[4]^2-30*a[1]*a[4]^2*a[5]-15*a[1]*a[5]*a[6]^2+25*a[1]^2*a[6]-a[6]-15*a[1]*a[2]^2*a[3]-30*a[1]*a[2]^2*a[5]-15*a[1]*a[2]*a[3]^2-30*a[1]*a[3]^2*a[5]-15*a[1]*a[5]^3+5*a[1]*a[5]';
strings{5}='-9*a[1]*a[3]*a[4]^2-18*a[1]*a[2]*a[4]*a[6]+9*a[1]^2*a[4]-a[4]-18*a[1]*a[3]*a[6]^2-3*a[1]*a[2]^3-18*a[1]*a[2]^2*a[3]-9*a[1]*a[2]^2*a[5]-18*a[1]*a[2]*a[3]*a[5]-9*a[1]*a[3]^3-18*a[1]*a[3]*a[5]^2+3*a[1]*a[3]';
strings{6}='6*a[2]*a[4]^2-3*a[4]^2*a[5]+6*a[2]*a[4]*a[6]+6*a[3]*a[4]*a[6]+6*a[2]*a[6]^2+3*a[2]^3+3*a[2]^2*a[3]+6*a[2]*a[3]^2+6*a[2]*a[3]*a[5]+6*a[2]*a[5]^2-a[2]+3*a[3]^2*a[5]';
polysys=lti2polysys(strings,[],[],'na',6);

% Ngai's Van der Pol oscillator
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

%% Numerical Problems = check doesn't check out

% Cassou
strings{1}='15*a[1]^4*a[2]*a[3]^2+6*a[1]^4*a[2]^3+21*a[1]^4*a[2]^2*a[3]-144*a[1]^2*a[2]-8*a[1]^2*a[2]^2*a[4]-28*a[1]^2*a[2]*a[3]*a[4]-648*a[1]^2*a[3]+36*a[1]^2*a[3]^2*a[4]+9*a[1]^4*a[3]^3-120';
strings{2}='30*a[2]^3*a[1]^4*a[3]-32*a[3]*a[4]^2*a[2]-720*a[3]*a[1]^2*a[2]-24*a[2]^3*a[1]^2*a[4]-432*a[2]^2*a[1]^2+576*a[4]*a[2]-576*a[3]*a[4]+16*a[2]*a[1]^2*a[3]^2*a[4]+16*a[3]^2*a[4]^2+16*a[4]^2*a[2]^2+9*a[2]^4*a[1]^4+5184+39*a[3]^2*a[1]^4*a[2]^2+18*a[3]^3*a[1]^4*a[2]-432*a[3]^2*a[1]^2+24*a[3]^3*a[1]^2*a[4]-16*a[2]^2*a[1]^2*a[3]*a[4]-240*a[2]';
strings{3}='216*a[3]*a[1]^2*a[2]-162*a[3]^2*a[1]^2-81*a[2]^2*a[1]^2+5184+1008*a[4]*a[2]-1008*a[3]*a[4]+15*a[2]^2*a[1]^2*a[3]*a[4]-15*a[2]^3*a[1]^2*a[4]-80*a[3]*a[4]^2*a[2]+40*a[3]^2*a[4]^2+40*a[4]^2*a[2]^2';
strings{4}='261+4*a[3]*a[1]^2*a[2]-3*a[3]^2*a[1]^2-4*a[2]^2*a[1]^2+22*a[4]*a[2]-22*a[3]*a[4]';
polysys=lti2polysys(strings,[],[],'na',4);

% STLS van Philippe, dreg = 21, 18 affine solutions op d=23
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
