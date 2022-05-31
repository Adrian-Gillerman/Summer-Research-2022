clear
clf
clc

%%%%%
%% convergence study
%%%%%
m         = 6;
N0        = 20;
lambda0   = 0.9;
tf        = 1.;
cOption   = 2;
iOption   = 3;
plottype  = ["nothing", "error"];
errarray  = zeros(1, m);
errmatrix = zeros(2, m);
h         = zeros(1, m);

for Ord = [2,4]
  for j = 0:m-1
    N             = N0*2^j;
    [x,u,bc]      = wave1D( Ord,N,lambda0,tf,cOption,iOption,plottype );
    ex            = Exact(cOption, iOption);
    uexact        = ex.uex( x,tf );
    err           = max(abs(u-uexact));
    errarray(j+1) = err;
    h(j+1)        = x(2) - x(1);
    % The following line prints the errors for the given number of grid
    % points.
    fprintf( 'N=%i, error=%e\n', N, err );
  end
  errmatrix(Ord/2, :) = errarray;
  convergencestudy( Ord, h, errarray)
end

legend( '2nd','h^2 ref','4th','h^4 ref','Location','SouthEast' )
title( sprintf( 'Convergence Study' ) )
plotName = sprintf( 'images/oneConv.eps' );
fprintf('Saving file=[%s]\n',plotName)
print('-depsc2',plotName)

sg = Spatial_grid( Ord,N );
iO = InputOutput( sg,bc,plottype,iOption );
%% Make a new error matrix file
% iO.output(errmatrix)

%% Check whether the errors are the correct ones
% iO.checkcorrect(errmatrix)
