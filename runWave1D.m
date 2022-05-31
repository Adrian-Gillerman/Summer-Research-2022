clear
clc

%%%%%
%% runs the wave1D.m file with the following parameters
%%%%%
N         = 320;
% Note that N = 640 will make 1/dx^4 to big for the code to be accurate.
lambda0   = 0.9;
tf        = 1.;
cOption   = 2;
iOption   = 3;
plottype  = ["nothing", "error"];
errarray  = zeros(1, 2);

for Ord = [2,4]
  [x,u,bc]        = wave1D( Ord,N,lambda0,tf,cOption,iOption,plottype );
  ex              = Exact(cOption, iOption);
  uexact          = ex.uex( x,tf );
  err             = max(abs(u-uexact));
  % The following line prints the errors for the given order.
  % fprintf( 'Ord=%i, N=%i, error=%e\n', Ord,N, err );
  errarray(Ord/2) = err;
end

sg = Spatial_grid( Ord,N );
iO = InputOutput( sg,bc,plottype,iOption );
%% The following code makes a new error array file
% iO.output(errarray)

%% Check whether the errors are the correct ones
if (iOption == 1 && cOption == 1)
    iO.checkcorrect(errarray)
end
