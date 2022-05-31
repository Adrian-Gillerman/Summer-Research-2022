classdef InputOutput
  %% InputOutput class
  % lw: line width
  % fs: font size
  % fr: frame rate
  % sg: spatial grid object
  % p: "nopause" if no pauses and "pause" if pauses
  % plottype: an array where the first element should be
  %            "IC" for just the initial condition
  %            "startfinish" for the initial conditon and final approx
  %            "movie" for every time step
  %            "nothing" for no plotting
  %            and the second element is either
  %            "solution" for the approximation or
  %            "error" for the error to the approximation
  % iOption:  determines solution, 
  %           1 = sinusoid, 2 = TZ quadratic in space and time
  
   
  properties
    lw = 3;
    fs = 16;
    fr = 10;
    sg
    bc
    p = "pause"
    plottype1
    plottype2
    iOption
  end
    
  methods
    function obj = InputOutput( sg,bc,plottype,iOption )
      obj.sg = sg;
      obj.bc = bc;
      if ismember(plottype(1), ["IC", "startfinish", "movie", "nothing"])
        obj.plottype1 = plottype(1);
      else
        error('unsupported plottype : in InputOutput')
      end
      if ismember(plottype(2), ["solution", "error"])
        obj.plottype2 = plottype(2);
      else
        error('unsupported plottype : in InputOutput')
      end 
      obj.iOption = iOption;
    end
    function plotIC( obj,u,uex )
      if ismember( obj.plottype1, ["IC", "startfinish", "movie"] )
        clf
        cplot( obj,u,uex )
      end
    end
    function movie( obj,n,u,uex )
      if (obj.plottype1 == "movie" && mod(n, obj.fr) == 0) 
        cplot( obj,u,uex )
      end
    end
    function plotfinal( obj,u,uex )
      if obj.plottype1 == "startfinish"
        cplot(obj,u,uex )
      end
    end
    function cplot( obj,u,uex )
      figure(1)
      xlabel( 'x' )
      if obj.plottype2 == "solution"
        plot( obj.sg.ig,u,'kx','LineWidth',obj.lw );
        ylabel( 'u' )
        if obj.iOption == 1
          axis( [obj.sg.xa,obj.sg.xb,-1.1,1.1] );
        end
      elseif obj.plottype2 == "error"
        plot( obj.sg.ig,u-uex,'kx','LineWidth',obj.lw );
        ylabel( 'error' )
      end
      set(gca,'FontSize',obj.fs);
      drawnow
      if obj.p == "pause"
        pause
      end
    end
    function output( obj,errarray )
      csvname = "errormatrix" + obj.bc.BCstr + ".csv";
      writematrix(errarray,csvname) 
      fprintf('Saving file=[%s]\n', csvname);
    end
    function checkcorrect( obj,errarray )
      csvname = "errormatrix" + obj.bc.BCstr + ".csv";
      correcterrors = readmatrix(csvname);
      fprintf('Checking whether we have the same errors...\n')
      errordif = max(abs(errarray - correcterrors));
      if errordif < 1e-13
        fprintf('Great job! You have the correct errors!\n')
      else
        fprintf('Unfortunately the errors are not correct.\n')
        fprintf('The error is off by [%e]. \n', errordif )
      end
    end
  end
end
