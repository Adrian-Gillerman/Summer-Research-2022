classdef Exact
    %Exact Provides exact solutions and its derivatives 
    %to the 1D wave equation
    %   Can calculate u,ux,uxx,uxxx,uxxxx,ut,utt,uttt,utttt,uxxt,uxtt,uxxtt
    %   as well as the forcing functions f,fx,fxx,ft,ftt from the exact
    %   solutions
    
    properties
        c
        iOption
        p
        pder
        p2der
        p3der
        p4der
    end
    
    methods
        function obj = Exact( c,iOption )
          %Exact Construct an instance of this class
          obj.c = c;
          obj.iOption = iOption;
          if floor(iOption) == iOption
            obj.p     = ones(1, iOption + 1);
            obj.pder  = polyder(obj.p);
            obj.p2der = polyder(obj.pder);
            obj.p3der = polyder(obj.p2der);
            obj.p4der = polyder(obj.p3der);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function uex = uex( obj,x,t )
        %   Gives exact solution of the 1D wave equation
          if obj.iOption == 1
            uex = sin(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uex = polyval(obj.p, x)*polyval(obj.p, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function f = f( obj,x,t )
          f = utt( obj,x,t )- obj.c^2*uxx( obj,x,t );
        end
        
        function fx = fx( obj,x,t )
          fx = uxtt( obj,x,t )- obj.c^2*uxxx( obj,x,t );
        end
        
        function fxx = fxx( obj,x,t )
          fxx = uxxtt( obj,x,t )- obj.c^2*uxxxx( obj,x,t );
        end
        
        function ft = ft( obj,x,t )
          ft = uttt( obj,x,t )- obj.c^2*uxxt( obj,x,t );
        end

        function ftt = ftt( obj,x,t )
          ftt = utttt( obj,x,t )- obj.c^2*uxxtt( obj,x,t );
        end
        
        function ux = ux( obj,x,t )
          if obj.iOption == 1
            ux = 5*cos(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            ux = polyval(obj.pder, x)*polyval(obj.p, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function uxx = uxx( obj,x,t )
          if obj.iOption == 1
            uxx = -(5)^2*sin(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uxx = polyval(obj.p2der, x)*polyval(obj.p, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function uxxx = uxxx( obj,x,t )
          if obj.iOption == 1
            uxxx = -(5)^3*cos(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uxxx = polyval(obj.p3der, x)*polyval(obj.p, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
 
        function uxxxx = uxxxx( obj,x,t )
          if obj.iOption == 1
            uxxxx = (5)^4*sin(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uxxxx = polyval(obj.p4der, x)*polyval(obj.p, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function ut = ut( obj,x,t )
          if(     obj.iOption == 1 )
            ut = -5*obj.c*cos(5*(x - obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            ut = polyval(obj.p, x)*polyval(obj.pder, t);
          else
            error('unsupported iOption : in Exact')
          end 
        end
 
        function utt = utt( obj,x,t )
          if obj.iOption == 1
            utt = -(-5*obj.c)^2*sin(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            utt = polyval(obj.p, x)*polyval(obj.p2der, t);
          else
            error('unsupported iOption : in Exact')
          end  
        end
        
        function uttt = uttt( obj,x,t )
          if obj.iOption == 1
            uttt = -(-5*obj.c)^3*cos(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uttt = polyval(obj.p, x)*polyval(obj.p3der, t);
          else
            error('unsupported iOption : in Exact')
          end
        end

        function utttt = utttt( obj,x,t )
          if obj.iOption == 1
            utttt = (-5*obj.c)^4*sin(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            utttt = polyval(obj.p, x)*polyval(obj.p4der, t);
          else
            error('unsupported iOption : in Exact')
          end 
        end
        
        function uxxt = uxxt( obj,x,t )
          if obj.iOption == 1
            uxxt = -(5)^2*(-5*obj.c)*cos(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uxxt = polyval(obj.p2der, x)*polyval(obj.pder, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function uxtt = uxtt( obj,x,t )
          if obj.iOption == 1
            uxtt = -(5)*(-5*obj.c)^2*cos(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uxtt = polyval(obj.pder, x)*polyval(obj.p2der, t);
          else
            error('unsupported iOption : in Exact')
          end
        end
        
        function uxxtt = uxxtt( obj,x,t )
          if obj.iOption == 1
            uxxtt = (5)^2*(-5*obj.c)^2*sin(5*(x'-obj.c*t));
          elseif ismember(obj.iOption, 2:5)
            uxxtt = polyval(obj.p2der, x)*polyval(obj.p2der, t);
          else
            error('unsupported iOption : in Exact')
          end
        end   
    end 
end

