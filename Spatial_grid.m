classdef Spatial_grid
    % spatial grid class
    % Ord:    order of the scheme
    % xa:     x-value of the left boundary
    % xb:     x-value of the right boundary
    % Nx:     number of interior cells
    % dx:     size of the spatial step
    % ng:     number of ghost points on each side
    % NXT:    total number of spatial points in the one-dim grid
    % ja:     index of the first non-ghost point
    % jb:     index of the last non-ghost point
    % jrange: array of indices that are not ghost points
    % xam:    x-value of the left-most point
    % xbp:    x-value of the right-most point
    % x:      array of grid points
    % ig:     array of interior (non-ghost) grid points
    
    properties
        Ord
        xa = 0
        xb = 1
        Nx
        dx
        ng  
        NXT    
        ja    
        jb
        jrange
        xam
        xbp
        x
        ig
    end
    
    methods
        function obj = Spatial_grid( Ord,Nx )
            % Spatial_grid: Construct an instance of this class
            % We set the following parameters
            obj.Ord    = Ord;
            obj.Nx     = Nx;
            obj.dx     = (obj.xb-obj.xa)/obj.Nx;
            if ismember(Ord, [2, 4])
              obj.ng   = Ord/2;
            else
              error('unsupported order : in Spatial_grid');
            end
            obj.NXT    = obj.Nx + 1 + 2 * obj.ng;
            obj.ja     = obj.ng + 1;
            obj.jb     = obj.NXT - obj.ng;
            obj.jrange = obj.ja:obj.jb;
            obj.xam    = obj.xa - obj.ng * obj.dx;
            obj.xbp    = obj.xb + obj.ng * obj.dx;
            obj.x      = linspace( obj.xam, obj.xbp, obj.NXT );
            obj.ig     = obj.x(obj.jrange);
        end
    end
end

