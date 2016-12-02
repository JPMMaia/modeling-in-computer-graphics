classdef UVPlot < handle
    properties
        figure
        mesh
        axes
        
        cline_plot
        cline_spacing
        unit_square_plot
        uv_plot
        selection_plot
        request_plot
        fixed_plot
        
        is_panning
        panning_point
        
        selection
        
        point_requested
        request_fun
    end
    
    methods
        function obj = UVPlot(figure, position)
            obj.figure = figure;
            obj.mesh = [];
            obj.cline_spacing = 0.05;
            
            obj.axes = axes('Parent',obj.figure,'Units','pixels',...
                'Position',position,'Units','normalized',...
                'DataAspectRatio',[1 1 1]);
            obj.cline_plot = patch('Vertices',[],'Faces',[],'Parent',obj.axes,...
                'EdgeColor','flat');
            obj.unit_square_plot = patch('Vertices',[0 0;1 0;1 1;0 1],...
                'Faces',[1 2;2 3;3 4;4 1],'EdgeColor',[0.4 0.4 0.4],...
                'LineWidth', 1.5, 'Parent',obj.axes);
            obj.uv_plot = patch('Vertices',[],'Faces',[],...
                'Parent',obj.axes,'EdgeColor',[0.2 0.2 0.2]);
            obj.request_plot = patch('Vertices',[],'Faces',[],...
                'Parent',obj.axes,'MarkerSize',6,'Marker','o',...
                'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.2 0.2]);
            obj.fixed_plot = patch('Vertices',[],'Faces',[],...
                'Parent',obj.axes,'MarkerSize',6,'Marker','o',...
                'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.2 0.2]);
            obj.selection_plot = patch('Vertices',[],'Faces',[],...
                'Parent',obj.axes,'MarkerSize',6,'Marker','o',...
                'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.5 0]);
            obj.selection = [];
            
            xlabel('u');
            ylabel('v');
            hold on;
            axis([-0.1 1.1 -0.1 1.1]);
            
            pan_obj = pan;
            setAllowAxesPan(pan_obj, obj.axes, false);
            rot3d_obj = rotate3d;
            setAllowAxesRotate(rot3d_obj, obj.axes, false);
            zoom_obj = zoom;
            setAllowAxesZoom(zoom_obj, obj.axes, false);
            
            obj.updateClinePlot();
            
            obj.point_requested = false;
        end
        
        function setFixedPoints(obj, vpos)
            obj.fixed_plot.Vertices = vpos;
            obj.fixed_plot.Faces = (1:size(vpos,1))';
        end
        
        function ret = getPosition(obj)
            obj.axes.Units = 'pixels';
            ret = obj.axes.Position;
            obj.axes.Units = 'normalized';
        end
        
        function requestPoint(obj, fcn)
            obj.point_requested = true;
            obj.request_fun = fcn;
        end
        
        function setClineSpacing(obj, val)
            obj.cline_spacing = val;
            obj.updateClinePlot();
        end
        
        function updateUVCoordinates(obj)
            return_early = false;
            if isempty(obj.mesh)
                return_early = true;
            end
            if ~isfield(obj.mesh.V_traits, 'uv')
                return_early = true;
            end
            if return_early
                obj.uv_plot.Vertices = [];
                obj.uv_plot.Faces = [];
                return;
            end
            he = obj.mesh.getAllEdges().halfedge();
            E = [he.from().index he.to().index];
            obj.uv_plot.Vertices = obj.mesh.getAllVertices().getTrait('uv');
            obj.uv_plot.Faces = E;
            obj.updateSelectionPlot();
        end
        
        function setMesh(obj, mesh)
            obj.mesh = mesh;
            obj.clearSelection();
            obj.updateUVCoordinates();
        end
        
        function buttonDown(obj, source, data)
            switch source.SelectionType
                case 'normal'
                    obj.is_panning = true;
                    obj.panning_point = obj.getCurrentPoint();
                case 'alt'
                    if obj.point_requested
                        obj.request_fun(obj.getCurrentPoint());
                        obj.point_requested = false;
                        obj.request_plot.Vertices = [];
                        obj.request_plot.Faces = [];
                    elseif ~obj.is_panning
                        obj.selectAt(obj.getCurrentPoint());
                    end
                otherwise
            end
        end
        
        function clearSelection(obj)
            obj.selection = [];
            obj.selection_plot.Vertices = [];
            obj.selection_plot.Faces = [];
        end
        
        function selectAt(obj, p)
            if isempty(obj.mesh)
                return;
            end
            if ~isfield(obj.mesh.V_traits, 'uv')
                return;
            end
            
            uv = obj.mesh.getAllVertices().getTrait('uv');
            dif = bsxfun(@minus, uv, p);
            [~,seli] = min(sum(dif.^2,2));
            obj.selection = seli;
            obj.updateSelectionPlot();
        end
        
        function select(obj, i)
            obj.selection = i;
            obj.updateSelectionPlot();
        end
        
        function updateSelectionPlot(obj)
            obj.selection_plot.Vertices = obj.mesh.getVertex(obj.selection).getTrait('uv');
            obj.selection_plot.Faces = 1;
        end
        
        function buttonUp(obj, source, data)
            obj.is_panning = false;
        end
        
        function buttonMotion(obj, source, data)
            if obj.is_panning
                p = obj.getCurrentPoint();
                delta = obj.panning_point - p;
                obj.axes.XLim = obj.axes.XLim+delta(1);
                obj.axes.YLim = obj.axes.YLim+delta(2);
                obj.updateClinePlot();
            end
            
            if obj.point_requested
                p = obj.getCurrentPoint();
                obj.request_plot.Vertices = p;
                obj.request_plot.Faces = 1;
            end
        end
        
        function scrollWheel(obj, source, data)
            p = obj.getCurrentPoint();
            if data.VerticalScrollCount > 0
                factor = 1/0.9;
            else
                factor = 0.9;
            end
            
            xlim = obj.axes.XLim;
            ylim = obj.axes.YLim;
            
            obj.axes.XLim = [p(1) - (p(1)-xlim(1))*factor ...
                p(1) + (xlim(2)-p(1))*factor];
            obj.axes.YLim = [p(2) - (p(2)-ylim(1))*factor ...
                p(2) + (ylim(2)-p(2))*factor];
            
            obj.updateClinePlot();
        end
        
        function ret = getCurrentPoint(obj)
            ret = obj.axes.CurrentPoint(1,[1 2]);
        end
        
        function updateClinePlot(obj)
            xlim = obj.axes.XLim;
            ylim = obj.axes.YLim;
            
            xstart = ceil(xlim(1)/obj.cline_spacing);
            xend = floor(xlim(2)/obj.cline_spacing);
            if xend-xstart < 200
                xcoord = (xstart:xend)' * obj.cline_spacing;
            else
                xcoord = zeros(0,2);
            end
            
            ystart = ceil(ylim(1)/obj.cline_spacing);
            yend = floor(ylim(2)/obj.cline_spacing);
            if yend-ystart < 200
                ycoord = (ystart:yend)' * obj.cline_spacing;
            else
                ycoord = zeros(0,2);
            end
            
            V1 = [kron(xcoord,[1;1]) repmat(ylim',size(xcoord,1),1)];
            V2 = [repmat(xlim',size(ycoord,1),1) kron(ycoord,[1;1])];
            V = [V1;V2];
            E = reshape(1:size(V,1),2,[])';
            C = [repmat([0.8 0.5 0.5],size(V1,1),1); repmat([0.5 0.5 0.8],size(V2,1),1)];
            
            obj.cline_plot.Vertices = V;
            obj.cline_plot.Faces = E;
            obj.cline_plot.FaceVertexCData = C;
        end
    end
end