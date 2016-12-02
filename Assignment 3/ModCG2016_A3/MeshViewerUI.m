classdef MeshViewerUI < handle
    properties(Access=private)
        fig
        axes
        lights
        uv_plot
        
        model_stats_text
        preserve_volume_ui
        wl_ui
        wp_ui
        eigenfunctions_ui
        shading_type_group
        fixed_list
        fixed_vi
        fixed_vpos
        add_fixed_button
        
        model_patch
        model_vertex_patch
        model_vn_patch
        model_fn_patch
        model_cline_patch_x
        model_cline_patch_y
        grid_patch
        boundary_patch
        bounding_box_patch
        faces_visible
        
        normal_scale
        vertex_normal_weighting
        mesh
        last_file
        
        L_uniform
        L_Cotangent
        L_MeanValue
        
        L_uniform_dirty
        L_Cotangent_dirty
        L_MeanValue_dirty
        
        basic_lambda
        recompute_every_iteration
        laplacian_weighting
        initial_mesh_volume
        initial_vertex_positions
        preserve_volume
        num_iterations
        lsq_weighting
        num_eigenfunctions
        laplacian_normalized
        quantile_range
        
        bdry_shape
        adaptive_spacing
    end
    
    properties(Access=private, Constant)
        edge_color = [0 0 0]
        vertex_color = [0 0 0]
        face_color = [0.6 0.6 1]
        bb_color = [0.4 0.4 0.4];
        grid_color = [0.7 0.7 0.7];
        boundary_color = [1 0.2 0.2];
        face_normal_color = [1 0 0.5];
        vertex_normal_color = [1 0.5 0];
        
        align_view_text = {'+X','-X','+Y','-Y','+Z','-Z'};
        align_view_angles = [-90 0;90 0;0 0;180 0;0 -90;0 90];
    end
    
    methods
        function obj = MeshViewerUI()
            obj.fig = figure('Name','ModCG Mesh Viewer',...
                'Visible','off','Position',[360,500,1600,800]);
            obj.fig.WindowButtonDownFcn = @obj.fig_button_down_cb;
            obj.fig.WindowButtonMotionFcn = @obj.fig_button_motion_cb;
            obj.fig.WindowButtonUpFcn = @obj.fig_button_up_cb;
            %obj.fig.WindowKeyPressFcn = @obj.fig_window_key_press_cb;
            %obj.fig.WindowKeyReleaseFcn = @obj.fig_window_key_release_cb;
            obj.fig.WindowScrollWheelFcn = @obj.fig_window_scroll_wheel_cb;
            
            % Model viewer axes
            obj.axes = axes('Parent',obj.fig,'Units','pixels',...
                'Position',[5 130 540 540],'Units','normalized',...
                'DataAspectRatio',[1 1 1]);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
            hold on;
            
            light_angles = [-45 45; -135 45; 90 -30];
            for i=1:size(light_angles,1)
                obj.lights{i} = lightangle(light_angles(i,1), light_angles(i,2));
            end
            
            obj.model_patch = patch('Vertices',[],'Faces',[],...
                'EdgeColor','none','FaceColor',obj.face_color,...
                'FaceLighting','none','SpecularStrength',0.5,'SpecularExponent',20,...
                'Parent',obj.axes);
            obj.model_vertex_patch = patch('Vertices',[],'Faces',[],...
                'MarkerFaceColor',obj.vertex_color,'MarkerSize',4,'Marker','o',...
                'MarkerEdgeColor','none','Visible','off','Parent',obj.axes);
            obj.grid_patch = patch('Vertices',[],'Faces',[],...
                'EdgeColor',obj.grid_color,'FaceColor','none','Parent',obj.axes);
            obj.bounding_box_patch = patch('Vertices',[],'Faces',[],...
                'EdgeColor',obj.bb_color,'FaceColor','none','LineWidth',1,...
                'Parent',obj.axes);
            obj.boundary_patch = patch('Vertices',[],'Faces',[],...
                'EdgeColor',obj.boundary_color,'FaceColor','none','LineWidth',2,...
                'Visible','off','Parent',obj.axes);
            obj.model_fn_patch = patch('Vertices',[],'Faces',[],...
                'EdgeColor',obj.face_normal_color,'FaceColor','none',...
                'LineWidth',1.5,'Visible','off','Parent',obj.axes);
            obj.model_vn_patch = patch('Vertices',[],'Faces',[],...
                'EdgeColor',obj.vertex_normal_color,'FaceColor','none',...
                'LineWidth',1.5,'Visible','off','Parent',obj.axes);
            obj.model_cline_patch_x = patch('Vertices',[],'Faces',[],...
                'EdgeColor',[0.5 0.2 0.2],'FaceColor','none',...
                'LineWidth',1.5,'Visible','on','Parent',obj.axes);
            obj.model_cline_patch_y = patch('Vertices',[],'Faces',[],...
                'EdgeColor',[0.2 0.2 0.5],'FaceColor','none',...
                'LineWidth',1.5,'Visible','on','Parent',obj.axes);
            
            % Parameter domain axes
            obj.uv_plot = UVPlot(obj.fig, [560 150 500 500]);
            
            % Model Stats
            obj.model_stats_text = uicontrol('Style','text','String','',...
                'Position',[0 782 200 20],'Units','normalized','HorizontalAlignment','left');
            
            % File Panel
            file_panel = uipanel(obj.fig,'Title','File','Units','pixels',...
                'Position',[1100 740 490 56],'Units','normalized');
            uicontrol(file_panel,'Style','pushbutton','String','Load Model...',...
                'Position',[5 20 80 20],'Units','normalized','Callback',@obj.loadModelPressed);
            
            % View Panel
            view_panel = uipanel(obj.fig,'Title','View','Units','pixels',...
                'Position',[1100 640 490 96],'Units','normalized');
            
            % V,E,F
            uicontrol(view_panel,'Style','checkbox','String','Vertices',...
                'Position',[5 64 60 20],'Units','normalized','Value',0,'Callback',@obj.viewOptionChanged);
            uicontrol(view_panel,'Style','checkbox','String','Edges',...
                'Position',[75 64 50 20],'Units','normalized','Value',0,'Callback',@obj.viewOptionChanged);
            uicontrol(view_panel,'Style','checkbox','String','Faces',...
                'Position',[135 64 60 20],'Units','normalized','Value',1,'Callback',@obj.viewOptionChanged);
            obj.faces_visible = 1;
            
            % Grid, BB, Boundary, UVs
            uicontrol(view_panel,'Style','checkbox','String','Grid',...
                'Position',[5 44 50 20],'Units','normalized','Value',1,'Callback',@obj.viewOptionChanged);
            uicontrol(view_panel,'Style','checkbox','String','Bounding Box',...
                'Position',[60 44 90 20],'Units','normalized','Value',1,'Callback',@obj.viewOptionChanged);
            uicontrol(view_panel,'Style','checkbox','String','Boundary',...
                'Position',[5 24 90 20],'Units','normalized','Value',0,'Callback',@obj.viewOptionChanged);
            uicontrol(view_panel, 'Style', 'checkbox','String','UV Lines',...
                'Position',[85 24 90 20],'Units','normalized','Value',1,'Callback',@obj.viewOptionChanged);
            
            % Normals
            uicontrol(view_panel,'Style','checkbox','String','Vertex Normals',...
                'Position',[240 64 100 20],'Units','normalized','Value',0,'Callback',@obj.viewOptionChanged);
            vn_type_group = uibuttongroup(view_panel,'BorderType','none',...
                'Units','pixels','Position',[350 64 200 20],'Units','normalized',...
                'SelectionChanged',@obj.vertexNormalWeightingChanged);
            uicontrol(vn_type_group, 'Style','radiobutton','String','area',...
                'Position',[0 0 50 20],'Units','normalized');
            uicontrol(vn_type_group, 'Style','radiobutton','String','angle',...
                'Position',[50 0 50 20],'units','normalized');
            obj.vertex_normal_weighting = 'area';
            
            uicontrol(view_panel,'Style','checkbox','String','Face Normals',...
                'Position',[240 44 100 20],'Units','normalized','Value',0,'Callback',@obj.viewOptionChanged);
            uicontrol(view_panel,'Style','text','String','Scale:',...
                'Position',[240 20 40 20],'Units','normalized','HorizontalAlignment','left');
            uicontrol(view_panel,'Style','slider','Min',-2,'Max',2,'SliderStep',[0.2/4 0.5/4],...
                'Position',[280 24 150 20],'Units','normalized','Value',0,'Callback',@obj.normalScaleChanged);
            obj.normal_scale = 0.04;
            
            % Align View
            uicontrol(view_panel,'Style','text','String','Align View',...
                'Position',[5 0 60 20],'Units','normalized','HorizontalAlignment','left');
            for i=1:6
                uicontrol(view_panel,'Style','pushbutton','String',obj.align_view_text{i},...
                    'Position',[65+(i-1)*30 2 25 20],'Units','normalized','UserData',i,...
                    'Callback',@obj.alignViewPressed);
            end
            
            % Shading Panel
            shading_panel = uipanel(obj.fig,'Title','Shading','Units','pixels',...
                'Position',[1100 580 490 56],'Units','normalized');
            
            uicontrol(shading_panel,'Style','text','String','Type:',...
                'Position',[5 20 60 20],'Units','normalized','HorizontalAlignment','left');
            obj.shading_type_group = uibuttongroup(shading_panel,'BorderType','none',...
                'Units','pixels','Position',[50 4 380 40],'Units','normalized',...
                'SelectionChanged',@obj.shadingTypeChanged);
            uicontrol(obj.shading_type_group,'Style','radiobutton','String','none',...
                'Position',[0 20 50 20],'Units','normalized');
            uicontrol(obj.shading_type_group,'Style','radiobutton','String','flat',...
                'Position',[50 20 50 20],'Units','normalized');
            uicontrol(obj.shading_type_group,'Style','radiobutton','String','Gouraud',...
                'Position',[90 20 70 20],'Units','normalized');
            uicontrol(obj.shading_type_group,'Style','radiobutton','String','mean curvature', ...
                'Position',[0 0 100 20],'Units','normalized');
            uicontrol(obj.shading_type_group,'Style','radiobutton','String','Gaussian curvature', ...
                'Position',[100 0 120 20],'Units','normalized');
            uicontrol(shading_panel,'Style','text', 'String', 'Quantile Range', ...
                'Position',[285 -1 80 20],'Units','normalized','HorizontalAlignment','left');
            uicontrol(shading_panel,'Style','edit','String','0.05',...
                'Position',[365 2 40 20],'Units','normalized','HorizontalAlignment','left',...
                'Callback',@obj.quantileRangeChanged,'UserData',{1,0.05});
            uicontrol(shading_panel,'Style','text', 'String', '-', ...
                'Position',[405 1 10 20],'Units','normalized');
            uicontrol(shading_panel,'Style','edit','String','0.95',...
                'Position',[415 2 40 20],'Units','normalized','HorizontalAlignment','left',...
                'Callback',@obj.quantileRangeChanged,'UserData',{2,0.95});
            obj.quantile_range = [0.05 0.95];
            
            % Laplacian Panel
            laplacian_panel = uipanel(obj.fig,'Title','Laplacian','Units','pixels',...
                'Position',[1100 520 490 56],'Units','normalized');
            
            % Laplacian Weighting
            uicontrol(laplacian_panel,'Style','text','String','Weighting:',...
                'Position',[5 20 60 20],'Units','normalized','HorizontalAlignment','left');
            laplacian_weighting_group = uibuttongroup(laplacian_panel,'BorderType','none',...
                'Units','pixels','Position',[65 24 250 20],'Units','normalized',...
                'SelectionChanged',@obj.laplacianWeightingChanged);
            uicontrol(laplacian_weighting_group,'Style','radiobutton','String','Uniform',...
                'Position',[0 0 65 20],'Units','normalized');
            uicontrol(laplacian_weighting_group,'Style','radiobutton','String','Cotangent',...
                'Position',[65 0 80 20],'Units','normalized');
            uicontrol(laplacian_weighting_group,'Style','radiobutton','String','Mean Value',...
                'Position',[145 0 80 20],'Units','normalized');
            uicontrol(laplacian_panel, 'Style','checkbox','String','normalized',...
                'Position',[300 24 100 20],'Units','normalized','Callback',@obj.laplacianNormalizedChanged,...
                'Value',1);
            uicontrol(laplacian_panel,'Style','checkbox','String','Recompute every iteration',...
                'Position',[5 4 150 20],'Units','normalized','Value',1,'Callback',@obj.recomputeCheckboxChanged);
            uicontrol(laplacian_panel,'Style','pushbutton','String','Recompute once',...
                'Position',[160 2 100 20],'Units','normalized','Callback',@obj.recomputeOnceCallback);
            
            obj.laplacian_weighting = 'Uniform';
            obj.laplacian_normalized = 1;
            obj.recompute_every_iteration = 1;
            obj.L_uniform_dirty = 1;
            obj.L_Cotangent_dirty = 1;
            obj.L_MeanValue_dirty = 1;
            
            % Smoothing Panel
            smoothing_panel = uipanel(obj.fig,'Title','Smoothing','Units','pixels',...
                'Position',[1100 300 490 216],'Units','normalized');
            
            % Num Iterations
            uicontrol(smoothing_panel, 'Style','text','String','Iterations',...
                'Position',[5 180 60 20],'Units','normalized','HorizontalAlignment','left');
            iterations_ui = NumericSliderEdit(smoothing_panel,[80 183 150 20],1,20,@obj.iterationsChanged);
            iterations_ui.setInteger(1);
            iterations_ui.setValue(1);
            iterations_ui.setSliderStep([1/19 5/19]);
            % Preserve Volume
            obj.preserve_volume_ui = uicontrol(smoothing_panel, 'Style','checkbox','String','Preserve volume',...
                'Position',[240 183 120 20],'Units','normalized','Callback',@obj.preserveVolumeChanged);
            
            % Lambda
            uicontrol(smoothing_panel, 'Style','text','String','Lambda',...
                'Position',[5 155 45 20],'Units','normalized','HorizontalAlignment','left');
            basic_lambda_ui = NumericSliderEdit(smoothing_panel,[80 158 150 20],0,5,@obj.basicLambdaChanged);
            
            % Explicit & Implicit
            uicontrol(smoothing_panel, 'Style','pushbutton','String','Explicit Smoothing',...
                'Position',[5 132 100 20],'Units','normalized','Callback',@obj.basicSmoothingCallback,...
                'UserData','explicit');
            uicontrol(smoothing_panel, 'Style','pushbutton','String','Implicit Smoothing',...
                'Position',[110 132 100 20],'Units','normalized','Callback',@obj.basicSmoothingCallback,...
                'UserData','implicit');
            % Reset
            uicontrol(smoothing_panel, 'Style','pushbutton','String','Reset',...
                'Position',[220 132 60 20],'Units','normalized','Callback',@obj.resetMeshCallback);
            
            obj.basic_lambda = 1;
            obj.preserve_volume = 0;
            obj.num_iterations = 1;
            
            % WL - WP Slider
            obj.wl_ui = uicontrol(smoothing_panel, 'Style','text','String','0.500 = WL',...
                'Position',[5 104 65 20],'Units','normalized','HorizontalAlignment','left');
            uicontrol(smoothing_panel, 'Style','slider','Min',0.001,'Max',0.999,'SliderStep',[0.1 0.2],...
                'Position',[70 106 100 20],'Value',0.5,'Units','normalized','Callback',@obj.wlwpSliderChanged);
            obj.wp_ui = uicontrol(smoothing_panel, 'Style','text','String','WP = 0.500',...
                'Position',[175 104 60 20],'Units','normalized','HorizontalAlignment','left');
            uicontrol(smoothing_panel, 'Style','pushbutton','String','LSQ Smoothing',...
                'Position',[5 80 100 20],'Units','normalized','Callback',@obj.basicSmoothingCallback,...
                'UserData','lsq');
            uicontrol(smoothing_panel, 'Style','pushbutton','String','Triangle Optimization',...
                'Position',[110 80 130 20],'Units','normalized','Callback',@obj.basicSmoothingCallback,...
                'UserData','triangle');
            
            obj.lsq_weighting = 0.5;
            
            % Spectral Smoothing
            uicontrol(smoothing_panel, 'Style','text','String','Eigenfunctions', ...
                'Position',[5 54 80 20],'Units','normalized','HorizontalAlignment','left');
            obj.eigenfunctions_ui = NumericSliderEdit(smoothing_panel,[90 56 150 20],1,20,@obj.eigenfunctionsChanged);
            obj.eigenfunctions_ui.setValue(1);
            obj.eigenfunctions_ui.setSliderStep([1/19 5/19]);
            obj.eigenfunctions_ui.setInteger(1);
            uicontrol(smoothing_panel, 'Style','pushbutton','String','Spectral Smoothing',...
                'Position',[5 30 120 20],'Units','normalized','Callback',@obj.basicSmoothingCallback,...
                'UserData','spectral');
            
            % Parameterization Panel
            param_panel = uipanel(obj.fig,'Title','Parameterization','Units','pixels',...
                'Position',[1100 50 490 246],'Units','normalized');
            
            % UV Line Spacing
            uicontrol(param_panel, 'Style','text','String','UV Spacing:',...
                'Position',[5 204 70 20],'Units','normalized','HorizontalAlignment','left');
            uvSpacingSlider = NumericSliderEdit(param_panel,[75 206 200 20],5,100,@obj.uvLineSpacingChanged);
            uvSpacingSlider.setSliderStep([1 5]/95);
            uvSpacingSlider.setInteger(1);
            uvSpacingSlider.setValue(20);
            
            % Boundary Shape
            uicontrol(param_panel, 'Style','text','String','Boundary:',...
                'Position',[5 180 80 20],'Units','normalized',...
                'HorizontalAlignment','left');
            boundary_shape_group = uibuttongroup(param_panel,'BorderType','none',...
                'Units','pixels','Position',[65 184 200 20],'Units','normalized',...
                'SelectionChanged',@obj.boundaryShapeChanged);
            uicontrol(boundary_shape_group,'Style','radiobutton','String','Circle',...
                'Position',[0 0 60 20],'Units','normalized','UserData','circle');
            uicontrol(boundary_shape_group,'Style','radiobutton','String','Square',...
                'Position',[60 0 80 20],'Units','normalized','UserData','square');
            obj.bdry_shape = 'circle';
            
            % Boundary Spacing
            uicontrol(param_panel, 'Style','text','String','Spacing:',...
                'Position',[220 180 80 20],'Units','normalized',...
                'HorizontalAlignment','left');
            boundary_spacing_group = uibuttongroup(param_panel,'BorderType','none',...
                'Units','pixels','Position',[270 184 200 20],'Units','normalized',...
                'SelectionChanged',@obj.boundarySpacingChanged);
            uicontrol(boundary_spacing_group,'Style','radiobutton','String','Uniform',...
                'Position',[0 0 60 20],'Units','normalized','UserData',false);
            uicontrol(boundary_spacing_group,'Style','radiobutton','String','Adaptive',...
                'Position',[65 0 80 20],'Units','normalized','UserData',true);
            
            obj.adaptive_spacing = false;
            
            uicontrol(param_panel, 'Style','pushbutton','String','Linear Barycentric Mapping',...
                'Position',[5 160 160 20],'Units','normalized','Callback',@obj.computeMappingCallback,...
                'UserData','barycentric');
            
            % Fixed vertex list
            str = cell(0,1);
            uicontrol(param_panel, 'Style','text','String','Fixed points:',...
                'Position',[5 132 100 20],'Units','normalized',...
                'HorizontalAlignment','left');
            obj.fixed_list = uicontrol(param_panel, 'Style','listbox','String',str,...
                'Position',[5 55 150 80],'Units','normalized','Callback',@obj.fixedListSelectionChanged);
            obj.add_fixed_button = uicontrol(param_panel, 'Style','pushbutton','String','Add',...
                'Position',[160 115 60 20],'Units','normalized','Callback',@obj.addFixedVertex);
            uicontrol(param_panel, 'Style','pushbutton','String','Remove',...
                'Position',[160 95 60 20],'Units','normalized','Callback',@obj.removeFixedVertex);
            
            uicontrol(param_panel, 'Style','pushbutton','String','Least-Squares Conformal Mapping',...
                'Position',[5 25 200 20],'Units','normalized','Callback',@obj.computeMappingCallback,...
                'UserData','lscm');
            uicontrol(param_panel, 'Style','pushbutton','String','Discrete Conformal Mapping',...
                'Position',[5 5 200 20],'Units','normalized','Callback',@obj.computeMappingCallback,...
                'UserData','dcp');
            
            movegui(obj.fig,'center');
            obj.fig.Visible = 'on';
            
            basic_lambda_ui.setValue(1);
        end
        
        function fixedListSelectionChanged(obj, source, data)
            v = source.Value;
            if v>=1 && v <= length(obj.fixed_vi)
                i = obj.fixed_vi(v);
                obj.uv_plot.select(i);
            end
        end
        
        function resetFixedVertices(obj)
            obj.fixed_vi = zeros(0,1);
            obj.fixed_vpos = zeros(0,2);
            obj.updateFixedList();
        end
        
        function addFixedVertex(obj, source, data)
            if ~isempty(obj.uv_plot.selection)
                obj.uv_plot.requestPoint(@obj.requestPointCallback);
                obj.add_fixed_button.Enable = 'off';
            end
        end
        
        function requestPointCallback(obj, p)
            vi = obj.uv_plot.selection;
            res = find(obj.fixed_vi == vi);
            if isempty(res)
                if isempty(obj.fixed_vi)
                    obj.fixed_vi = vi;
                    obj.fixed_vpos = p;
                else
                    obj.fixed_vi = [obj.fixed_vi; vi];
                    obj.fixed_vpos = [obj.fixed_vpos; p];
                end
            else
                obj.fixed_vpos(res,:) = p;
            end
            obj.updateFixedList();
            obj.add_fixed_button.Enable = 'on';
        end
        
        function updateFixedList(obj)
            n = length(obj.fixed_vi);
            str = cell(n,1);
            for i=1:n
                str{i} = sprintf('v %i: (%.3f, %.3f)', obj.fixed_vi(i), ...
                    obj.fixed_vpos(i,1), obj.fixed_vpos(i,2));
            end
            obj.fixed_list.String = str;
            obj.uv_plot.setFixedPoints(obj.fixed_vpos);
            obj.fixed_list.Value = 1;
        end
        
        function removeFixedVertex(obj, source, data)
            if isempty(source.Value)
                return;
            end
            if source.Value >= 1 && source.Value <= length(source.String)
                v = obj.fixed_list.Value;
                obj.fixed_vpos = [obj.fixed_vpos(1:(v-1),:); obj.fixed_vpos((v+1):end,:)];
                obj.fixed_vi = [obj.fixed_vi(1:(v-1)); obj.fixed_vi((v+1):end)];
                obj.updateFixedList();
            end
        end
        
        function uvLineSpacingChanged(obj, val)
            obj.uv_plot.setClineSpacing(1/val);
            obj.updateClinePlot();
        end
        
        function computeMappingCallback(obj, source, data)
            switch source.UserData
                case 'barycentric'
                    uv = MeshParameterization.barycentricParameterization(obj.mesh, ...
                        obj.getLaplacian(), obj.bdry_shape, obj.adaptive_spacing);
                case 'lscm'
                    if length(obj.fixed_vi) < 2
                        error('Fix at least 2 points!');
                    end
                    uv = MeshParameterization.lscmParameterization(obj.mesh, obj.fixed_vi, obj.fixed_vpos);
                case 'dcp'
                    if length(obj.fixed_vi) < 2
                        error('Fix at least 2 points!');
                    end
                    uv = MeshParameterization.dcpParameterization(obj.mesh, obj.fixed_vi, obj.fixed_vpos);
                otherwise
                    error('Unknown parameterization technique!');
            end
            obj.mesh.getAllVertices().setTrait('uv', uv);
            obj.uv_plot.updateUVCoordinates();
            
            obj.updateClinePlot();
        end
        
        function updateClinePlot(obj)
            [v1,e1,v2,e2] = MeshHelper.computeCoordinateLines(obj.mesh, obj.uv_plot.cline_spacing);
            obj.model_cline_patch_x.Vertices = v1;
            obj.model_cline_patch_x.Faces = e1;
            obj.model_cline_patch_y.Vertices = v2;
            obj.model_cline_patch_y.Faces = e2;
        end
        
        function resetClinePlot(obj)
            obj.model_cline_patch_x.Vertices = [];
            obj.model_cline_patch_x.Faces = [];
            obj.model_cline_patch_y.Vertices = [];
            obj.model_cline_patch_y.Faces = [];
        end
        
        function boundaryShapeChanged(obj, source, data)
            obj.bdry_shape = source.SelectedObject.UserData;
        end
        
        function boundarySpacingChanged(obj, source, data)
            obj.adaptive_spacing = source.SelectedObject.UserData;
        end
        
        function quantileRangeChanged(obj, source, data)
            val = str2double(source.String);
            if isnan(val)
                source.String = source.UserData{2};
            else
                if val < 0 || val > 1
                    val = min(max(val,0),1);
                    source.String = val;
                end
                source.UserData{2} = val;
                obj.quantile_range(source.UserData{1}) = val;
                obj.updateShading(obj.shading_type_group.SelectedObject.String);
            end
        end
        
        function laplacianNormalizedChanged(obj, source, data)
            obj.laplacian_normalized = source.Value;
            obj.L_uniform_dirty = 1;
            obj.L_Cotangent_dirty = 1;
            obj.L_MeanValue_dirty = 1;
        end
        
        function eigenfunctionsChanged(obj, val)
            obj.num_eigenfunctions = val;
        end
        
        function wlwpSliderChanged(obj,source,data)
            obj.lsq_weighting = source.Value;
            obj.wl_ui.String = sprintf('%0.3f = WL', 1-obj.lsq_weighting);
            obj.wp_ui.String = sprintf('WP = %0.3f', obj.lsq_weighting);
        end
        
        function resetMeshCallback(obj,source,data)
            obj.mesh.getAllVertices().setTrait('position',obj.initial_vertex_positions);
            obj.updateModel();
            obj.L_uniform_dirty = 1;
            obj.L_Cotangent_dirty = 1;
            obj.L_MeanValue_dirty = 1;
        end
        
        function iterationsChanged(obj,val)
            obj.num_iterations = val;
        end
        
        function basicSmoothingCallback(obj,source,data)
            for i=1:obj.num_iterations
                switch source.UserData
                    case 'implicit'
                        V_smooth = MeshSmoothing.implicitSmoothing(obj.mesh,obj.getLaplacian(),obj.basic_lambda);
                    case 'explicit'
                        V_smooth = MeshSmoothing.explicitSmoothing(obj.mesh,obj.getLaplacian(),obj.basic_lambda);
                    case 'lsq'
                        V_smooth = MeshSmoothing.lsqSmoothing(obj.mesh,obj.getLaplacian(),1-obj.lsq_weighting, obj.lsq_weighting);
                    case 'triangle'
                        V_smooth = MeshSmoothing.triangleSmoothing(obj.mesh, obj.getLaplacian('Uniform'), obj.getLaplacian('Cotangent'), ...
                            1-obj.lsq_weighting, obj.lsq_weighting);
                    case 'spectral'
                        V_smooth = MeshSmoothing.spectralSmoothing(obj.mesh, obj.getLaplacian(), obj.num_eigenfunctions);
                    otherwise
                        V_smooth = obj.mesh.toFaceVertexMesh();
                end
                obj.mesh.getAllVertices().setTrait('position',V_smooth);
                if obj.recompute_every_iteration
                    obj.L_uniform_dirty = 1;
                    obj.L_Cotangent_dirty = 1;
                    obj.L_MeanValue_dirty = 1;
                    if i~=obj.num_iterations
                        MeshHelper.calculateHalfedgeTraits(obj.mesh);
                    end
                end
            end
            
            if obj.preserve_volume
                new_vol = MeshHelper.computeVolume(obj.mesh);
                if abs(new_vol) > 10e-10
                    MeshHelper.scaleMesh(obj.mesh, (obj.initial_mesh_volume / new_vol)^(1/3));
                end
            end
            obj.updateModel(0);
            
        end
        
        function L = getLaplacian(obj,type)
            if nargin < 2
                type = obj.laplacian_weighting;
            end
            switch type
                case 'Uniform'
                    if obj.L_uniform_dirty
                        obj.L_uniform = MeshLaplacian.computeUniformLaplacian(obj.mesh,obj.laplacian_normalized);
                        obj.L_uniform_dirty = 0;
                    end
                    L = obj.L_uniform;
                    return;
                case 'Cotangent'
                    if obj.L_Cotangent_dirty
                        obj.L_Cotangent = MeshLaplacian.computeCotangentLaplacian(obj.mesh,obj.laplacian_normalized);
                        obj.L_Cotangent_dirty = 0;
                    end
                    L = obj.L_Cotangent;
                    return;
                case 'Mean Value'
                    if obj.L_MeanValue_dirty
                        obj.L_MeanValue = MeshLaplacian.computeMeanValueLaplacian(obj.mesh,obj.laplacian_normalized);
                        obj.L_MeanValue_dirty = 0;
                    end
                    L = obj.L_MeanValue;
                    return;
            end
        end
        
        function preserveVolumeChanged(obj,source,data)
            obj.preserve_volume = source.Value;
        end
        
        function laplacianWeightingChanged(obj,source,data)
            obj.laplacian_weighting = data.NewValue.String;
        end
        
        function recomputeOnceCallback(obj, source, data)
            obj.L_uniform_dirty = 1;
            obj.L_Cotangent_dirty = 1;
            obj.L_MeanValue_dirty = 1;
        end
        
        function recomputeCheckboxChanged(obj,source,data)
            obj.recompute_every_iteration = source.Value;
            if obj.recompute_every_iteration
                obj.L_uniform_dirty = 1;
                obj.L_Cotangent_dirty = 1;
                obj.L_MeanValue_dirty = 1;
            end
        end
        
        function basicLambdaChanged(obj,val)
            obj.basic_lambda = val;
        end
        
        function alignViewPressed(obj,source,data)
            i = source.UserData;
            view(obj.align_view_angles(i,:));
        end
        
        function vertexNormalWeightingChanged(obj,source,data)
            obj.vertex_normal_weighting = data.NewValue.String;
            MeshHelper.calculateVertexNormals(obj.mesh, obj.vertex_normal_weighting);
            obj.updateVertexNormalModel();
        end
        
        function loadModelPressed(obj,source,data)
            filters = {'*.obj', 'OBJ Files (*.mat)';...
                '*.*', 'All Files (*.*)'};
            if isempty(obj.last_file)
                [filename, pathname] = uigetfile(filters, 'Load Model...');
            else
                [filename, pathname] = uigetfile(filters, 'Load Model...', obj.last_file);
            end
            if ~(isnumeric(filename) && isnumeric(pathname))
                obj.last_file = [pathname filename];
                obj.loadModel(obj.last_file);
            end
        end
        
        function shadingTypeChanged(obj,source,data)
            obj.updateShading(data.NewValue.String);
        end
        
        function updateShading(obj, str)
            axes(obj.axes);
            
            obj.model_patch.FaceVertexCData = [];
            
            if ~obj.faces_visible
                obj.model_patch.FaceColor = 'none';
                obj.model_patch.FaceVertexCData = [];
            else
                obj.model_patch.FaceColor = obj.face_color;
                v_colors = [];
                switch str
                    case 'none'
                        obj.model_patch.FaceLighting = str;
                    case 'flat'
                        obj.model_patch.FaceLighting = str;
                    case 'Gouraud'
                        obj.model_patch.FaceLighting = str;
                    case 'mean curvature'
                        v_colors = obj.mesh.getAllVertices().getTrait('mean_curv');
                    case 'Gaussian curvature'
                        v_colors = obj.mesh.getAllVertices().getTrait('gauss_curv');
                end
                if ~isempty(v_colors)
                    obj.model_patch.FaceLighting = 'none';
                    obj.model_patch.FaceColor = 'interp';
                    c_min = quantile(v_colors, obj.quantile_range(1));
                    c_max = quantile(v_colors, obj.quantile_range(2));
                    v_colors(v_colors < c_min) = c_min;
                    v_colors(v_colors > c_max) = c_max;
                    obj.model_patch.FaceVertexCData = v_colors;
                    colorbar
                else
                    colorbar('off');
                end
            end
        end
        
        function normalScaleChanged(obj,source,data)
            obj.normal_scale = 0.04*(3^source.Value);
            obj.updateFaceNormalModel();
            obj.updateVertexNormalModel();
        end
        
        function viewOptionChanged(obj, source, ~)
            switch(source.String)
                case 'Vertices'
                    obj.model_vertex_patch.Visible = val2vis(source.Value);
                case 'Edges'
                    if source.Value==1
                        obj.model_patch.EdgeColor = obj.edge_color;
                    else
                        obj.model_patch.EdgeColor = 'none';
                    end
                case 'Faces'
                    obj.faces_visible = source.Value;
                    obj.updateShading(obj.shading_type_group.SelectedObject.String);
                case 'Vertex Normals'
                    obj.model_vn_patch.Visible = val2vis(source.Value);
                case 'Face Normals'
                    obj.model_fn_patch.Visible = val2vis(source.Value);
                case 'Grid'
                    obj.grid_patch.Visible = val2vis(source.Value);
                case 'Bounding Box'
                    obj.bounding_box_patch.Visible = val2vis(source.Value);
                case 'Boundary'
                    obj.boundary_patch.Visible = val2vis(source.Value);
                case 'UV Lines'
                    obj.model_cline_patch_x.Visible = val2vis(source.Value);
                    obj.model_cline_patch_y.Visible = val2vis(source.Value);
                otherwise
            end
        end
        
        function loadModel(obj, filename)
            obj.mesh = ModelLoader.loadOBJ(filename);
            obj.uv_plot.setMesh(obj.mesh);
            
            obj.updateModelStats();
            
            obj.updateModel();
            obj.resetFixedVertices();
            
            obj.initial_mesh_volume = MeshHelper.computeVolume(obj.mesh);
            obj.initial_vertex_positions = obj.mesh.toFaceVertexMesh();
            obj.L_uniform_dirty = 1;
            obj.L_Cotangent_dirty = 1;
            obj.L_MeanValue_dirty = 1;
            hasBoundary = MeshHelper.hasBoundary(obj.mesh);
            obj.preserve_volume_ui.Enable = val2vis(~hasBoundary);
            if hasBoundary
                obj.preserve_volume_ui.Value = 0;
            end
            obj.eigenfunctions_ui.setMin(1);
            obj.eigenfunctions_ui.setMax(obj.mesh.num_vertices-2);
            obj.eigenfunctions_ui.setSliderStep([1 5]/(obj.mesh.num_vertices-3));
        end
        
        function updateModel(obj,reset_camera)
            axes(obj.axes);
            
            if nargin < 2
                reset_camera = 1;
            end
            [V,F] = obj.mesh.toFaceVertexMesh();
            
            obj.model_patch.Vertices = V;
            obj.model_patch.Faces = F;
            
            obj.model_vertex_patch.Vertices = V;
            obj.model_vertex_patch.Faces = (1:obj.mesh.num_vertices)';
            
            obj.resetClinePlot();
            
            obj.updateAxes(reset_camera);
            
            MeshHelper.calculateFaceTraits(obj.mesh);
            MeshHelper.calculateVertexTraits(obj.mesh);
            MeshHelper.calculateHalfedgeTraits(obj.mesh);
            MeshHelper.calculateVertexNormals(obj.mesh, obj.vertex_normal_weighting);
            MeshHelper.calculateDiscreteCurvatures(obj.mesh);
            obj.updateShading(obj.shading_type_group.SelectedObject.String);
            
            obj.updateGridModel();
            obj.updateBoundingBoxModel();
            obj.updateBoundaryModel();
            obj.updateFaceNormalModel();
            obj.updateVertexNormalModel();
        end
        
        function updateModelStats(obj)
            nv = obj.mesh.num_vertices;
            ne = obj.mesh.num_edges;
            nf = obj.mesh.num_faces;
            
            obj.model_stats_text.String = sprintf(...
                ' v %i e %i f %i',nv,ne,nf);
        end
        
        function updateAxes(obj,reset_camera)
            axes(obj.axes);
            
            if nargin < 2
                reset_camera = 1;
            end
            V = obj.mesh.toFaceVertexMesh();
            
            pmin = min(V,[],1);
            pmax = max(V,[],1);
            
            view_offset = max(max(0.1*(pmax-pmin)), 1e-3);
            view_min = pmin - view_offset;
            view_max = pmax + view_offset;
            
            axis vis3d;
            if reset_camera
                view(-60,30);
                zoom out;
                zoom(0.7);
            end
            
            xlim([view_min(1) view_max(1)]);
            ylim([view_min(2) view_max(2)]);
            zlim([view_min(3) view_max(3)]);
        end
        
        function updateGridModel(obj)
            axes(obj.axes);
            
            r = max([obj.axes.XLim(2)-obj.axes.XLim(1)...
                obj.axes.YLim(2)-obj.axes.YLim(1)]);
            step = 10^round(log(r/10)/log(10));
            x_ticks = ((ceil(obj.axes.XLim(1)/step)*step):step:(floor(obj.axes.XLim(2)/step)*step))';
            y_ticks = ((ceil(obj.axes.YLim(1)/step)*step):step:(floor(obj.axes.YLim(2)/step)*step))';
            v1 = [kron(x_ticks,[1;1]) repmat(obj.axes.YLim(:),length(x_ticks),1)];
            v2 = [repmat(obj.axes.XLim(:),length(y_ticks),1) kron(y_ticks,[1;1])];
            
            height = 0;
            zl = zlim;
            if zl(1)>0
                height = zl(1)*0.99+zl(2)*0.01;
            end
            obj.grid_patch.Vertices = [[v1;v2] height*ones(2*(length(x_ticks)+length(y_ticks)),1)];
            obj.grid_patch.Faces = reshape(1:size(obj.grid_patch.Vertices,1),2,[])';
        end
        
        function updateBoundingBoxModel(obj)
            if isempty(obj.mesh)
                return;
            end
            [p_min, p_max] = MeshHelper.getBoundingBox(obj.mesh);
            [V,E] = GeometryHelper.buildBoxEdges(p_min,p_max);
            
            obj.bounding_box_patch.Vertices = V;
            obj.bounding_box_patch.Faces = E;
        end
        
        function updateBoundaryModel(obj)
            [V_start, V_end] = MeshHelper.getBoundaryEdges(obj.mesh);
            
            obj.boundary_patch.Vertices = reshape([V_start'; V_end'],3,[])';
            obj.boundary_patch.Faces = reshape(1:size(obj.boundary_patch.Vertices,1),2,[])';
        end
        
        function updateFaceNormalModel(obj)
            if ~isfield(obj.mesh.F_traits, 'normal') || ...
                    ~isfield(obj.mesh.F_traits, 'centroid')
                obj.model_fn_patch.Vertices = [];
                obj.model_fn_patch.Faces = [];
                return;
            end
            
            f = obj.mesh.getAllFaces();
            v1 = f.getTrait('centroid');
            v2 = v1 + f.getTrait('normal')*obj.normal_scale;
            
            obj.model_fn_patch.Vertices = reshape([v1'; v2'],3,[])';
            obj.model_fn_patch.Faces = reshape(1:size(obj.model_fn_patch.Vertices,1),2,[])';
        end
        
        function updateVertexNormalModel(obj)
            if ~isfield(obj.mesh.V_traits, 'normal')
                obj.model_vn_patch.Vertices = [];
                obj.model_vn_patch.Faces = [];
                return;
            end
            
            v = obj.mesh.getAllVertices();
            v1 = v.getTrait('position');
            v2 = v1 + v.getTrait('normal')*obj.normal_scale;
            
            obj.model_vn_patch.Vertices = reshape([v1'; v2'],3,[])';
            obj.model_vn_patch.Faces = reshape(1:size(obj.model_vn_patch.Vertices,1),2,[])';
        end
        
        function fig_button_down_cb(obj, source, data)
            if obj.currentPointInUVPlot()
                obj.uv_plot.buttonDown(source, data);
            end
        end
        
        function fig_button_motion_cb(obj, source, data)
            obj.uv_plot.buttonMotion(source, data);
        end
        
        function fig_button_up_cb(obj, source, data)
            obj.uv_plot.buttonUp(source, data);
        end
        
        function fig_window_scroll_wheel_cb(obj, source, data)
            if obj.currentPointInUVPlot()
                obj.uv_plot.scrollWheel(source, data);
            end
        end
        
        function ret = currentPointInUVPlot(obj)
            pos = obj.uv_plot.getPosition();
            r1 = pos([1 2]);
            r2 = r1 + pos([3 4]);
            cp = obj.fig.CurrentPoint;
            ret = all(r1 < cp) && all(cp < r2);
        end
    end
end