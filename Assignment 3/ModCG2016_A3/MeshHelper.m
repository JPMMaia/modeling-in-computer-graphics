classdef MeshHelper < handle
    methods(Static)
        function l = getBoundaryLoop(mesh)
            % Returns a p-by-1 array of all halfedge indices that belong to
            % the boundary loop of the mesh. The first index is
            % arbitrary, but from there on they are listed in the CW order
            % given by the next()-pointer of the boundary halfedges.
            % If the mesh is a closed manifold, returns [].
            
            % TODO_A3 Task 1a
            % Extract the boundary loop of the loaded mesh in the order
            % given by the halfedge next()-method.
            
            % Challenge: Can you find a parallel solution to this task,
            % without iterating over the elements of the boundary loop?
            % I could not :(

            % Get all halfedges:
            halfedges = mesh.getAllHalfedges();
            
            % Each halfedges references a face:
            faces = halfedges.face();
            facesIndices = faces.index;
            
            % Find all halfedge indices which belong to the boundary:
            boundaryHalfedgesIndices = find(facesIndices == 0);
            
            % If the mesh has no boundary, then return empty list:
            if(size(boundaryHalfedgesIndices, 1) == 0)
               l = [];
               return;
            end
            
            % Allocate memory for the boundary loop array:
            l = zeros(size(boundaryHalfedgesIndices, 1), 1);
            
            % Get the first boundary halfedge of the list and add it to the
            % boundary loop:
            firstHalfedge = mesh.getHalfedge(boundaryHalfedgesIndices(1));
            l(1, 1) = firstHalfedge.index;
            
            % Iterate through the boundary loop and add each one to the
            % boundary loop list:
            index = 2;
            currentHalfedge = firstHalfedge.next();
            while(currentHalfedge.index ~= firstHalfedge.index)
               
                % Add halfedge index to the boundary loop:
                l(index, 1) = currentHalfedge.index;
                index = index + 1;
                
                % Go for the next halfedge of the bounrary:
                currentHalfedge = currentHalfedge.next();
                
            end
            
        end
        
        function b = hasBoundary(mesh)
            % Returns 1 if the mesh has a boundary, and 0 otherwise.
            b = any(mesh.getAllHalfedges().face().index==0);
        end
        
        function vol = computeVolume(mesh)
            % Returns the signed volume of the mesh.

            vol = 1;
        end
        
        function scaleMesh(mesh, factor)
            % Scales the mesh by a gives factor relative to the midpoint of
            % the bounding box.
            V = mesh.toFaceVertexMesh();
            p_min = min(V,[],1);
            p_max = max(V,[],1);
            p_center = 0.5*(p_max+p_min);
            V_new = bsxfun(@plus, bsxfun(@minus, V, p_center)*factor, p_center);
            mesh.getAllVertices().setTrait('position', V_new);
        end
        
        function [p_min, p_max] = getBoundingBox(mesh)
            % Returns the points with minimal and maximal coordinates of the
            % smallest axis-aligned bounding box of the mesh.

            V = mesh.getAllVertices().getTrait('position');
            p_min = min(V,[],1);
            p_max = max(V,[],1);
        end
        
        function [V_start, V_end] = getBoundaryEdges(mesh)
            % Returns a list of line segments describing the boundary of
            % the mesh. Returns two nbe-by-3 arrays (nbe=number of
            % boundary edges), such that the i-th row of V_start and the
            % i-th row of V_end contain the coordinates of the vertices of
            % a boundary edge. The order of the boundary edges is
            % arbitrary.
            
            he = mesh.getAllHalfedges();
            he_bdry = mesh.getHalfedge(he.face().index == 0);

            V_start = he_bdry.from().getTrait('position');
            V_end = he_bdry.to().getTrait('position');
        end
        
        
        function calculateFaceTraits(mesh)
            % Fills in a number of face traits in the TriangleMesh mesh.
            % Each face stores its surface area (trait 'area'), its
            % centroid, which is the arithmetic mean of its corner vertices
            % (trait 'centroid), and its normal (trait 'normal').

            f = mesh.getAllFaces();
            he1 = f.halfedge();
            v1 = he1.from().getTrait('position');
            v2 = he1.to().getTrait('position');
            v3 = he1.next().to.getTrait('position');
            fn_weighted = cross(v2-v1,v3-v1);
            areas = 0.5*sqrt(sum(fn_weighted .* fn_weighted, 2));
            f.setTrait('area',areas);
            f.setTrait('centroid',(v1 + v2 + v3) / 3);
            f.setTrait('normal',normr(fn_weighted));
        end
        
        function calculateVertexTraits(mesh)
            % Computes the degree of each vertex and stores it in the
            % vertex trait 'degree'.
            v = mesh.getAllVertices();
            he1 = v.halfedge();
            he_current = he1.twin().next();
            degs = zeros(mesh.num_vertices,1);
            i=1;
            while any(degs==0)
                degs(degs==0 & he_current.index == he1.index) = i;
                he_current = he_current.twin().next();
                i = i+1;
            end
            v.setTrait('degree',degs);
        end
        
        
        function calculateHalfedgeTraits(mesh)
            % Computes the 'angle' halfedge trait, which gives the angle
            % between the halfedge and its previous halfedge in radians,
            % the 'cot_angle' trait which gives its Cotangent, and 'tan_mv'
            % which gives the Tangent of half the angle.

            he = mesh.getAllHalfedges();
            v = he.from().getTrait('position');
            v_next = he.to().getTrait('position');
            v_prev = he.prev().from().getTrait('position');
            d = v_next-v;
            d_len = sqrt(sum(d.*d,2));
            d1 = bsxfun(@times, d, 1./d_len);
            d2 = normr(v_prev-v);
            angles = acos(sum(d1.*d2,2));
            he.setTrait('angle',angles);

            temp = cot(angles);
            temp(he.face().index==0) = 0;
            he.setTrait('cot_angle',temp);
        end
        
        function calculateDiscreteCurvatures(mesh)
            % Computes vertex traits 'mean_curv' for the discrete mean
            % curvature and 'gauss_curv' for the discrete Gaussian
            % curvature.

            mesh.getAllVertices().setTrait('gauss_curv', 1);
            mesh.getAllVertices().setTrait('mean_curv', 1);
        end
        
        
        function calculateVertexNormals(mesh, weighting)
            % Computes vertex normals as a weighted mean of face normals.
            % The parameter 'weighting' can be one of the following:
            % 'area': The face normal weights equal the face surface areas.
            % 'angle': The face normal weights equal the opening angle
            %    of the face at the vertex.
            if nargin<2
                weighting='area';
            end
            switch weighting
                case 'area'
                    f = mesh.getAllFaces();
                    fn_weighted = bsxfun(@times, f.getTrait('normal'), f.getTrait('area'));

                    vi1 = f.halfedge().from().index;
                    vi2 = f.halfedge().to().index;
                    vi3 = f.halfedge().next().to().index;

                    vi = repmat(reshape([vi1 vi2 vi3]',[],1),3,1);
                    vj = kron([1;2;3],ones(mesh.num_faces*3,1));
                    vals = kron(fn_weighted(:),[1;1;1]);

                    vn = normr(full(sparse(vi,vj,vals)));
                    mesh.getAllVertices().setTrait('normal',vn);

                case 'angle'
                    he = mesh.getAllHalfedges();
                    he_inner = mesh.getHalfedge(he.index(he.face().index ~= 0));
                    vi = he_inner.from().index;
                    fn_weighted = bsxfun(@times, he_inner.face().getTrait('normal'), he_inner.getTrait('angle'));
                    vn = zeros(mesh.num_vertices,3);
                    for i=1:3
                        vn(:,i) = full(sparse(vi, ones(size(vi)), fn_weighted(:,i)));
                    end
                    mesh.getAllVertices().setTrait('normal',normr(vn));
            end
        end
        
        function [v1,e1,v2,e2] = computeCoordinateLines(mesh, cline_spacing)
            % Computes the uv coordinate lines on the mesh. [v1,e1] are
            % u-lines, and [v2,e2] are v-lines. cline_spacing is the
            % distance between coordinate lines in uv-space.
            
            if isempty(mesh)
                v1 = [];
                e1 = [];
                v2 = [];
                e2 = [];
                return;
            end
            if ~isfield(mesh.V_traits, 'uv')
                v1 = [];
                e1 = [];
                v2 = [];
                e2 = [];
                return;
            end
            
            % get uv coords of all face corner vertices
            uv = mesh.getAllVertices().getTrait('uv') * (1/cline_spacing);
            [~,F] = mesh.toFaceVertexMesh();
            uv_f = cell(1,3);
            for i=1:3
                uv_f{i} = uv(F(:,i),:);
            end
            temp_x = [uv_f{1}(:,1) uv_f{2}(:,1) uv_f{3}(:,1)];
            temp_y = [uv_f{1}(:,2) uv_f{2}(:,2) uv_f{3}(:,2)];
            f_bbox_min = [min(temp_x,[],2) min(temp_y,[],2)];
            f_bbox_max = [max(temp_x,[],2) max(temp_y,[],2)];
            
            % get line segment numbers for all faces
            num_x_is = floor(f_bbox_max(:,1)) - ceil(f_bbox_min(:,1)) + 1;
            num_y_is = floor(f_bbox_max(:,2)) - ceil(f_bbox_min(:,2)) + 1;
            
            % coordinates of all line segments in uv coords
            all_x_line_coords = cell2mat(arrayfun(@(a,b)(a:b)', ceil(f_bbox_min(:,1)), floor(f_bbox_max(:,1)), 'UniformOutput', 0));
            all_y_line_coords = cell2mat(arrayfun(@(a,b)(a:b)', ceil(f_bbox_min(:,2)), floor(f_bbox_max(:,2)), 'UniformOutput', 0));
            
            % Generate uv coords of face vertices for each line segment
            uv_mult_x = cell(1,3);
            uv_mult_y = cell(1,3);
            for i=1:3
                uv_mult_x{i} = cell2mat(arrayfun(@(a,b,c) repmat([a b],c,1), uv_f{i}(:,1), uv_f{i}(:,2), num_x_is, 'UniformOutput', 0));
                uv_mult_y{i} = cell2mat(arrayfun(@(a,b,c) repmat([a b],c,1), uv_f{i}(:,1), uv_f{i}(:,2), num_y_is, 'UniformOutput', 0));
            end
            
            % Clip lines with triangles to yield final line segment
            % coordinates in uv space
            x_lines = repmat([-Inf Inf], size(all_x_line_coords,1),1);
            y_lines = repmat([-Inf Inf], size(all_y_line_coords,1),1);
            for i=1:3
                i1 = i;
                i2 = mod(i,3)+1;
                
                from_below = uv_mult_x{i2}(:,1) >= uv_mult_x{i1}(:,1);
                inters = uv_mult_x{i1}(:,2) + ((all_x_line_coords - uv_mult_x{i1}(:,1)) ./ (uv_mult_x{i2}(:,1) - uv_mult_x{i1}(:,1))) .* ...
                    (uv_mult_x{i2}(:,2) - uv_mult_x{i1}(:,2));
                inters(isinf(inters)) = nan;
                x_lines(from_below,1) = max(x_lines(from_below,1), inters(from_below));
                x_lines(~from_below,2) = min(x_lines(~from_below,2), inters(~from_below));
                
                from_left = uv_mult_y{i2}(:,2) <= uv_mult_y{i1}(:,2);
                inters = uv_mult_y{i1}(:,1) + ((all_y_line_coords - uv_mult_y{i1}(:,2)) ./ (uv_mult_y{i2}(:,2) - uv_mult_y{i1}(:,2))) .* ...
                    (uv_mult_y{i2}(:,1) - uv_mult_y{i1}(:,1));
                inters(isinf(inters)) = nan;
                y_lines(from_left,1) = max(y_lines(from_left,1), inters(from_left));
                y_lines(~from_left,2) = min(y_lines(~from_left,2), inters(~from_left));
            end
            
            % Compute barycentric coordinates of line segments wrt faces
            % vertical lines (x-lines)
            y2_min_y3 = uv_mult_x{2}(:,2) - uv_mult_x{3}(:,2);
            x3_min_x2 = uv_mult_x{3}(:,1) - uv_mult_x{2}(:,1);
            x1_min_x3 = uv_mult_x{1}(:,1) - uv_mult_x{3}(:,1);
            y1_min_y3 = uv_mult_x{1}(:,2) - uv_mult_x{3}(:,2);
            denom = y2_min_y3 .* x1_min_x3 + x3_min_x2 .* y1_min_y3;
            
            x_temp = (all_x_line_coords - uv_mult_x{3}(:,1)) ./ denom;
            y1_temp = (x_lines(:,1) - uv_mult_x{3}(:,2)) ./ denom;
            y2_temp = (x_lines(:,2) - uv_mult_x{3}(:,2)) ./ denom;
            
            lambda_1_x = y2_min_y3 .* x_temp;
            lambda_1_y1 = x3_min_x2 .* y1_temp;
            lambda_1_y2 = x3_min_x2 .* y2_temp;
            l1_x = bsxfun(@plus, [lambda_1_y1 lambda_1_y2], lambda_1_x);
            
            lambda_2_x = -y1_min_y3 .* x_temp;
            lambda_2_y1 = x1_min_x3 .* y1_temp;
            lambda_2_y2 = x1_min_x3 .* y2_temp;
            l2_x = bsxfun(@plus, [lambda_2_y1 lambda_2_y2], lambda_2_x);
            
            l3_x = 1-l1_x-l2_x;
            
            % horizonal lines (y-lines)
            y2_min_y3 = uv_mult_y{2}(:,2) - uv_mult_y{3}(:,2);
            x3_min_x2 = uv_mult_y{3}(:,1) - uv_mult_y{2}(:,1);
            x1_min_x3 = uv_mult_y{1}(:,1) - uv_mult_y{3}(:,1);
            y1_min_y3 = uv_mult_y{1}(:,2) - uv_mult_y{3}(:,2);
            denom = y2_min_y3 .* x1_min_x3 + x3_min_x2 .* y1_min_y3;
            
            x1_temp = (y_lines(:,1) - uv_mult_y{3}(:,1)) ./ denom;
            x2_temp = (y_lines(:,2) - uv_mult_y{3}(:,1)) ./ denom;
            y_temp = (all_y_line_coords - uv_mult_y{3}(:,2)) ./ denom;
            
            lambda_1_x1 = y2_min_y3 .* x1_temp;
            lambda_1_x2 = y2_min_y3 .* x2_temp;
            lambda_1_y = x3_min_x2 .* y_temp;
            l1_y = bsxfun(@plus, [lambda_1_x1 lambda_1_x2], lambda_1_y);
            
            lambda_2_x1 = -y1_min_y3 .* x1_temp;
            lambda_2_x2 = -y1_min_y3 .* x2_temp;
            lambda_2_y = x1_min_x3 .* y_temp;
            l2_y = bsxfun(@plus, [lambda_2_x1 lambda_2_x2], lambda_2_y);
            
            l3_y = 1-l1_y-l2_y;
            
            % Get face vertex positions
            V = mesh.getAllVertices().getTrait('position');
            FV = [V(F(:,1),:) V(F(:,2),:) V(F(:,3),:)];
            FV_cell = mat2cell(FV,size(FV,1),ones(1,9));
            rep_fun = @(n,varargin)(repmat(cell2mat(varargin),n,1));
            FV_rep_x = cell2mat(arrayfun(rep_fun, num_x_is,FV_cell{:}, 'UniformOutput', 0));
            FV_rep_y = cell2mat(arrayfun(rep_fun, num_y_is,FV_cell{:}, 'UniformOutput', 0));
            
            px_1 = bsxfun(@times, FV_rep_x(:,1:3), l1_x(:,1)) + bsxfun(@times, FV_rep_x(:,4:6), l2_x(:,1)) + bsxfun(@times, FV_rep_x(:,7:9), l3_x(:,1));
            px_2 = bsxfun(@times, FV_rep_x(:,1:3), l1_x(:,2)) + bsxfun(@times, FV_rep_x(:,4:6), l2_x(:,2)) + bsxfun(@times, FV_rep_x(:,7:9), l3_x(:,2));
            
            py_1 = bsxfun(@times, FV_rep_y(:,1:3), l1_y(:,1)) + bsxfun(@times, FV_rep_y(:,4:6), l2_y(:,1)) + bsxfun(@times, FV_rep_y(:,7:9), l3_y(:,1));
            py_2 = bsxfun(@times, FV_rep_y(:,1:3), l1_y(:,2)) + bsxfun(@times, FV_rep_y(:,4:6), l2_y(:,2)) + bsxfun(@times, FV_rep_y(:,7:9), l3_y(:,2));
            
            v1 = reshape([px_1 px_2]',3,[])';
            e1 = reshape(1:(2*size(px_1,1)),2,[])';
            v2 = reshape([py_1 py_2]',3,[])';
            e2 = reshape(1:(2*size(py_1,1)),2,[])';
        end
    end
end