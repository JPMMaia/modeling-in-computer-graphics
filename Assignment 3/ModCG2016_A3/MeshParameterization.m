classdef MeshParameterization < handle
    properties
    end
    
    methods(Static)
        function uv = barycentricParameterization(mesh, laplacian, bdry_shape, adaptive_spacing)
            % Computes a linear barycentric parameterization of the mesh,
            % given a mesh Laplacian.
            % bdry_shape: either 'circle' or 'square'
            % adaptive_spacing: either 0 or 1
            % uv: an nv-by-2 array, where the i-th row contains the uv
            % coordinates of the i-th vertex in mesh.
            
            bdry_hei = MeshHelper.getBoundaryLoop(mesh);
            if isempty(bdry_hei)
                error('Surface has no boundary!');
            end
            
            if strcmp(bdry_shape, 'circle')
                uv_fixed = MeshParameterization.generateCircleBoundary(mesh, bdry_hei, adaptive_spacing);
            elseif strcmp(bdry_shape, 'square')
                uv_fixed = MeshParameterization.generateSquareBoundary(mesh, bdry_hei, adaptive_spacing);
            else
                error('Bdry Shape unknown!');
            end

            % TODO_A3 Task 1c
            % Compute a mesh parameterization as a linear barycentric
            % mapping. Given a fixed uv-boundary and a Laplacian,
            % substitute the rows associated with known uv values by
            % identity rows, and put the known uv values into the
            % right-hand side vector.

            % Get all boundary halfedges:
            boundaryHalfedges = mesh.getHalfedge(bdry_hei);
            
            % Get all the boundary vertices indices. These indices
            % correspond to the rows associated with known uv values.
            knownUVIndices = boundaryHalfedges.from().index;
            
            % Find all non zero values of the laplacian matrix:
            [nonZeroLaplacianRows, nonZeroLaplacianColumns, nonZeroLaplacianValues] = find(laplacian);
            
            % Remove all the rows associated with the known UV values:
            [indicesToKeep, ~, ~] = find(~ismember(nonZeroLaplacianRows, knownUVIndices));
            nonZeroLaplacianRows = nonZeroLaplacianRows(indicesToKeep);
            nonZeroLaplacianColumns = nonZeroLaplacianColumns(indicesToKeep);
            nonZeroLaplacianValues = nonZeroLaplacianValues(indicesToKeep);
            
            % Substitute the rows associated with known uv values by
            % identity rows:
            nonZeroLaplacianRows = [ nonZeroLaplacianRows ; knownUVIndices ];
            nonZeroLaplacianColumns = [ nonZeroLaplacianColumns ; knownUVIndices ];
            nonZeroLaplacianValues = [ nonZeroLaplacianValues ; ones(size(knownUVIndices, 1), 1) ];
            
            % Create a new laplacian where the rows associated with known
            % uv values are identity rows:
            laplacian = sparse(nonZeroLaplacianRows, nonZeroLaplacianColumns, nonZeroLaplacianValues, size(laplacian, 1), size(laplacian, 2));
            
            % Put the known uv values into the right-hand side vector b:
            b = zeros(mesh.num_vertices, 2);
            b(knownUVIndices, :) = uv_fixed;
            
            % Solve the linear system laplacian . uv = b:
            uv = laplacian \ b;
            
        end
        
        function uvpos = generateCircleBoundary(mesh, bdry_hei, adaptive_spacing)
            % Distributes the boundary vertices of a mesh along the circle
            % with center (0.5,0.5) and radius 0.5.
            % If adaptive_spacing is 0, the spacing of the points is
            % equidistant. If adaptive_spacing is 1, the spacing of the
            % points is proportial to the original edge lengths between
            % boundary vertices in the mesh.
            % bdry_hei: a p-by-1 array as returned by
            % MeshHelper.getBoundaryLoop(...)
            % uvpos: a p-by-2 array, where the i-th row contains the
            % position of the i-th boundary vertex on the circle.

            % TODO_A3 Task 1d
            % Distribute boundary points adaptively along a circle.
            
            % TODO_A3 Task 1b
            % Distribute boundary points uniformly along a circle.
            
            % Define the radius and center of the sphere:
            radius = 0.5;
            center = [ 0.5, 0.5 ];

            % Get the number of halfedges:
            halfedgeCount = length(bdry_hei);

            if adaptive_spacing
                
                % Get the boundary halfedges:
                halfedges = mesh.getHalfedge(bdry_hei);
                
                % Calculate the length of each edge:
                edgeLength = sqrt(sum((halfedges.to().getTrait('position') - halfedges.from().getTrait('position')).^2, 2));
                
                % Calculate the length of the boundary:
                totalEdgeLength = sum(edgeLength, 1);
                
                % Calculate a weight for each point between 0 and 1:
                edgeWeight = edgeLength ./ totalEdgeLength;
                
                % Map weights from [0, 1] to [0, 2.0 * pi]:
                deltaAngle = edgeWeight .* 2.0 .* pi;
                
            else
                
                % Give an equal angle to all points in the range[0,
                % 2.0*pi]:
                deltaAngle = repmat(2.0 * pi / halfedgeCount, halfedgeCount, 1);
                
            end
            
            % Calculate the angle associated with each vertex:
            angles = cumsum(deltaAngle);

            % The position of each vertex is calculated using polar
            % coordinates:
            % [centerX + radius * cos(angle), centerY + radius * sin(angle)]
            uvpos = repmat(center, halfedgeCount, 1)  + radius .* [ cos(angles) , sin(angles) ];
            
        end
        
        function uvpos = generateSquareBoundary(mesh, bdry_hei, adaptive_spacing)
            % Distributes the boundary vertices of a mesh along the sides
            % of the unit square.
            % If adaptive_spacing is 0, the spacing of the points is
            % equidistant. If adaptive_spacing is 1, the spacing of the
            % points is proportial to the original edge lengths between
            % boundary vertices in the mesh.
            % bdry_hei: a p-by-1 array as returned by
            % MeshHelper.getBoundaryLoop(...)
            % uvpos: a p-by-2 array, where the i-th row contains the
            % position of the i-th boundary vertex on the square.

            % TODO_A3 Task 1d
            % Distribute boundary points adaptively along a square.
            
            % TODO_A3 Task 1b
            % Distribute boundary points uniformly along a square.
            
            % Define the limits of the square boundary:
            squareBoundaryLimits = [ 0.0, 0.0 ; 1.0, 0.0 ; 1.0, 1.0 ; 0.0, 1.0 ; 0.0, 0.0 ];

            % Get the number of halfedges:
            halfedgeCount = length(bdry_hei);

            if adaptive_spacing
                
                % Get the boundary halfedges:
                halfedges = mesh.getHalfedge(bdry_hei);
                
                % Calculate the length of each edge:
                edgeLength = sqrt(sum((halfedges.to().getTrait('position') - halfedges.from().getTrait('position')).^2, 2));
                
                % Calculate the length of the boundary:
                totalEdgeLength = sum(edgeLength, 1);
                
                % Calculate a weight for each point between 0 and 1:
                edgeWeight = edgeLength ./ totalEdgeLength;
                
                % Map from [0, 1] to [0, 4]:
                boundaryWeights = 4.0 .* edgeWeight;
                
            else
                
                % Give an equal weight for each point in the range[0, 4]:
                boundaryWeights = ones(halfedgeCount, 1) .* (4 / halfedgeCount);
                
            end
            
            % Create a row vector with values in range [1, 5]. This
            % represents the position of the vertex in a stretched
            % boundary of length 4 (the perimeter of the square):
            boundaryPosition = 1.0 + cumsum(boundaryWeights);

            % Calculate the floor and ceil of the boundary position.
            % For instance, if boundaryPosition = 1.4, the floor will
            % be 1 and the ceil will be 2:
            floorLimit = floor(boundaryPosition);
            ceilLimit = ceil(boundaryPosition);

            % Calculate the position along an edge. For instance, if
            % boundaryPosition = 1.4, the positionAlongEdge = 0.4.
            positionAlongEdge = boundaryPosition - floorLimit;

            % Get the vertices belonging to each edge:
            vertexEdge0 = squareBoundaryLimits(int64(floorLimit), :);
            vertexEdge1 = squareBoundaryLimits(int64(ceilLimit), :);

            % The position is calculated by interpolation between the
            % two vertices of the edge using the positionAlongEdge
            % (which has a value between 0 and 1):
            uvpos = (1.0 - positionAlongEdge) .* vertexEdge0 + positionAlongEdge .* vertexEdge1;
            
        end
        
        function normalizedVector = normalize(vector)
            
            normalizedVector = vector ./ sqrt(sum(vector .^ 2, 2));
            
        end
        
        function [orthonormalBasisX, orthonormalBasisY]  = computeOrthonormalBasis(halfedge)
            
            % Get each point of the triangle:
            xi = halfedge.from().getTrait('position');
            xj = halfedge.to().getTrait('position');
            xk = halfedge.next().to().getTrait('position');
            
            % Create the orthonormal basis vector X:
            orthonormalBasisX = (xj - xi);
            orthonormalBasisX = MeshParameterization.normalize(orthonormalBasisX);
            
            % Compute the normal to the basis:
            basisNormal = cross(orthonormalBasisX, (xk - xi), 2);
            basisNormal = MeshParameterization.normalize(basisNormal);
            
            % Create the orthonormal basis vector Y:
            orthonormalBasisY = cross(basisNormal, orthonormalBasisX, 2);
            
        end
        
        function M = computeM(mesh)
            
            % Get the halfedge associated with each face:
            halfedge = mesh.getAllFaces().halfedge();
            
            % Compute the orthonormal basis for the given triangle:
            [basisX, basisY] = MeshParameterization.computeOrthonormalBasis(halfedge);
            
            % Get the area of the triangle associated with the halfedge:
            triangleArea = halfedge.face().getTrait('area');
            
            % Compute the scalar which is applied to each element of the
            % matrix:
            mScalar = 1.0 ./ (2.0 .* triangleArea);
            
            % Repeat each element into a 2-by-1 block of a new matrix:
            mScalar = repelem(mScalar, 2, 1);
            
            % Compute the first row:
            yjMinusYk = basisY(:, 2) - basisY(:, 3);
            ykMinusYi = basisY(:, 3) - basisY(:, 1);
            yiMinusYj = basisY(:, 1) - basisY(:, 2);
            
            % Compute the second row:
            xkMinusXj = basisX(:, 3) - basisX(:, 2);
            xiMinusXk = basisX(:, 1) - basisX(:, 3);
            xjMinusXi = basisX(:, 2) - basisX(:, 1);
            
            % Create the sparse matrix Mt:
            rowsX = repmat(2 .* halfedge.face().index', 1, 3);
            rowsY = rowsX - 1;
            columnsX = [ halfedge.from().index', halfedge.to().index', halfedge.next().to().index' ];
            columnsY = columnsX;
            values = [ yjMinusYk', ykMinusYi', yiMinusYj', xkMinusXj', xiMinusXk', xjMinusXi' ];
            Mt = mScalar .* sparse([rowsY, rowsX], [columnsY, columnsX], values, 2 * mesh.num_faces, mesh.num_vertices);
            
            % Create the sparse matrix Mt rotated:
            rotatedValues = [ -xkMinusXj', -xiMinusXk', -xjMinusXi', yjMinusYk', ykMinusYi', yiMinusYj' ];
            MtRotated = mScalar .* sparse([rowsY, rowsX], [columnsY, columnsX], rotatedValues, 2 * mesh.num_faces, mesh.num_vertices);
            
            % Compute the matrix M:
            triangleArea = repelem(triangleArea, 2, 1);
            M = [ -triangleArea .* MtRotated, triangleArea .* Mt ];
            
        end
        
        function uv = lscmParameterization(mesh, bdry_vi, bdry_uvpos)
            % Computes a mesh parameterization according to the
            % least-squares conformal map (LSCM) algorithm.
            % mesh: a triangle mesh with 1 boundary loop
            % bdry_vi: a p-by-1 vector with p>=2. Contains the indices of
            % fixed vertices.
            % bdry_uvpos: a p-by-2 matrix. The i-th row contains the uv
            % coordinates of the fixed vertex with index bdry_vi(i).
            % uv: an nv-by-2 matrix. The i-th row contains the uv
            % coordinates of the i-th vertex.
            
            % TODO_A3 Task 2a
            % Compute a mesh parameterization with the LSCM algorithm.
            % Guidelines: your first goal is to assemble a
            % (2*nf)-by-(2*nv) matrix M, where nf is the number of
            % faces, and nv is the number of vertices in the mesh.
            % Every pair of rows corresponds to the conformality energy
            % E_LCSM of a triangle. Every pair of columns corresponds
            % to the unknown u- and v-coordinates of one vertex.
            % You should assume that the unknowns of the
            % linear system defined by M have the order (u_1, v_1, u_2,
            % v_2, ..., u_nv, v_nv).
            % During the computation of M, the fixed vertices are not
            % taken into account yet! The fixed vertex boundary
            % conditions are handled by the method
            % MeshParameterization.setBCAndSolve(...), which you should
            % implement next.
            
            M = MeshParameterization.computeM(mesh); 
            
            uv = MeshParameterization.setBCAndSolve(M, bdry_vi, bdry_uvpos);
        end
        
        function uv = dcpParameterization(mesh, bdry_vi, bdry_uvpos)
            % Computes a mesh parameterization according to the
            % discrete conformal parameterization (DCP) algorithm.
            % mesh: a triangle mesh with 1 boundary loop
            % bdry_vi: a p-by-1 vector with p>=2. Contains the indices of
            % fixed vertices.
            % bdry_uvpos: a p-by-2 matrix. The i-th row contains the uv
            % coordinates of the fixed vertex with index bdry_vi(i).
            % uv: an nv-by-2 matrix. The i-th row contains the uv
            % coordinates of the i-th vertex.

            % TODO_A3 Task 3
            % Compute a mesh parameterization with the DCP algorithm.
            % The assembly of the DCP system matrix M proceeds in a way
            % similar to that of LSCM. One difference is that M is
            % square this time, such that both the number of rows and
            % the number of columns equal 2*nv. Each pair of rows
            % corresponds to the derivative of the Dirichlet energy
            % around the 1-ring of a vertex, see Eq. 7 in [Desbrun et
            % al. 2002] and the lecture slides.

            M = sparse(2*mesh.num_vertices, 2*mesh.num_vertices);
            uv = MeshParameterization.setBCAndSolve(M, bdry_vi, bdry_uvpos);
        end
        
        function uv = setBCAndSolve(M, bdry_vi, bdry_uvpos)
            % Substitutes boundary conditions into a linear system and
            % solves Mu = 0.
            % M is a linear system matrix in (2*nv) unknowns that
            % correspond to the u- and v-coordinates of all nv vertices in
            % a mesh. Therefore M has 2*nv columns, and number of rows >=
            % 2*nv. The unknowns are assumed to have the order (u_1, v_1,
            % u_2, v_2, ..., u_nv, v_nv).
            % bdry_vi is a list of p indices corresponding to fixed
            % vertices, with p >= 2. bdry_uvpos is a p-by-2 matrix, where
            % the i-th row contains the fixed uv coordinates of the vertex
            % with index bdry_vi(i).
            % This puts the columns of M that correspond to fixed u- and
            % v-coordinates to the right-hand side by multiplying them with
            % the corresponding fixed coordinate values. This reduced the
            % number of columns in the linear system by 2*p.
            % The system is then solved, which yields the uv-coordinates of
            % the free (that is, non-fixed) vertices. The uv-coordinates of
            % free and fixed vertices are then recombined in an nv-by-2
            % matrix, such that the i-th row contains the uv-coordinates of
            % the i-th vertex in the mesh, and returned.
            
            % TODO_A3 Task 2b
            % Given a linear system, substitute fixed vertex boundary
            % conditions and solve.
            
            % Calculate vertex count, fixed uv count and unfixed uv count:
            vertexCount = size(M,2)*0.5;
            fixedUVCount = size(bdry_vi, 1);
            unfixedUVCount = vertexCount - fixedUVCount;
            
            % Set which columns from the M matrix should be selected (u and
            % v):
            columnIndices = [ bdry_vi; bdry_vi + 5 ];
            
            % Select the columns of the known uv positions:
            mUvColumns = M(:, columnIndices);
            
            % Remove the columns of the known uv positions from the M
            % matrix:
            M(:, columnIndices) = [];
            
            % Compute the right hand side of the equation -A_l * UV_fixed:
            fixedUV = reshape(bdry_uvpos, 2 * fixedUVCount, 1);
            b = -mUvColumns * fixedUV;
            
            % Solve the linear system A_f * UV_nonFixed = b:
            nonFixedUV = M \ b;
           
            % Recombine fixed and non-fixed UV values into a matrix nv-by-2 matrix:
            rows = [ repmat(setdiff(1:vertexCount, bdry_vi'), 1, 2), repmat(bdry_vi', 1, 2) ]; 
            columns = [ ones(1, unfixedUVCount), 2 .* ones(1, unfixedUVCount), ones(1, fixedUVCount), 2 .* ones(1, fixedUVCount) ];
            values = [ nonFixedUV', fixedUV' ];
            uv = full(sparse(rows, columns, values, vertexCount, 2));
            
        end
    end
end