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

            uv = zeros(mesh.num_vertices, 2);
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

            if adaptive_spacing
                % TODO_A3 Task 1d
                % Distribute boundary points adaptively along a circle.
                uvpos = zeros(length(bdr_hei), 2);
            else
                % TODO_A3 Task 1b
                % Distribute boundary points uniformly along a circle.
                uvpos = zeros(length(bdr_hei), 2);
            end
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

            if adaptive_spacing
                % TODO_A3 Task 1d
                % Distribute boundary points adaptively along a square.
                uvpos = zeros(length(bdr_hei), 2);
            else
                % TODO_A3 Task 1b
                % Distribute boundary points uniformly along a square.
                uvpos = zeros(length(bdr_hei), 2);
            end
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

            M = sparse(2*mesh.num_faces, 2*mesh.num_vertices);
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
            uv = zeros(size(M,2)*0.5,2);
        end
    end
end