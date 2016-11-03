classdef MeshSmoothing < handle
    methods(Static)
        function V_smooth = explicitSmoothing(mesh, L, lambda)
            % Computes a forward Euler step to smooth the mesh with a
            % Laplacian matrix L, and a factor of lambda.

            % TODO_A2 Task 2a
            %
            % Perform explicit mesh smoothing, as described in the
            % slides and in [Desbrun1999], Section 2.2.
            
            % Compute (I + lambda*L):
            numberOfVertices = mesh.num_vertices;
            iPlusLamdaL = (speye(numberOfVertices) + lambda .* L);
            
            % Get all vertices positions:
            verticesPositions = mesh.getAllVertices().getTrait('position');
            
            % Compute (I + lambda*L) * x:
            V_smooth = iPlusLamdaL * verticesPositions;
            
            % TODO FIX SMOOTHING
            
        end
        
        function V_smooth = implicitSmoothing(mesh, L, lambda)
            % Computes a backward Euler step to smooth the mesh with a
            % Laplacian matrix L, and a factor of lambda.

            % TODO_A2 Task 2b
            %
            % Perform implicit mesh smoothing, as described in the
            % slides and in [Desbrun1999], Section 2.3.
            
            % Compute (I - lambda*dt*L) (ignoring dt):
            numberOfVertices = mesh.num_vertices;
            iMinusLambdaL = (speye(numberOfVertices) - lambda .* L);
            
            % Get all vertices positions:
            verticesPositions = mesh.getAllVertices().getTrait('position');
            
            % Solve the linear system: (I - lambda*dt*L) * X^(n+1) = X^(n)
            V_smooth = iMinusLambdaL \ verticesPositions;
            
        end
        
        function V_smooth = lsqSmoothing(mesh,L,wl,wp)
            % Performs least-squares mesh smoothing as described in
            % [Nealen2006]. wl are weights for the smoothing rows, and wp
            % the weights for the shape-preserving rows.
            
            % TODO_A2 Task 3
            %
            % Implement least-squares mesh smoothing as described in
            % the slides and in [Nealen2006].
            
            vertexCount = mesh.num_vertices;
            m = 10;
            
            % Get all vertices positions:
            vertexPositions = mesh.getAllVertices().getTrait('position');
            
            mVertices = vertexPositions(1:m, :);
            Wp1 = eye(m, vertexCount) .* wp;
            Wp2 = eye(m, vertexCount + m) .* wp;
            Wl1 = eye(vertexCount, vertexCount) .* wl;
            Wl2 = eye(vertexCount, vertexCount + m) .* wl;
            
            A = [ L ; Wp1 ];
            b = [ zeros(vertexCount, 3) ; mVertices ];
            
            V_smooth = inv(A' * A) * A' * b;
            
            % TODO
            
        end
        
        function V_smooth = triangleSmoothing(mesh, L_uniform, L_Cotangent, wl, wp)
            % Performs detail preserving triangle shape optimization as
            % described in [Nealen2006]. wl are the weights for the
            % triangle shape optimization rows, and wp the weights for the
            % detail-preserving rows.

            % TODO_A2 Task 4
            %
            % Implement detail preserving triangle shape optimization
            % mesh smoothing as described in the slides and in
            % [Nealen2006].
            
            

            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = spectralSmoothing(mesh, L, k)
            % Performs spectral mesh smoothing through a low-pass filtering
            % of the Laplacian eigenvectors, in which only the k lowest
            % frequencies are preserved.

            % TODO_A2 Task 5b
            %
            % Perform spectral smoothing. In order to do that, perform
            % a sparse eigendecomposition of the Laplacian L that
            % computes only the eigenvectors associated with the k
            % smallest-magnitude eigenvalues. Then project the vertex
            % positions onto the basis spanned by these eigenvectors
            % and reconstruct a filtered version of the mesh.

            V_smooth = mesh.toFaceVertexMesh();
        end
    end
end