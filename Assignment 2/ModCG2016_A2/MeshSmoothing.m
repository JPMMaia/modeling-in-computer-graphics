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
            
            % Compute X^(n+1) = (I + lambda*L) * X^(n):
            V_smooth = iPlusLamdaL * verticesPositions;
            
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
            
            % Get vertex count:
            vertexCount = mesh.num_vertices;
            
            % Get all vertices positions:
            vertexPositions = mesh.getAllVertices().getTrait('position');
            
            % Create Wl matrix:
            Wl = wl .* eye(vertexCount, vertexCount);
            
            % Create Wp matrix:
            Wp = wp .* eye(vertexCount, vertexCount);
            
            % Solve the system, with f = 0
            % [ W_L * L ]        [  W_L * f  ]
            % [ ------- ] V'_d = [ --------- ]
            % [   W_p   ]        [ W_p * V_d ]
            A = [ Wl * L ; Wp ];
            b = [ zeros(vertexCount, 3) ; Wp * vertexPositions ];
            V_smooth = inv(A' * A) * A' * b;
            
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
            
            % Get vertex count:
            vertexCount = mesh.num_vertices;
            
            % Get all vertices positions:
            vertexPositions = mesh.getAllVertices().getTrait('position');
            
            % Create Wl matrix:
            Wl = wl .* eye(vertexCount, vertexCount);
            
            % Create Wp matrix:
            Wp = wp .* eye(vertexCount, vertexCount);
            
            % For detail preserving triangle shape optimization, f =
            % delta_(d,c):
            f = L_Cotangent * vertexPositions;
            
            % Solve the system, with f = delta_(d,c)
            % [ W_L * L ]        [  W_L * f  ]
            % [ ------- ] V'_d = [ --------- ]
            % [   W_p   ]        [ W_p * V_d ]
            A = [ Wl * L_uniform ; Wp ];
            b = [ Wl * f ; Wp * vertexPositions ];
            V_smooth = inv(A' * A) * A' * b;
            
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

            % Spectral theorem:
            % A = Q * D * inv(Q)
            % where
            % - D is as diagonal matrix and the entries of D are the
            % eigenvalues of A;
            % - Q is a unitary matrix in which the column vectors are the
            % eigenvectors of A
            
            % Calculate the eigenvectors associated with the k smallest
            % magnitude eigenvalues of the matrix L:
            [eigenvectors, ~] = eigs(L, k, 'sm');
            
            % In some cases (i. e. it can happen if the the laplacian
            % matrix is not symmetric), the eigenvectors become imaginary
            % numbers. In this case, return the original mesh:
            if(~isreal(eigenvectors))
                V_smooth = mesh.getAllVertices().getTrait('position');
                return;
            end
            
            % Get all vertices:
            vertexPositions = mesh.getAllVertices().getTrait('position');
            
            % Project the vertex positions onto the basis spanned by the
            % calculated eigenvectors:
            projectedVertexPositions = eigenvectors' * vertexPositions;
            
            % Reconstruct the filtered version of the mesh:
            V_smooth = eigenvectors * projectedVertexPositions;
            
        end
    end
end