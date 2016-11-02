classdef MeshSmoothing < handle
    methods(Static)
        function V_smooth = explicitSmoothing(mesh, L, lambda)
            % Computes a forward Euler step to smooth the mesh with a
            % Laplacian matrix L, and a factor of lambda.

            % TODO_A2 Task 2a
            %
            % Perform explicit mesh smoothing, as described in the
            % slides and in [Desbrun1999], Section 2.2.

            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = implicitSmoothing(mesh, L, lambda)
            % Computes a backward Euler step to smooth the mesh with a
            % Laplacian matrix L, and a factor of lambda.

            % TODO_A2 Task 2b
            %
            % Perform implicit mesh smoothing, as described in the
            % slides and in [Desbrun1999], Section 2.3.

            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = lsqSmoothing(mesh,L,wl,wp)
            % Performs least-squares mesh smoothing as described in
            % [Nealen2006]. wl are weights for the smoothing rows, and wp
            % the weights for the shape-preserving rows.
            
            % TODO_A2 Task 3
            %
            % Implement least-squares mesh smoothing as described in
            % the slides and in [Nealen2006].

            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = triangleSmoothing(mesh,L_uniform,L_Cotangent,wl,wp)
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