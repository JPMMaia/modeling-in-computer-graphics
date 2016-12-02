classdef MeshSmoothing < handle
    methods(Static)
        function V_smooth = explicitSmoothing(mesh, L, lambda)
            % Computes a forward Euler step to smooth the mesh with a
            % Laplacian matrix L, and a factor of lambda.

            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = implicitSmoothing(mesh, L, lambda)
            % Computes a backward Euler step to smooth the mesh with a
            % Laplacian matrix L, and a factor of lambda.
            
            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = lsqSmoothing(mesh,L,wl,wp)
            % Performs least-squares mesh smoothing as described in
            % [Nealen2006]. wl are weights for the smoothing rows, and wp
            % the weights for the shape-preserving rows.
            
            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = triangleSmoothing(mesh,L_uniform,L_Cotangent,wl,wp)
            % Performs detail preserving triangle shape optimization as
            % described in [Nealen2006]. wl are the weights for the
            % triangle shape optimization rows, and wp the weights for the
            % detail-preserving rows.
            
            V_smooth = mesh.toFaceVertexMesh();
        end
        
        function V_smooth = spectralSmoothing(mesh, L, k)
            % Performs spectral mesh smoothing through a low-pass filtering
            % of the Laplacian eigenvectors, in which only the k lowest
            % frequencies are preserved.
            
            V_smooth = mesh.toFaceVertexMesh();
        end
    end
end