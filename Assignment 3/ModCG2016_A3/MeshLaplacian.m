classdef MeshLaplacian < handle
    methods(Static)
        function L = computeUniformLaplacian(mesh, normalized)
            % Returns the mesh Laplacian with uniform weights as a sparse
            % nv-by-nv matrix, where nv is the number of vertices in mesh.
            % The sign convention is chosen such that diagonal entries are
            % negative and off-diagonal entries are positive.
            % If the 'normalized' flag is set, each row is normalized such
            % that the diagonal entry is -1, and the sum of each row is 0.

            nv = mesh.num_vertices;
            nhe = mesh.num_halfedges;
            he = mesh.getAllHalfedges();
            if normalized
                L = sparse([(1:nv)'; he.from().index], ...
                    [(1:nv)'; he.to().index], ...
                    [-ones(nv,1); 1 ./ he.from().getTrait('degree')], ...
                    nv,nv);
            else
                L = sparse([(1:nv)'; he.from().index], ...
                    [(1:nv)'; he.to().index], ...
                    [-mesh.getAllVertices().getTrait('degree'); ones(nhe,1)]);
            end
        end
        
        function L = computeCotangentLaplacian(mesh, normalized)
            % Returns the mesh Laplacian with Cotangent weights as a sparse
            % nv-by-nv matrix, where nv is the number of vertices in mesh.
            % The sign convention is chosen such that diagonal entries are
            % negative and off-diagonal entries are positive. For details,
            % see the lecture slides and [Nealen2006].
            % If the 'normalized' flag is set, each row is normalized such
            % that the diagonal entry is -1, and the sum of each row is 0.
            
            nv = mesh.num_vertices;
            he = mesh.getAllHalfedges();
            val_od = he.prev().getTrait('cot_angle') .* (he.face().index > 0) + ...
                he.twin().prev().getTrait('cot_angle') .* (he.twin().face().index > 0);
            L_od = sparse(he.from().index, he.to().index, val_od);
            if normalized
                L = bsxfun(@rdivide, L_od, full(sum(L_od,2))) - speye(nv);
            else
                L = L_od - sparse((1:nv)',(1:nv)',full(sum(L_od,2)));
            end
        end
        
        function L = computeMeanValueLaplacian(mesh, normalized)
            % Returns the mesh Laplacian with mean value weights as a sparse
            % nv-by-nv matrix, where nv is the number of vertices in mesh.
            % The sign convention is chosen such that diagonal entries are
            % negative and off-diagonal entries are positive. For details,
            % see the lecture slides.
            % If the 'normalized' flag is set, each row is normalized such
            % that the diagonal entry is -1, and the sum of each row is 0.

            % TODO_A3 Task 1e
            % Implement the Laplacian with mean value weights.

            nv = mesh.num_vertices;
            if normalized
                L = sparse(nv, nv);
            else
                L = sparse(nv, nv);
            end
        end
    end
end