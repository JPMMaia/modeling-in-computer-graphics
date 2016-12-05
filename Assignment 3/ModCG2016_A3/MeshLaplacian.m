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

            % Get the number of vertices:
            vertexCount = mesh.num_vertices;
            
            % Get all halfedges:
            halfedges = mesh.getAllHalfedges();
            
            % Get the alpha_ij and beta_ji angles (the angles of boundary
            % faces are set to 0):
            alphaiMinus1 = halfedges.twin().next().getTrait('angle') .* (halfedges.twin().face().index > 0);
            alphai = halfedges.getTrait('angle') .* (halfedges.face().index > 0);
            
            % Calculate the vector from vi to v0:
            viv0 = halfedges.to().getTrait('position') - halfedges.from().getTrait('position');
            
            % Compute the length of the vector from vi to v0:
            lengthOfViv0 = sqrt(sum(viv0.^2, 2));
            
            % Compute the weight values according to Floater, Michael S.
            % "Mean value coordinates." Computer aided geometric design 20,
            % no. 1 (2003): 19-27.:
            weightValues = (tan(alphaiMinus1 ./ 2.0) + tan(alphai / 2.0)) ./ lengthOfViv0;
            
            % Create a laplacian matrix with the weight values:
            weightValuesMatrix = sparse(halfedges.from().index, halfedges.to().index, weightValues);
            
            if normalized
               
                % Normalize laplacian matrix by diving each row by the sum
                % of all elements in that row. Then add -1 to the diagonal
                % of the matrix:
                L = bsxfun(@rdivide, weightValuesMatrix, full(sum(weightValuesMatrix, 2))) - speye(vertexCount);
                
            else
                
                % Set the diagonal values to the sum of each row (so that
                % the sum of each row is equal to 0):
                L = weightValuesMatrix - sparse((1:vertexCount)', (1:vertexCount)', full(sum(weightValuesMatrix, 2)));
                
            end
            
        end
    end
end