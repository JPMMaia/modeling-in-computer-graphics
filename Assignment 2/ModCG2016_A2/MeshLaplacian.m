classdef MeshLaplacian < handle
    methods(Static)
        
        function rowVectorsLength = lengthOfRowVectors(rowVectors)
           
            rowVectorsLength = sqrt(sum(rowVectors.^2, 2));
            
        end
        
        function [angle] = computeAngleBetweenVectors(vector1, vector2)
            
            numerator = dot(vector1, vector2, 2);
            denominator = MeshLaplacian.lengthOfRowVectors(vector1) .* MeshLaplacian.lengthOfRowVectors(vector2);
            angle = acos(numerator ./ denominator);
            
        end
        
        function [alpha, beta] = computeAngles(halfedge)
           
            % Get all needed vertices (vi, vj, vj-1, vj+1):
            vi = halfedge.from().getTrait('position');
            vj = halfedge.to().getTrait('position');
            vjMinus1 = halfedge.prev().from().getTrait('position');
            vjPlus1 = halfedge.twin().next().to().getTrait('position');
            
            % Calculate alpha which is the angle between the vectors [vj-1
            % to vi] and [vj-1 to vj]:
            vjMinus1ToVi = vi - vjMinus1;
            vjMinus1ToVj = vj - vjMinus1;
            alpha = MeshLaplacian.computeAngleBetweenVectors(vjMinus1ToVi, vjMinus1ToVj);
            
            % Calculate beta which is the angle between the vectors [vj+1
            % to vi] and [vj+1 to vj]:
            vjPlus1ToVi = vi - vjPlus1;
            vjPlus1ToVj = vj - vjPlus1;
            beta = MeshLaplacian.computeAngleBetweenVectors(vjPlus1ToVi, vjPlus1ToVj);
            
        end
        
        function L = computeUniformLaplacian(mesh, normalized)
            % Returns the mesh Laplacian with uniform weights as a sparse
            % nv-by-nv matrix, where nv is the number of vertices in the mesh.
            % The sign convention is chosen such that diagonal entries are
            % negative and off-diagonal entries are positive.
            % If the 'normalized' flag is set, each row is normalized such
            % that the diagonal entry is -1, and the sum of each row is 0.
            
            % TODO_A2 Task 1a
            %
            % Compute the normalized mesh Laplacian with uniform 
            % weights. Use the sparse constructor to specify the row
            % indices, column indices, and values of non-zeros matrix
            % entries. For a detailed defintion, see the lecture slides
            % and [Nealen2006].
            
            % The laplacian matrix is a 'n' x 'n' matrix, in which 'n'
            % denotes the number of vertices:
            numberOfVertices = mesh.num_vertices;
            
            % L(i, j) = -1, if i = j:
            diagonalRows = [ 1:numberOfVertices ];
            diagonalColumns = [ 1:numberOfVertices ];
            diagonalValues = -ones(1, numberOfVertices);
            
            % L(i, j) = w[i,j], if i and j belong to the same edge.
            % For uniform laplacian (non-normalized), w[i, j] = 1:
            halfedges = mesh.getAllHalfedges(); 
            weightRows = halfedges.from().index';
            weightColumns = halfedges.to().index';
            weightValues = ones(1, mesh.num_halfedges);
            
            % TODO_A2 Task 5a
            %
            % Extend this method to compute a non-normalized version of
            % the mesh Laplacian with uniform weights.
            
            % Calculate the non-normalized version of the mesh Laplacian:
            L = sparse([diagonalRows weightRows], [diagonalColumns weightColumns], [diagonalValues weightValues], numberOfVertices, numberOfVertices);
            
            % If the 'normalized' flag is set, then normalize the matrix:
            if normalized
                
                % Normalize matrix. The sum of each row should be 0. The
                % sum of each row of the non-normalized matrix is equal to
                % the sum of all weights plus the diagonal element which
                % has a value of -1 (therefore it must be offsetted by + 1).
                L = L ./ (sum(L, 2) + 1);
                
            end
            
            % TODO something must be wrong...
            
        end
        
        function L = computeCotangentLaplacian(mesh, normalized)
            % Returns the mesh Laplacian with Cotangent weights as a sparse
            % nv-by-nv matrix, where nv is the number of vertices in the mesh.
            % The sign convention is chosen such that diagonal entries are
            % negative and off-diagonal entries are positive. For details,
            % see the lecture slides and [Nealen2006].
            % If the 'normalized' flag is set, each row is normalized such
            % that the diagonal entry is -1, and the sum of each row is 0.
            
            % TODO_A2 Task 1b
            %
            % Compute the normalized mesh Laplacian with Cotangent 
            % weights. Use the sparse constructor to specify the row
            % indices, column indices, and values of non-zeros matrix
            % entries. For a detailed defintion, see the lecture slides
            % and [Nealen2006]. Note that a halfedge trait 'cot_angle'
            % is already added by MeshHelper.calculateHalfedgeTraits().

            % The laplacian matrix is a 'n' x 'n' matrix, in which 'n'
            % denotes the number of vertices:
            numberOfVertices = mesh.num_vertices;
            
            % L(i, j) = -1, if i = j:
            diagonalValues = -ones(1, numberOfVertices);
            diagonalRows = [ 1:numberOfVertices ];
            diagonalColumns = [ 1:numberOfVertices ];
            
            % L(i, j) = w[i,j], if i and j belong to the same edge.
            % For cotangent laplacian (non-normalized), w[i, j] = cot(alpha) + cot(beta):
            halfedges = mesh.getAllHalfedges();
            [alphas, betas] = MeshLaplacian.computeAngles(halfedges);
            weightValues = (alphas + betas)';
            weightRows = halfedges.from().index';
            weightColumns = halfedges.to().index';
            
            % TODO_A2 Task 5a
            %
            % Extend this method to compute a non-normalized version of
            % the mesh Laplacian with Cotangent weights.

            % Calculate the non-normalized version of the mesh Laplacian:
            L = sparse([diagonalRows weightRows], [diagonalColumns weightColumns], [diagonalValues weightValues], numberOfVertices, numberOfVertices);
            
            % If the 'normalized' flag is set, then normalize the matrix:
            if normalized
                
                % Normalize matrix. The sum of each row should be 0. The
                % sum of each row of the non-normalized matrix is equal to
                % the sum of all weights plus the diagonal element which
                % has a value of -1 (therefore it must be offsetted by + 1).
                L = L ./ (sum(L, 2) + 1);
                
            end
            
        end
    end
end