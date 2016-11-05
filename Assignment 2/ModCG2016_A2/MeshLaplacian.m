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
        
        function L = computeLaplacian(mesh, normalized, weightValues)
           
            % The laplacian matrix is a 'n' x 'n' matrix, in which 'n'
            % denotes the number of vertices:
            numberOfVertices = mesh.num_vertices;
            
            % L(i, j) = -1, if i = j:
            diagonalRows = 1:numberOfVertices;
            diagonalColumns = 1:numberOfVertices;
            
            % L(i, j) = w[i,j], if i and j belong to the same edge:
            halfedges = mesh.getAllHalfedges(); 
            weightRows = halfedges.from().index';
            weightColumns = halfedges.to().index';
            
            if normalized
               
                % Set diagonal values to -1:
                diagonalValues = -ones(1, numberOfVertices);
                
                % Calculate a sparse matrix with the weights only:
                LWeights = sparse(weightRows, weightColumns, weightValues, numberOfVertices, numberOfVertices);
                
                % Calculate the sum all weights by row:
                sumsByRow = sum(LWeights, 2);
                
                % Normalize weights:
                LWeights = LWeights ./ sumsByRow;
                
                % Compute the normalized version of the mesh Laplacian, by
                % adding the diagonals to the normalized weights:
                L = LWeights + sparse(diagonalRows, diagonalColumns, diagonalValues, numberOfVertices, numberOfVertices);
                
            else
                
                % Build an adjacency matrix, where A(i, j) = 1 if i and j
                % belong to the same edge:
                adjacencyMatrix = sparse(weightRows, weightColumns, ones(1, mesh.num_halfedges), numberOfVertices, numberOfVertices);
                
                % Each diagonal value is equals to degree of that vertex:
                diagonalValues = -sum(adjacencyMatrix ~= 0, 1);
                
                % Calculate the non-normalized version of the mesh Laplacian, as described in [Sorkine2006]:
                L = sparse([diagonalRows weightRows], [diagonalColumns weightColumns], [diagonalValues weightValues], numberOfVertices, numberOfVertices);
                
            end
            
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
            
            % TODO_A2 Task 5a
            %
            % Extend this method to compute a non-normalized version of
            % the mesh Laplacian with uniform weights.
            
            % For uniform laplacian, w[i, j] = 1:
            weightValues = ones(1, mesh.num_halfedges);
            
            L = MeshLaplacian.computeLaplacian(mesh, normalized, weightValues);
            
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
            
            % TODO_A2 Task 5a
            %
            % Extend this method to compute a non-normalized version of
            % the mesh Laplacian with Cotangent weights.

            % For cotangent laplacian, w[i, j] = cot(alpha) + cot(beta):
            halfedges = mesh.getAllHalfedges();
            [alphas, betas] = MeshLaplacian.computeAngles(halfedges);
            weightValues = (cot(alphas) + cot(betas))';
            
            L = MeshLaplacian.computeLaplacian(mesh, normalized, weightValues);
            
        end
    end
end