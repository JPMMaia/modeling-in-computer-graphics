classdef MeshHelper < handle
    methods(Static)
        
        function rowVectors = normalizeRowVectors(rowVectors)
            
            rowVectors = rowVectors ./ MeshHelper.lengthOfRowVectors(rowVectors);
            
        end
        
        function rowVectorsLength = lengthOfRowVectors(rowVectors)
           
            rowVectorsLength = sqrt(sum(rowVectors.^2, 2));
            
        end
        
        function sum = sumRowVectorsByIndex(indices, rowVectors)
            
            % Replicate the row and column indices:
            numColumns = size(rowVectors, 2);
            subscripts = [ 
                repmat(indices(:), numColumns, 1) ...
                kron(1:numColumns, ones(1, numel(indices))).'
                ];
            
            % Sum all row vectors by index:
            sum = accumarray(subscripts, rowVectors(:));
            
        end
        
        function array = removeValuesFromArray(array, values)
            
            array(intersect(array, values)) = [];
            array = setdiff(array, values);
            
        end
        
        function boundaryHalfedges = getAllBoundaryHalfedges(mesh)
            
            % Get all halfedges:
            halfedges = mesh.getAllHalfedges();
            
            % For boundary halfedges, the faces returned have index 0:
            incidentFaces = halfedges.face();
            
            % Find all faces' indices which have index 0:
            boundaryFacesIndices = incidentFaces.index(:, 1) == 0;
            
            % Get all boundary halfedges indices using the indices of the faces:
            boundaryHalfedgesIndices = halfedges.index(1, boundaryFacesIndices');
            
            % If there are no boundary halfedges, then the model has no
            % boundaries:
            if(isempty(boundaryHalfedgesIndices))
                boundaryHalfedges = [];
                return;
            end
            
            % Get all boundary halfedges using the previous calculated indices:
            boundaryHalfedges = mesh.getHalfedge(boundaryHalfedgesIndices);
            
        end
        
        function nonBoundaryHalfedges = getAllNonBoundaryHalfedges(mesh)
           
            % Get all halfedges:
            halfedges = mesh.getAllHalfedges();
            
            % Get all boundary halfedges:
            boundaryHalfedges = MeshHelper.getAllBoundaryHalfedges(mesh);
            
            if(isempty(boundaryHalfedges))
                nonBoundaryHalfedges = halfedges;
                return;
            end
            
            % Remove all boundary halfedges:
            nonBoundaryHalfedgesIndices = MeshHelper.removeValuesFromArray(halfedges.index, boundaryHalfedges.index);
            nonBoundaryHalfedges = mesh.getHalfedge(nonBoundaryHalfedgesIndices);
            
        end
        
        
        function [p_min, p_max] = getBoundingBox(mesh)
            % Returns the points with minimal and maximal coordinates of the
            % smallest axis-aligned bounding box of the mesh.

            % TODO_A1 Task 1
            % 
            % Find the axis-aligned bounding box of the mesh and return its
            % minimal and maximal corner vertices. Use the vertex trait
            % 'position' to find them.
            
            % Get all vertices' positions:
            vertices = mesh.getAllVertices();
            positions = vertices.getTrait('position');
            
            % As 'positions' is a matrix, p_min and p_max will be row vectors
            % containing the min and max of each collumn:
            p_min = min(positions);
            p_max = max(positions);
            
        end
        
        
        function [V_start, V_end] = getBoundaryEdges(mesh)
            % Returns a list of line segments describing the boundary of
            % the mesh. Returns two nbe-by-3 arrays (nbe=number of
            % boundary edges), such that the i-th row of V_start and the
            % i-th row of V_end describe the two end points of the ith boundary
            % edge. The order of boundary edges is arbitrary.

            % TODO_A1 Task 2
            % 
            % Find all boundary edges of the mesh. You can achieve this by
            % finding all halfedges that do not have an incident face (i.e. its
            % face index equals zero). Make sure to test your
            % implementation for meshes with and without boundary.

            % Get all boundary halfedges:
            boundaryHalfedges = MeshHelper.getAllBoundaryHalfedges(mesh);
            
            % If there are no boundary halfedges, then the model has no
            % boundaries:
            if(isempty(boundaryHalfedges))
                V_start = [];
                V_end = [];
                return;
            end
            
            % Get all starting and ending vertices of the boundary halfedges:
            startVertices = boundaryHalfedges.from();
            endVertices = boundaryHalfedges.to();
            
            % Output the results:
            V_start = startVertices.getTrait('position');
            V_end = endVertices.getTrait('position');
            
        end
        
        
        function calculateFaceTraits(mesh)
            % Fills in a number of face traits in the TriangleMesh mesh.
            % Each face stores its surface area (trait 'area'), its
            % centroid, which is the arithmetic mean of its corner vertices
            % (trait 'centroid), and its normal (trait 'normal')

            % TODO_A1 Task 3
            % 
            % Fill in the face traits a) 'area', b) 'centroid', and  c) 'normal'.
            % 'area' is the surface area of a triangular face. 'centroid' is
            % the mean of the three corner vertices. 'normal' is the uniquely
            % defined outwards-facing normal of the face, given CCW winding of
            % the three corner vertices.
            
            % Get all faces of the mesh:
            faces = mesh.getAllFaces();
            
            % Get the halfedge of each face:
            halfedges = faces.halfedge();
            
            % Get the vertices positions of each face (triangle):
            vertexPosition1 = halfedges.from().getTrait('position');
            vertexPosition2 = halfedges.to().getTrait('position');
            vertexPosition3 = halfedges.next().to().getTrait('position');
            
            % Get the edges of the triangle:
            edge1 = vertexPosition2 - vertexPosition1;
            edge2 = vertexPosition3 - vertexPosition2;
            edge3 = vertexPosition1 - vertexPosition3;
            
            % Calculate the length of each length of the triangle:
            edge1Length = sqrt(sum((edge1).^2, 2));
            edge2Length = sqrt(sum((edge2).^2, 2));
            edge3Length = sqrt(sum((edge3).^2, 2));
            
            % Calculate the area of the triangle, using the Heron's
            % Formula:
            halfPerimeter = (edge1Length + edge2Length + edge3Length) ./ 2;
            area = sqrt(halfPerimeter .* (halfPerimeter - edge1Length) .* (halfPerimeter - edge2Length) .* (halfPerimeter - edge3Length));
            
            % Calculate the centroid of each triangle, which is just the
            % arithmetic mean of its vertices:
            centroid = (vertexPosition1 + vertexPosition2 + vertexPosition3) ./ 3;
            
            % Calculate the normal of the face, using the cross product and
            % then normalize them:
            normals = cross(edge1, edge2, 2);
            normals = normals ./ sqrt(sum(normals.^2, 2));
            
            % Write traits:
            faces.setTrait('area', area);
            faces.setTrait('centroid', centroid);
            faces.setTrait('normal', normals);
                        
        end
        
        function calculateVertexTraits(mesh)
            % Computes the degree of each vertex and stores it in the
            % vertex trait 'degree'.
            v = mesh.getAllVertices();
            he1 = v.halfedge();
            he_current = he1.twin().next();
            degs = zeros(mesh.num_vertices,1);
            i=1;
            while any(degs==0)
                degs(degs==0 & he_current.index == he1.index) = i;
                he_current = he_current.twin().next();
                i = i+1;
            end
            v.setTrait('degree',degs);
        end
        
        
        function calculateHalfedgeTraits(mesh)
            % Computes the 'angle' halfedge trait, which gives the angle
            % between the halfedge and its previous halfedge in radians.
            
            % TODO_A1 Task 5a
            %
            % Calculate the angle between a halfedge and its previous halfedge.
            % Store the resulting angle in the halfedge trait 'angle'.
            
            % Get current halfedges and its corresponding previous ones:
            currentHalfedges = mesh.getAllHalfedges();
            previousHalfedges = currentHalfedges.prev();
            
            % Get the three vertices positions that form the two edges (the
            % two edges share a vertex):
            vertexPosition1 = previousHalfedges.from().getTrait('position');
            vertexPosition2 = currentHalfedges.from().getTrait('position');
            vertexPosition3 = currentHalfedges.to().getTrait('position');
            
            % Calculate the two edge vectors:
            edgeVector1 = vertexPosition2 - vertexPosition1;
            edgeVector2 = vertexPosition3 - vertexPosition2;
            
            % Calculate the angle between two edges using the geometric
            % definition of the dot product:
            % (e1 . e2) = cos(e1^e2) * ||e1|| * ||e2|| <=>
            % cos(e1^e2) = (e1 . e2) ./ (||e1|| * ||e2||) <=>
            % angle(e1^2) = acos((e1 . e2) ./ (||e1|| * ||e2||))
            angles = acos(sum(edgeVector1 .* edgeVector2, 2) ./ (MeshHelper.lengthOfRowVectors(edgeVector1) .* MeshHelper.lengthOfRowVectors(edgeVector2)));
            
            % Output result:
            currentHalfedges.setTrait('angle', angles);
            
        end
        
        
        function calculateVertexNormals(mesh, weighting)
            % Computes vertex normals as a weighted mean of face normals.
            % The parameter 'weighting' can be one of the following:
            % 'area': The face normal weights equal the face surface areas.
            % 'angle': The face normal weights equal the opening angle
            %    of the face at the vertex.
            % Store the results in the vertex trait 'normal'.
            if nargin<2
                weighting='area';
            end
            switch weighting
                case 'area'
                    % TODO_A1 Task 4
                    % 
                    % Fill in the 'area' branch of this function. Calculate
                    % the vertex normals as weighted averages of the
                    % adjacent face normals, where the weight is given by
                    % the surface area of the face. Don't forget to
                    % normalize!                
                    
                    % Get all non-boundary halfedges:
                    halfedges = MeshHelper.getAllNonBoundaryHalfedges(mesh);
                                        
                    % Each non-boundary halfedge has an incident face:
                    faces = halfedges.face();

                    % Calculate the normal associated with each halfedge by
                    % multiplying the normal and the area of its incident
                    % face:
                    normals = faces.getTrait('normal') .* faces.getTrait('area');
                    
                    % Each halfedge has an origin vertex:
                    vertices = halfedges.from();
                
                    % Sum all normals associated with each vertex and normalize them:
                    normalsByVertex = MeshHelper.sumRowVectorsByIndex(vertices.index, normals);
                    normalsByVertex = MeshHelper.normalizeRowVectors(normalsByVertex);
                    
                    % Output result:
                    mesh.getAllVertices().setTrait('normal', normalsByVertex);
                    
                case 'angle'
                    % TODO_A1 Task 5b
                    %
                    % Fill in the 'angle' branch of this function.
                    % Calculate the vertex normals as weighted averages of
                    % the adjacent face normals, where the weight is given
                    % by the angle that the face confines at the vertex.
                    % Use the 'angle' halfedge trait computed in Task 5a for this.
                    
                    % Get all non-boundary halfedges:
                    halfedges = MeshHelper.getAllNonBoundaryHalfedges(mesh);
                                        
                    % Each non-boundary halfedge has an incident face:
                    faces = halfedges.face();

                    % Calculate the normal associated with each halfedge by
                    % multiplying the normal of its incident face by the
                    % angle of the halfedge:
                    normals = faces.getTrait('normal') .* halfedges.getTrait('angle');
                    
                    % Each halfedge has an origin vertex:
                    vertices = halfedges.from();
                
                    % Sum all normals associated with each vertex and normalize them:
                    normalsByVertex = MeshHelper.sumRowVectorsByIndex(vertices.index, normals);
                    normalsByVertex = MeshHelper.normalizeRowVectors(normalsByVertex);
                    
                    % Output result:
                    mesh.getAllVertices().setTrait('normal', normalsByVertex);

            end
        end
    end
end