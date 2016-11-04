classdef MeshHelper < handle
    methods(Static)
        
        function b = hasBoundary(mesh)
            % Returns 1 if the mesh has a boundary, and 0 otherwise.
            b = any(mesh.getAllHalfedges().face().index==0);
        end
        
        function magnitude = computeMagnitude(vector)
            
            magnitude = sqrt(sum(vector.^2, 2));
            
        end
        
        function volume = computeSignedVolumeOfTriangle(face)
            
            % Each face has a halfedge associated with it:
            halfedge = face.halfedge();
            
            % Get the position of each vertex of the face:
            vertexPosition1 = halfedge.from().getTrait('position');
            vertexPosition2 = halfedge.to().getTrait('position');
            vertexPosition3 = halfedge.next().to().getTrait('position');
            
            % Calculate the volume of the triangle as described in
            % [zhang2001efficient]
            volume = ...
            - vertexPosition3(:, 1) .* vertexPosition2(:, 2) .* vertexPosition1(:, 3) ...
            + vertexPosition2(:, 1) .* vertexPosition3(:, 2) .* vertexPosition1(:, 3) ...
            + vertexPosition3(:, 1) .* vertexPosition1(:, 2) .* vertexPosition2(:, 3) ...
            - vertexPosition1(:, 1) .* vertexPosition3(:, 2) .* vertexPosition2(:, 3) ...
            - vertexPosition2(:, 1) .* vertexPosition1(:, 2) .* vertexPosition3(:, 3) ...
            + vertexPosition1(:, 1) .* vertexPosition2(:, 2) .* vertexPosition3(:, 3);
            
        end
        
        function vol = computeVolume(mesh)
            % Returns the signed volume of the mesh.

            % TODO_A2 Task 2c
            %
            % Compute the signed volume of the triangle mesh for the
            % volume preservation code to work. See [Desbrun1999],
            % Section 3 for details on how volume preservation is
            % implemented.
            % Do your own research on how to compute the volume
            % of a triangle mesh.
            
            % Get all triangled faces:
            faces = mesh.getAllFaces();
            
            % Calculating the volume of a triangle mesh, as described in 
            % [Zhang, Cha, and Tsuhan Chen. "Efficient feature extraction
            % for 2D/3D objects in mesh representation." Image Processing,
            % 2001. Proceedings. 2001 International Conference on. Vol. 3.
            % IEEE, 2001]([zhang2001efficient])
            volumes = MeshHelper.computeSignedVolumeOfTriangle(faces);
            
            vol = sum(volumes, 1);
        end
        
        function scaleMesh(mesh, factor)
            % Scales the mesh by a gives factor relative to the midpoint of
            % the bounding box.
            V = mesh.toFaceVertexMesh();
            p_min = min(V,[],1);
            p_max = max(V,[],1);
            p_center = 0.5*(p_max+p_min);
            V_new = bsxfun(@plus, bsxfun(@minus, V, p_center)*factor, p_center);
            mesh.getAllVertices().setTrait('position', V_new);
        end
        
        function [p_min, p_max] = getBoundingBox(mesh)
            % Returns the points with minimal and maximal coordinates of the
            % smallest axis-aligned bounding box of the mesh.
            
            V = mesh.getAllVertices().getTrait('position');
            p_min = min(V,[],1);
            p_max = max(V,[],1);
        end
        
        
        function [V_start, V_end] = getBoundaryEdges(mesh)
            % Returns a list of line segments describing the boundary of
            % the mesh. Returns two nbe-by-3 arrays (nbe=number of
            % boundary edges), such that the i-th row of V_start and the
            % i-th row of V_end describe the two end points of a boundary
            % edge.
            
            he = mesh.getAllHalfedges();
            he_bdry = mesh.getHalfedge(he.face().index == 0);

            V_start = he_bdry.from().getTrait('position');
            V_end = he_bdry.to().getTrait('position');
        end
        
        
        function calculateFaceTraits(mesh)
            % Fills in a number of face traits in the TriangleMesh mesh.
            % Each face stores its surface area (trait 'area'), its
            % centroid, which is the arithmetic mean of its corner vertices
            % (trait 'centroid), and its normal (trait 'normal').
            
            f = mesh.getAllFaces();
            he1 = f.halfedge();
            v1 = he1.from().getTrait('position');
            v2 = he1.to().getTrait('position');
            v3 = he1.next().to.getTrait('position');
            fn_weighted = cross(v2-v1,v3-v1);
            areas = 0.5*sqrt(sum(fn_weighted .* fn_weighted, 2));
            f.setTrait('area',areas);
            f.setTrait('centroid',(v1 + v2 + v3) / 3);
            f.setTrait('normal',normr(fn_weighted));
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
            % between the halfedge and its previous halfedge in radians, as
            % well as the 'cot_angle' trait which gives its Cotangent.
            
            he = mesh.getAllHalfedges();
            v = he.from().getTrait('position');
            v_next = he.to().getTrait('position');
            v_prev = he.prev().from().getTrait('position');
            d1 = normr(v_next-v);
            d2 = normr(v_prev-v);
            angles = acos(sum(d1.*d2,2));
            he.setTrait('angle',angles);
            he.setTrait('cot_angle',cot(angles));

        end
        
        function nonBoundaryHalfedges = getAllNonBoundaryHalfedges(mesh)
            
            % Get all halfedges:
            halfedges = mesh.getAllHalfedges();
            
            % For boundary halfedges, the faces returned have index 0:
            incidentFaces = halfedges.face();
            
            % Find all faces' indices which have index != 0:
            nonBoundaryFacesIndices = incidentFaces.index(:, 1) ~= 0;
            
            % Get all non-boundary halfedges indices using the indices of the faces:
            boundaryHalfedgesIndices = halfedges.index(1, nonBoundaryFacesIndices');
            
            % Get all non-boundary halfedges using the previous calculated indices:
            nonBoundaryHalfedges = mesh.getHalfedge(boundaryHalfedgesIndices);
            
        end
        
        function isObtuse = isTriangleObtuse(angleP, angleQ, angleR)
            
            isObtuse = angleP > 90 | angleQ > 90 | angleR > 90;
            
        end
        
        function voronoiArea = computeVoronoiArea(vertexPositionP, vertexPositionQ, vertexPositionR, cotagentOfAngleQ, cotagentOfAngleR)
           
            vectorPQ = vertexPositionQ - vertexPositionP;
            vectorPR = vertexPositionR - vertexPositionP;
            
            vectorPQDistanceSquared = sum(vectorPQ .^ 2, 2);
            vectorPRDistanceSquared = sum(vectorPR .^ 2, 2);
            
            voronoiArea = ( vectorPRDistanceSquared .* cotagentOfAngleQ + vectorPQDistanceSquared .* cotagentOfAngleR) ./ 8;
            
        end
        
        function triangleArea = computeAreaOfTriangle(vertexPositionP, vertexPositionQ, vertexPositionR)
            
            % Get the edges of the triangle:
            vectorPQ = vertexPositionQ - vertexPositionP;
            vectorPR = vertexPositionR - vertexPositionP;
            vectorQR = vertexPositionR - vertexPositionQ;
            
            % Calculate the length of each length of the triangle:
            vectorPQLength = sqrt(sum((vectorPQ).^2, 2));
            vectorPRLength = sqrt(sum((vectorPR).^2, 2));
            vectorQRLength = sqrt(sum((vectorQR).^2, 2));
            
            % Calculate the area of the triangle, using the Heron's
            % Formula:
            halfPerimeter = (vectorPQLength + vectorPRLength + vectorQRLength) ./ 2;
            triangleArea = sqrt(halfPerimeter .* (halfPerimeter - vectorPQLength) .* (halfPerimeter - vectorPRLength) .* (halfPerimeter - vectorQRLength));
            
        end
        
        function aMixed = computeAMixed(halfedgeFromP)
            
            % Get remainder halfedges of the triangle:
            halfedgeFromQ = halfedgeFromP.next();
            halfedgeFromR = halfedgeFromQ.next();
            
            % Get vertex position of each vertex of the face:
            vertexPositionP = halfedgeFromP.from().getTrait('position');
            vertexPositionQ = halfedgeFromQ.from().getTrait('position');
            vertexPositionR = halfedgeFromR.from().getTrait('position');
            
            % Mark halfedges which have an associated obtuse triangle:
            isObtuse = MeshHelper.isTriangleObtuse(halfedgeFromP.getTrait('angle'), halfedgeFromQ.getTrait('angle'), halfedgeFromR.getTrait('angle'));
            
            % Compute voronoi area for non-obtuse triangles:
            aMixed = ~isObtuse .* MeshHelper.computeVoronoiArea(vertexPositionP, vertexPositionQ, vertexPositionR, halfedgeFromQ.getTrait('cot_angle'), halfedgeFromR.getTrait('cot_angle'));
            
            % Calculate the triangle areas of each triangle.
            % 'triangleAreas' will contain only the areas of the obtuse
            % triangles.
            triangleAreas = isObtuse .* MeshHelper.computeAreaOfTriangle(vertexPositionP, vertexPositionQ, vertexPositionR);
            
            % Mark halfedges which have an angle of triangle at vertex
            % obtuse:
            isAngleObtuseAtVertex = halfedgeFromP.getTrait('angle') > 90;
            
            % If the angle of triangle at vertex is obtuse, add area(T) /
            % 2:
            aMixed = aMixed + isAngleObtuseAtVertex .* (triangleAreas ./ 2);
            
            % Else, add area(T) / 4:
            aMixed = aMixed + ~isAngleObtuseAtVertex .* (triangleAreas ./ 4);
            
            % If the triangle is obtuse:
%             if (isTriangleObtuse(vertexPositionP, vertexPositionQ, vertexPositionR))
%                
%                 % If the angle of triangle at vertex is obtuse:
%                 if (halfedgeFromP.getTrait('angle') > 90)
%                     
%                     % Add area(T)/2
%                     % TODO
%                     
%                 else
%                     
%                     % Add area(T)/4
%                     % TODO
%                     
%                 end
%                 
%             % If the triangle is non-obtuse:
%             else 
%                 
%                 % Use the Voronoi formula:
%                 aMixed = computeVoronoiArea(vertexPosition1, vertexPosition2, vertexPosition3, halfedgeFromQ.getTrait('cot_angle'), halfedgeFromR.getTrait('cot_angle'));
%                 
%             end
            
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
        
        function meanCurvature = computeMeanCurvatureNormal(halfedges, aMixedPerVertex)
            
            % Get origin of each halfedge:
            origins = halfedges.from().getTrait('position');
            
            % Get destination of each halfedge
            destinations = halfedges.to().getTrait('position');
            
            % Get alpha and beta angles associated with each halfedge, as
            % described in [Meyer2002, Section 3.1]:
            cotangentOfAlpha = halfedges.prev().getTrait('cot_angle');
            cotangentOfBeta = halfedges.twin().prev().getTrait('cot_angle');
            
            % Calculate the contribution of each halfedge to the mean
            % curvature of its origin vertex:
            meanCurvaturePerHalfedge = (cotangentOfAlpha + cotangentOfBeta) .* (origins - destinations);
            
            % Each halfedge has an origin vertex:
            vertices = halfedges.from();
            
            % Sum the contributions of each halfedge towards the mean
            % curvature of its origin vertex, as described in [Meyer2002,
            % Section 3.5]:
            scalar = 1 ./ ( 2 .* aMixedPerVertex);
            meanCurvature = scalar .* MeshHelper.sumRowVectorsByIndex(vertices.index, meanCurvaturePerHalfedge);
            
        end
        
        function gaussianCurvature = computeGaussianCurvature(halfedges, aMixedPerVertex)
            
            % Get theta which is the angle of a neighbor face at the vertex
            % which each halfedge origins, as described in [Meyer2002,
            % Section 4.1]:
            thetaAngles = halfedges.getTrait('angle');
            
            % Each halfedge has an origin vertex:
            vertices = halfedges.from();
            
            % Sum all theta angles and store the result per vertex. This
            % way, each row (corresponding to a vertex) will have the sum
            % of all the theta angles associated with that vertex:
            sumOfThetaAnglesPerVertex = MeshHelper.sumRowVectorsByIndex(vertices.index, thetaAngles);
            
            % Calculate Gaussian Curvature as described in [Meyer2002,
            % Section 4.2]:
            gaussianCurvature = (2 .* pi - sumOfThetaAnglesPerVertex) ./ aMixedPerVertex;
            
        end
        
        function calculateDiscreteCurvatures(mesh)
            % Computes vertex traits 'mean_curv' for the discrete mean
            % curvature and 'gauss_curv' for the discrete Gaussian
            % curvature.
            
            % TODO_A2 Task 6
            % Compute the discrete Gauss curvature and discrete mean
            % curvature for every vertex. See [Meyer2002], Sections
            % 3.3-3.5, 4 for a detailed description. The first step
            % should be to compute A_mixed for every adjacent face of a
            % vertex. You can store it as a halfedge trait, because it
            % describes a property between a halfedge and its
            % previous halfedge (similar to the 'angle' trait). From
            % there, you can sum up the mixed area for each vertex, and
            % then compute the discrete curvatures.
            % Note that a 'cot_angle' halfedge trait has already been
            % added (see calculateHalfedgeTrait()).

            % Get all non-boundary halfedges:
            halfedges = MeshHelper.getAllNonBoundaryHalfedges(mesh);
            
            % Compute the A_mixed of every adjacent face of the
            % vertex associated with origin of the halfedge:
            aMixed = MeshHelper.computeAMixed(halfedges);
            halfedges.setTrait('a_mixed', aMixed);
            
            % Each halfedge has an origin vertex:
            vertices = halfedges.from();
            
            % Sum the mix area for each vertex:
            aMixedPerVertex = MeshHelper.sumRowVectorsByIndex(vertices.index, aMixed);
            
            % Calculate the mean curvature per vertex:
            meanCurvatureNormal = MeshHelper.computeMeanCurvatureNormal(halfedges, aMixedPerVertex);
            
            % The mean curvature value K_H is half the magnitude of the
            % mean curvature normal, as described in [Meyer2002, Section
            % 3.5]:
            meanCurvatureValues = 0.5 .* MeshHelper.computeMagnitude(meanCurvatureNormal);
            mesh.getAllVertices().setTrait('mean_curv', meanCurvatureValues);
            
            % Compute the gaussian curvature, as described in [Meyer2002,
            % Section 4]:
            gaussianCurvatureValues = MeshHelper.computeGaussianCurvature(halfedges, aMixedPerVertex);
            mesh.getAllVertices().setTrait('gauss_curv', gaussianCurvatureValues);
            
        end
        
        
        function calculateVertexNormals(mesh, weighting)
            % Computes vertex normals as a weighted mean of face normals.
            % The parameter 'weighting' can be one of the following:
            % 'area': The face normal weights equal the face surface areas.
            % 'angle': The face normal weights equal the opening angle
            %    of the face at the vertex.
            if nargin<2
                weighting='area';
            end
            switch weighting
                case 'area'
                    f = mesh.getAllFaces();
                    fn_weighted = bsxfun(@times, f.getTrait('normal'), f.getTrait('area'));

                    vi1 = f.halfedge().from().index;
                    vi2 = f.halfedge().to().index;
                    vi3 = f.halfedge().next().to().index;

                    vi = repmat(reshape([vi1 vi2 vi3]',[],1),3,1);
                    vj = kron([1;2;3],ones(mesh.num_faces*3,1));
                    vals = kron(fn_weighted(:),[1;1;1]);

                    vn = normr(full(sparse(vi,vj,vals)));
                    mesh.getAllVertices().setTrait('normal',vn);

                case 'angle'
                    he = mesh.getAllHalfedges();
                    he_inner = mesh.getHalfedge(he.index(he.face().index ~= 0));
                    vi = he_inner.from().index;
                    fn_weighted = bsxfun(@times, he_inner.face().getTrait('normal'), he_inner.getTrait('angle'));
                    vn = zeros(mesh.num_vertices,3);
                    for i=1:3
                        vn(:,i) = full(sparse(vi, ones(size(vi)), fn_weighted(:,i)));
                    end
                    mesh.getAllVertices().setTrait('normal',normr(vn));
            end
        end
    end
end