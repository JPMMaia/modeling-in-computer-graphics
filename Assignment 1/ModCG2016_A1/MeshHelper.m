classdef MeshHelper < handle
    methods(Static)
        
        function [p_min, p_max] = getBoundingBox(mesh)
            % Returns the points with minimal and maximal coordinates of the
            % smallest axis-aligned bounding box of the mesh.

            % TODO_A1 Task 1
            % 
            % Find the axis-aligned bounding box of the mesh and return its
            % minimal and maximal corner vertices. Use the vertex trait
            % 'position' to find them.

            p_min = [0 0 0];
            p_max = [0.2 0.2 0.2];
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

            V_start = zeros(0,3);
            V_end = zeros(0,3);
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
                case 'angle'
                    % TODO_A1 Task 5b
                    %
                    % Fill in the 'angle' branch of this function.
                    % Calculate the vertex normals as weighted averages of
                    % the adjacent face normals, where the weight is given
                    % by the angle that the face confines at the vertex.
                    % Use the 'angle' halfedge trait computed in Task 5a for this.
            end
        end
    end
end