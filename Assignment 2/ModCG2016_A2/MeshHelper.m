classdef MeshHelper < handle
    methods(Static)
        
        function b = hasBoundary(mesh)
            % Returns 1 if the mesh has a boundary, and 0 otherwise.
            b = any(mesh.getAllHalfedges().face().index==0);
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
            
            vol = 1;
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

            mesh.getAllVertices().setTrait('gauss_curv', 1);
            mesh.getAllVertices().setTrait('mean_curv', 1);
            
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