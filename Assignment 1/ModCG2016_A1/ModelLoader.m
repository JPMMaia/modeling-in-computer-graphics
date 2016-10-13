classdef ModelLoader < handle
    properties
    end
    
    methods(Static)
        function mesh = loadOBJ(filename)
            [V,F] = readOBJ_mex(filename);
            % -Z forward, Y up
            V = V*[0 -1 0; 0 0 1; -1 0 0];
            % uniformly resize to fit into unit cube
            r = max(range(V,1));
            V = V/r;
            % put in the center of the ground plane
            V = bsxfun(@minus, V, [mean(V(:,1:2),1) min(V(:,3),[],1)]);
            mesh = TriangleMesh(V,F);
        end
    end
end