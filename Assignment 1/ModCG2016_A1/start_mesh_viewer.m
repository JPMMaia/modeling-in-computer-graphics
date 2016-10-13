load_paths();

mesh_viewer = MeshViewerUI();

baseModelPath = 'data/';

models = [
  {'nefertiti.obj'}
  {'armadillo_6k.obj'} 
  {'bunny.obj'}
  {'cube_1.obj'}
  {'cube_2.obj'}
  {'cube_3.obj'}
  {'cube3.obj'}
  {'fan1.obj'}
  {'fan2.obj'}
  {'fertilty.obj'}
  {'patch1.obj'}
  {'venus.obj'}
];

mesh_viewer.loadModel(strcat(baseModelPath, char(models(4))));
while waitforbuttonpress ~= 1
end;

for i=1:size(models)
    
    mesh_viewer.loadModel(strcat(baseModelPath, char(models(i))));
    
    while waitforbuttonpress ~= 1
    end;
    
end

clear mesh_viewer;