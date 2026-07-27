Avec `if(std::abs(covered - faceArea) > tol_ * 0.5 * perimeter)`
```
[0]      =============== Mesh Info ===============

             Cell-numbering              = Global
             Number of 'Node'            = 81
             Number of 'Poly2D'          = 64
             Number of 'Poly3D'          = 16
             -- No parents --
             Cell attributes of 'Node'   = ['GlobalNumber']
             Cell attributes of 'Poly2D' = ['CADGroupID', 'GlobalNumber']
             Datasets:
               'Coordinates' spanning across ( 'Node' ):
                 'CoordinateX', range = [0 <---> 0.001]
                 'CoordinateY', range = [0 <---> 0.001]
                 'CoordinateZ', range = [0.001 <---> 0.001]
             
[0]      FSMeshSelection::ApplyInto() 0.00061 [s] (wall clock time)

[0]      =============== Mesh Info ===============

             Cell-numbering              = Global
             Number of 'Node'            = 79
             Number of 'Poly2D'          = 62
             Number of 'Poly3D'          = 23
             -- No parents --
             Cell attributes of 'Node'   = ['GlobalNumber']
             Cell attributes of 'Poly2D' = ['CADGroupID', 'GlobalNumber']
             Datasets:
               'Coordinates' spanning across ( 'Node' ):
                 'CoordinateX', range = [0 <---> 0.001]
                 'CoordinateY', range = [0 <---> 0.001]
                 'CoordinateZ', range = [0.001 <---> 0.001]
```

Avec `if(std::abs(covered - faceArea) > tol_)`

```
[0]      =============== Mesh Info ===============

             Cell-numbering              = Global
             Number of 'Node'            = 81
             Number of 'Poly2D'          = 64
             Number of 'Poly3D'          = 16
             -- No parents --
             Cell attributes of 'Node'   = ['GlobalNumber']
             Cell attributes of 'Poly2D' = ['CADGroupID', 'GlobalNumber']
             Datasets:
               'Coordinates' spanning across ( 'Node' ):
                 'CoordinateX', range = [0 <---> 0.001]
                 'CoordinateY', range = [0 <---> 0.001]
                 'CoordinateZ', range = [0.001 <---> 0.001]
             
[0]      FSMeshSelection::ApplyInto() 0.00061 [s] (wall clock time)

[0]      =============== Mesh Info ===============

             Cell-numbering              = Global
             Number of 'Node'            = 79
             Number of 'Poly2D'          = 62
             Number of 'Poly3D'          = 23
             -- No parents --
             Cell attributes of 'Node'   = ['GlobalNumber']
             Cell attributes of 'Poly2D' = ['CADGroupID', 'GlobalNumber']
             Datasets:
               'Coordinates' spanning across ( 'Node' ):
                 'CoordinateX', range = [0 <---> 0.001]
                 'CoordinateY', range = [0 <---> 0.001]
                 'CoordinateZ', range = [0.001 <---> 0.001]
```