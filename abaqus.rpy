# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-2 replay file
# Internal Version: 2014_08_22-16.00.46 134497
# Run by gcastellazzi on Thu Dec 05 10:00:40 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.17187, 1.14583), width=172.5, 
    height=113.667)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('extract_example_with_functions.py', __main__.__dict__)
#: Directory already exists: X:/GitHub/SAPT/results/
#: SAPT_functions successfully imported.
#: Warning: An output database lock file X:\GitHub\SAPT\example.lck has been detected. This may indicate that the output database is opened for writing by another application.
#: X:/GitHub/SAPT/example.odb will be opened read-only. 
#: Model: X:/GitHub/SAPT/example.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       10
#: Number of Node Sets:          10
#: Number of Steps:              2
#: Successfully opened ODB file: X:/GitHub/SAPT/example.odb
#: #### Entering the function: scan_odb_by_instance ####
#: Element types by instance in ODB:
#: Instance 'PART-1-1':
#:   C3D8: 1000
#: Scanned ODB by instance to identify element types.
#: #### Entering the function: get_step_names ####
#: Steps in ODB: ['Step-1', 'Step-2']
#: Running in non-GUI mode, skipping viewport setup.
#: Setup complete. Ready to perform operations on the ODB.
#: #### Entering the function: extract_node_coordinates ####
#: #### Entering the function: extract_node_coordinates ####
#: #### Entering the function: extract_node_coordinates ####
#: Running in non-GUI mode, skipping viewport setup.
#: #### Entering the function: identify_element_types ####
#: Element types and counts for instance 'PART-1-1':
#:   C3D8: 1000
#: ############################################################################################
#: ('Instance: PART-1-1Count of Hex8 elements:', 1000)
#: ('Instance: PART-1-1Count of Tet4 elements:', 0)
#: ############################################################################################
#: #
#: ############################################
#: ('Instance: PART-1-1 #### Step:', 'Step-1')
#: ############################################
#: #### Entering the function: write_node_ids_and_coordinates ####
#: Node IDs and coordinates have been written to : X:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_node_coordinates.txt
#: #### Entering the function: transform_coordinates_dict ####
#: Original Basis:
#: [[-1.  0.  0.]
#:  [ 0.  1.  0.]
#:  [ 0.  0. -1.]]
#: New Basis:
#: [[ 1.  0.  0.]
#:  [ 0.  0.  1.]
#:  [ 0. -1.  0.]]
#: Rotation Matrix:
#: [[-1.  0.  0.]
#:  [ 0.  0.  1.]
#:  [ 0.  1.  0.]]
#: Translation Vector:
#: [ 50. -50.   0.]
#: #### End transform_coordinates_dict ############################
#: #### Entering the function: write_boundary_nodes_and_connectivity ####
#: Boundary nodes and connectivity have been written to : X:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_C3D8_boundary_nodes_connectivity.txt
#: #### Entering the function: extract_hex_elements ####
#: Extracted 1000 Hex8 elements from instancePART-1-1
#: #### Entering the function: write_triangular_faces ####
#: Triangular face connectivity has been written toX:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_C3D8_triangular_faces.txt
#: #### Entering the function: write_boundary_displacements ####
#: Boundary displacements written to X:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_C3D8_boundary_displacements.txt
#: VTK tri start
#: #### Entering the function: write_vtk_file_warp ####
#: ##   Writing VTK file...
#: ##   VTK file written to:X:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_C3D8_triangular_mesh_U.vtk
#: VTK tri done
#: VTK Hexa start
#: #### Entering the function: write_vtk_hex8_warp ####
#: VTK file for Hex8 mesh with vector displacements written to:X:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_C3D8_Hexa_mesh_U.vtk
#: VTK Hexa end
#: #### Entering the function: write_ply_with_displacement_colors ####
#: PLY file written to X:/GitHub/SAPT/results/test_example_PART-1-1_Step-1_C3D8_nodes_with_displacement_colors.ply
#: ############################################
#: ('Instance: PART-1-1 #### Step:', 'Step-2')
#: ############################################
#: #### Entering the function: write_node_ids_and_coordinates ####
#: Node IDs and coordinates have been written to : X:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_node_coordinates.txt
#: #### Entering the function: transform_coordinates_dict ####
#: Original Basis:
#: [[-1.  0.  0.]
#:  [ 0.  1.  0.]
#:  [ 0.  0. -1.]]
#: New Basis:
#: [[ 1.  0.  0.]
#:  [ 0.  0.  1.]
#:  [ 0. -1.  0.]]
#: Rotation Matrix:
#: [[-1.  0.  0.]
#:  [ 0.  0.  1.]
#:  [ 0.  1.  0.]]
#: Translation Vector:
#: [ 50. -50.   0.]
#: #### End transform_coordinates_dict ############################
#: #### Entering the function: write_boundary_nodes_and_connectivity ####
#: Boundary nodes and connectivity have been written to : X:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_C3D8_boundary_nodes_connectivity.txt
#: #### Entering the function: extract_hex_elements ####
#: Extracted 1000 Hex8 elements from instancePART-1-1
#: #### Entering the function: write_triangular_faces ####
#: Triangular face connectivity has been written toX:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_C3D8_triangular_faces.txt
#: #### Entering the function: write_boundary_displacements ####
#: Boundary displacements written to X:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_C3D8_boundary_displacements.txt
#: VTK tri start
#: #### Entering the function: write_vtk_file_warp ####
#: ##   Writing VTK file...
#: ##   VTK file written to:X:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_C3D8_triangular_mesh_U.vtk
#: VTK tri done
#: VTK Hexa start
#: #### Entering the function: write_vtk_hex8_warp ####
#: VTK file for Hex8 mesh with vector displacements written to:X:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_C3D8_Hexa_mesh_U.vtk
#: VTK Hexa end
#: #### Entering the function: write_ply_with_displacement_colors ####
#: PLY file written to X:/GitHub/SAPT/results/test_example_PART-1-1_Step-2_C3D8_nodes_with_displacement_colors.ply
print 'RT script done'
#: RT script done
