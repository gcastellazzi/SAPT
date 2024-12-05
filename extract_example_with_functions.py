#######################################################
#### Simulia Abaqus Python Tool (SAPT) - Functions ####
#######################################################

# Clear all user-defined variables from the global namespace
# This prevents any pre-existing variables from interfering with the script
for name in dir():
    if not name.startswith("_"):
        del globals()[name]

# Import necessary modules
import sys
import os
from abaqus import *  # Core Abaqus module
from abaqusConstants import *  # Abaqus constants
from odbAccess import *  # For accessing ODB files
import numpy as np  # For numerical operations

#### INITIAL SETUP ####
# Define the ODB file name and paths
odb_name = 'example.odb'  # Name of the .odb file
odb_folder_path = 'X:/GitHub/SAPT/'  # Folder containing the .odb file
sapt_functions_folder = 'X:/GitHub/SAPT/'  # Folder containing SAPT_functions.py
results_folder = os.path.join(odb_folder_path, 'results/')  # Folder for saving output files
additional_suffix = 'test_example_'  # Prefix for output file names

# Construct the full path to the ODB file
odb_path = os.path.join(odb_folder_path, odb_name)

# Ensure the results folder exists
if not os.path.exists(results_folder):
    os.makedirs(results_folder)  # Create the folder if it does not exist
    print("Directory created: %s" % results_folder)
else:
    print("Directory already exists: %s" % results_folder)

# Add the folder containing SAPT_functions.py to the Python path
if sapt_functions_folder not in sys.path:
    sys.path.append(sapt_functions_folder)

# Import functions from SAPT_functions.py
try:
    from SAPT_functions import *
    print("SAPT_functions successfully imported.")
except ImportError as e:
    print("Error importing SAPT_functions: %s" % str(e))
    sys.exit("Ensure SAPT_functions.py is in the specified folder.")

#### START THE ANALYSIS OF THE ODB ####
try:
    # Open the output database
    odb = openOdb(path=odb_path)
    print("Successfully opened ODB file: %s" % odb_path)
except Exception as e:
    sys.exit("Failed to open ODB file %s. Error: %s" % (odb_path, str(e)))

# Scan ODB by instance to get element types
try:
    instance_element_types = scan_odb_by_instance(odb)
    print("Scanned ODB by instance to identify element types.")
except Exception as e:
    sys.exit("Error scanning ODB by instance. Error: %s" % str(e))

# Get the list of steps
try:
    list_of_steps = get_step_names(odb)
    print("Steps in ODB: %s" % str(list_of_steps))
except Exception as e:
    sys.exit("Error retrieving steps from ODB. Error: %s" % str(e))

#### SETUP RESULTS VIEW ####
# Optional: Setup for results visualization if in GUI mode
if hasattr(session, 'getInputs'):  # True only in Abaqus/CAE GUI mode
    print("Configuring viewport for results in CAE mode.")
    try:
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
    except KeyError:
        print("Viewport: 1 not found. Ensure the CAE GUI is set up correctly.")
else:
    print("Running in non-GUI mode, skipping viewport setup.")

print("Setup complete. Ready to perform operations on the ODB.")

#### COORDINATE SYSTEM UPDATE
change_coordinates = True
# Define the nodes that form the current coordinate system
# Origin_Node_id: The node representing the original origin of the coordinate system
# Node_1_id: A node used to define the x-axis of the original coordinate system
# Node_2_id: A node used to define the plane (and implicitly the y-axis and z-axis) of the original system
Origin_Node_id = 11
Node_1_id = 1221
Node_2_id = 1211

# Extract the coordinates of the three defining nodes from the ODB
# This defines the current (original) coordinate system in terms of the model's global coordinates.
original_origin = tuple(extract_node_coordinates(odb, instance_element_types.keys()[0], Origin_Node_id))
original_point1 = tuple(extract_node_coordinates(odb, instance_element_types.keys()[0], Node_1_id))
original_point2 = tuple(extract_node_coordinates(odb, instance_element_types.keys()[0], Node_2_id))

# Define the new coordinate system
# The new_origin is where the original_origin will be translated to.
# new_point1 and new_point2 define the orientation of the new x-axis and the plane for the new coordinate system.
# This system will rotate the original points to align with these definitions.
new_origin = (0.0, 0.0, 0.0)  # New origin location (e.g., global origin)
new_point1 = (1.0, 0.0, 0.0)  # Defines the new x-axis direction
new_point2 = (1.0, 0.0, 1.0)  # Defines the plane for the new coordinate system (new z-axis is perpendicular to this plane)

# Summary:
# 1. Translation: Shifts the original_origin to the new_origin.
# 2. Rotation: Aligns the original coordinate system axes with the new system, ensuring the original points align as follows:
#     - original_origin -> new_origin
#     - original_point1 -> new_point1
#     - original_point2 -> new_point2
# The transformation ensures no scaling or deformation of the mesh-only rigid body translation and rotation.

# Conditional check to ensure the script only tries to access `session` in CAE mode
if hasattr(session, 'getInputs'):  # True only in GUI mode
    # This will only run in CAE mode with GUI
    print("Viewport set in CAE mode.")
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)
else:
    print("Running in non-GUI mode, skipping viewport setup.")


for instance_name in instance_element_types:
    # Identify element types
    element_types = identify_element_types(odb, instance_name)
    count_of_hex8 = element_types.get('C3D8', 0)
    print("############################################################################################")
    print("Instance: " + instance_name + "Count of Hex8 elements:", count_of_hex8)
    count_of_tet4 = element_types.get('C3D4', 0)
    print("Instance: " + instance_name + "Count of Tet4 elements:", count_of_tet4)
    print("############################################################################################")
    print("#")
    for step_name in list_of_steps:
        # COORDINATES
        print("############################################")
        print("Instance: " + instance_name + " #### Step:", step_name)
        print("############################################")
        output_filename = results_folder +  additional_suffix + instance_name + '_' + step_name + '_' + 'node_coordinates.txt'
        output_filename_t = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'transformed_node_coordinates.txt'
        nodes_in_model, nodes_in_model_dic = write_node_ids_and_coordinates(odb, step_name, instance_name, output_filename)
        # Transform coordinates
        if change_coordinates:
            nodes_in_model = transform_coordinates_dict(original_origin, original_point1, original_point2, new_origin, new_point1, new_point2, nodes_in_model_dic, output_filename_t)
        #
        if instance_element_types[instance_name].get('C3D8',0)>0:            
            # BOUNDARY NODES AND CONNECTIVITY (ORIGINAL)
            output_filename = results_folder +  additional_suffix + instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'boundary_nodes_connectivity.txt'
            boundary_nodes, q4_faces_in_model = write_boundary_nodes_and_connectivity(odb, instance_name, output_filename)
            #
            # Extract Hex8 elements
            hex_elements = extract_hex_elements(odb, instance_name)
            #
            # TRIANGULARE FACES
            output_filename = results_folder +  additional_suffix + instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'triangular_faces.txt'
            faces_triangular_in_model = write_triangular_faces(odb, instance_name, boundary_nodes, output_filename)
            #
            # DISPLACEMENTS
            output_filename = results_folder +  additional_suffix + instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'boundary_displacements.txt'
            displacements_in_model = write_boundary_displacements(odb, step_name, instance_name, boundary_nodes, output_filename)
            #
            # VTK file
            print('VTK tri start')
            output_filename = results_folder +  additional_suffix + instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'triangular_mesh_U.vtk'
            write_vtk_file_warp(output_filename, nodes_in_model, faces_triangular_in_model, displacements_in_model)
            print('VTK tri done')
            #
            # VTK file
            print('VTK Hexa start')
            output_filename = results_folder +  additional_suffix + instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'Hexa_mesh_U.vtk'
            write_vtk_hex8_warp(output_filename, nodes_in_model, hex_elements, displacements_in_model)
            print('VTK Hexa end')
            #
            # PLY FILE
            output_ply_file = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + "nodes_with_displacement_colors.ply"
            write_ply_with_displacement_colors(output_ply_file, nodes_in_model, displacements_in_model)
            #
        if instance_element_types[instance_name].get('C3D4',0)>0:
            # Nodes and Coordinates
            nodes_output = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_nodes.txt'
            nodes_in_model, nodes_in_model_dic = write_node_ids_and_coordinates_2(odb, instance_name, nodes_output)
            #
            # Tet4 Connectivity
            connectivity_output = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_connectivity.txt'
            tet4_elements = extract_tet4_elements(odb, instance_name, connectivity_output)
            #
            # Boundary Nodes
            boundary_nodes = identify_boundary_nodes(tet4_elements)
            #
            # Boundary Displacements
            displacements_output = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_boundary_displacements.txt'
            displacements_in_model = extract_displacements(odb, step_name, instance_name, boundary_nodes, displacements_output)
            #
            # Write VTK file
            vtk_output = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_mesh_U.vtk'
            write_vtk_file_warp(vtk_output, nodes_in_model, tet4_elements, displacements_in_model)
            # PLY FILE
            output_ply_file = results_folder +  additional_suffix +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + "nodes_with_displacement_colors.ply"
            write_ply_with_displacement_colors(output_ply_file, nodes_in_model, displacements_in_model)
            #
        if instance_element_types[instance_name].get('C3D4',0)>0 and instance_element_types[instance_name].get('C3D8',0)>0:
            print("Element type not found") 

close_odb = True #False
if close_odb:
    odb.close()
