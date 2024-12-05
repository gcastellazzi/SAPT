# Abaqus functions
import sys
import os
from abaqus import *
from abaqusConstants import *
from odbAccess import *
import numpy as np
import inspect

# Utility function to print messages based on verbose mode
def verbose_print(message, verbose):
    if verbose == 'On':
        print(message)

# Function to write node IDs and coordinates to a file
def write_node_ids_and_coordinates(odb, step_name, instance_name, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    # Access the specified instance
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance: " + instance_name + " not found in ODB.", verbose)
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    nodes_in_model = []
    # Prepare to write node IDs and coordinates
    with open(output_filename, 'w') as f:
        f.write("Node ID, X, Y, Z\n")
        
        # Loop through nodes in the instance and write their coordinates
        for node in instance.nodes:
            node_id = node.label
            x, y, z = node.coordinates
            f.write("{}, {}, {}, {}\n".format(node_id, x, y, z))
            # Append node information to nodes_in_model list
            nodes_in_model.append((node_id, x, y, z))
    
    nodes_in_model_dic = {node_id: (x, y, z) for node_id, x, y, z in nodes_in_model}
    verbose_print("Node IDs and coordinates have been written to : " + str(output_filename), verbose)
    
    return nodes_in_model, nodes_in_model_dic

# Function to write boundary nodes and their connectivity to a file
def write_boundary_nodes_and_connectivity(odb, instance_name, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    # Check if instance exists
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance: " + instance_name + " not found in ODB.", verbose)
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    q4_faces_in_model = []
    
    # Initialize a dictionary to count occurrences of each node
    node_count = {}
    
    # Loop through elements in the instance to count node occurrences
    for element in instance.elements:
        for node in element.connectivity:
            if node in node_count:
                node_count[node] += 1
            else:
                node_count[node] = 1
    
    # Identify boundary nodes (connected to fewer than a certain number of elements)
    boundary_nodes = {node for node, count in node_count.items() if count <= 4}
    
    # Define 8-node hexahedral element face definitions (each face is a quadrilateral)
    face_definitions = [
        (0, 1, 2, 3),  # Bottom face
        (4, 5, 6, 7),  # Top face
        (0, 1, 5, 4),  # Side face 1
        (1, 2, 6, 5),  # Side face 2
        (2, 3, 7, 6),  # Side face 3
        (3, 0, 4, 7)   # Side face 4
    ]
    
    # Prepare to write boundary connectivity
    with open(output_filename, 'w') as f:
        f.write("Boundary Node ID, Connectivity\n")
        
        # Loop through elements and extract boundary faces
        for element in instance.elements:
            element_connectivity = list(element.connectivity)
            
            # Skip elements with fewer than 8 nodes
            if len(element_connectivity) < 8:
                continue
            
            # Check each face of the element to see if it lies on the boundary
            for face_nodes_indices in face_definitions:
                face_nodes = [element_connectivity[i] for i in face_nodes_indices]
                
                # If all face nodes are in the boundary node set, save this face
                if all(node in boundary_nodes for node in face_nodes):
                    # Write face connectivity as a quadrilateral
                    f.write("Face: {}, {}, {}, {}\n".format(face_nodes[0], face_nodes[1], face_nodes[2], face_nodes[3]))
                    # Append node information to nodes_in_model list
                    q4_faces_in_model.append([face_nodes[0], face_nodes[1], face_nodes[2], face_nodes[3]])
    
    #nodes_in_model_dic = {node_id: (f1, f2, f3, f4) for f1, f2, f3, f4 in q4_faces_in_model}
    verbose_print("Boundary nodes and connectivity have been written to : " + output_filename, verbose)
    return boundary_nodes, q4_faces_in_model

# Function to write external face connectivity as triangles
def write_triangular_faces(odb, instance_name, boundary_nodes, filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    # Check if instance exists
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance: " + instance_name + " not found in ODB.", verbose)
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    results = []
    faces_triangular_in_model = []
    
    # Define face connectivity for an 8-node hexahedral element (C3D8)
    face_definitions = [
        (0, 1, 2, 3),  # Bottom face
        (4, 5, 6, 7),  # Top face
        (0, 1, 5, 4),  # Side face 1
        (1, 2, 6, 5),  # Side face 2
        (2, 3, 7, 6),  # Side face 3
        (3, 0, 4, 7)   # Side face 4
    ]
    
    # Convert boundary nodes list to a set for faster lookup
    boundary_nodes_set = set(boundary_nodes)
    
    # Loop through elements and check for quadrilateral boundary faces
    for element in instance.elements:
        element_connectivity = list(element.connectivity)
        
        # Skip elements that do not have exactly 8 nodes
        if len(element_connectivity) != 8:
            continue
        
        # Check each face of the element
        for face_indices in face_definitions:
            face_nodes = [element_connectivity[i] for i in face_indices]
            
            # Check if all nodes in this face are in the boundary node set
            if all(node in boundary_nodes_set for node in face_nodes):
                # Divide the quadrilateral face into two triangles
                triangle1 = [face_nodes[0], face_nodes[1], face_nodes[2]]
                triangle2 = [face_nodes[0], face_nodes[2], face_nodes[3]]
                results.extend([triangle1, triangle2])
    
    # Write triangular faces to a file
    with open(filename, 'w') as f:
        f.write("Triangle Connectivity\n")
        for triangle in results:
            f.write("{}, {}, {}\n".format(triangle[0], triangle[1], triangle[2]))
            faces_triangular_in_model.append((triangle[0], triangle[1], triangle[2]))
    
    verbose_print("Triangular face connectivity has been written to" + filename, verbose)
    return faces_triangular_in_model

def write_boundary_displacements(odb, step_name, instance_name, boundary_nodes, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Extracts displacements at nodal positions for boundary nodes and writes them to a file.
    """
    # Check if the instance exists
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance " + str(instance_name) + " not found in ODB.", verbose)
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    displacements_in_model = []
    
    # Check if the step exists
    if step_name not in odb.steps:
        verbose_print("Step " + str(step_name) + " not found in ODB.", verbose)
        odb.close()
        return
    
    step = odb.steps[step_name]
    frame = step.frames[-1]  # Use the last frame
    
    # Retrieve displacement field output
    displacement_field = frame.fieldOutputs['U']
    
    # Filter displacements for boundary nodes
    with open(output_filename, 'w') as f:
        f.write("Node ID, Ux, Uy, Uz\n")
        for value in displacement_field.values:
            node_label = value.nodeLabel
            if node_label in boundary_nodes:  # Filter by boundary node IDs
                ux, uy, uz = value.data
                f.write("{}, {}, {}, {}\n".format(node_label, ux, uy, uz))
                displacements_in_model.append((node_label, ux, uy, uz))
    
    verbose_print("Boundary displacements written to " + str(output_filename), verbose)
    return displacements_in_model

def write_boundary_displacements_old2(odb, step_name, instance_name, boundary_nodes, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    # Check if the instance exists
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance " + str(instance_name) + " not found in ODB.", verbose)
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    displacements_in_model = []
    NODAL = 'NODAL'
    # Check if the step exists
    if step_name not in odb.steps:
        verbose_print("Step" +str(step_name) + " not found in ODB.", verbose)
        odb.close()
        return
    
    step = odb.steps[step_name]
    frame = step.frames[-1]  # Use the last frame
    
    # Retrieve displacement field output for the specified frame
    displacement_field = frame.fieldOutputs['U']
    
    # Filter displacements at nodal positions for boundary nodes
    nodal_displacement = displacement_field.getSubset(position=NODAL)
    
    # Write boundary displacements to file
    with open(output_filename, 'w') as f:
        f.write("Node ID, Ux, Uy, Uz\n")
        for value in nodal_displacement.values:
            node_label = value.nodeLabel
            if node_label in boundary_nodes:
                ux, uy, uz = value.data
                f.write("{}, {}, {}, {}\n".format(node_label, ux, uy, uz))
                displacements_in_model.append((node_label, ux, uy, uz))
    
    verbose_print("Boundary displacements written to " + str(output_filename), verbose)
    return displacements_in_model

# Function to extract Ux, Uy, Uz displacements at boundary nodes
def write_boundary_displacements_old(odb, step_name, instance_name, boundary_nodes, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    # Check if instance exists
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance: " + instance_name + " not found in ODB.", verbose)
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    displacements_in_model = []
    
    # Check if the step exists
    if step_name not in odb.steps:
        verbose_print("Step: " + step_name + " not found in ODB.", verbose)
        odb.close()
        return
    
    step = odb.steps[step_name]
    frame = step.frames[-1]  # Last frame, or change to a specific frame if needed
    
    # Convert boundary nodes to a set for faster lookup
    boundary_nodes_set = set(boundary_nodes)
    
    # Retrieve displacement field output for the specified frame
    displacement_field = frame.fieldOutputs['U']
    
    # Filter displacements at nodal positions for boundary nodes
    nodal_displacement = displacement_field.getSubset(position=NODAL)
    
    # Prepare to write boundary displacements
    with open(output_filename, 'w') as f:
        f.write("Node ID, Ux, Uy, Uz\n")
        
        # Loop through the nodal displacement values and filter for boundary nodes
        for value in nodal_displacement.values:
            node_label = value.nodeLabel
            if node_label in boundary_nodes_set:
                ux, uy, uz = value.data  # Ux, Uy, Uz components
                f.write("{}, {}, {}, {}\n".format(node_label, ux, uy, uz))
                displacements_in_model.append((node_label, ux, uy, uz))
    
    verbose_print("Boundary displacements (Ux, Uy, Uz) have been written to: " + output_filename, verbose)
    return displacements_in_model

# Function to extract Hex8 elements from an ODB file
def extract_hex_elements(odb, instance_name, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Extracts Hex8 (8-node hexahedral) elements from an ODB instance.
    
    :param odb: The ODB object.
    :param instance_name: Name of the instance from which to extract elements.
    :return: List of Hex8 elements with node connectivity [[node1, node2, ..., node8], ...].
    """
    # Check if the instance exists
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance: " + instance_name + " not found in ODB.", verbose)
        odb.close()
        return []
    
    instance = odb.rootAssembly.instances[instance_name]
    hex_elements = []
    
    # Loop through elements in the instance
    for element in instance.elements:
        # Only consider elements with 8 nodes (Hex8 elements)
        if len(element.connectivity) == 8:
            # Append the connectivity as a list of node IDs
            hex_elements.append(list(element.connectivity))
    
    verbose_print("Extracted " + str(len(hex_elements)) + " Hex8 elements from instance" + str(instance_name), verbose)
    return hex_elements


def write_vtk_file(filename, nodes, triangles, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for triangular mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        for _, x, y, z in nodes:
            f.write("{} {} {}\n".format(x, y, z))
        
        # Write triangle connectivity using remapped indices
        f.write("CELLS {} {}\n".format(len(triangles), len(triangles) * 4))
        for triangle in triangles:
            f.write("3 {} {} {}\n".format(
                node_index_map[triangle[0]], node_index_map[triangle[1]], node_index_map[triangle[2]]
            ))
        
        # Cell types
        f.write("CELL_TYPES {}\n".format(len(triangles)))
        f.write("5\n" * len(triangles))  # Triangle type in VTK
        
        # Write displacements as scalar fields
        f.write("POINT_DATA {}\n".format(len(nodes)))
        
        # Ux displacements
        f.write("SCALARS Ux float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            ux = displacements.get(node_id, (0.0, 0.0, 0.0))[0]
            f.write("{}\n".format(ux))
        
        # Uy displacements
        f.write("SCALARS Uy float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            uy = displacements.get(node_id, (0.0, 0.0, 0.0))[1]
            f.write("{}\n".format(uy))
        
        # Uz displacements
        f.write("SCALARS Uz float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            uz = displacements.get(node_id, (0.0, 0.0, 0.0))[2]
            f.write("{}\n".format(uz))

def write_vtk_file_warp(filename, nodes, triangles, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    # Convert displacements_list to a dictionary for easier access
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Handle both list and dictionary structures for nodes
    if isinstance(nodes, dict):
        node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    elif isinstance(nodes, list):
        node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    else:
        raise TypeError("Unsupported nodes data structure. Expected dict or list.")
    
    verbose_print('##   Writing VTK file...', verbose)
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for triangular mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        if isinstance(nodes, dict):
            for node_id, (x, y, z) in nodes.items():
                f.write("{} {} {}\n".format(x, y, z))
        elif isinstance(nodes, list):
            for node_id, x, y, z in nodes:
                f.write("{} {} {}\n".format(x, y, z))
        
        # Write triangle connectivity using remapped indices
        f.write("CELLS {} {}\n".format(len(triangles), len(triangles) * 4))
        for triangle in triangles:
            f.write("3 {} {} {}\n".format(
                node_index_map[triangle[0]], node_index_map[triangle[1]], node_index_map[triangle[2]]
            ))
            
        # Cell types
        f.write("CELL_TYPES {}\n".format(len(triangles)))
        f.write("5\n" * len(triangles))  # Triangle type in VTK
        
        # Write displacements as a vector field
        f.write("POINT_DATA {}\n".format(len(nodes)))
        f.write("VECTORS displacement float\n")
        if isinstance(nodes, dict):
            for node_id in nodes.keys():
                ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
                f.write("{} {} {}\n".format(ux, uy, uz))
        elif isinstance(nodes, list):
            for node_id, x, y, z in nodes:
                ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
                f.write("{} {} {}\n".format(ux, uy, uz))
    
    verbose_print("##   VTK file written to:" + filename, verbose)

def write_vtk_file_warp_old(filename, nodes, triangles, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for triangular mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        for _, x, y, z in nodes:
            f.write("{} {} {}\n".format(x, y, z))
        
        # Write triangle connectivity using remapped indices
        f.write("CELLS {} {}\n".format(len(triangles), len(triangles) * 4))
        for triangle in triangles:
            f.write("3 {} {} {}\n".format(
                node_index_map[triangle[0]], node_index_map[triangle[1]], node_index_map[triangle[2]]
            ))
        
        # Cell types
        f.write("CELL_TYPES {}\n".format(len(triangles)))
        f.write("5\n" * len(triangles))  # Triangle type in VTK
        
        # Write displacements as a vector field
        f.write("POINT_DATA {}\n".format(len(nodes)))
        f.write("VECTORS displacement float\n")
        for node_id, _, _, _ in nodes:
            ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to (0.0, 0.0, 0.0) if no data
            f.write("{} {} {}\n".format(ux, uy, uz))
    
    verbose_print("##   VTK file for triangular mesh with vector displacements written to:" + filename, verbose)

def write_vtk_hex8(filename, nodes, hex_elements, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for Hex8 mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        for _, x, y, z in nodes:
            f.write("{} {} {}\n".format(x, y, z))
        
        # Write Hex8 connectivity using remapped indices
        f.write("CELLS {} {}\n".format(len(hex_elements), len(hex_elements) * 9))
        for hex_element in hex_elements:
            f.write("8 {} {} {} {} {} {} {} {}\n".format(
                node_index_map[hex_element[0]], node_index_map[hex_element[1]], node_index_map[hex_element[2]], node_index_map[hex_element[3]],
                node_index_map[hex_element[4]], node_index_map[hex_element[5]], node_index_map[hex_element[6]], node_index_map[hex_element[7]]
            ))
        
        # Cell types
        f.write("CELL_TYPES {}\n".format(len(hex_elements)))
        f.write("12\n" * len(hex_elements))  # Hexahedron type in VTK
        
        # Write displacements as scalar fields
        f.write("POINT_DATA {}\n".format(len(nodes)))
        
        # Ux displacements
        f.write("SCALARS Ux float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            ux = displacements.get(node_id, (0.0, 0.0, 0.0))[0]
            f.write("{}\n".format(ux))
        
        # Uy displacements
        f.write("SCALARS Uy float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            uy = displacements.get(node_id, (0.0, 0.0, 0.0))[1]
            f.write("{}\n".format(uy))
        
        # Uz displacements
        f.write("SCALARS Uz float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            uz = displacements.get(node_id, (0.0, 0.0, 0.0))[2]
            f.write("{}\n".format(uz))

def write_vtk_hex8_warp(filename, nodes, hex_elements, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Writes a VTK file for Hex8 elements with displacement fields.
    
    :param filename: Output file name for the VTK file.
    :param nodes: Node coordinates (dictionary {node_id: (x, y, z)} or list [(node_id, x, y, z)]).
    :param hex_elements: Hexahedral elements connectivity.
    :param displacements_list: List of displacements [(node_id, ux, uy, uz), ...].
    """
    # Convert displacements_list to a dictionary for easier access
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Handle both dictionary and list structures for nodes
    if isinstance(nodes, dict):
        node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    elif isinstance(nodes, list):
        node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    else:
        raise TypeError("Unsupported nodes data structure. Expected dict or list.")
    
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for Hex8 mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        if isinstance(nodes, dict):
            for node_id, (x, y, z) in nodes.items():
                f.write("{} {} {}\n".format(x, y, z))
        elif isinstance(nodes, list):
            for node_id, x, y, z in nodes:
                f.write("{} {} {}\n".format(x, y, z))
        
        # Write Hex8 connectivity using remapped indices
        f.write("CELLS {} {}\n".format(len(hex_elements), len(hex_elements) * 9))
        for hex_element in hex_elements:
            f.write("8 {} {} {} {} {} {} {} {}\n".format(
                node_index_map[hex_element[0]], node_index_map[hex_element[1]], 
                node_index_map[hex_element[2]], node_index_map[hex_element[3]],
                node_index_map[hex_element[4]], node_index_map[hex_element[5]], 
                node_index_map[hex_element[6]], node_index_map[hex_element[7]]
            ))
            
        # Cell types
        f.write("CELL_TYPES {}\n".format(len(hex_elements)))
        f.write("12\n" * len(hex_elements))  # Hexahedron type in VTK
        
        # Write displacements as a vector field
        f.write("POINT_DATA {}\n".format(len(nodes)))
        f.write("VECTORS displacement float\n")
        if isinstance(nodes, dict):
            for node_id in nodes.keys():
                ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
                f.write("{} {} {}\n".format(ux, uy, uz))
        elif isinstance(nodes, list):
            for node_id, x, y, z in nodes:
                ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
                f.write("{} {} {}\n".format(ux, uy, uz))
        
    verbose_print("VTK file for Hex8 mesh with vector displacements written to:" + filename, verbose)

def write_vtk_hex8_warp_old(filename, nodes, hex_elements, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for Hex8 mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        for _, x, y, z in nodes:
            f.write("{} {} {}\n".format(x, y, z))
        
        # Write Hex8 connectivity using remapped indices
        f.write("CELLS {} {}\n".format(len(hex_elements), len(hex_elements) * 9))
        for hex_element in hex_elements:
            f.write("8 {} {} {} {} {} {} {} {}\n".format(
                node_index_map[hex_element[0]], node_index_map[hex_element[1]], node_index_map[hex_element[2]], node_index_map[hex_element[3]],
                node_index_map[hex_element[4]], node_index_map[hex_element[5]], node_index_map[hex_element[6]], node_index_map[hex_element[7]]
            ))
        
        # Cell types
        f.write("CELL_TYPES {}\n".format(len(hex_elements)))
        f.write("12\n" * len(hex_elements))  # Hexahedron type in VTK
        
        # Write displacements as a vector field
        f.write("POINT_DATA {}\n".format(len(nodes)))
        f.write("VECTORS displacement float\n")
        for node_id, _, _, _ in nodes:
            ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to (0.0, 0.0, 0.0) if no data
            f.write("{} {} {}\n".format(ux, uy, uz))
    
    verbose_print("VTK file for Hex8 mesh with vector displacements written to:" + filename, verbose)

# Function to write node IDs and coordinates to a file and return the node list and dictionary
def write_node_ids_and_coordinates_2(odb, instance_name, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance " + str(instance_name) + " not found in ODB.", verbose)
        odb.close()
        return [], {}
    
    instance = odb.rootAssembly.instances[instance_name]
    nodes_in_model = []
    
    # Write node coordinates to file
    with open(output_filename, 'w') as f:
        f.write("Node ID, X, Y, Z\n")
        for node in instance.nodes:
            node_id = node.label
            x, y, z = node.coordinates
            f.write("{}, {}, {}, {}\n".format(node_id, x, y, z))
            nodes_in_model.append((node_id, x, y, z))
    
    nodes_in_model_dic = {node_id: (x, y, z) for node_id, x, y, z in nodes_in_model}
    verbose_print("Node IDs and coordinates written to:" + output_filename, verbose)
    return nodes_in_model, nodes_in_model_dic

# Function to extract Tet4 elements and write connectivity to a file
def extract_tet4_elements(odb, instance_name, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance " + str(instance_name) + " not found in ODB.", verbose)
        odb.close()
        return []
    
    instance = odb.rootAssembly.instances[instance_name]
    tet4_elements = []
    
    # Write connectivity to file
    with open(output_filename, 'w') as f:
        f.write("Element ID, Node1, Node2, Node3, Node4\n")
        for element in instance.elements:
            if len(element.connectivity) == 4:  # Tet4 elements
                f.write("{}, {}, {}, {}, {}\n".format(element.label, *element.connectivity))
                tet4_elements.append(list(element.connectivity))
    
    verbose_print("Tet4 element connectivity written to:" + output_filename, verbose)
    return tet4_elements

# Function to identify boundary nodes based on element connectivity
def identify_boundary_nodes(tet4_elements, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    node_count = {}
    
    # Count occurrences of each node
    for element in tet4_elements:
        for node in element:
            node_count[node] = node_count.get(node, 0) + 1
    
    # Boundary nodes are those connected to fewer than 4 elements
    boundary_nodes = {node for node, count in node_count.items() if count < 4}
    return boundary_nodes

# Function to extract displacements for boundary nodes
def extract_displacements(odb, step_name, instance_name, boundary_nodes, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance " + str(instance_name) + " not found in ODB.", verbose)
        odb.close()
        return []
    
    instance = odb.rootAssembly.instances[instance_name]
    displacements_in_model = []
    
    if step_name not in odb.steps:
        verbose_print("Step '{}' not found in ODB.".format(step_name), verbose)
        odb.close()
        return []
    
    step = odb.steps[step_name]
    frame = step.frames[-1]  # Use the last frame
    
    displacement_field = frame.fieldOutputs['U']
    nodal_displacement = displacement_field.getSubset(position=NODAL)
    
    # Write boundary displacements to file
    with open(output_filename, 'w') as f:
        f.write("Node ID, Ux, Uy, Uz\n")
        for value in nodal_displacement.values:
            node_label = value.nodeLabel
            if node_label in boundary_nodes:
                ux, uy, uz = value.data
                f.write("{}, {}, {}, {}\n".format(node_label, ux, uy, uz))
                displacements_in_model.append((node_label, ux, uy, uz))
    
    verbose_print("Boundary displacements written to:" + output_filename, verbose)
    return displacements_in_model

# Function to write VTK file for Tet4 elements with displacements
def write_vtk_tet4(filename, nodes, tet4_elements, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Map node IDs to zero-based indices
    node_index_map = {node_id: idx for idx, node_id in enumerate(nodes.keys())}
    
    with open(filename, 'w') as f:
        # VTK Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK output for Tet4 mesh with displacement fields\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Write node coordinates
        f.write("POINTS {} float\n".format(len(nodes)))
        for _, x, y, z in nodes:
            f.write("{} {} {}\n".format(x, y, z))
        
        # Write Tet4 connectivity using remapped indices
        num_elements = len(tet4_elements)
        f.write("CELLS {} {}\n".format(num_elements, num_elements * 5))  # 4 nodes + 1 count per element
        for element in tet4_elements:
            f.write("4 {} {} {} {}\n".format(
                node_index_map[element[0]], node_index_map[element[1]],
                node_index_map[element[2]], node_index_map[element[3]]
            ))
        
        # Cell types (10 for tetrahedral cells in VTK)
        f.write("CELL_TYPES {}\n".format(num_elements))
        f.write("10\n" * num_elements)  # Tet4 cell type in VTK
        
        # Write displacements as scalar fields
        f.write("POINT_DATA {}\n".format(len(nodes)))
        
        # Ux displacements
        f.write("SCALARS Ux float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            ux = displacements.get(node_id, (0.0, 0.0, 0.0))[0]
            f.write("{}\n".format(ux))
        
        # Uy displacements
        f.write("SCALARS Uy float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            uy = displacements.get(node_id, (0.0, 0.0, 0.0))[1]
            f.write("{}\n".format(uy))
        
        # Uz displacements
        f.write("SCALARS Uz float 1\nLOOKUP_TABLE default\n")
        for node_id, _, _, _ in nodes:
            uz = displacements.get(node_id, (0.0, 0.0, 0.0))[2]
            f.write("{}\n".format(uz))
    
    verbose_print("VTK file for Tet4 mesh written to:" + filename, verbose)

def identify_element_types(odb, instance_name, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Identifies and counts element types within a specified instance.
    
    :param odb: The ODB object.
    :param instance_name: Name of the instance to inspect.
    :return: Dictionary with element types and their counts.
    """
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance '{}' not found in ODB.".format(instance_name), verbose)
        odb.close()
        return {}
    
    instance = odb.rootAssembly.instances[instance_name]
    element_type_counts = {}
    
    # Loop over elements in the instance and count each type
    for element in instance.elements:
        element_type = element.type
        if element_type in element_type_counts:
            element_type_counts[element_type] += 1
        else:
            element_type_counts[element_type] = 1
    
    verbose_print("Element types and counts for instance '{}':".format(instance_name), verbose)
    for element_type, count in element_type_counts.items():
        verbose_print("  {}: {}".format(element_type, count), verbose)
    
    return element_type_counts

def scan_odb_by_instance(odb, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Scans the ODB by instance and constructs a dictionary with element types for each instance.
    
    :param odb: The ODB object.
    :return: Dictionary with instance names as keys, each containing a dictionary of element types and counts.
    """
    instance_element_types = {}
    
    # Loop through each instance in the root assembly
    for instance_name, instance in odb.rootAssembly.instances.items():
        element_type_counts = {}
        
        # Loop over elements in the instance and count each type
        for element in instance.elements:
            element_type = element.type
            if element_type in element_type_counts:
                element_type_counts[element_type] += 1
            else:
                element_type_counts[element_type] = 1
        
        # Store the element type counts for the instance
        instance_element_types[instance_name] = element_type_counts
    
    # Print summary
    verbose_print("Element types by instance in ODB:", verbose)
    for instance_name, element_types in instance_element_types.items():
        verbose_print("Instance '{}':".format(instance_name), verbose)
        for element_type, count in element_types.items():
            verbose_print("  {}: {}".format(element_type, count), verbose)
    
    return instance_element_types

def get_step_names(odb, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Returns a list of all step names in the ODB.

    :param odb: The ODB object.
    :return: List of step names.
    """
    return list(odb.steps.keys())

def transform_coordinates_dict(origin, point1, point2, new_origin, new_point1, new_point2, node_coordinates, output_filename, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Transforms node coordinates by translating and rotating to a new reference system.
    """
    
    def compute_basis(p0, p1, p2):
        """Compute orthonormal basis vectors from three points."""
        v1 = np.array(p1) - np.array(p0)
        v2 = np.array(p2) - np.array(p0)
        x_axis = v1 / np.linalg.norm(v1)
        v2_proj = v2 - np.dot(v2, x_axis) * x_axis
        y_axis = v2_proj / np.linalg.norm(v2_proj)
        z_axis = np.cross(x_axis, y_axis)
        return np.array([x_axis, y_axis, z_axis])
    
    
    # Compute orthonormal bases
    original_basis = compute_basis(origin, point1, point2)
    new_basis = compute_basis(new_origin, new_point1, new_point2)
    verbose_print("Original Basis:\n" + str(original_basis), verbose)
    verbose_print("New Basis:\n" + str(new_basis), verbose)
    
    # Compute the rotation matrix
    rotation_matrix = np.dot(new_basis.T, original_basis)
    verbose_print("Rotation Matrix:\n" + str(rotation_matrix), verbose)
    
    # Compute the translation vector (to align the origins)
    translation_vector = np.array(new_origin) - np.dot(rotation_matrix, np.array(origin))
    verbose_print("Translation Vector:\n" + str(translation_vector), verbose)
    
    # Transform node coordinates
    transformed_nodes = {}
    with open(output_filename, 'w') as f:
        f.write("Node ID, X, Y, Z\n")
        for node_id, coord in node_coordinates.items():
            # Step 1: Translate coordinates to the new origin
            translated = np.array(coord) - np.array(origin)
            # Step 2: Apply rotation
            rotated = np.dot(rotation_matrix, translated)
            # Step 3: Translate to the final position
            transformed = rotated + np.array(new_origin)
            transformed_nodes[node_id] = tuple(transformed)
            f.write("{}, {}, {}, {}\n".format(node_id, transformed[0], transformed[1], transformed[2]))
    
    verbose_print("#### End transform_coordinates_dict ############################", verbose)
    return transformed_nodes

def extract_node_coordinates(odb, instance_name, node_number, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Extracts the coordinates of a specific node by its number from an Abaqus ODB file.
    
    :param odb: The ODB object.
    :param instance_name: Name of the instance containing the node.
    :param node_number: The node number whose coordinates need to be extracted.
    :return: Tuple of coordinates (x, y, z) or None if the node is not found.
    """
    if instance_name not in odb.rootAssembly.instances:
        verbose_print("Instance: " + instance_name + " not found in ODB.", verbose)
        odb.close()
        return None
    
    instance = odb.rootAssembly.instances[instance_name]
    
    # Find the target node in the instance
    for node in instance.nodes:
        if node.label == node_number:
            return node.coordinates
    
    verbose_print("Node number" + str(node_number) + " not found in instance " + str(instance_name) + ".", verbose)
    return None


def write_ply_with_displacement_colors(filename, nodes, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Writes a .ply file where the nodes are represented as points, and their color is based on displacement magnitude.
    
    :param filename: The output .ply file name.
    :param nodes: Dictionary of nodes {node_id: (x, y, z)} or list [(node_id, x, y, z)].
    :param displacements_list: List of displacements [(node_id, ux, uy, uz), ...].
    """
    # Convert displacements_list to a dictionary for easier access
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Determine node structure
    if isinstance(nodes, dict):
        node_items = nodes.items()
    elif isinstance(nodes, list):
        node_items = [(node_id, (x, y, z)) for node_id, x, y, z in nodes]
    else:
        raise TypeError("Unsupported nodes data structure. Expected dict or list.")
    
    # Compute displacement magnitudes and map to colors
    points_with_colors = []
    for node_id, (x, y, z) in node_items:
        ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
        displacement_magnitude = np.sqrt(ux**2 + uy**2 + uz**2)
        
        # Map displacement magnitude to RGB color (simple normalization example)
        max_displacement = 1.0  # You can adjust this based on your data
        normalized_magnitude = min(displacement_magnitude / max_displacement, 1.0)
        red = int(normalized_magnitude * 255)
        green = int((1 - normalized_magnitude) * 255)
        blue = 0  # Fixed blue component
        points_with_colors.append((x, y, z, red, green, blue))

# Write .ply file
def write_ply_with_displacement_colors(filename, nodes, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Writes a .ply file where the nodes are represented as points, and their color is based on displacement magnitude.
    
    :param filename: The output .ply file name.
    :param nodes: Dictionary of nodes {node_id: (x, y, z)} or list [(node_id, x, y, z)].
    :param displacements_list: List of displacements [(node_id, ux, uy, uz), ...].
    """
    # Convert displacements_list to a dictionary for easier access
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Determine node structure
    if isinstance(nodes, dict):
        node_items = nodes.items()
    elif isinstance(nodes, list):
        node_items = [(node_id, (x, y, z)) for node_id, x, y, z in nodes]
    else:
        raise TypeError("Unsupported nodes data structure. Expected dict or list.")
    
    # Compute displacement magnitudes and map to colors
    points_with_colors = []
    for node_id, (x, y, z) in node_items:
        ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
        # if ux == 0.0:
        #     verbose_print('ux is zero')
        displacement_magnitude = np.sqrt(ux**2 + uy**2 + uz**2)
        
        # Map displacement magnitude to RGB color (simple normalization example)
        max_displacement = 1.0  # You can adjust this based on your data
        normalized_magnitude = min(displacement_magnitude / max_displacement, 1.0)
        red = int(normalized_magnitude * 255)
        green = int((1 - normalized_magnitude) * 255)
        blue = 0  # Fixed blue component
        points_with_colors.append((x, y, z, red, green, blue))
    
    # Write .ply file
    with open(filename, 'w') as f:
        # PLY Header
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points_with_colors)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        
        # Write points with colors
        for x, y, z, r, g, b in points_with_colors:
            f.write("{} {} {} {} {} {}\n".format(x, y, z, r, g, b))
    
    
    verbose_print("PLY file written to " + str(filename), verbose)

def write_ply_with_displacement_and_colors(filename, nodes, displacements_list, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Writes a .ply file where the nodes are represented as deformed points, and their color is based on displacement magnitude.
    
    :param filename: The output .ply file name.
    :param nodes: Dictionary of nodes {node_id: (x, y, z)} or list [(node_id, x, y, z)].
    :param displacements_list: List of displacements [(node_id, ux, uy, uz), ...].
    :param verbose: 'On' to enable verbose output, 'Off' to disable.
    """
    def verbose_print(message, verbose):
        if verbose == 'On':
            print(message)
    
    # Convert displacements_list to a dictionary for easier access
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Determine node structure
    if isinstance(nodes, dict):
        node_items = nodes.items()
    elif isinstance(nodes, list):
        node_items = [(node_id, (x, y, z)) for node_id, x, y, z in nodes]
    else:
        raise TypeError("Unsupported nodes data structure. Expected dict or list.")
    
    # Compute deformed points and map displacement magnitudes to colors
    points_with_colors = []
    for node_id, (x, y, z) in node_items:
        ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
        deformed_x, deformed_y, deformed_z = x + ux, y + uy, z + uz  # Deformed coordinates
        displacement_magnitude = np.sqrt(ux**2 + uy**2 + uz**2)
        
        # Map displacement magnitude to RGB color (simple normalization example)
        max_displacement = 1.0  # You can adjust this based on your data
        normalized_magnitude = min(displacement_magnitude / max_displacement, 1.0)
        red = int(normalized_magnitude * 255)
        green = int((1 - normalized_magnitude) * 255)
        blue = 0  # Fixed blue component
        points_with_colors.append((deformed_x, deformed_y, deformed_z, red, green, blue))
    
    # Write .ply file
    with open(filename, 'w') as f:
        # PLY Header
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points_with_colors)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        
        # Write deformed points with colors
        for x, y, z, r, g, b in points_with_colors:
            f.write("{} {} {} {} {} {}\n".format(x, y, z, r, g, b))
    
    verbose_print("PLY file written to " + str(filename), verbose)

def write_ply_boundary_with_displacement(filename, nodes, displacements_list, boundary_nodes, verbose='On'):
    verbose_print("#### Entering the function: " + inspect.currentframe().f_code.co_name + " ####", verbose)
    """
    Writes a .ply file where the boundary nodes are represented as deformed points, 
    and their color is based on displacement magnitude.
    
    :param filename: The output .ply file name.
    :param nodes: Dictionary of nodes {node_id: (x, y, z)} or list [(node_id, x, y, z)].
    :param displacements_list: List of displacements [(node_id, ux, uy, uz), ...].
    :param boundary_nodes: List or set of boundary node IDs to include.
    :param verbose: 'On' to enable verbose output, 'Off' to disable.
    """
    def verbose_print(message, verbose):
        if verbose == 'On':
            print(message)
    
    # Convert displacements_list to a dictionary for easier access
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Ensure boundary_nodes is a set for efficient filtering
    boundary_nodes_set = set(boundary_nodes)
    
    # Determine node structure and filter by boundary nodes
    if isinstance(nodes, dict):
        node_items = {node_id: coords for node_id, coords in nodes.items() if node_id in boundary_nodes_set}
    elif isinstance(nodes, list):
        node_items = [(node_id, (x, y, z)) for node_id, x, y, z in nodes if node_id in boundary_nodes_set]
    else:
        raise TypeError("Unsupported nodes data structure. Expected dict or list.")
    
    # Compute deformed points and map displacement magnitudes to colors
    points_with_colors = []
    for node_id, (x, y, z) in node_items.items() if isinstance(node_items, dict) else node_items:
        ux, uy, uz = displacements.get(node_id, (0.0, 0.0, 0.0))  # Default to zero displacement
        deformed_x, deformed_y, deformed_z = x + ux, y + uy, z + uz  # Deformed coordinates
        displacement_magnitude = np.sqrt(ux**2 + uy**2 + uz**2)
        
        # Map displacement magnitude to RGB color (simple normalization example)
        max_displacement = 1.0  # You can adjust this based on your data
        normalized_magnitude = min(displacement_magnitude / max_displacement, 1.0)
        red = int(normalized_magnitude * 255)
        green = int((1 - normalized_magnitude) * 255)
        blue = 0  # Fixed blue component
        points_with_colors.append((deformed_x, deformed_y, deformed_z, red, green, blue))
    
    # Write .ply file
    with open(filename, 'w') as f:
        # PLY Header
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points_with_colors)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        
        # Write deformed boundary points with colors
        for x, y, z, r, g, b in points_with_colors:
            f.write("{} {} {} {} {} {}\n".format(x, y, z, r, g, b))
    
    verbose_print("PLY file written to " + str(filename), verbose)