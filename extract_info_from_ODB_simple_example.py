# Clear all user-defined variables from the global namespace
for name in dir():
    if not name.startswith("_"):
        del globals()[name]

from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os


# Function to write node IDs and coordinates to a file
def write_node_ids_and_coordinates(odb, step_name, instance_name, output_filename):
    # Access the specified instance
    if instance_name not in odb.rootAssembly.instances:
        print("Instance: " + instance_name + " not found in ODB.")
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
    print("Node IDs and coordinates have been written to : " + str(output_filename))
    return nodes_in_model, nodes_in_model_dic

# Function to write boundary nodes and their connectivity to a file
def write_boundary_nodes_and_connectivity(odb, instance_name, output_filename):
    # Check if instance exists
    if instance_name not in odb.rootAssembly.instances:
        print("Instance: " + instance_name + " not found in ODB.")
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
    print("Boundary nodes and connectivity have been written to : " + output_filename)
    return boundary_nodes, q4_faces_in_model

# Function to write external face connectivity as triangles
def write_triangular_faces(odb, instance_name, boundary_nodes, filename):
    # Check if instance exists
    if instance_name not in odb.rootAssembly.instances:
        print("Instance: " + instance_name + " not found in ODB.")
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
    
    print("Triangular face connectivity has been written to", filename)
    return faces_triangular_in_model

# Function to extract Ux, Uy, Uz displacements at boundary nodes
def write_boundary_displacements(odb, step_name, instance_name, boundary_nodes, output_filename):
    # Check if instance exists
    if instance_name not in odb.rootAssembly.instances:
        print("Instance: " + instance_name + " not found in ODB.")
        odb.close()
        return
    
    instance = odb.rootAssembly.instances[instance_name]
    displacements_in_model = []
    
    # Check if the step exists
    if step_name not in odb.steps:
        print("Step: " + step_name + " not found in ODB.")
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
    
    print("Boundary displacements (Ux, Uy, Uz) have been written to: " + output_filename)
    return displacements_in_model

# Function to extract Hex8 elements from an ODB file
def extract_hex_elements(odb, instance_name):
    """
    Extracts Hex8 (8-node hexahedral) elements from an ODB instance.
    
    :param odb: The ODB object.
    :param instance_name: Name of the instance from which to extract elements.
    :return: List of Hex8 elements with node connectivity [[node1, node2, ..., node8], ...].
    """
    # Check if the instance exists
    if instance_name not in odb.rootAssembly.instances:
        print("Instance: " + instance_name + " not found in ODB.")
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
    
    print("Extracted {} Hex8 elements from instance '{}'.".format(len(hex_elements), instance_name))
    return hex_elements


def write_vtk_file(filename, nodes, triangles, displacements_list):
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    
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

def write_vtk_file_warp(filename, nodes, triangles, displacements_list):
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    
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
    
    print("VTK file for triangular mesh with vector displacements written to:", filename)

def write_vtk_hex8(filename, nodes, hex_elements, displacements_list):
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    
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

def write_vtk_hex8_warp(filename, nodes, hex_elements, displacements_list):
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Create zero-based index mapping
    node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    
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
    
    print("VTK file for Hex8 mesh with vector displacements written to:", filename)

# Function to write node IDs and coordinates to a file and return the node list and dictionary
def write_node_ids_and_coordinates_2(odb, instance_name, output_filename):
    if instance_name not in odb.rootAssembly.instances:
        print("Instance '{}' not found in ODB.".format(instance_name))
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
    print("Node IDs and coordinates written to:", output_filename)
    return nodes_in_model, nodes_in_model_dic

# Function to extract Tet4 elements and write connectivity to a file
def extract_tet4_elements(odb, instance_name, output_filename):
    if instance_name not in odb.rootAssembly.instances:
        print("Instance '{}' not found in ODB.".format(instance_name))
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
    
    print("Tet4 element connectivity written to:", output_filename)
    return tet4_elements

# Function to identify boundary nodes based on element connectivity
def identify_boundary_nodes(tet4_elements):
    node_count = {}
    
    # Count occurrences of each node
    for element in tet4_elements:
        for node in element:
            node_count[node] = node_count.get(node, 0) + 1
    
    # Boundary nodes are those connected to fewer than 4 elements
    boundary_nodes = {node for node, count in node_count.items() if count < 4}
    return boundary_nodes

# Function to extract displacements for boundary nodes
def extract_displacements(odb, step_name, instance_name, boundary_nodes, output_filename):
    if instance_name not in odb.rootAssembly.instances:
        print("Instance '{}' not found in ODB.".format(instance_name))
        odb.close()
        return []
    
    instance = odb.rootAssembly.instances[instance_name]
    displacements_in_model = []
    
    if step_name not in odb.steps:
        print("Step '{}' not found in ODB.".format(step_name))
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
    
    print("Boundary displacements written to:", output_filename)
    return displacements_in_model

# Function to write VTK file for Tet4 elements with displacements
def write_vtk_tet4(filename, nodes, tet4_elements, displacements_list):
    displacements = {node_id: (ux, uy, uz) for node_id, ux, uy, uz in displacements_list}
    
    # Map node IDs to zero-based indices
    node_index_map = {node_id: idx for idx, (node_id, x, y, z) in enumerate(nodes)}
    
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
    
    print("VTK file for Tet4 mesh written to:", filename)

def identify_element_types(odb, instance_name):
    """
    Identifies and counts element types within a specified instance.
    
    :param odb: The ODB object.
    :param instance_name: Name of the instance to inspect.
    :return: Dictionary with element types and their counts.
    """
    if instance_name not in odb.rootAssembly.instances:
        print("Instance '{}' not found in ODB.".format(instance_name))
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
    
    print("Element types and counts for instance '{}':".format(instance_name))
    for element_type, count in element_type_counts.items():
        print("  {}: {}".format(element_type, count))
    
    return element_type_counts

def scan_odb_by_instance(odb):
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
    print("Element types by instance in ODB:")
    for instance_name, element_types in instance_element_types.items():
        print("Instance '{}':".format(instance_name))
        for element_type, count in element_types.items():
            print("  {}: {}".format(element_type, count))
    
    return instance_element_types

def get_step_names(odb):
    """
    Returns a list of all step names in the ODB.

    :param odb: The ODB object.
    :return: List of step names.
    """
    return list(odb.steps.keys())


#####################################################

# Usage
odb_name = 'example.odb'
odb_folder_path = 'C:/temp\Extract_info_from_ODB/'
odb_path = odb_folder_path + odb_name
results_folder = odb_folder_path + 'results/'
# Check if the folder exists
if not os.path.exists(results_folder):
    # Create the folder if it does not exist
    os.makedirs(results_folder)
    print("Directory created:", results_folder)
else:
    print("Directory already exists:", results_folder)

# Open the output database
odb = openOdb(path=odb_path)
# Scan ODB by instance to get element types
instance_element_types = scan_odb_by_instance(odb)
# Get list of steps
list_of_steps = get_step_names(odb)
#set the Results view
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
    print("Instance: " + instance_name + "Count of Hex8 elements:", count_of_hex8)
    count_of_tet4 = element_types.get('C3D4', 0)
    print("Instance: " + instance_name + "Count of Tet4 elements:", count_of_tet4)
    for step_name in list_of_steps:
        # COORDINATES
        output_filename = results_folder +  instance_name + '_' + step_name + '_' + 'node_coordinates.txt'
        nodes_in_model, nodes_in_model_dic = write_node_ids_and_coordinates(odb, step_name, instance_name, output_filename)
        #
        if instance_element_types[instance_name].get('C3D8',0)>0:            
            # BOUNDARY NODES AND CONNECTIVITY (ORIGINAL)
            output_filename = results_folder +  instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'boundary_nodes_connectivity.txt'
            boundary_nodes, q4_faces_in_model = write_boundary_nodes_and_connectivity(odb, instance_name, output_filename)
            #
            # Extract Hex8 elements
            hex_elements = extract_hex_elements(odb, instance_name)
            #
            # TRIANGULARE FACES
            output_filename = results_folder +  instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'triangular_faces.txt'
            faces_triangular_in_model = write_triangular_faces(odb, instance_name, boundary_nodes, output_filename)
            #
            # DISPLACEMENTS
            output_filename = results_folder +  instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'boundary_displacements.txt'
            displacements_in_model = write_boundary_displacements(odb, step_name, instance_name, boundary_nodes, output_filename)
            #
            # VTK file
            output_filename = results_folder +  instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'triangular_mesh_U.vtk'
            write_vtk_file_warp(output_filename, nodes_in_model, faces_triangular_in_model, displacements_in_model)
            #
            # VTK file
            output_filename = results_folder +  instance_name + '_' + step_name + '_' + 'C3D8' + '_'  + 'Hexa_mesh_U.vtk'
            write_vtk_hex8_warp(output_filename, nodes_in_model, hex_elements, displacements_in_model)
        if instance_element_types[instance_name].get('C3D4',0)>0:
            # Nodes and Coordinates
            nodes_output = results_folder +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_nodes.txt'
            nodes_in_model, nodes_in_model_dic = write_node_ids_and_coordinates_2(odb, instance_name, nodes_output)
            #
            # Tet4 Connectivity
            connectivity_output = results_folder +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_connectivity.txt'
            tet4_elements = extract_tet4_elements(odb, instance_name, connectivity_output)
            #
            # Boundary Nodes
            boundary_nodes = identify_boundary_nodes(tet4_elements)
            #
            # Boundary Displacements
            displacements_output = results_folder +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_boundary_displacements.txt'
            displacements_in_model = extract_displacements(odb, step_name, instance_name, boundary_nodes, displacements_output)
            #
            # Write VTK file
            vtk_output = results_folder +  instance_name + '_' + step_name + '_' + 'C3D4' + '_'  + 'tet4_mesh_U.vtk'
            write_vtk_file_warp(vtk_output, nodes_in_model, tet4_elements, displacements_in_model)
        if instance_element_types[instance_name].get('C3D4',0)>0 and instance_element_types[instance_name].get('C3D8',0)>0:
            print("Element type not found") 

close_odb = True #False
if close_odb:
    odb.close()

