# SAPT - Simulia Abaqus Python Tools

**SAPT** is a simple set of Python tools designed to inspect Abaqus files (.cae and .odb). This initial version focuses on extracting data from the `.odb` file.

## Overview

SAPT offers a range of functionalities to extract information from Abaqus ODB files, including:
- **Node and element data**: Retrieves node coordinates, boundary node connectivity, and element connectivity.
- **Displacement fields**: Extracts Ux, Uy, and Uz displacements at boundary nodes.
- **Element type detection**: Identifies the type of elements in each instance (e.g., `C3D8`, `C3D4`).
- **VTK file generation**: Produces VTK files with either scalar or vector fields for displacement visualization in compatible software.

The code also supports `Hex8` and `Tet4` element types, handling boundary nodes, connectivity, and displacements for each element type.

## How to Use SAPT

1. **Place the `extract_info_from_ODB.py` script** in your desired working directory, with access to Abaqus and a `.odb` file for testing.

2. **Run the script in Abaqus**:
   - In `CAE` mode:
     ```bash
     abaqus cae script=extract_info_from_ODB.py
     ```
   - In `nogui` mode (without graphical interface):
     ```bash
     abaqus cae nogui=extract_info_from_ODB.py
     ```

## Requirements

- **Abaqus** with Python support (the script imports Abaqus libraries and will run only within Abaqus).
- **An Abaqus ODB file** for analysis (replace `odb_name` and `odb_folder_path` variables as needed).

## Files Generated

For each instance and step in the ODB file, SAPT generates:

1. **Node Coordinates**: Files listing node coordinates in each instance.
2. **Boundary Nodes and Connectivity**: Lists boundary nodes and their connectivity.
3. **Element Connectivity**: Separate files for `Hex8` and `Tet4` elements, detailing connectivity.
4. **Displacements**: Ux, Uy, and Uz displacements at boundary nodes.
5. **VTK Files**: Triangular and Hexahedral VTK files for mesh and displacement visualization.

Files are saved in a `results` folder created in the output directory if it doesn’t already exist.

## Sample Script Explanation

Here’s a breakdown of the main functions:

- **`write_node_ids_and_coordinates()`**: Writes node coordinates for each instance in a step.
- **`write_boundary_nodes_and_connectivity()`**: Determines boundary nodes and writes their connectivity.
- **`write_triangular_faces()`**: Converts quadrilateral boundary faces into triangles for VTK file generation.
- **`extract_displacements()`**: Extracts Ux, Uy, and Uz displacement fields for nodes.
- **`write_vtk_file()` / `write_vtk_hex8()`**: Generates VTK files for both `Hex8` and `Tet4` meshes, with scalar or vector displacement fields for post-processing.

## Example Output

Each function writes the extracted data to a `.txt` file or a `.vtk` file in the `results` folder. Here is a sample VTK output for visualizing the displacement field:

```plaintext
# vtk DataFile Version 3.0
VTK output for Hex8 mesh with displacement fields
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 100 float
0.0 0.0 0.0
...
CELL_TYPES 50
12
...
POINT_DATA 100
VECTORS displacement float
0.0 0.0 0.0
...