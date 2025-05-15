# Topological Indices Calculation Web App (Streamlit version of original script)
# Authors: Joudy Allam, Nour Saad
# Date: April 2025

# to run streamlit, 
# open Anaconda prompt > cd filepath > streamlit run app.py

# ================= Importing necessary libraries =================

import streamlit as st   # Web app framework to create interactive UI
import os                # For interacting with the operating system (e.g., file paths)
import tempfile          # For creating temporary files/folders during processing
import networkx as nx    # NetworkX for graph representations and topological index calculations
import plotly.graph_objs as go # Plotly for interactive graph visualizations
import zipfile           # For creating downloadable ZIP files
from rdkit import Chem   # RDKit for molecule parsing and chemical structure handling
from rdkit.Chem import AllChem     # RDKit's 3D coordinate generator and sanitization tools
import numpy as np                 # NumPy for fast numerical operations (used for distance matrix)
from collections import deque      # Deque for efficient BFS traversal in shortest path calculation
import io                # for text download

# ================= Original calculation functions preserved =================

# ----------------- Function that calculates edge density of a molecule -----------------
def edge_density(graph):
    """
                                       e
    Formula: edge density (ED) = -------------
                                 (n * (n-1))/2
                with e = actual number of bonds e, 
                     n = actual number of atoms (nodes)
                     (n * (n-1))/2 = maximum number of theoretical bonds possible

    Parameters: graph representation of molecule

    Returns: float: Edge density value between 0 and 1
    """
    n = graph.number_of_nodes()
    e = graph.number_of_edges()
    max_edges = n * (n - 1) / 2
    ED = e / max_edges if max_edges > 0 else 0
    return ED

# ----------------- Function that calculates Wiener Index of a molecule -----------------
def wiener_index(graph):
    """
    Formula: Wiener Index (W) = (1/2) * ∑∑ d(ni, nj) 
                                 for all pairs of nodes n i ≠ j

                where d(ni, nj) is the shortest path length (number of edges) 
                between nodes ni and nj in the graph.

    Parameters: graph representation of molecule

    Returns: float: Wiener index (always ≥ 0)
    """
    W = 0
    for node in graph.nodes():
        # Use BFS to calculate shortest path lengths from this node
        lengths = bfs_shortest_path(graph, node)
        W += sum(lengths.values())
    return W / 2

## Function to calculate Shortest Path for the Wiener Index
def bfs_shortest_path(graph, start_node):
    distances = {node: float('inf') for node in graph.nodes()}  # Initialize distances to infinity
    distances[start_node] = 0  # Distance to itself is zero
    queue = deque([start_node])  # Use a deque for efficient pop from the front

    while queue:
        current_node = queue.popleft()
        current_distance = distances[current_node]

        for neighbor in graph.neighbors(current_node):
            if distances[neighbor] == float('inf'):  # Not visited yet
                distances[neighbor] = current_distance + 1
                queue.append(neighbor)

    return distances


# ----------------- Function that calculates PetitJean Index of a molecule -----------------
def petitjean_index(graph):
    """
    Formula: Petitjean Index (PJ) = (D - R) / R

    where:
        D = diameter of the graph = max eccentricity (longest shortest path from any node)
        R = radius of the graph = min eccentricity (shortest longest path from any node)

    Parameters: graph representation of molecule

    Returns: float: Petitjean index (≥ 0)
    """
    try:
        ecc = calculate_eccentricity(graph)
        R = min(ecc.values())
        D = max(ecc.values())
        return (D - R) / R if R > 0 else 0
    except Exception as e:
        return f"Error: {str(e)}"

## Function that calculates eccentricities for the PetitJean Index
def calculate_eccentricity(graph):
    """Way to calculate eccentricity without using nx.eccentricity (which is used in the verification step)"""
    eccentricity = {}
    for node in graph.nodes():
        # Calculate shortest path lengths from this node to all others
        lengths = nx.single_source_shortest_path_length(graph, node)
        # The eccentricity is the maximum path length from the node
        if lengths:
            eccentricity[node] = max(lengths.values())
        else:
            eccentricity[node] = float('inf')  # If no paths are found (disconnected graph)
    return eccentricity

# ================= Validates custom topological index calculations =================
# ================= by comparing them with RDKit's built-in methods =================

# Calculation of topological indices using RDKit's built-in methods
def compute_rdkit_wiener(mol):
    dists = Chem.GetDistanceMatrix(mol)
    return np.sum(dists) / 2

def compute_rdkit_edge_density(mol):
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    max_edges = num_atoms * (num_atoms - 1) / 2
    return num_bonds / max_edges if max_edges > 0 else 0

# Validates custom topological index calculations (Wiener, Edge Density) 
# Note: Topological Petitjean Index is computed only using our custom implementation, 
#       and has been manually verified on paper for small molecules 
#       due to the lack of a built-in RDKit function.

def validate_with_builtins(mol_block, graph, exclude_h=True):
        """
        Note: If hydrogen atoms are excluded by user choice (`exclude_h=True`), they will also
              be removed from the RDKit molecule before computing built-in values. 
              This ensures consistency with the heavy-atom graph used in custom index calculations.
              
        Parameters:
            mol_block (str): Text content of a MOL file
            graph (networkx.Graph): Graph representation of the molecule
    
        Returns:
            dict: Comparison of (custom value, RDKit value) for each index
        """
        # Parse MOL block into RDKit molecule object 
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False)
        Chem.SanitizeMol(mol)  # Ensure valency and connectivity correctness
        
        # Respect user's choice: Remove hydrogens from the RDKit molecule only if exclude_h is True
        if exclude_h:
            mol = Chem.RemoveHs(mol)  # Exclude H atoms to match the heavy-atom graph
        
        # Generate 3D coordinates if needed (required for some descriptors)
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    
        # Generate 3D coordinates if needed
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        # RDKit built-in calculations
        builtin_wiener = compute_rdkit_wiener(mol)
        builtin_edge_density = compute_rdkit_edge_density(mol)

        # Custom graph-based calculations
        custom_wiener = wiener_index(graph)
        custom_edge_density = edge_density(graph)
        custom_petitjean = petitjean_index(graph)

        return {
            "Edge Density": (custom_edge_density, builtin_edge_density),
            "Wiener Index": (custom_wiener, builtin_wiener),
            "Petitjean Index": (custom_petitjean, "No built-in function available for the topological PetitJean.")
        }

# ================= Interactive visualization using Plotly with atom numbering =================

# ----------------- Define Atom colors for visualization with Plotly -----------------

# This dictionary maps chemical atom types to specific colors 
atom_color_map = {
    "C": "black",
    "H": "white",
    "O": "red",
    "N": "blue",
    "S": "yellow",
    "P": "purple",
    "Cl": "green",
    "F": "cyan",
    "Br": "brown",
    "I": "orange"
    }



# ----------------- Interactive visualization -----------------

def visualize_molecule(graph, title):
    
    """
    Automatically adapts to whether the molecule has 3D coordinates or not.
    If 3D, uses real spatial coordinates for plotting; otherwise, falls back to a 2D spring layout.
    Atoms are colored by element type, labeled with their IDs, and display hover information
    including atom type and coordinates. Bonds are drawn as edges between atoms.
    """

    # Step 1: Check if the graph contains 3D coordinates
    has_3d = False
    for node in graph.nodes():
        coordinates = graph.nodes[node].get('coordinates', (0, 0, 0))
        z_coord = coordinates[2]  # Extract the Z-coordinate
        if z_coord != 0:
            has_3d = True
            break  # No need to check further if one 3D coordinate is found

    # Step 2: Generate the layout based on the dimensionality
    if has_3d:
        # Uses 3D coordinates stored in the graph
        pos = {node: graph.nodes[node]['coordinates'] for node in graph.nodes()}
    else:
        # Back to 2D layout (spring layout)
        pos = nx.spring_layout(graph)

    # Step 3: Prepare node and edge data
    node_x, node_y, node_z = [], [], []
    node_color = []
    node_text = []  # Atom IDs to display as text on nodes
    hover_text = []

    for node in graph.nodes():
        # Extract coordinates, falling back to (0, 0, 0) if missing
        x, y, z = pos[node] if has_3d else (*pos[node], 0)
        node_x.append(x)
        node_y.append(y)
        node_z.append(z)
        
        # Atom type and color assignment
        atom_type = graph.nodes[node].get('type', 'Unknown')
        coordinates = graph.nodes[node].get('coordinates', (0, 0, 0))
        x_coord, y_coord, z_coord = coordinates
        color = atom_color_map.get(atom_type, 'lightgray') # Default gray for unknowns
        node_color.append(color)

        # Atom label on the node (numeric ID)
        node_text.append(str(node))

        # Hover tooltip showing atom type and coordinates
        hover_text.append(
            f"Atom: {atom_type} (ID: {node})<br>X: {x_coord:.3f}<br>Y: {y_coord:.3f}<br>Z: {z_coord:.3f}"
        )


    # Extract edge positions for drawing bonds
    edge_x, edge_y, edge_z = [], [], []
    for edge in graph.edges():
        x0, y0, z0 = pos[edge[0]] if has_3d else (*pos[edge[0]], 0)
        x1, y1, z1 = pos[edge[1]] if has_3d else (*pos[edge[1]], 0)
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])

    # Step 4: Plot in 3D if Z-coordinates are present
    
    # Plotly 3D scatter objects are used for both atoms (nodes) and bonds (edges). 
    
    if has_3d:
        edge_trace = go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode="lines",
            line=dict(width=2, color="#888"),
            hoverinfo="none"
        )

        node_trace = go.Scatter3d(
            x=node_x, y=node_y, z=node_z,
            mode="markers+text",
            marker=dict(
                size=8,
                color=node_color,
                line=dict(width=2, color="white")
            ),
            text=node_text,
            textposition="top center",
            textfont=dict(size=10, color="white"),
            hovertext=hover_text,
            hoverinfo="text"
        )
        
        # We hid axis lines and labels for a cleaner look.
        
        fig = go.Figure(data=[edge_trace, node_trace],
               layout=go.Layout(
                   title=title,
                   showlegend=False,
                   margin=dict(b=0, l=0, r=0, t=40),
                   scene=dict(
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title="", showbackground=False, visible=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title="", showbackground=False, visible=False),
                       zaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title="", showbackground=False, visible=False),
                   )
               ))


    # Step 5: Plot in 2D if no Z-coordinates are present
    
    # 2D plots use Plotly's 2D scatter plots. Coordinates are from a layout algorithm (spring_layout).
    else:
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=2, color="#888"),
            hoverinfo="none",
            mode="lines"
        )

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode="markers+text",
            marker=dict(
                size=20,
                color=node_color,
                line=dict(width=2, color="white")
            ),
            text=node_text,
            textposition="middle center",
            textfont=dict(size=12, color="white"),
            hovertext=hover_text,
            hoverinfo="text"
        )

        fig = go.Figure(data=[edge_trace, node_trace],
                       layout=go.Layout(
                           title=title,
                           showlegend=False,
                           hovermode="closest",
                           margin=dict(b=0, l=0, r=0, t=40),
                           xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                       ))

    # Step 6: Display the plot on Streamlit
    st.plotly_chart(fig, use_container_width=True)

# ================= Reading Input Data =================

# ----------------- Reading molecular files (.mol) -----------------
def read_mol_file(filepath, fallback_name=None, exclude_h=True):

    """
    Reads a .mol file and converts it into a NetworkX graph representation,
    while also extracting the molecule name from the first line.
    If no name exists in the file, uses the filename (without extension) as fallback.
    
    The "exclude_h" option can optionally exclude hydrogen atoms from the graph.
    When set to True, only non-hydrogen atoms are considered, which is useful 
    for the heavy-atom-based topological index calculations.
    
    Parameters:
        filepath (str): Path to the .mol file.
        fallback_name (str, optional): Alternative name if the molecule name is not found in the file.
        exclude_h (bool): If True, hydrogen atoms (type 'H') are ignored in the graph construction.
    
    Returns:
        tuple: (molecule_name, networkx.Graph)
    """

    G = nx.Graph()
    with open(filepath, 'r') as file:
        lines = file.readlines()

        # Use first line for molecule name or fallback to filename
        if lines and lines[0].strip():
            molecule_name = lines[0].strip()
        else:
            # Use fallback_name (e.g., original uploaded name) without extension
            if fallback_name:
                molecule_name = os.path.splitext(os.path.basename(fallback_name))[0]
            else:
                molecule_name = os.path.splitext(os.path.basename(filepath))[0]


        try:
            atom_count = int(lines[3][:3].strip())
            bond_count = int(lines[3][3:6].strip())
        except (IndexError, ValueError):
            st.error("Error: Invalid atom or bond count.")
            return molecule_name, G

        atoms = lines[4:4 + atom_count]
        bonds = lines[4 + atom_count:4 + atom_count + bond_count]

        if exclude_h:
            # When excluding hydrogens, skip atoms of type 'H' and reindex the graph
            atom_index_map = {}
            new_index = 1
            
            for i, atom in enumerate(atoms):
                atom_type = atom[31:34].strip()
                if atom_type == 'H':
                    continue  # Skip hydrogen
                try:
                    x_coord = float(atom[0:10].strip())
                    y_coord = float(atom[10:20].strip())
                    z_coord = float(atom[20:30].strip())
                    G.add_node(new_index, type=atom_type, coordinates=(x_coord, y_coord, z_coord))
                    atom_index_map[i + 1] = new_index  # original index → new index
                    new_index += 1
                except ValueError:
                    st.warning(f"Warning: Invalid atom coordinates in line {i+1}")
            # Add edges only between non-hydrogen atoms
            for bond in bonds:
                try:
                    start = int(bond[:3].strip())
                    end = int(bond[3:6].strip())
                    if start in atom_index_map and end in atom_index_map:
                        G.add_edge(atom_index_map[start], atom_index_map[end])
                except ValueError:
                    st.warning("Warning: Invalid bond data.")
        else:
            # if chosen not to exclude the hydrogen atoms, 
            # include all atoms (with hydrogens) with original indices
            for i, atom in enumerate(atoms):
                try:
                    atom_type = atom[31:34].strip()
                    x = float(atom[0:10].strip())
                    y = float(atom[10:20].strip())
                    z = float(atom[20:30].strip())
                    G.add_node(i + 1, type=atom_type, coordinates=(x, y, z))
                except ValueError:
                    st.warning(f"Warning: Invalid atom coordinates in line {i + 1}")
            for bond in bonds:
                try:
                    start = int(bond[:3].strip())
                    end = int(bond[3:6].strip())
                    G.add_edge(start, end)
                except ValueError:
                    st.warning("Warning: Invalid bond data.")
    return molecule_name, G

# ----------------- Read folder of MOL files as zip file -----------------
def read_mol_folder_zip(folderpath, exclude_h):
    """
    Reads all .mol files from a folder extracted from a ZIP archive, and processes
    each using `read_mol_file()` to build molecular graphs.
    
    Parameters:
        folderpath (str): Path to the extracted ZIP folder containing .mol files.
        exclude_h (bool): If True, hydrogen atoms are excluded from the graphs.
    
    Returns:
        list of tuples: (molecule_name, graph, mol_block) for each molecule.
    """

    graphs = []
    
    # Loop through all files in the extracted zip folder
    for filename in os.listdir(folderpath):
        if filename.endswith('.mol'):
            # Construct full path to the file
            filepath = os.path.join(folderpath, filename)
            
            # Read and convert the mol file into a graph, using filename as fallback name
            molecule_name, graph = read_mol_file(filepath, fallback_name=filename, exclude_h=exclude_h)

            
            with open(filepath, 'r') as f:
                mol_block = f.read()
            
            # Append molecule name, graph, and raw mol content to the list
            graphs.append((molecule_name, graph, mol_block))
            
    return graphs

# ----------------- Read folder of MOL files from Directory path  -----------------
def read_mol_folder_dir(folderpath, exclude_h):
    """
    Reads all .mol files from a specified directory path.
    Each .mol file is read using read_mol_file(), which converts the structure into a graph.
    Also reads and stores the raw content of each .mol file (mol_block) for later use.
    
    Parameters:
        folderpath (str): Path to the folder with .mol files
        exclude_h (bool): If True, hydrogen atoms are excluded from the graphs.
    
    Returns:
        list of tuples: (molecule_name, graph, mol_block) for each molecule
    """
    
    graphs = []

    # Loop through every file in the folder
    for filename in os.listdir(folderpath):
        if filename.endswith('.mol'):
            filepath = os.path.join(folderpath, filename)

            # Read molecule structure, using filename as fallback if needed
            molecule_name, graph = read_mol_file(filepath, fallback_name=filename, exclude_h=exclude_h)

            with open(filepath, 'r') as f:
                mol_block = f.read()

            # Append molecule name, graph, and raw mol content to the list
            graphs.append((molecule_name, graph, mol_block))

    return graphs

# ----------------- Read SDF file -----------------
def read_sdf_file(filepath, exclude_h):
    """
    Reads an .sdf file containing multiple molecules and converts each into a graph.
    
    Each molecule block is extracted and saved temporarily as a .mol file.
    The molecule name is extracted from the first line of the block.
    If the first line is empty, the name defaults to the temporary filename. 

    
    Parameters:
        filepath (str): Path to the .sdf file
        exclude_h (bool): If True, hydrogen atoms are excluded from the graphs.
    
    Returns:
        list of tuples: (molecule_name, graph, mol_block) for each molecule
    """
    molecules = []
    with open(filepath, 'r') as file:
        content = file.read()
        # Since each molecule is separated by "$$$$"
        # Split molecules using '$$$$' followed by a newline
        sdf_blocks = [m.strip() for m in content.split('$$$$\n') if m.strip()]
        
        # Process each molecule block individually
        for i, block in enumerate(sdf_blocks):
            # Try to extract a molecule name from the first line of the block
            molecule_name = block.splitlines()[0].strip() if block.splitlines() and block.splitlines()[0].strip() else f"molecule_{i+1}"

            # Create a temporary MOL file
            with tempfile.NamedTemporaryFile(delete=False, suffix=".mol", mode='w') as temp_file:
                temp_file.write(block + '\n')  # # Write molecule block and add a newline at the end
                temp_filepath = temp_file.name

            try:
                # Convert the temporary MOL file into a graph
                name, graph = read_mol_file(temp_filepath, fallback_name=molecule_name, exclude_h=exclude_h)
                
                mol_block = block + '\n'
                
                # Store tuple (name, graph, block) for later processing
                molecules.append((name, graph, mol_block))
            finally:
                # Clean up temporary file 
                os.remove(temp_filepath)

    return molecules   # Return list of all parsed molecules as graphs

# ============================================================================================================
# NOTE on Input Source:
# The input files used in this application may come from different chemical databases: PubChem, ChemSpider...
#
# These databases assign different IDs to the same molecule.
# Because of this, we cannot reliably link the input file to its original database entry.
# In most cases, even the user might not remember the exact source of the file.
#
# Furthermore, .mol and .sdf formats do not explicitly store the source database name.
# Therefore, no automatic linking to PubChem or ChemSpider web pages is implemented from within this tool.
# ============================================================================================================


# ================= Interface =================

# ----------------- Streamlit UI -----------------

# Set Streamlit page configuration: title and layout
st.set_page_config(page_title="Topological Indices Calculator", layout="wide")
st.markdown("""
    <style>
    .block-container {
        padding-top: 2rem;
    }

    .top-right-credit {
        position: absolute;
        top: 1rem;
        right: 2rem;
        font-size: 1rem;
        font-family: 'Merriweather', serif;
        color: #999;
    }

    .title-text {
        font-size: 3rem;
        font-family: 'Orbitron', serif;
        font-weight: 600;
        color: white;
        text-align: center;
        letter-spacing: 1px;
        text-transform: uppercase;
        margin-top: 40px; 
        padding: 1rem 0;
        border-bottom: 2px solid white;
        margin-bottom: 20px;
    }
    </style>
""", unsafe_allow_html=True)

# Display the app title and authors
st.markdown("""
<div class='top-right-credit'>Developed by Joudy Allam and Nour Saad</div>
<div class='title-text'>🧪 Topological Indices Calculator</div>
""", unsafe_allow_html=True)


# ----------------- Upload Input Data -----------------

# User can choose whether to include or exclude hydrogen atoms from the molecular graph.
# Topological indices are typically computed using only heavy atoms (non-H),
# so by default (value=True), hydrogen atoms are excluded unless the user unchecks the box.
exclude_h = st.checkbox("Exclude Hydrogen Atoms from Calculations", value = True)

# User selects whether they want to upload a single file (whether .mol or .sdf) or a folder of files
upload_type = st.radio("Select Input Type:", ["Single File (.mol or .sdf)", "Folder of .mol Files"])
# User selects which index (or all) to calculate
selected_index = st.radio("Select Index to Calculate:", ["Edge Density", "Wiener Index", "PetitJean Index", "All"])

# Initialize molecule storage list
molecules = []

# Case 1: Upload a single .mol or .sdf file 
if upload_type == "Single File (.mol or .sdf)":
    file = st.file_uploader("Upload a .mol or .sdf file", type=["mol", "sdf"])
    if file is not None:
        # Temporarily save the uploaded file
        suffix = os.path.splitext(file.name)[-1]
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            tmp.write(file.read())
            tmp_path = tmp.name
            
        # Read MOL files directly
        if file.name.endswith(".mol"):
            with open(tmp_path, 'r') as f:
                mol_block = f.read()
            name, graph = read_mol_file(tmp_path, fallback_name=file.name, exclude_h = exclude_h)
            molecules = [(name, graph, mol_block)]
    
        # For SDF files, use custom reader that extracts each molecule
        elif file.name.endswith(".sdf"):
             molecules = read_sdf_file(tmp_path, exclude_h = exclude_h)  

# Case 2: Upload a folder of .mol files 
if upload_type == "Folder of .mol Files":
    # User selects whether to upload a ZIP folder or enter path manually
    folder_input_type = st.radio("Select Folder Input Method:", ["Upload ZIP file", "Specify Folder Path"])
    
    # --- Upload ZIP file containing MOL files ---
    if folder_input_type == "Upload ZIP file":
        folder_zip = st.file_uploader("Upload a ZIP file containing .mol files", type=["zip"])
        if folder_zip is not None:
            with tempfile.TemporaryDirectory() as tmpdir:
                # Extract the uploaded ZIP to a temporary directory
                with zipfile.ZipFile(folder_zip, 'r') as zip_ref:
                    zip_ref.extractall(tmpdir)
                
                # Use your helper function to read all .mol files inside
                for name, graph, mol_block in read_mol_folder_zip(tmpdir, exclude_h=exclude_h):
                    molecules.append((name, graph, mol_block))

                            
    # --- Manually specify a local directory path containing .mol files ---
    elif folder_input_type == "Specify Folder Path":
        folder_path = st.text_input("Enter the folder path containing .mol files:")
        if folder_path:
            for name, graph, mol_block in read_mol_folder_dir(folder_path, exclude_h=exclude_h):
                molecules.append((name, graph, mol_block))


# ----------------- Display -----------------

# Once molecules are loaded, display calculated indices and visualizations
if molecules:
    st.success(f"Loaded {len(molecules)} molecule(s).")
    export_lines = []  # collect lines for exportable text summary

    for name, graph, mol_block in molecules:
        # Section header for each molecule
        st.markdown(f"### 🔬 Molecule: {name}")
        cols = st.columns([1, 2])

        # Left column: Show the selected topological indices
        with cols[0]:
            ed, wi, pj = None, None, None
            if selected_index == "Edge Density" or selected_index == "All":
                ed = edge_density(graph)
                st.metric("Edge Density", f"{ed:.4f}")
            if selected_index == "Wiener Index" or selected_index == "All":
                wi = wiener_index(graph)
                st.metric("Wiener Index", f"{wi:.4f}")
            if selected_index == "PetitJean Index" or selected_index == "All":
                pj = petitjean_index(graph)
                st.metric("PetitJean Index", f"{pj:.4f}")

            # --- Add Button for Calculated Indices ---
            export_text = f"{name}\n"
            if ed is not None:
                export_text += f"Edge Density: {ed:.4f}\n"
            if wi is not None:
                export_text += f"Wiener Index: {wi:.4f}\n"
            if pj is not None:
                export_text += f"PetitJean Index: {pj:.4f}\n"
            export_lines.append(export_text)

        # Right column: Show interactive molecule visualization
        with cols[1]:
            visualize_molecule(graph, name)

        # -- Check Match Button --
        if mol_block and st.button(f"🔍 Check Match for {name}"):
            result = validate_with_builtins(mol_block, graph, exclude_h)
            if "Error" in result:
                st.error(result["Error"])
            else:
                st.markdown("#### 🔎 Custom vs Built-in Comparison using RDKit Library")
                for key, (custom_val, builtin_val) in result.items():
                    st.write(f"{key}")
                    st.write(f"• Our code: {custom_val:.4f}")
                    if isinstance(builtin_val, (float, int)):
                        st.write(f"• Built-in : {builtin_val:.4f}")
                        diff = abs(custom_val - builtin_val)
                        color = "green" if diff < 0.01 else "red"
                        st.markdown(f"<span style='color:{color}'>Difference: {diff:.4f}</span>", unsafe_allow_html=True)
                    else:
                        st.write(f"• Built-in : {builtin_val}")
                st.markdown("---")

    # Download button for all results combined
    combined_result = "\n\n".join(export_lines)
    buffer = io.BytesIO(combined_result.encode("utf-8"))  # Convert to BytesIO for Streamlit compatibility
    st.download_button("📄 Download All Indices as TXT", buffer, file_name="molecule_indices.txt", mime="text/plain")

# If no file is uploaded yet
else:
    if upload_type == "Single File (.mol or .sdf)":
        st.info("Please upload a .mol or .sdf file.")
    elif upload_type == "Folder of .mol Files":
        st.info("Please upload a ZIP archive or specify a valid folder path containing .mol files.")