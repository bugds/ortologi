# Instructions

## Viewing the Results

After the COG analysis, a Pyvis HTML file is generated in the **Results** folder. Upon opening the file, a loading bar will appear. If the loading bar gets stuck at 0%, there may be an issue with your file.

## Interface

Once the file has loaded, a COG graph will appear at the top of the screen, and various submenus will be displayed at the bottom.

### COG Graph

The COG graph consists of **nodes** representing genes, connected by **branches** that indicate orthology between two genes.

- To select nodes, click on them with the **left mouse button**.
- To select multiple genes, hold the **Control** (or **Command** on macOS) key while clicking.
- Nodes can be dragged to rearrange the graph for better visualization.

In standard mode:
- Genes from the organism used in the query are highlighted in **red**.
- Nodes in the largest maximal cliques containing the red nodes are colored **yellow**.
- All other nodes are displayed in **dark blue (navy)**.

### Select Submenu

The **Select** submenu includes the following features:

- **Search Box**: Use the upper input box to search for specific nodes. Enter a fragment of the desired string (e.g., an organism name like "Homo sapiens") and click to select all nodes matching the fragment.
- **Group Selection**: Use the dropdown menu at the bottom of the submenu to add nodes to the current selection. Options include:
  - All nodes
  - Only visible nodes
  - Only hidden nodes
  - Nodes connected to currently selected nodes
  - Markov clusters containing currently selected nodes
- Click the **Select group** button to apply your selection. User-defined groups can also be added to the current selection using this dropdown menu.

### Editing Submenu

The **Editing** submenu allows you to switch between two modes:

1. **Painting Nodes**: Choose between coloring nodes based on the largest maximal cliques or Markov clustering. Select the desired mode and click the **Change colors** button to apply it.
2. **Physics Toggle**: Use the **Physics** checkbox to enable or disable node physics.

### Textbox

- **Selected Genes**: Names of selected genes will appear in the textbox upon selection.
- **Hide/Reveal Nodes**: Use the **Hide** button to hide selected nodes or the **Reveal** button to make them visible again.
- **Color Customization**: Use the **Paint** option to color selected nodes. Enter the desired color in the **rgb(r,g,b)** input box (e.g., `rgb(255,0,0)` for red).
- **User-Defined Groups**: Assign nodes to user-defined groups by typing the group name into the input box at the bottom. Click the **Assign group** button to save the group.

### Inspect Textbox Content

The textbox includes three buttons for further analysis:

1. **Get FASTA**: Retrieves sequences of representative isoforms for the selected genes.
2. **Description (visible)**: Provides information on the connections and composition of Markov clusters for the visible selected nodes.
3. **Build MSA**: Uses the EMBL-EBI Clustal Omega service to build a multiple sequence alignment (MSA) of the protein sequences of the selected nodes.