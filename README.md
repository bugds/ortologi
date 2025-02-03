# Pavlov's COGs

A program for building clusters of orthologous groups (COG) graphs.

![image](https://github.com/user-attachments/assets/d297e72d-f18b-4839-b810-1f6e52797b2c)

## Building the Database for Search

The database installation process is outlined in the **cogcreatedb.ipynb** file.

### Dependencies for Database Creation

- [Pandas](https://pandas.pydata.org/)
- [NCBI Entrez Programming Utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [NCBI Datasets Command-Line Tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/)
- [BLAST Command-Line Applications](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

Carefully review, modify, and run the **cogcreatedb.ipynb** file to create your BLAST database. By default, the script will create a database of eukaryotic RefSeq protein sequences from three organisms of each taxonomic class.

## Running the COG Analysis

### Dependencies for the COG Analysis

A **dependencies.yml** file is provided for installing the required dependencies. It is recommended to use this file with specified versions of packages:

- [NetworkX](https://networkx.org/)
- [Markov Clustering](https://github.com/guyallard/markov_clustering)
- [Biopython](https://biopython.org/)
- [Pyvis](https://pyvis.readthedocs.io/) (new versions of Pyvis might be laggy)
- [BLAST Command-Line Applications](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

### Running the Analysis

The **cog.py** script performs the COG analysis. It requires two mandatory arguments:

1. **Input file**: Contains identifiers of the queried sequences.
2. **Output folder**: A desired empty folder for storing results.

#### Additional Arguments

- `-b/--init`, `-e/--eval`, and `-q/--qcov`: Specify the number of hits, E-value, and query coverage limits for the initial BLAST search.
- `-t/--threadNum`: Number of threads to use.
- `-o/--orthology`: Minimum fraction of best BLAST hits among gene isoforms required to establish orthology between two genes.
- `-s/--stage`: Run a specific stage of the analysis:
  - `1`: Initial BLAST search.
  - `2`: Orthology analysis.
  - `3`: Building the COG graph.
- `-m/--merge`: Merge the resulting graphs into one (useful when studying a family of genes).
- `-r/--removeXML`: Keep the (usually large) XML files of BLAST search results.
- `-a/--algorithm`: Select the algorithm:
  - `strict`: Uses all protein isoforms of a gene.
  - `best`: Uses only the representative isoform (as in [OMA](https://doi.org/10.1093/nar/gkaa1007)).
- `-g/--gravity`: Adjust the gravity constant, which affects the graph visualization.
- `-c/--config`: Specify a custom configuration file containing:
  - Paths to the table linking gene and organism identifiers.
  - Taxonomy table.
  - Name of the BLAST database to use.
  - Paths to `blastp` and `blastdbcmd` utilities.

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
