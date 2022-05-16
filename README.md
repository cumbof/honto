# honto

**HO**mophily **N**etwork **TO**ol

<img src="https://anaconda.org/conda-forge/honto/badges/version.svg"> <img src="https://anaconda.org/conda-forge/honto/badges/downloads.svg">

`honto` is a tool designed for assessing and measuring homophily in networks whose nodes have categorical attributes, namely when the nodes of networks come partitioned into classes.

Homophily evaluation is performed through the comparison between the relative edge density of the subgraphs, induced by each class, and the corresponding expected relative edge density under a null model.

The novelty of our approach consists in prescribing an endogenous null model, namely, the sample space of the null model is built on the input network itself. This allows us to give exact explicit expressions for the z-scores of the relative edge density of each class as well as other related statistics.

## Install

`honto` is available through `pip` and `conda`.
Please, use one of the following commands to start working with `honto`:

```
# Install with pip
pip install honto

# Install with conda
conda install -c conda-forge honto
```

Please note that `honto` is also available as a Galaxy tool. It's wrapper is available under the official Galaxy ToolShed at [https://toolshed.g2.bx.psu.edu/view/fabio/honto](https://toolshed.g2.bx.psu.edu/view/fabio/honto)

## Usage

Once installed, you can start running `honto` by specifying some arguments:

```
honto --input_edges ~/edges.txt \
      --input_nodes ~/nodes.txt \
      --verbose
```

List of standard arguments:
```
--input_edges           -- Path to the file with the list of edges
--input_nodes           -- Path to the file with nodes and colors
--weight_threshold      -- Threshold for considering edges based in their weight
--isolated              -- Insert isolated nodes
--nproc                 -- Make the computation of the z-scores parallel for singletons
--overwrite             -- Overwrite results if already exist
--verbose               -- Print results in real time
-v, --version           -- Print current honto version and exit
```

List of arguments for log-transforming z-scores:
```
--log_transform     -- Log-transform z-scores
--scale_factor      -- Rescale z-scores with this constant before log-transforming values
--scale_from_one    -- Set z-scores to 1 if lower than 1 before log-transforming values
```

List of arguments for customizing the heatmap
```
--cmap      -- Heatmap colormap
--vmin      -- Min value to anchor the colormap
--vmax      -- Max value to anchor the colormap
--center    -- The value at which to center the colormap when plotting divergant data
--cbar      -- Whether to draw a colorbar
```

## Contributing

Long-term discussion and bug reports are maintained via GitHub Issues, while code review is managed via GitHub Pull Requests.

Please, (i) be sure that there are no existing issues/PR concerning the same bug or improvement before opening a new issue/PR; (ii) write a clear and concise description of what the bug/PR is about; (iii) specifying the list of steps to reproduce the behavior in addition to versions and other technical details is highly recommended.
