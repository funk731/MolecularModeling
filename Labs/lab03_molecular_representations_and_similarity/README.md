# Molecular Representations and Similarity

Author(s): Enes Kelestemur

**Notes:** Generative AI tools, such as ChatGPT and Gemini, were used while creating this lecture notebook.

## Overview

This tutorial explores various molecular representations and similarity metrics used in cheminformatics. It covers key metrics for comparing molecular structures and illustrates their applications using Python. The notebook provides practical examples, visualizations, and comparisons of various molecular similarity measures.

## Lab Structure

### 1. Distance and Similarity Metrics

- **Euclidean Distance**: Understanding the straight-line distance between two points in a multi-dimensional space.
- **Manhattan Distance**: Measures the distance between points by summing the absolute differences of their coordinates.
- **Chebyshev Distance**: Uses the maximum absolute difference between the coordinates of the points.
- **Minkowski Distance**: A generalized distance metric that includes both Euclidean and Manhattan distances.
- **Hamming Distance**: Counts the number of differing positions between two binary vectors.
- **Jaccard Distance and Tanimoto Similarity**: Evaluates the similarity between sets, particularly useful in comparing molecular fingerprints.
- **Dice Similarity**: Focuses on the overlap between sets, emphasizing the intersection over their total size.
- **Cosine Similarity**: Measures the cosine of the angle between two vectors, often used in high-dimensional spaces.

### 2. Descriptor Representations

- **Chemical Descriptors**: Exploration of numerical values that represent physical, chemical, and structural properties of molecules.
- **Fingerprints**: Methods like Morgan (ECFP) and RDKit fingerprints to encode molecular structures.
  - **Visualizing Fingerprint Bits**: Techniques to visualize and interpret the bits in molecular fingerprints.

### 3. Language Representations

- **SMILES**: Simplified Molecular Input Line Entry System for encoding molecules as strings.
- **SMARTS**: A language for specifying molecular patterns and querying substructures.

### 4. Graph Representations

- **Intro to Graphs**: Introduction to graph theory concepts relevant to molecular structures.
- **Molecular Graphs**: Representing molecules as graphs, with atoms as nodes and bonds as edges.
- **Advanced Topics for Molecular Graphs**: Delving deeper into graph-based representations.
  - **Invariant and Equivariant Graph Representations**: Exploring advanced graph representations that are invariant to permutations and equivariant under group actions.

## Learning Outcomes

By the end of this tutorial, participants will be able to:

- Apply a variety of distance and similarity metrics for molecular comparison.
- Use Python libraries to calculate molecular similarities and explore different representation techniques.
- Differentiate between descriptor, language, and graph-based molecular representations.
- Utilize appropriate similarity measures in cheminformatics applications.

## Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [Scikit-learn Documentation](https://scikit-learn.org/stable/)
- [Seaborn Documentation](https://seaborn.pydata.org/)
- [Networkx Documentation](https://networkx.org/)
- [Distance metric explanation and images](https://medium.com/@jodancker/a-brief-introduction-to-distance-measures-ac89cbd2298)
- [TeachOpenCADD T036 E(3) Equivariant GNN](https://github.com/volkamerlab/teachopencadd/tree/master/teachopencadd/talktorials/T036_e3_equivariant_gnn)
- [SE(3)-Transformers: 3D Roto-Translation Equivariant Attention Networks](https://proceedings.neurips.cc/paper/2020/file/15231a7ce4ba789d13b722cc5c955834-Paper.pdf)
- [E(n) Equivariant Graph Neural Networks](https://proceedings.mlr.press/v139/satorras21a.html)

