# BME205 Asgn4: Longest Path in a Directed Acyclic Graph (DAG)

This Python script is designed to solve Problem 15 from Rosalind, which involves finding the longest path in a directed acyclic graph (DAG) between two specified nodes. It utilizes a modified adjacency list format to represent the graph, where edges are denoted as "source->target:weight".

## Usage

To use this script, you need to provide an input file containing the graph data in the specified format. The input file should contain three lines:

1. An integer representing the source node of the graph.
2. An integer representing the sink node of the graph.
3. Edge-weighted graph represented in the modified adjacency list format.

Here's an example of how to run the script:

```bash
python longest_path_in_dag.py input.txt
```

Replace `input.txt` with the path to your input file.

## Input Format

The input file should follow the format described above. For example:

```makefile
0
4
0->1:7
0->2:4
2->3:2
1->4:1
3->4:3
```

## Output Format

The script will output the length of the longest path in the graph followed by the longest path itself. The output format is as follows:

```rust
11
0->2->3->4
```

This indicates that the longest path has a length of 11 and consists of the nodes 0, 2, 3, and 4 in sequence.