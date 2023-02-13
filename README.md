# PM

## Steiner tree problem

The input of the classic STP consists of an undirected edge-weighted  graph $G$ with vertex set $V(G)$, edge set $E(G)$, a non-negative weight function on its edges $c_G:E(G)\rightarrow \mathbb{R}^{+}$, and a set of \emph{terminals} $A\subseteq V(G)$.
The problem is to determine a minimum Steiner tree, i.e., a tree $T$ spanning all the terminals in $A$ (and possibly other vertices) and minimizing $\sum_{e\in E(T)}c_G(e)$.

## Algorithm in large graphs

We introduced a new partition-and-merge (PM) algorithm for effectively solving large-scale STP instances. To ensure the effectiveness of the algorithm, we investigated different merging heuristics, local search algorithms, and new solution storing techniques.

The algorithm breaks the input graph into small subgraphs and then merges the subgraphs in a bottom-up manner. During the merging procedure, partial Steiner trees in the subgraphs are also created and optimized by an efficient local optimization. When the merging procedure ends, the algorithm terminates and reports the final solution for the input graph.

We compared the algorithm with the state-of-the-art STP heuristics, including the winners of the 11th DIMACS and 3rd PACE competitions. Experiments on well-known large benchmark instances validated the competitiveness of our algorithm, for solving these large STP instances.

## Acknowledgements

We thank the authors of Staynerd, PUW, CIMAT and other solvers in the paper for providing their codes. We also thank the organizers of 11th DIMACS challenge and 3rd PACE competition for hosting the data.
