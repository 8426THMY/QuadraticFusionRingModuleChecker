# Haagerup-Izumi Module Checker
This program looks for special modules (so-called "NIM-reps") over Haagerup-Izumi-like fusion rings. Let $G$ be a finite Abelian group, $Q \leq G$ a subgroup and $\psi : G \to G$ an order-two group automorphism. From this data, we construct a ring with generators $\{g\}_{g \in G} \sqcup \{gX\}_{g \in G}$ and fusion rules

$$g \otimes X = gX, \quad X \otimes g = \psi(g)X, \quad X \otimes X = 𝟙 \oplus \bigoplus_{q \in Q}{qX}.$$

Here, we write $𝟙$ for the group identity. These rings show up as decategorifications of Haagerup-Izumi fusion categories ($Q = G$, $\psi(g) = g^{-1}$) and their crossed products.


# Usage
This program accepts one or more ring files as input,

`hi_algebra_checker ring1 ring2 ... ringn`.

The ring files themselves should contain the following information.

```
num_pgroups
pgroup1
...
pgroupn

num_subgroups
subgroup1
...
subgroupm

Q_index
psi_g1 ... psi_gk
```

Here, a p-group is specified by a prime, the number of cyclic factors and the prime powers of each factor. That is, `prime num_factors power1 ... powern`.

Currently, the program does not compute the subgroups of $G$ by itself, so they have to be specified manually. A subgroup is specified by its order followed by the indices of its elements, so `order g1_id ... gn_id`. For $G = \mathbb{Z}/p_1^{e_1}\mathbb{Z} \times \dots \times \mathbb{Z}/p_n^{e_n}\mathbb{Z}$, the index $g_{\textnormal{id}}$ of the group element $g = (g_1, \dots, g_n)$ is

$$g_{\textnormal{id}} = \sum_{i=1}^{n}{g_i\prod_{j=1}^{i-1}{p_j^{e_j}}}.$$

The order of these cyclic factors is determined by the order of declaration in the ring file. Finally, `Q_index` is simply the (zero-based) index of $Q$ in the list of subgroups and `psi_gi` is the index of $\psi(g_i)$, where $g_i$ is the $i$-th group element (ordered as above).


# Output

A module over these rings has exactly two orbits of generators under the action of $G$. For each module, we currently only output the action matrices for a single member of each orbit. Let $M$ be a generator of the module and $F_M$ its action matrix (as output by the program). Then $(F_M)_{i,j}$ is the multiplicity of the $j$-th generator of the module in the product of $M$ with the $i$-th generator of the ring.