# Quadratic Fusion Ring Module Checker
This program looks for special modules (so-called "NIM-reps") over quadratic fusion rings. Currently, we only support quadratic fusion rings over finite Abelian groups and whose non-invertible generators have irrational Frobenius-Perron dimension. Such rings are determined by a finite Abelian group $G$, a subgroup $Q \leq G$, an order-two group automorphism $\psi : G/Q \to G/Q$, an fixed point $r(Q) \in G/Q$ of $\psi$ and non-negative integers $\\{m\_{gQ}\\}\_{gQ \in G/Q}$ satisfying $m\_{g^{-1}r(Q)} = m\_{gQ}$ and $m\_{g\psi(hQ)} = m\_{ghQ}$ for all $gQ, hQ \in G/Q$.

Given the data described above, we obtain a quadratic fusion ring $R$ with generators $\\{g\\}\_{g \in G} \sqcup \\{X\_{gQ}\\}\_{gQ \in G/Q}$ and fusion rules

$$g \otimes X\_Q = X\_{gQ}, \quad X\_Q \otimes g = X\_{\psi(gQ)}, \quad X \otimes X = \bigoplus\_{q \in Q}{q} \oplus \bigoplus\_{gQ \in G/Q}{X\_{gQ}^{\oplus m\_{gQ}}}.$$

The duals are given by $X\_{gQ}^{\*} = X\_{r(Q)\psi(g^{-1}Q)}$. These rings show up everywhere in the theory of fusion categories, with famous examples being decategorifications of near-group ($Q = G$) and Haagerup-Izumi ($Q = 𝟙$, $\psi(g) = g^{-1}$, $m_g = 1$) fusion categories. They also show up as crossed products of the latter family, which was the original inspiration for this program.


## Usage
This program accepts one or more ring files as input,

`quadratic_module_checker ring1 ring2 ... ringn`.

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
r(Q)
m_g1 ... m_gk
```

Here, a p-group is specified by a prime, the number of cyclic factors and the prime powers of each factor. That is, `prime num_factors power1 ... powern`.

Currently, the program does not compute the subgroups of $G$ by itself, so they have to be specified manually. A subgroup is specified by its order followed by the indices of its elements, so `order g1_id ... gn_id`. For $G = \mathbb{Z}/p_1^{e_1}\mathbb{Z} \times \dots \times \mathbb{Z}/p_n^{e_n}\mathbb{Z}$, the index $g_{\textnormal{id}}$ of the group element $g = (g_1, \dots, g_n)$ is

$$g_{\textnormal{id}} = \sum_{i=1}^{n}{g_i\prod_{j=1}^{i-1}{p_j^{e_j}}}.$$

The order of these cyclic factors is determined by the order of declaration in the ring file. Finally, `Q_index` is simply the (zero-based) index of $Q$ in the list of subgroups, `psi_gi` is the index of any element of $\psi(g_i Q)$ where $g_i$ is the $i$-th group element (ordered as above), `r(Q)` is the index of any element of $r(Q)$ and `m_gi` is the value of $m\_{g_i Q} \in \mathbb{Z}\_{\geq 0}$.


## Output

Since we assume the non-invertible generators have irrational Frobenius-Perron dimension, a module over these rings has exactly two orbits of generators under the action of $G$. For each module, we compute the fusion matrices for a single member of each orbit. Let $M$ be a generator of the module and $F_M$ its action matrix (as output by the program). Then $(F_M)_{i,j}$ is the multiplicity of the $j$-th generator of the module in the product of $M$ with the $i$-th generator of the ring.

Despite computing the fusion matrices internally, we currently only output the corresponding algebra objects in the fusion category. I'll probably get around to updating this soon.