# Finite Element Analysis with Adaptive Mesh Refinement (AMR) 
**Author:** Anurag Bhattacharyya

**Dependencies:**
The dependencies are listed in [`.github/workflows/Python_env.yml`](https://github.com/anurag-bha/AdaptiveFiniteElements/blob/main/.github/workflows/Python_env.yml):
* [MeshPy](https://pypi.org/project/MeshPy/)
* SciPy
 
**Problem Formulation:**
The main problem investigated is the 2d plane stress elasticity problem for a L-shaped domain defined mathematically:
```math
 
```
```math
\begin{split}
  \Delta^2 u(x, y) & = -F_b \hspace{5 mm} on \hspace{5 mm} \Omega  \\
  u & = 0 \hspace{5 mm} on \hspace{5 mm} \Omega_b
\end{split}
```
Here $F_b$ is the body force and $\Omega_b$ represents the boundary.

# Goals:
In this elasticity problem, there is also a body force acting on the entire structure
in the downward direction. From elasticity studies we know that the region
connecting the overhanging portion to the top  xed boundary will be highly
stressed and as a result the error in these elements will be high as well. I expect
that the adaptive mesh re nement process will capture these regions of high
internal stresses, consequently more error prone, and the mesh will be re ned
in this region to lower the peak in the internal stress values and subsequenly a
reduction in the error norm.


**Initial stress distribution over the domain:**
![](https://github.com/anurag-bha/AdaptiveFiniteElements/blob/main/Figs/Internal%20stress%20distribution%20over%20domain.png)

# Implementation:
The adaptive re nement algorithm follows a cycle of
_solve_->_estimate_->_mark_->_refine_.
We start with solve. A structural  nite element solver was fully developed
in-house coupled with the [MeshPy](https://pypi.org/project/MeshPy/) grid generator. Linear triangular elements
was used for meshing the design domain having 2 degrees of freedom at each
vertex. A body force to simulate loading due to gravity was impelmented. The
undeformed mesh (blue) and the deformed mesh (green) is shown in Figure.
After calculating the required values on the initial mesh, we estimate the errors
in the elements. For this we use a posteriori error estimator.
![](https://github.com/anurag-bha/AdaptiveFiniteElements/blob/main/Figs/Undeformed%20and%20deformed%20FEA%20meshes.png)


After the workflow is executed the adaptively refined mesh file is updated in the repo and is visualized below along with the deformed file
![](https://github.com/anurag-bha/AdaptiveFiniteElements/blob/main/Figs/Adaptive%20mesh%20refinement.png)


**Final stress distribution over the refined domain:**
![](https://github.com/anurag-bha/AdaptiveFiniteElements/blob/main/Figs/Internal%20stress%20distribution%20over%20refined%20mesh.png)

# Comparing AMR with full domain refinement


**Adaptive mesh refinement performance**
| D.O.Fs       | time(s)          | $$L_{inf}$$  |
| ------------- |:-------------:| -----:|
| 96            | 0.02997       |  0.34533|
| 508           | 0.15364       |   0.1145 |

**Whole domain refinement performance**
| D.O.Fs       | time(s)          | $$L_{inf}$$  |
| ------------- |:-------------:| -----:|
| 452            | 0.13805       |  0.3198|
| 564          | 0.15363       |   0.2338 |
| 674          | 0.16926       |   0.2139 |
| 812          | 0.2161       |   0.2149 |

