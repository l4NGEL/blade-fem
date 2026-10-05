# Formulation notes

Short notes on what the code discretises and why each choice was made. Notation: N_a are
the linear shape functions of a triangle with area A, and b_a = ∂N_a/∂x, c_a = ∂N_a/∂y
are constant over the element.

## Steady conduction

Strong form: −∇·(k∇T) = Q in Ω, with three boundary conditions:

- T = T_D on Γ_D;
- −k ∂T/∂n = h(T − T∞) on Γ_h;
- k ∂T/∂n = q_in on Γ_q.

Weak form: find T with T = T_D on Γ_D such that, for all v vanishing on Γ_D,

    ∫Ω k ∇T·∇v dA + ∫Γh h T v ds = ∫Ω Q v dA + ∫Γh h T∞ v ds + ∫Γq q_in v ds

Element and edge matrices for linear triangles and 2-node edges of length l:

- Conduction: K_ab = k A (b_a b_b + c_a c_b).
- Convection: (h l / 6)·[[2, 1], [1, 2]] on each convective edge, with load h T∞ l / 2 per node.
- Source: integrated with the 7-point Dunavant rule (degree 5), so smooth manufactured
  sources add no quadrature error at the rates being measured.

Dirichlet values are eliminated rather than penalised. The free-free block stays symmetric
positive definite, and the "reaction" r = K T − f at a fixed node is the consistent nodal
heat flow through that boundary. The reported heat balance therefore closes to round-off:
Σ(Dirichlet reactions) + Σ(convective heat) + Σ(prescribed flux) + ∫Q = 0.

## Linear elasticity (CST)

With ε = (εxx, εyy, γxy) = B u_e (B is constant per element), the stiffness is
K_e = A Bᵀ D B and the stress is σ = D (B u_e − ε₀).

| Model | D | Free thermal strain ε₀ | σzz |
|---|---|---|---|
| Plane stress | E/(1−ν²)·[[1, ν, 0], [ν, 1, 0], [0, 0, (1−ν)/2]] | αΔT (1, 1, 0) | 0 |
| Plane strain | E/((1+ν)(1−2ν))·[[1−ν, ν, 0], [ν, 1−ν, 0], [0, 0, (1−2ν)/2]] | (1+ν) αΔT (1, 1, 0) | ν(σxx+σyy) − EαΔT |

ΔT is the element mean of the nodal temperatures minus T_ref. The thermal load vector is
f_e = A Bᵀ D ε₀. Surface tractions are integrated with 3-point Gauss on each edge, so the
non-uniform exact tractions of the reference problems are applied consistently.

## Generalized plane strain

The axial strain is a plane over the section, measured from its area centroid (x₀, y₀):

    εzz(x, y) = e₀ + κx (x − x₀) + κy (y − y₀)

The three numbers (e₀, κx, κy) are global unknowns appended after the 2N displacements.
The strain vector becomes (εxx, εyy, γxy, εzz) and the element operator gains a fourth row
[0 … 0 | 1, x − x₀, y − y₀]. D is the full 3D isotropic matrix restricted to those four
components, and ε₀ = αΔT (1, 1, 0, 1).

No load is applied to the three global unknowns, so their equations state

    ∫ σzz dA = 0,   ∫ σzz (x − x₀) dA = 0,   ∫ σzz (y − y₀) dA = 0,

that is, zero axial force and zero bending moments: the ends of a long body are free.
Because εzz is linear over each element, the stiffness integrand is quadratic and is
integrated exactly with a 3-point rule.

Checks:

- A uniform ΔT gives e₀ = αΔT and zero stress.
- A linear T gives the matching curvature.
- A heated thick cylinder matches the free-ended closed form in all four components.

## Rigid-body modes and singular systems

A free body (the blade) is held by a pin and a roller, which is statically determinate,
so the self-equilibrated thermal load produces zero reactions up to round-off. That is
asserted in the tests.

If constraints are missing, the reduced matrix is singular. `LinearSystem` detects this
from the LDLᵀ pivots (smallest < 10⁻¹² × largest) and throws, instead of returning a
solution dominated by round-off.

## Error measures

- L² and H¹ errors use the 7-point rule on each element against the exact field.
- Stress errors compare the element (constant) stress with the exact stress at the
  centroid, weighted by area.
- The observed order between two meshes is log(e₁/e₂) / log(h₁/h₂).

Expected orders for linear triangles on smooth problems:

- displacement or temperature in L²: 2;
- gradient, strain and stress: 1;
- the Kirsch SCF read from node-averaged stresses happened to converge at about 2 here.
