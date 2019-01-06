var documenterSearchIndex = {"docs": [

{
    "location": "basics/hartreefockbogoliubov/#",
    "page": "Hartree-Fock-Bogoliubov Theory",
    "title": "Hartree-Fock-Bogoliubov Theory",
    "category": "page",
    "text": ""
},

{
    "location": "basics/hartreefockbogoliubov/#Hartree-Fock-Bogoliubov-Theory-1",
    "page": "Hartree-Fock-Bogoliubov Theory",
    "title": "Hartree-Fock-Bogoliubov Theory",
    "category": "section",
    "text": "This document contains basics of Hartree-Fock-Bogoliubov Theory. We will follow [Goodman][Goodman80].  H =  sum_ij T_ij c_i^* c_j\n     + frac14 sum_ijkl V_ijkl c_i^* c_j^* c_l c_kThe finite-temperature properties of the system described by the above Hamiltonian can be computed using the density matrixbeginalign\nD = frace^-beta HZ \nZ = mathrmTr e^-beta H\nendalignIn Hartree-Fock-Bogoliubov theory, we approximate the Hamiltonian H by a non-interacting Hamiltonian of \"quasiparticles\"H_textHFB = E_0 + sum_n E_n a_n^* a_nThe quasiparticle operator a is related to c in the following way:c_i = sum_n left( U_in a_n + V_in a_n^* right)The HFB density matrix is given bybeginalign\n  D_textHFB = frac1Z_textHFB e^-beta sum_n E_n hatn_n  \n  Z_textHFB = mathrmTr e^-beta sum_n E_n hatn_n \nendalignwherehatn_n = a_n^* a_nis the quasiparticle number operator.The HFB partition function can be calculated:Z_textHFB = prod_n left( 1 + e^-beta E_n right)and the density matrix is given in terms of the quasiparticle number operatorD_textHFB = Z_textHFB^-1\nprod_n\nleft\n  e^-beta E_n hatn_n + (1 - hatn_n)\nright\n= prod_n left f_n hatn_n + (1-f_n) (1 - hatn_n) rightwheref_n = frac11 + e^beta E_nis the Fermi-Dirac distribution function for the nth quasiparticleSingle-quasparticle density matrix overlinerho  and pairing tensor overlinetbeginalign\noverlinerho_ij\n  = leftlangle a_j^* a_i rightrangle\n   = mathrmTr left( D a_j^* a_i right) \noverlinet_ij\n  = leftlangle a_j a_i rightrangle\n   = mathrmTr left( D a_j a_i right)\nendalignWithin HFB,beginalign\n  overlinerho_ij = delta_ij f_i \n  overlinet_ij = 0\nendalignThe single-particle density matrix and pairing tensor arebeginalign\nrho_ij\n  = leftlangle c_j^* c_i rightrangle\n  = mathrmTr left( D c_j^* c_i right) \nt_ij\n  = leftlangle c_j c_i rightrangle\n  = mathrmTr left( D c_j c_i right)\nendalignbeginalign\n  rho = U f U^dagger + V (1-f) V^dagger \n  t    = U f V^intercal + V (1-f) U^intercal\nendalignwhere f_ij = delta_ij f_i."
},

{
    "location": "basics/hartreefockbogoliubov/#Expectation-Values-1",
    "page": "Hartree-Fock-Bogoliubov Theory",
    "title": "Expectation Values",
    "category": "section",
    "text": "Wick\'s theoremleftlangle c_i^* c_j^* c_l c_k rightrangle\n=\nleftlangle c_i^* c_k rightrangle\nleftlangle c_j^* c_l rightrangle\n- leftlangle c_i^* c_l rightrangle\nleftlangle c_j^* c_k rightrangle\n+ leftlangle c_i^* c_j^* rightrangle\nleftlangle c_l c_k rightranglebeginalign\nE = mathrmtr\nleft\n  left( T + frac12 Gamma right) rho + frac12 Delta t^dagger\nright \nS = - k_B sum_i left f_i ln f_i + (1-f_i) ln (1-f_i)right \nN = mathrmtr rho\nendalignbeginalign\nGamma_ij = sum_kl V_ikjl rho_lk \nDelta_ij = frac12 sum_kl V_ijkl t_kl \nendalignThe grand potentialOmega\n= sum_ij (T - mu)_ij rho_ji\n+ frac12 sum_ijkl V_ijkl rho_lj rho_ki\n+ frac14 sum_ijkl V_ijkl t_ij^* t_kl\n+ k_B T sum_i left f_i ln f_i + (1-f_i) ln (1-f_i)right[Goodman80]: http://escholarship.org/uc/item/3xm630cr \"Alan Goodman, Finite-Temperature Hartree-Fock-Bogoliubov Theory, LBNL Paper LBL-11151 (1980)\""
},

{
    "location": "basics/hfbdecomposition/#",
    "page": "Hartree-Fock-Bogoliubov Decomposition",
    "title": "Hartree-Fock-Bogoliubov Decomposition",
    "category": "page",
    "text": ""
},

{
    "location": "basics/hfbdecomposition/#Hartree-Fock-Bogoliubov-Decomposition-1",
    "page": "Hartree-Fock-Bogoliubov Decomposition",
    "title": "Hartree-Fock-Bogoliubov Decomposition",
    "category": "section",
    "text": "We can write the Hamiltonian in the way that makes the Hermiticity manifest, taking into account duplicates properlybeginalign\nH = sum_ij T_ij c_i^* c_j\n + frac14 sum_ijkl c_i^* c_j^* c_l c_k\n\n  = sum_i T_ii c_i^* c_i\n    + sum_i lt j (T_ij c_i^* c_j + T_ij^* c_j^* c_i) \n  quad\n  + frac14 sum_i  j\n  left\n      V_ijij c_i^* c_j^* c_j c_i\n    - V_ijij c_j^* c_i^* c_j c_i\n    - V_ijij c_i^* c_j^* c_i c_j\n    + V_ijij c_j^* c_i^* c_i c_j\n  right\n  \n  quad\n  + frac14 sum_substack(ij) (kl)  text and   (ij)  (kl) \n  left\n      V_ijkl     c_i^* c_j^* c_l c_k\n    - V_ijkl     c_j^* c_i^* c_l c_k\n    - V_ijkl     c_i^* c_j^* c_k c_l\n    + V_ijkl     c_j^* c_i^* c_k c_l right\n  qquad left\n    + V_ijkl^* c_k^* c_l^* c_j c_i\n    - V_ijkl^* c_l^* c_k^* c_j c_i\n    - V_ijkl^* c_k^* c_l^* c_i c_j\n    + V_ijkl^* c_l^* c_k^* c_i c_j\n  right\n\n= sum_i T_ii c_i^* c_i\n  + sum_i lt j (T_ij c_i^* c_j + T_ij^* c_j^* c_i) \nquad\n  + sum_i  j\n  left\n      V_ijij c_i^* c_j^* c_j c_i\n  right\n  + sum_substack(ij) (kl)  text and   (ij)  (kl) \n  left\n      V_ijkl     c_i^* c_j^* c_l c_k\n    + V_ijkl^* c_k^* c_l^* c_j c_i\n  right\nendalignbeginalign\n  Gamma_ij = sum_kl V_ikjl rho_lk \n  Gamma_ik = sum_jl V_ijkl rho_lj \n  Delta_ij = frac12 sum_kl V_ijkl t_kl\nendalign"
},

{
    "location": "basics/hfbdecomposition/#Diagonal-Interaction-1",
    "page": "Hartree-Fock-Bogoliubov Decomposition",
    "title": "Diagonal Interaction",
    "category": "section",
    "text": "Let us first consider the diagonal interaction of the form  V c_1^* c_2^* c_2 c_1 which is the Hermitian conjugate of itself. Following [Goodman][Goodman80], this interaction decomposes into the following mean fields in the particle-hole and particle-particle channels:Interaction PH Channel PP Channel\nV(1212) =  V Γ(11) =  V ρ(22) Δ(12) =  ½ V t(12)\nV(1221) = -V Γ(12) = -V ρ(12) Δ(12) = -½ V t(21)\nV(2112) = -V Γ(21) = -V ρ(21) Δ(21) = -½ V t(12)\nV(2121) =  V Γ(22) =  V ρ(11) Δ(21) =  ½ V t(21)Here, the equal sign \'=\' does not indicate equality. Rather, it refers to the contribution of the right hand side to the mean field on the left hand side. Some of the mean fields are redundant, as required by the Hermiticity of the Hamiltonian. Overall, we end up withbeginalign\nGamma_11 =  V rho_22 \nGamma_12 = -V rho_12 \nGamma_22 =  V rho_11\nendalignin the particle-hole channel, andbeginalign\nDelta_12 =  V t_12\nendalignin the particle-particle channel."
},

{
    "location": "basics/hfbdecomposition/#Offdiagonal-Interaction-1",
    "page": "Hartree-Fock-Bogoliubov Decomposition",
    "title": "Offdiagonal Interaction",
    "category": "section",
    "text": "Now let us consider the offdiagonal interaction term of the form  V c_1^* c_2^* c_4 c_3Given this term in the interaction, it is implied that its Hermitian conjugate  V^* c_3^* c_4^* c_2 c_1is also included in the Hamiltonian. The two terms (the explicit term and its Hermitian conjugate) decomposes intoInteraction PH Channel PP-Channel\nV(1234) =  V Γ(13) =  V ρ(42) Δ(12) =  ½ V t(34)\nV(1243) = -V Γ(14) = -V ρ(32) Δ(12) = -½ V t(43)\nV(2134) = -V Γ(23) = -V ρ(41) Δ(21) = -½ V t(34)\nV(2143) =  V Γ(24) =  V ρ(31) Δ(21) =  ½ V t(43)\n  \nV(3412) =  V* Γ(31) =  V ρ(24) Δ(34) =  ½ V t(12)\nV(3421) = -V* Γ(32) = -V ρ(14) Δ(34) = -½ V t(21)\nV(4312) = -V* Γ(41) = -V ρ(23) Δ(43) = -½ V t(12)\nV(4321) =  V* Γ(42) =  V ρ(13) Δ(43) =  ½ V t(21)As we have mentioned above for the diagonal interaction terms, this table contains redundant mean fields. Overwall. we havebeginalign\nGamma_13 =  V  rho_24^* \nGamma_14 = -V  rho_23^* \nGamma_23 = -V  rho_14^* \nGamma_24 =  V  rho_13^*\nendalignin the particle-hole channel andbeginalign\nDelta_12 =  V   t_34 \nDelta_34 =  V^* t_12\nendalignin the particle-particle channel."
},

{
    "location": "basics/hfbequation/#",
    "page": "HFB Equations",
    "title": "HFB Equations",
    "category": "page",
    "text": ""
},

{
    "location": "basics/hfbequation/#HFB-Equations-1",
    "page": "HFB Equations",
    "title": "HFB Equations",
    "category": "section",
    "text": "delta Omega = 0"
},

{
    "location": "basics/momentumspace/#",
    "page": "Momentum Space Formulation",
    "title": "Momentum Space Formulation",
    "category": "page",
    "text": ""
},

{
    "location": "basics/momentumspace/#Momentum-Space-Formulation-1",
    "page": "Momentum Space Formulation",
    "title": "Momentum Space Formulation",
    "category": "section",
    "text": ""
},

{
    "location": "basics/momentumspace/#Hopping-Elements-1",
    "page": "Momentum Space Formulation",
    "title": "Hopping Elements",
    "category": "section",
    "text": "sum_mathbfR\n    c_alpha^*(mathbfR + boldsymbolrho_alpha)\n    c_beta(mathbfR + boldsymbolrho_beta)\n  =\n    sum_mathbfk\n    e^i mathbfk cdot (boldsymbolrho_beta - boldsymbolrho_alpha)\n    c_alpha^*(mathbfk)\n    c_beta(mathbfk)Pairing Elementssum_mathbfR\n    c_alpha^*(mathbfR + boldsymbolrho_alpha)\n    c_beta^*(mathbfR + boldsymbolrho_beta)\n  =\n    sum_mathbfk\n    e^i mathbfk cdot (boldsymbolrho_beta - boldsymbolrho_alpha)\n    c_alpha^*(mathbfk)\n    c_beta^*(-mathbfk)sum_mathbfR\n    c_alpha(mathbfR + boldsymbolrho_alpha)\n    c_beta(mathbfR + boldsymbolrho_beta)\n  =\n    sum_mathbfk\n    e^i mathbfk cdot (boldsymbolrho_beta - boldsymbolrho_alpha)\n    c_alpha(-mathbfk)\n    c_beta(mathbfk)"
},

{
    "location": "basics/momentumspace/#Expectation-Values-1",
    "page": "Momentum Space Formulation",
    "title": "Expectation Values",
    "category": "section",
    "text": "rho_alpha beta\n  =\n      frac1N\n      sum_mathbfk\n      e^-i mathbfk cdot ( boldsymbolrho_beta - boldsymbolrho_alpha ) \n      sum_n\n       f(epsilon_n) U_alpha n U_beta n^*t_alphabeta\n  =\n      frac1N\n      sum_mathbfk_alpha mathbfk_beta\n      e^-i mathbfk cdot ( boldsymbolrho_beta\n                        -boldsymbolrho_alpha ) \n      sum_n\n        f(epsilon_n )\n        U_alpha n\n        V_beta n^*"
},

{
    "location": "basics/momentumspace/#Grand-Potential-1",
    "page": "Momentum Space Formulation",
    "title": "Grand Potential",
    "category": "section",
    "text": "E = mathrmtr left left( T + frac12 Gamma right) rho + frac12 Delta t^dagger right\n  = E_T + E_Gamma + E_Deltawherebeginalign\nE_T = mathrmtr left( T rho right) \nE_Gamma = frac12 mathrmtr left( Gamma rho        right) \nE_Delta = frac12 mathrmtr left( Delta t^dagger right)\nendalignIn momentum spacebeginalign\nT_ij^mathbfk = T_ij e^i mathbfk cdot boldsymbolrho_ij\nrho_ij = frac1N sum_mathbfk rho_ij^mathbfk e^-i mathbfk cdot boldsymbolrho_ij\nendalignThusbeginalign\nE_T\n  = sum_ij T_ij rho_ji\n   = sum_ij T_ij\n     frac1N sum_mathbfk rho_ji^mathbfk e^-i mathbfk cdot boldsymbolrho_ji \n   = frac1N sum_mathbfk sum_ij T_ij e^i mathbfk cdot boldsymbolrho_ij \n      rho_ji^mathbfk\n   = frac1N sum_ij\n      sum_mathbfk\n      T_ij^mathbfk\n      rho_ji^mathbfk\nendalignSimilarly for Γ and Δ:beginalign\nE_Gamma\n  = frac12 frac1N sum_ij\n    sum_mathbfk\n    Gamma_ij^mathbfk\n    rho_ji^mathbfk \nE_Delta\n  = frac12 frac1N sum_ij\n    sum_mathbfk\n    Delta_ij^mathbfk\n    t_ij^mathbfk ^* \nendalignHere however, Γ and Δ need to be consistent with ρ and t."
},

{
    "location": "basics/topology/#",
    "page": "Topology",
    "title": "Topology",
    "category": "page",
    "text": ""
},

{
    "location": "basics/topology/#Topology-1",
    "page": "Topology",
    "title": "Topology",
    "category": "section",
    "text": "Topology blah."
},

{
    "location": "basics/topology/#Chern-Number-1",
    "page": "Topology",
    "title": "Chern Number",
    "category": "section",
    "text": "[Kohmoto][Kohmoto85] Haldane [Fukui, Hatsugai, and Suzuki][Fukui05]"
},

{
    "location": "basics/topology/#Z2-Topological-Index-of-a-Time-Reversal-Invariant-System-1",
    "page": "Topology",
    "title": "Z2 Topological Index of a Time-Reversal-Invariant System",
    "category": "section",
    "text": "[Fu and Kane][Fu06]. [Fukui and Hatsugai][Fukui07].[Kohmoto85]: https://doi.org/10.1016/0003-4916%2885%2990148-4 \"Mahito Kohmoto, Topological invariant and the quantization of the Hall conductance, Ann. Phys. 160, 343 (1985)\"[Fukui05]: https://doi.org/10.1143/JPSJ.74.1674 \"Takahiro Fukui, Yasuhiro Hatsugai, and Hiroshi Suzuki, Chern Numbers in Discretized Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances, J. Phys. Soc. Jpn. 74, 1674 (2005).\"[Fu06]: https://doi.org/10.1103/PhysRevB.74.195312 \"Liang Fu and C. L. KAne, Time reversal polarization and a Z_2 adiabatic spin pump, Phys. Rev. B 74, 195312 (2006).\"[Fukui07]: https://doi.org/10.1143/JPSJ.76.053702 \"Takahiro Fukui and Yasuhiro Hatsugai, Quantum Spin Hall Effect in Three Dimensional Materials: Lattice Computation of Z2 Topological Invariants and Its Application to Bi and Sb, J. Phys. Soc. Jpn. 76, 053702 (2007).\""
},

{
    "location": "guide/example/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": ""
},

{
    "location": "guide/hamiltonian/#",
    "page": "Hamiltonian: Full Interacting Hamiltonian",
    "title": "Hamiltonian: Full Interacting Hamiltonian",
    "category": "page",
    "text": ""
},

{
    "location": "guide/hamiltonian/#Hamiltonian:-Full-Interacting-Hamiltonian-1",
    "page": "Hamiltonian: Full Interacting Hamiltonian",
    "title": "Hamiltonian: Full Interacting Hamiltonian",
    "category": "section",
    "text": ""
},

{
    "location": "guide/hamiltonian/#Terms-of-Hamiltonian-1",
    "page": "Hamiltonian: Full Interacting Hamiltonian",
    "title": "Terms of Hamiltonian",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl currently supports Hamiltonian with hopping and (quartic) interaction terms. The hopping terms can be grouped into \"diagonal\" terms and \"offdiagonal\" terms:beginalign\nH_textdiagonal-hopping = sum_i t_ii c_i^* c_i \nH_textoffdiagonal-hopping = sum_i neq j t_ij c_i^* c_j\nendalignThe hermiticity of the Hamiltonian requires that t_ii be real, and t_ij = t_ji^*. Thus we can rewrite the offdiagonal term asH_textoffdiagonal-hopping = sum_i lt j t_ij c_i^* c_j + t_ij^* c_j^* c_iTo incorporate this we define two structs: HoppingDiagonal and HoppingOffdiagonal:struct HoppingDiagonal{R<:Real}\n    amplitude ::R\n    i ::Int\n    Ri ::Vector{Int}\nend\n\nstruct HoppingOffdiagonal{C<:Number}\n    amplitude ::C\n    i ::Int\n    j ::Int\n    Ri ::Vector{Int}\n    Rj ::Vector{Int}\nendi and j are integers representing the \"index\" of the orbital, and Ri and Rj represents which unitcell they are in. The \"unitcell coordinates\" Ri and Rj are necessary for the representation of the Hamiltonian in the momentum space. A single HoppingOffdiagonal represents both t_ij c_i^* c_j and its hermitian conjugate.Similarly for the interaction, two structs are defined:struct InteractionDiagonal{R<:Real}\n    amplitude ::R\n    i ::Int\n    j ::Int\n    Ri ::Vector{Int}\n    Rj ::Vector{Int}\nend\n\nstruct InteractionOffdiagonal{C<:Number}\n    amplitude ::C\n    i ::Int\n    j ::Int\n    k ::Int\n    l ::Int\n    Ri ::Vector{Int}\n    Rj ::Vector{Int}\n    Rk ::Vector{Int}\n    Rl ::Vector{Int}\nendwhich representU c_i^* c_j^* c_j c_i qquad i lt jandU     c_i^* c_j^* c_l c_k +\nU^* c_k^* c_l^* c_j c_i qquad\ni lt j wedge k lt l wedge  i lt k vee ( i = k wedge j lt l )In other words, (i,j) and (k,l) each need to be in lexicographical order, as well as (ij)  (kl).The following unions are definedconst Hopping = Union{HoppingDiagonal, HoppingOffdiagonal}\nconst Interaction = Union{InteractionDiagonal, InteractionOffdiagonal}"
},

{
    "location": "guide/hamiltonian/#Functions-Creating-Terms-1",
    "page": "Hamiltonian: Full Interacting Hamiltonian",
    "title": "Functions Creating Terms",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl provides functions that allows generation of hopping terms and interaction terms from the orbital names and their CarteCoord: hoppingbycarte and interactionbycarte."
},

{
    "location": "guide/hamiltonian/#FullHamiltonian-1",
    "page": "Hamiltonian: Full Interacting Hamiltonian",
    "title": "FullHamiltonian",
    "category": "section",
    "text": "The full interacting Hamiltonian is represented by the FullHamiltonian classmutable struct FullHamiltonian{O}\n    unitcell ::UnitCell{O}\n    hoppings ::Vector{Hopping}\n    interactions ::Vector{Interaction}\nendaddhopping! and addinteraction!."
},

{
    "location": "guide/hfb/#",
    "page": "HFB: Hartree-Fock-Bogoliubov",
    "title": "HFB: Hartree-Fock-Bogoliubov",
    "category": "page",
    "text": ""
},

{
    "location": "guide/hfb/#HFB:-Hartree-Fock-Bogoliubov-1",
    "page": "HFB: Hartree-Fock-Bogoliubov",
    "title": "HFB: Hartree-Fock-Bogoliubov",
    "category": "section",
    "text": ""
},

{
    "location": "guide/hfb/#HFBComputer-1",
    "page": "HFB: Hartree-Fock-Bogoliubov",
    "title": "HFBComputer",
    "category": "section",
    "text": "HFBComputer is a struct which is used to calculate the Hartree-Fock-Bogoliubov Hamiltonian from mean-field parameters, and the parameters from the eigenstates of Hartree-Fock-Bogoliubov Hamiltonian.mutable struct HFBComputer{O}\n    unitcell ::UnitCell{O}\n    hoppings ::Vector{Spec.Hopping}\n    temperature ::Float64\n\n    fermi ::Function\n    ρ_registry ::Vector{CollectRow}\n    t_registry ::Vector{CollectRow}\n    Γ_registry ::Vector{DeployRow}\n    Δ_registry ::Vector{DeployRow}\nendThe type CollectRow keeps record of the \"observables\" that need to be measured in order to calculate the mean-field parameters. They are the ρ(i,j) and t(i,j) as defined by Goodman.const CollectRow = Tuple{Bool, Int, Int, Vector{Float64}}DeployRow contains necessary information to calculate Γ(i,j) and Δ(i,j), using the measured ρ(i,j) and t(i,j)const DeployRow = Tuple{Bool, Int, Int, Vector{Float64}, Vector{Tuple{Int, ComplexF64, Bool}}}"
},

{
    "location": "guide/hfb/#HFBSolver-1",
    "page": "HFB: Hartree-Fock-Bogoliubov",
    "title": "HFBSolver",
    "category": "section",
    "text": "HFBSolver is blah.mutable struct HFBSolver{O}\n    # Originals\n    hamiltonian ::FullHamiltonian{O}\n    size ::Vector{Int}\n    temperature ::Float64\n\n    # Derivatives\n    momentumgrid ::Array{Vector{Float64}}\n    hfbhamiltonian ::HFBHamiltonian{O}\n    hfbcomputer ::HFBComputer{O}\n    greencollector ::Function\nendHFBSolver contains:The full interacting Hamiltonian,\nSize of the system (i.e. number of k-points in the Brillouin zone)\nTemperature (needed to calculate the ρ and t)And using these properties, BLAH."
},

{
    "location": "guide/hfb/#Using-HFBSolver-to-Find-Self-Consistent-Solution-1",
    "page": "HFB: Hartree-Fock-Bogoliubov",
    "title": "Using HFBSolver to Find Self-Consistent Solution",
    "category": "section",
    "text": "loop(solver, solution, 100)"
},

{
    "location": "guide/lattice/#",
    "page": "Lattice",
    "title": "Lattice",
    "category": "page",
    "text": ""
},

{
    "location": "guide/lattice/#Lattice-1",
    "page": "Lattice",
    "title": "Lattice",
    "category": "section",
    "text": ""
},

{
    "location": "guide/lattice/#CarteCoord-and-FractCoord-1",
    "page": "Lattice",
    "title": "CarteCoord and FractCoord",
    "category": "section",
    "text": "CarteCoord is simply an alias for Vector{Float64}, which represents the coordinates in the real-space. FractCoord on the other hand, represents a location in the units of the lattice vectors. For example, fc in the following code represents a location in two-dimensional space whose CarteCoord is 21 mathbfa_1 + 32 mathbfa_2 where mathbfa_1 and mathbfa_2 are lattice vectors.fc = FractCoord([2, 3], [0.1, 0.2])FractCoord can also be created in the following waysfc2 = FractCoord([2.3, 3.2])\nfc3 = FractCoord(2)fc2 is the same as fc, and fc3 represents origin in two dimensions.FractCoord supports the following arithmetic operationsfc = FractCoord([2, 3], [0.1, 0.2])\n\nfc + fc\nfc - fc\nfc + [1, 0]\nfc - [1, 0]Conversion between CarteCoord and FractCoord:latticevectors = [2.0 0.0; 0.0 1.0]\n\nfc = FractCoord([2, 3], [0.1, 0.2])\ncc = fract2carte(latticevectors, fc)\nfc2 = carte2fract(latticevectors, cc)cc is [4.2, 3.2] and fc ≈ fc2 is true."
},

{
    "location": "guide/lattice/#UnitCell-1",
    "page": "Lattice",
    "title": "UnitCell",
    "category": "section",
    "text": "A UnitCell is defined by lattice vectors and orbitals (orbitals here contain all possible degrees of freedom, such as spin, physical \"orbital\", sublattice, etc.).unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=Tuple{Symbol, Symbol})\naddorbital!(unitcell, (:A, :UP), FractCoord([0,0], [0.0, 0.0]))\naddorbital!(unitcell, (:A, :DN), FractCoord([0,0], [0.0, 0.0]))"
},

{
    "location": "guide/linearizedgap/#",
    "page": "LinearizedGap",
    "title": "LinearizedGap",
    "category": "page",
    "text": ""
},

{
    "location": "guide/linearizedgap/#LinearizedGap-1",
    "page": "LinearizedGap",
    "title": "LinearizedGap",
    "category": "section",
    "text": "Linearized gap equation. Refer to basics."
},

{
    "location": "guide/topology/#",
    "page": "Topology",
    "title": "Topology",
    "category": "page",
    "text": ""
},

{
    "location": "guide/topology/#Topology-1",
    "page": "Topology",
    "title": "Topology",
    "category": "section",
    "text": ""
},

{
    "location": "guide/topology/#Chern-1",
    "page": "Topology",
    "title": "Chern",
    "category": "section",
    "text": ""
},

{
    "location": "guide/topology/#Z2-1",
    "page": "Topology",
    "title": "Z2",
    "category": "section",
    "text": "Z2 index of time-reversal invariant topological insulator/superconductor.z2index(...)"
},

{
    "location": "#",
    "page": "HartreeFockBogoliubov.jl Documentation",
    "title": "HartreeFockBogoliubov.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#HartreeFockBogoliubov.jl-Documentation-1",
    "page": "HartreeFockBogoliubov.jl Documentation",
    "title": "HartreeFockBogoliubov.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "HartreeFockBogoliubov.jl Documentation",
    "title": "Overview",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl is a Julia package for solving interacting fermion Hamiltonian using Hartree-Fock-Bogoliubov (HFB) approach. This project aims to assist fast development of HFB solver for a generic Hamiltonian of the formmathcalH = sum_alphabeta c_alpha^dagger t_alphabeta c_beta\n+ frac12 sum_alphabetagammadelta V_alphabetagammadelta c_alpha^dagger c_beta^dagger c_delta c_gamma"
},

{
    "location": "#Installation-1",
    "page": "HartreeFockBogoliubov.jl Documentation",
    "title": "Installation",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl is not a registered package. Install it byHartreeFockBogoliubov.jl is currently not registered in the Julia Package Listing. To install, typeadd https://github.com/kyungminlee/HartreeFockBogoliubov.jlto Julia\'s package manager."
},

{
    "location": "#Usage-1",
    "page": "HartreeFockBogoliubov.jl Documentation",
    "title": "Usage",
    "category": "section",
    "text": "To construct HFB solverSpecify a unitcell, i.e. the lattice constants.\nSpecify the \"orbitals\" of the unitcell, and their locations. Here the word \"orbital\" refers to any fermionic degrees of freedom, and includes sites, orbitals in the conventional sense, and spins, etc.\nSpecify the hoppings between the orbitals, together with their real-space displacements. The displacement is required since, with the periodic boundary condition, hopping from orbital α to orbital β within the same unitcell is different from the hopping from orbital α to β in different unitcells.\nSpecify the interactions. As for the hoppings, you need to specify the orbitals and the displacement of the interactions.\nThe full interacting Hamiltonian constructed so far can be passed to a HFBSolver which decomposes it into different HFB channels."
},

{
    "location": "#License-1",
    "page": "HartreeFockBogoliubov.jl Documentation",
    "title": "License",
    "category": "section",
    "text": "MIT License\n\nCopyright (c) 2016 Kyungmin Lee\n\nPermission is hereby granted, free of charge, to any person obtaining a copy\nof this software and associated documentation files (the \"Software\"), to deal\nin the Software without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE."
},

{
    "location": "internals/#",
    "page": "Internal Documentation",
    "title": "Internal Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "internals/#Internal-Documentation-1",
    "page": "Internal Documentation",
    "title": "Internal Documentation",
    "category": "section",
    "text": "This page lists all the documented internals of the HartreeFockBogoliubov module and submodules."
},

{
    "location": "internals/#Contents-1",
    "page": "Internal Documentation",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]"
},

{
    "location": "internals/#Index-1",
    "page": "Internal Documentation",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module.Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]"
},

{
    "location": "internals/hartreefockbogoliubov/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals/hartreefockbogoliubov/#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": ""
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.CarteCoord",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.CarteCoord",
    "category": "type",
    "text": "CarteCoord\n\nCartesian coordinates. Vector{Float64}.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.FractCoord",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.FractCoord",
    "category": "type",
    "text": "FractCoord\n\nFractional coordinates.\n\nMembers\n\nwhole ::Vector{Int}: Integer part of fractional coordinates\nfraction ::Vector{Float64}: [0,1) part of fractional coordinates\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.UnitCell",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.UnitCell",
    "category": "type",
    "text": "UnitCell{O}\n\nParameters\n\nO: type of \"orbital\". Any type can be used, but we recommend using String or tuple of String and Int      for compatibility with JSON.\n\nMembers\n\nlatticevectors ::Array{Float64, 2}: Lattice vectors\nreducedreciprocallatticevectors ::Array{Float64, 2}: Reduced reciprocal lattice vectors (transpose of inverse of latticevectors)\nreciprocallatticevectors ::Array{Float64, 2}: Reciprocal lattice vectors. 2π * reducedreciprocallatticevectors\norbitals ::Vector{Tuple{T, FractCoord}}: List of orbitals within unit cell\norbitalindices ::Dict{T, Int}: Indices of orbitals\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.addorbital!-Union{Tuple{O}, Tuple{UnitCell{O},O,FractCoord}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.addorbital!",
    "category": "method",
    "text": "addorbital!\n\nAdd an orbital to the unit cell.\n\nArguments\n\nuc ::UnitCell{T}\norbitalname ::{T}\norbitalcoord ::FractCoord\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.carte2fract-Tuple{AbstractArray{#s23,2} where #s23<:AbstractFloat,Array{Float64,1}}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.carte2fract",
    "category": "method",
    "text": "carte2fract\n\nArguments\n\nlatticevectors ::AbstractArray{<:AbstractFloat, 2}: square matrix whose columns are lattice vectors.\ncc ::CarteCoord: cartesian coordinates\ntol ::Real=sqrt(eps(Float64)): tolerance\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.carte2fract-Tuple{UnitCell,Array{Float64,1}}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.carte2fract",
    "category": "method",
    "text": "carte2fract\n\nArguments\n\nlatticevectors ::Array{Float64, 2}\ncc ::CarteCoord\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.dimension-Tuple{FractCoord}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.dimension",
    "category": "method",
    "text": "dimension\n\nDimension of the fractional coordinates\n\nArguments\n\nfc ::FractCoord: Fractional coordinates.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.dimension-Tuple{UnitCell}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.dimension",
    "category": "method",
    "text": "dimension\n\nSpatial dimension of the unit cell.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.fract2carte-Tuple{AbstractArray{#s31,2} where #s31<:AbstractFloat,FractCoord}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.fract2carte",
    "category": "method",
    "text": "fract2carte\n\nArguments\n\nlatticevectors ::AbstractArray{<:AbstractFloat, 2}: square matrix whose columns are lattice vectors.\nfc ::FractCoord: fractional coordinates\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.fract2carte-Tuple{UnitCell,FractCoord}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.fract2carte",
    "category": "method",
    "text": "fract2carte\n\nArguments\n\nlatticevectors ::Array{Float64, 2}\nfc ::FractCoord\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbital-Tuple{UnitCell,Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbital",
    "category": "method",
    "text": "getorbital\n\nArguments\n\nuc ::UnitCell{T}\nindex ::Integer\n\nReturn\n\n(orbitalname, fractcoord)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbital-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbital",
    "category": "method",
    "text": "getorbital\n\nGet the orbital (its orbital name and its fractional coordinates) with the given name.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\nReturn\n\n(orbitalname, fractcoord)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbitalcoord-Tuple{UnitCell,Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbitalcoord",
    "category": "method",
    "text": "getorbitalcoord\n\nArguments\n\nuc ::UnitCell\nidx ::Integer\n\nReturn\n\nfractcoord\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbitalcoord-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbitalcoord",
    "category": "method",
    "text": "getorbitalcoord\n\nGet the fractional coordinates of the orbital with the given name.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\nReturn\n\nfractcoord\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbitalindex-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbitalindex",
    "category": "method",
    "text": "getorbitalindex\n\nGet index of the given orbital.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbitalindexcoord-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbitalindexcoord",
    "category": "method",
    "text": "getorbitalindexcoord\n\nArguments\n\nuc ::UnitCell{T}\nname ::T\n\nReturn\n\n(index, fractcoord)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.getorbitalname-Tuple{UnitCell,Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.getorbitalname",
    "category": "method",
    "text": "getorbitalname\n\nArguments\n\nuc ::UnitCell\nindex ::Integer\n\nReturn\n\norbitalname\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.hasorbital-Union{Tuple{O}, Tuple{UnitCell{O},O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.hasorbital",
    "category": "method",
    "text": "hasorbital{T}\n\nTest whether the unit cell contains the orbital of given name.\n\nArguments\n\nuc ::UnitCell{O}\nname ::O\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.make_unitcell-Tuple{AbstractArray{#s28,2} where #s28<:AbstractFloat}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.make_unitcell",
    "category": "method",
    "text": "UnitCell\n\nConstruct an n-dimensional lattice.\n\nArguments\n\nlatticevectors ::AbstractArray{<:AbstractFloat, 2}: Lattice vectors\nOrbitalType::DataType\n\nOptional Arguments\n\ntol=sqrt(eps(Float64)): Epsilon\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.make_unitcell-Tuple{AbstractFloat}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.make_unitcell",
    "category": "method",
    "text": "UnitCell\n\nConstruct a one-dimensional lattice.\n\nArguments\n\nlatticeconstant ::Float64: Lattice constant\nOrbitalType: List of orbitals\n\nOptional Arguments\n\ntol=sqrt(eps(Float64)): Tolerance\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.momentumgrid-Tuple{UnitCell,AbstractArray{#s33,1} where #s33<:Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.momentumgrid",
    "category": "method",
    "text": "momentumgrid\n\nGenerate an n-dimensional grid of momenta of given shape\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.numorbital-Tuple{UnitCell}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.numorbital",
    "category": "method",
    "text": "numorbital\n\nNumber of orbitals of the unit cell.\n\nArguments\n\nuc ::UnitCell\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.whichunitcell-Union{Tuple{O}, Tuple{UnitCell{O},O,Array{Float64,1}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.whichunitcell",
    "category": "method",
    "text": "whichunitcell\n\nReturn\n\nR ::Vector{Int}: which unit cell the specificied orbital/cartesian coordinates belongs to.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Lattice.whichunitcell-Union{Tuple{O}, Tuple{UnitCell{O},O,FractCoord}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Lattice.whichunitcell",
    "category": "method",
    "text": "whichunitcell\n\nReturn\n\nR ::Vector{Int}: which unit cell the specificied orbital/cartesian coordinates belongs to.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#Base.isapprox-Tuple{FractCoord,FractCoord}",
    "page": "Internals",
    "title": "Base.isapprox",
    "category": "method",
    "text": "isapprox(x, y; rtol::Real=atol>0 ? 0 : √eps, atol::Real=0, nans::Bool=false, norm::Function)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#Lattice-1",
    "page": "Internals",
    "title": "Lattice",
    "category": "section",
    "text": "Modules = [HartreeFockBogoliubov.Lattice]"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec",
    "category": "module",
    "text": "Submodule `Spec`\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.Hopping",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.Hopping",
    "category": "constant",
    "text": "Hopping\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.Interaction",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.Interaction",
    "category": "constant",
    "text": "Interaction\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.FullHamiltonian",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.FullHamiltonian",
    "category": "type",
    "text": "FullHamiltonian\n\nMembers\n\nunitcell ::UnitCell\nhoppings ::Vector{Hopping}\ninteractions ::Vector{Interaction}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.FullHamiltonian-Union{Tuple{UnitCell{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.FullHamiltonian",
    "category": "method",
    "text": "Hamiltonian\n\nCreate an empty Hamiltonian\n\nArguments\n\nunitcell ::UnitCell\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.HoppingDiagonal",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.HoppingDiagonal",
    "category": "type",
    "text": "HoppingDiagonal{R<:Real}\n\nRepresents\n\n  t c_i^* c_i\n\nMembers\n\namplitude ::R\ni ::Int: index of orbital\nRi ::Vector{Int}: which unit cell? (indexed by a1, and a2)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.HoppingDiagonal-Union{Tuple{R}, Tuple{R,Integer,AbstractArray{#s33,1} where #s33<:Integer}} where R<:Real",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.HoppingDiagonal",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.HoppingOffdiagonal",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.HoppingOffdiagonal",
    "category": "type",
    "text": "HoppingOffdiagonal{C<:Number}\n\nRepresents\n\n  t c_i^* c_j + t^* c_j^* c_i\n\nt, i, j and the unitcell-coordinates Ri and Rj are stored. Require that (i, Ri) <= (j, Rj)\n\nMembers\n\namplitude :: C\ni, j ::T\nRi, Rj ::Vector{Int}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.InteractionDiagonal",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.InteractionDiagonal",
    "category": "type",
    "text": "InteractionDiagonal{R<:Real}\n\nRepresents\n\n    U c_i^* c_j^* c_j c_i\n\nMembers\n\namplitude ::R\ni, j ::Int\nRi, Rj ::Vector{Int}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.InteractionOffdiagonal",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.InteractionOffdiagonal",
    "category": "type",
    "text": "InteractionOffdiagonal{C<:Number}\n\ni < j, k < l, i < k or (i == k and j < l)\n\nRepresents\n\n   U     c_i^* c_j^* c_l c_k\n + U^* c_k^* c_l^* c_j c_i\n\nOnly keep the first term (and require i < j, k < l, i <= k)\n\nOrdering of orbitals by (i, Ri)\n\nMembers\n\namplitude ::C\ni, j, k, l ::Int\nRi, Rj, Rk, Rl ::Vector{Int}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.addhopping!-Union{Tuple{C}, Tuple{O}, Tuple{FullHamiltonian{O},HoppingOffdiagonal{C}}} where C<:Number where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.addhopping!",
    "category": "method",
    "text": "addhopping!\n\nArguments\n\nhamiltonian ::Hamiltonian\nhopping ::HoppingOffdiagonal\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.addhopping!-Union{Tuple{R}, Tuple{O}, Tuple{FullHamiltonian{O},HoppingDiagonal{R}}} where R<:Real where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.addhopping!",
    "category": "method",
    "text": "addhopping!\n\nArguments\n\nhamiltonian ::Hamiltonian\nhopping ::HoppingDiagonal\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.addinteraction!-Union{Tuple{C}, Tuple{O}, Tuple{FullHamiltonian{O},InteractionOffdiagonal{C}}} where C<:Number where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.addinteraction!",
    "category": "method",
    "text": "addinteraction!\n\nArguments\n\nhamiltonian ::Hamiltonian\ninteraction ::InteractionOffdiagonal\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.addinteraction!-Union{Tuple{R}, Tuple{O}, Tuple{FullHamiltonian{O},InteractionDiagonal{R}}} where R<:Real where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.addinteraction!",
    "category": "method",
    "text": "addinteraction!\n\nArguments\n\nhamiltonian ::Hamiltonian\ninteraction ::InteractionDiagonal\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.hoppingbycarte-Union{Tuple{C}, Tuple{O}, Tuple{UnitCell{O},C,O,O,Array{Float64,1},Array{Float64,1}}} where C<:Number where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.hoppingbycarte",
    "category": "method",
    "text": "hoppingbycarte{T}\n\nMake a hopping element with cartesian coordinates.\n\nArguments\n\nuc ::UnitCell{T}\namplitude ::Number\ni ::T\nj ::T\nri ::CarteCoord\nrj ::CarteCoord\ntol ::Real : Optional. Defaults to sqrt(eps(Float64))\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.hoppingbycarte-Union{Tuple{R}, Tuple{O}, Tuple{UnitCell{O},R,O,Array{Float64,1}}} where R<:Real where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.hoppingbycarte",
    "category": "method",
    "text": "hoppingbycarte{T}\n\nMake a hopping element with cartesian coordinates.\n\nArguments\n\nuc ::UnitCell{T}\namplitude ::Real\ni ::T\nri ::CarteCoord\ntol ::Real : Optional. Defaults to sqrt(eps(Float64))\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.interactionbycarte-Union{Tuple{C}, Tuple{O}, Tuple{UnitCell{O},C,O,O,O,O,Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}} where C<:Number where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.interactionbycarte",
    "category": "method",
    "text": "interactionbycarte{T}\n\nMake an interaction element with cartesian coordinates.\n\nArguments\n\n* `uc ::UnitCell{T}`\n* `amplitude ::Number`\n* `i ::T`\n* `j ::T`\n* `k ::T`\n* `l ::T`\n* `ri ::CarteCoord`\n* `rj ::CarteCoord`\n* `rk ::CarteCoord`\n* `rl ::CarteCoord`\n* `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.interactionbycarte-Union{Tuple{R}, Tuple{O}, Tuple{UnitCell{O},R,O,O,Array{Float64,1},Array{Float64,1}}} where R<:Real where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.interactionbycarte",
    "category": "method",
    "text": "interactionbycarte{T}\n\nMake an interaction element with cartesian coordinates.\n\nArguments\n\n* `uc ::UnitCell{T}`\n* `amplitude ::Number`\n* `i ::T`\n* `j ::T`\n* `ri ::CarteCoord`\n* `rj ::CarteCoord`\n* `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.islocal-Union{Tuple{HoppingDiagonal{R}}, Tuple{R}} where R<:Real",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.islocal",
    "category": "method",
    "text": "islocal\n\nCheck if the hopping element is local (i.e. Ri is zero)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Spec.localized-Union{Tuple{HoppingDiagonal{R}}, Tuple{R}} where R<:Real",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Spec.localized",
    "category": "method",
    "text": "localized\n\nReturn a hopping element that is local (i.e. Ri is zero)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#Spec-1",
    "page": "Internals",
    "title": "Spec",
    "category": "section",
    "text": "Modules = [HartreeFockBogoliubov.Spec]"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator",
    "category": "module",
    "text": "Generator submodule\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{C}, Tuple{O}, Tuple{UnitCell{O},HoppingOffdiagonal{C}}} where C<:Number where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{FullHamiltonian{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{HoppingDiagonal,1},AbstractArray{HoppingOffdiagonal,1}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{HoppingDiagonal,1}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{HoppingOffdiagonal,1}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{Union{HoppingDiagonal, HoppingOffdiagonal},1}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Generator.hopping_inplace-Union{Tuple{R}, Tuple{O}, Tuple{UnitCell{O},HoppingDiagonal{R}}} where R<:Real where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Generator.hopping_inplace",
    "category": "method",
    "text": "hopping_inplace\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#Generator-1",
    "page": "Internals",
    "title": "Generator",
    "category": "section",
    "text": "Modules = [HartreeFockBogoliubov.Generator]"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.isvalidtimereversalmatrix-Tuple{AbstractArray{#s69,2} where #s69<:Number}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.isvalidtimereversalmatrix",
    "category": "method",
    "text": "isvalidtimereversalmatrix\n\nTest whether the given matrix is a valid unitary matrix for the time reversal operation.\n\n```math\nT = U ⋅ K\n```\n\n``U`` must satisfy the two conditions:\n1. ``U U^{\\dagger} = 1`` (from unitarity of ``U``)\n2. ``U = - U^{\\mathsf{T}}`` (from `T^2 = -1`)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.squarify-Union{Tuple{FullHamiltonian{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.squarify",
    "category": "method",
    "text": "squarify\n\nArguments\n\nuc::Spec.FullHamiltonian{O}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.squarify-Union{Tuple{HFBHamiltonian{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.squarify",
    "category": "method",
    "text": "squarify\n\nArguments\n\nuc::HFB.HFBHamiltonian{O}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.squarify-Union{Tuple{UnitCell{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.squarify",
    "category": "method",
    "text": "squarify\n\nIn order to make the Hamiltoinian a periodic function of momentum, bring all the sites to the origin. In addition, make the unitcell into a square.\n\nArguments\n\nuc::Lattice.UnitCell{O}\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.timereversalindexgrid-Tuple{Integer,Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.timereversalindexgrid",
    "category": "method",
    "text": "generate k-space grid (which has (2 n1, 2 n2) points TOTAL in the Brillouin zone)\n\nExample\n\nWhen n1 = 4, n2 = 3, this function returns an OrderedDict that represents the following structure\n\ni2|\n  |\n5 | -i -i -i -i -i -i -i -i\n4 | -i -i -i -i -i -i -i -i\n3 | 0h +h +h +h 0h -h -h -h\n2 | +i +i +i +i +i +i +i +i\n1 | +i +i +i +i +i +i +i +i\n0 | 0z +z +z +z 0z -z -z -z\n--+----------------------------\n  |  0  1  2  3  4  5  6  7  i1\n\nwhere 0z, +z, -z are represented respectively by :TRIZERO, :POSZERO, and :NEGZERO, and   0h, +h, -h by :TRIHALF, :POSHALF, and :NEGHALF, and   +i, -i by :POSINT, :NEGINT.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.chernnumber-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{Union{HoppingDiagonal, HoppingOffdiagonal},1},Integer,Integer,AbstractArray{#s30,1} where #s30<:Integer}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.chernnumber",
    "category": "method",
    "text": "chernnumber\n\nCompute chern number of the band structure defined by the hoppings and the selected bands.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.getnambuphase-Union{Tuple{O}, Tuple{UnitCell{O},Function,AbstractArray{#s33,2} where #s33<:Number,Integer,Integer}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.getnambuphase",
    "category": "method",
    "text": "Get the phase of hamiltonian function.     Returns nan if no phase\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.Topology.z2index-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{Union{HoppingDiagonal, HoppingOffdiagonal},1},AbstractArray{T,2} where T,Integer,Integer,AbstractArray{#s16,1} where #s16<:Integer}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.Topology.z2index",
    "category": "method",
    "text": "z2index\n\nCompute Z2 index of time-reversal-invariant Hamiltonian.\n\nArguments\n\nuc::UnitCell{O}\nhops::AbstractVector{Hopping}\ntimereversal::AbstractMatrix\nn1 ::Integer\nn2 ::Integer\nselectpairs::AbstractVector{<:Integer}\n\nOptional Arguments\n\ntol ::Real = sqrt(eps(Float64))\n\nReturns\n\n(The Z2 index,  max| Hₖ - T⁻¹HₖT | for k in TRIMs)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#Topology-1",
    "page": "Internals",
    "title": "Topology",
    "category": "section",
    "text": "Modules = [HartreeFockBogoliubov.Topology]"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBComputer",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBComputer",
    "category": "type",
    "text": "HFBConmputer is a type holding the ρ, t and Γ, Δ of a Hartree-Fock-Bogoliubov Hamiltonian.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBComputer-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{HoppingDiagonal,1},AbstractArray{HoppingOffdiagonal,1},Real}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBComputer",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBHamiltonian",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBHamiltonian",
    "category": "type",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBHamiltonian-Union{Tuple{FullHamiltonian{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBHamiltonian",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBHamiltonian-Union{Tuple{O}, Tuple{UnitCell{O},AbstractArray{HoppingDiagonal,1},AbstractArray{HoppingOffdiagonal,1},AbstractArray{HoppingMeanField,1},AbstractArray{PairingMeanField,1}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBHamiltonian",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBSolver",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBSolver",
    "category": "type",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HFBSolver-Union{Tuple{O}, Tuple{FullHamiltonian{O},AbstractArray{#s70,1} where #s70<:Integer,Real}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HFBSolver",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.addinteraction!-Union{Tuple{C}, Tuple{O}, Tuple{HFBHamiltonian{O},InteractionOffdiagonal{C}}} where C<:Number where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.addinteraction!",
    "category": "method",
    "text": "Add offdiagonal interaction\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.addinteraction!-Union{Tuple{R}, Tuple{O}, Tuple{HFBHamiltonian{O},InteractionDiagonal{R}}} where R<:Real where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.addinteraction!",
    "category": "method",
    "text": "addinteraction!\n\nAdd diagonal interaction\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.compute_hfbfield!-Tuple{HartreeFockBogoliubov.HFB.HFBField,HartreeFockBogoliubov.HFB.HFBComputer,HartreeFockBogoliubov.HFB.HFBAmplitude}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.compute_hfbfield!",
    "category": "method",
    "text": "Compute Γ and Δ from ρ and t.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.isvalid-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,HartreeFockBogoliubov.HFB.HFBAmplitude}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.isvalid",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.isvalid-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,HartreeFockBogoliubov.HFB.HFBField}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.isvalid",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.loop-Tuple{HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBField,Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.loop",
    "category": "method",
    "text": "loop\n\nPerform selfconsistency loop a number of times with the given precondition and given update functions.\n\nArguments\n\nsolver ::HFBSolver{T}\nsol::HFBAmplitude\nrun::Integer\n\nOptional Arguments\n\nupdate::Function=simpleupdate\nprecondition::Function=identity: on target field\ncallback::Function=_noop: Function called after every update as\n\ncallback(i, run)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.loop_threaded-Tuple{HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBAmplitude,Integer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.loop_threaded",
    "category": "method",
    "text": "loop\n\nPerform selfconsistency loop a number of times with the given precondition and given update functions.\n\nArguments\n\nsolver ::HFBSolver{T}\nsf::HFBAmplitude\nrun::Integer\n\nOptional Arguments\n\nupdate::Function=simpleupdate\nprecondition::Function=identity\ncallback::Function=_noop: Function called after every update as\n\ncallback(i, run)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_Deltamatrix-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,AbstractArray{#s33,1} where #s33<:Number}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_Deltamatrix",
    "category": "method",
    "text": "Return a generator of Δ matrix (which is a function of momentum)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_Gammamatrix-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,AbstractArray{#s33,1} where #s33<:Number}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_Gammamatrix",
    "category": "method",
    "text": "Return a generator of Γ matrix (which is a function of momentum)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_greencollector-Union{Tuple{HFBComputer{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_greencollector",
    "category": "method",
    "text": "make_greencollector\n\nReturns a function which has the following signature\n\ncollector(k, eigenvalues, eigenvectors, ρout, tout)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hamiltonian-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,AbstractArray{#s25,1} where #s25<:Number,AbstractArray{#s24,1} where #s24<:Number}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hamiltonian",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hfbamplitude-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,Function,Function}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hfbamplitude",
    "category": "method",
    "text": "func : (idx, i, j, r) -> val\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hfbamplitude-Tuple{HartreeFockBogoliubov.HFB.HFBComputer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hfbamplitude",
    "category": "method",
    "text": "func : (idx, i, j, r) -> 0\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hfbamplitude-Union{Tuple{O}, Tuple{HFBComputer{O},HFBAmplitudeHint{O}}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hfbamplitude",
    "category": "method",
    "text": "Check if hint contains ρ\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hfbfield-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,HartreeFockBogoliubov.HFB.HFBAmplitude}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hfbfield",
    "category": "method",
    "text": "Compute Γ and Δ from ρ and t.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hfbfield-Tuple{HartreeFockBogoliubov.HFB.HFBComputer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hfbfield",
    "category": "method",
    "text": "Compute Γ and Δ from ρ and t.\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hint-Union{Tuple{O}, Tuple{HFBComputer{O},HFBAmplitude}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hint",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_hoppingmatrix-Tuple{HartreeFockBogoliubov.HFB.HFBComputer}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_hoppingmatrix",
    "category": "method",
    "text": "Return a generator of hopping matrix (which is a function of momentum)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.next_hfbamplitude!-Tuple{HartreeFockBogoliubov.HFB.HFBAmplitude,HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBField}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.next_hfbamplitude!",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.next_hfbamplitude-Tuple{HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBAmplitude}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.next_hfbamplitude",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.next_hfbamplitude-Tuple{HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBField}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.next_hfbamplitude",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.next_hfbamplitude_threaded!-Tuple{HartreeFockBogoliubov.HFB.HFBAmplitude,HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBField}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.next_hfbamplitude_threaded!",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.next_hfbamplitude_threaded-Tuple{HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBField}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.next_hfbamplitude_threaded",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.randomize!-Tuple{HartreeFockBogoliubov.HFB.HFBComputer,HartreeFockBogoliubov.HFB.HFBAmplitude}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.randomize!",
    "category": "method",
    "text": "Randomize a hfbamplitude\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.CollectRow",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.CollectRow",
    "category": "type",
    "text": "CollectRow is holds info on how to compute ρ or t. Its elements are:\n\nIs diagonal? (only for rho)\nrow orbital\ncol orbital\ndisplacement r(col) - r(row)\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.DeployRow",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.DeployRow",
    "category": "type",
    "text": "DeployRow is holds info on how to compute Γ or Δ. Its elements are:\n\nIs diagonal\nrow orbital\ncol orbital\ndisplacement r(col) - r(row)\nlist of sources, each of which is a tuple of\nindex of ρ or t from which to compute this Γ or Δ.\namplitude (coefficient to multiply to ρ or t)\nboolean indicating whether   (1) conjugation is needed (for ρ/Γ) or   (2) minus sign is needed (for t/Δ).\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.HoppingMeanField",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.HoppingMeanField",
    "category": "type",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.PairingMeanField",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.PairingMeanField",
    "category": "type",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_particleholeregistry-Tuple{HartreeFockBogoliubov.HFB.HFBHamiltonian}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_particleholeregistry",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.make_particleparticleregistry-Tuple{HartreeFockBogoliubov.HFB.HFBHamiltonian}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.make_particleparticleregistry",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.HFB.next_hfbfield-Tuple{HartreeFockBogoliubov.HFB.HFBSolver,HartreeFockBogoliubov.HFB.HFBField}",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.HFB.next_hfbfield",
    "category": "method",
    "text": "\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#HFB-1",
    "page": "Internals",
    "title": "HFB",
    "category": "section",
    "text": "Modules = [HartreeFockBogoliubov.HFB]"
},

{
    "location": "internals/hartreefockbogoliubov/#HartreeFockBogoliubov.LinearizedGap.linearizedpairingkernel-Union{Tuple{HFBSolver{O}}, Tuple{O}} where O",
    "page": "Internals",
    "title": "HartreeFockBogoliubov.LinearizedGap.linearizedpairingkernel",
    "category": "method",
    "text": "linearizedpairingkernel\n\nCompute the kernel Γ of the linearized gap equation in the pairing channel which is written as Δ = Γ⋅Δ\n\n\n\n\n\n"
},

{
    "location": "internals/hartreefockbogoliubov/#LinearizedGap-1",
    "page": "Internals",
    "title": "LinearizedGap",
    "category": "section",
    "text": "Modules = [HartreeFockBogoliubov.LinearizedGap]"
},

]}
