var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#HartreeFockBogoliubov.jl-Documentation-1",
    "page": "Home",
    "title": "HartreeFockBogoliubov.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl is a Julia package for solving interacting fermion Hamiltonian using Hartree-Fock-Bogoliubov (HFB) approach. This project aims to assist fast development of HFB solver for a generic Hamiltonian of the formmathcalH = sum_alphabeta c_alpha^dagger t_alphabeta c_beta\n+ frac12 sum_alphabetagammadelta V_alphabetagammadelta c_alpha^dagger c_beta^dagger c_delta c_gamma"
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl is not a registered package. Install it byHartreeFockBogoliubov.jl is currently not registered in the Julia Package Listing. To install, typeadd https://github.com/kyungminlee/HartreeFockBogoliubov.jlto Julia\'s package manager."
},

{
    "location": "#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "To construct HFB solverSpecify a unitcell, i.e. the lattice constants.\nSpecify the \"orbitals\" of the unitcell, and their locations. Here the word \"orbital\" refers to any fermionic degrees of freedom, and includes sites, orbitals in the conventional sense, and spins, etc.\nSpecify the hoppings between the orbitals, together with their real-space displacements. The displacement is required since, with the periodic boundary condition, hopping from orbital α to orbital β within the same unitcell is different from the hopping from orbital α to β in different unitcells.\nSpecify the interactions. As for the hoppings, you need to specify the orbitals and the displacement of the interactions.\nThe full interacting Hamiltonian constructed so far can be passed to a HFBSolver which decomposes it into different HFB channels."
},

{
    "location": "#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "MIT License\n\nCopyright (c) 2016 Kyungmin Lee\n\nPermission is hereby granted, free of charge, to any person obtaining a copy\nof this software and associated documentation files (the \"Software\"), to deal\nin the Software without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE."
},

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
    "location": "basics/hfbequation/#",
    "page": "Hartree-Fock-Bogoliubov Equations",
    "title": "Hartree-Fock-Bogoliubov Equations",
    "category": "page",
    "text": ""
},

{
    "location": "basics/hfbequation/#HFB-Equations-1",
    "page": "Hartree-Fock-Bogoliubov Equations",
    "title": "HFB Equations",
    "category": "section",
    "text": "delta Omega = 0"
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
    "page": "Topological Invariants",
    "title": "Topological Invariants",
    "category": "page",
    "text": ""
},

{
    "location": "basics/topology/#Topology-1",
    "page": "Topological Invariants",
    "title": "Topology",
    "category": "section",
    "text": "Topology blah."
},

{
    "location": "basics/topology/#Chern-Number-1",
    "page": "Topological Invariants",
    "title": "Chern Number",
    "category": "section",
    "text": "[Kohmoto][Kohmoto85] Haldane [Fukui, Hatsugai, and Suzuki][Fukui05]"
},

{
    "location": "basics/topology/#Z2-Topological-Index-of-a-Time-Reversal-Invariant-System-1",
    "page": "Topological Invariants",
    "title": "Z2 Topological Index of a Time-Reversal-Invariant System",
    "category": "section",
    "text": "[Fu and Kane][Fu06]. [Fukui and Hatsugai][Fukui07].[Kohmoto85]: https://doi.org/10.1016/0003-4916%2885%2990148-4 \"Mahito Kohmoto, Topological invariant and the quantization of the Hall conductance, Ann. Phys. 160, 343 (1985)\"[Fukui05]: https://doi.org/10.1143/JPSJ.74.1674 \"Takahiro Fukui, Yasuhiro Hatsugai, and Hiroshi Suzuki, Chern Numbers in Discretized Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances, J. Phys. Soc. Jpn. 74, 1674 (2005).\"[Fu06]: https://doi.org/10.1103/PhysRevB.74.195312 \"Liang Fu and C. L. KAne, Time reversal polarization and a Z_2 adiabatic spin pump, Phys. Rev. B 74, 195312 (2006).\"[Fukui07]: https://doi.org/10.1143/JPSJ.76.053702 \"Takahiro Fukui and Yasuhiro Hatsugai, Quantum Spin Hall Effect in Three Dimensional Materials: Lattice Computation of Z2 Topological Invariants and Its Application to Bi and Sb, J. Phys. Soc. Jpn. 76, 053702 (2007).\""
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
    "location": "guide/hamiltonian/#",
    "page": "Hamiltonian",
    "title": "Hamiltonian",
    "category": "page",
    "text": ""
},

{
    "location": "guide/hamiltonian/#Hamiltonian:-Full-Interacting-Hamiltonian-1",
    "page": "Hamiltonian",
    "title": "Hamiltonian: Full Interacting Hamiltonian",
    "category": "section",
    "text": ""
},

{
    "location": "guide/hamiltonian/#Terms-of-Hamiltonian-1",
    "page": "Hamiltonian",
    "title": "Terms of Hamiltonian",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl currently supports Hamiltonian with hopping and (quartic) interaction terms. The hopping terms can be grouped into \"diagonal\" terms and \"offdiagonal\" terms:beginalign\nH_textdiagonal-hopping = sum_i t_ii c_i^* c_i \nH_textoffdiagonal-hopping = sum_i neq j t_ij c_i^* c_j\nendalignThe hermiticity of the Hamiltonian requires that t_ii be real, and t_ij = t_ji^*. Thus we can rewrite the offdiagonal term asH_textoffdiagonal-hopping = sum_i lt j t_ij c_i^* c_j + t_ij^* c_j^* c_iTo incorporate this we define two structs: HoppingDiagonal and HoppingOffdiagonal:struct HoppingDiagonal{R<:Real}\n    amplitude ::R\n    i ::Int\n    Ri ::Vector{Int}\nend\n\nstruct HoppingOffdiagonal{C<:Number}\n    amplitude ::C\n    i ::Int\n    j ::Int\n    Ri ::Vector{Int}\n    Rj ::Vector{Int}\nendi and j are integers representing the \"index\" of the orbital, and Ri and Rj represents which unitcell they are in. The \"unitcell coordinates\" Ri and Rj are necessary for the representation of the Hamiltonian in the momentum space. A single HoppingOffdiagonal represents both t_ij c_i^* c_j and its hermitian conjugate.Similarly for the interaction, two structs are defined:struct InteractionDiagonal{R<:Real}\n    amplitude ::R\n    i ::Int\n    j ::Int\n    Ri ::Vector{Int}\n    Rj ::Vector{Int}\nend\n\nstruct InteractionOffdiagonal{C<:Number}\n    amplitude ::C\n    i ::Int\n    j ::Int\n    k ::Int\n    l ::Int\n    Ri ::Vector{Int}\n    Rj ::Vector{Int}\n    Rk ::Vector{Int}\n    Rl ::Vector{Int}\nendwhich representU c_i^* c_j^* c_j c_i qquad i lt jandU     c_i^* c_j^* c_l c_k +\nU^* c_k^* c_l^* c_j c_i qquad\ni lt j wedge k lt l wedge  i lt k vee ( i = k wedge j lt l )In other words, (i,j) and (k,l) each need to be in lexicographical order, as well as (ij)  (kl).The following unions are definedconst Hopping = Union{HoppingDiagonal, HoppingOffdiagonal}\nconst Interaction = Union{InteractionDiagonal, InteractionOffdiagonal}"
},

{
    "location": "guide/hamiltonian/#Functions-Creating-Terms-1",
    "page": "Hamiltonian",
    "title": "Functions Creating Terms",
    "category": "section",
    "text": "HartreeFockBogoliubov.jl provides functions that allows generation of hopping terms and interaction terms from the orbital names and their CarteCoord: hoppingbycarte and interactionbycarte."
},

{
    "location": "guide/hamiltonian/#FullHamiltonian-1",
    "page": "Hamiltonian",
    "title": "FullHamiltonian",
    "category": "section",
    "text": "The full interacting Hamiltonian is represented by the FullHamiltonian classmutable struct FullHamiltonian{O}\n    unitcell ::UnitCell{O}\n    hoppings ::Vector{Hopping}\n    interactions ::Vector{Interaction}\nendaddhopping! and addinteraction!."
},

{
    "location": "guide/hfb/#",
    "page": "HFB",
    "title": "HFB",
    "category": "page",
    "text": ""
},

{
    "location": "guide/hfb/#HFB:-Hartree-Fock-Bogoliubov-1",
    "page": "HFB",
    "title": "HFB: Hartree-Fock-Bogoliubov",
    "category": "section",
    "text": ""
},

{
    "location": "guide/hfb/#HFBComputer-1",
    "page": "HFB",
    "title": "HFBComputer",
    "category": "section",
    "text": "HFBComputer is a struct which is used to calculate the Hartree-Fock-Bogoliubov Hamiltonian from mean-field parameters, and the parameters from the eigenstates of Hartree-Fock-Bogoliubov Hamiltonian.mutable struct HFBComputer{O}\n    unitcell ::UnitCell{O}\n    hoppings ::Vector{Spec.Hopping}\n    temperature ::Float64\n\n    fermi ::Function\n    ρ_registry ::Vector{CollectRow}\n    t_registry ::Vector{CollectRow}\n    Γ_registry ::Vector{DeployRow}\n    Δ_registry ::Vector{DeployRow}\nendThe type CollectRow keeps record of the \"observables\" that need to be measured in order to calculate the mean-field parameters. They are the ρ(i,j) and t(i,j) as defined by Goodman.const CollectRow = Tuple{Bool, Int, Int, Vector{Float64}}DeployRow contains necessary information to calculate Γ(i,j) and Δ(i,j), using the measured ρ(i,j) and t(i,j)const DeployRow = Tuple{Bool, Int, Int, Vector{Float64}, Vector{Tuple{Int, ComplexF64, Bool}}}"
},

{
    "location": "guide/hfb/#HFBSolver-1",
    "page": "HFB",
    "title": "HFBSolver",
    "category": "section",
    "text": "HFBSolver is blah.mutable struct HFBSolver{O}\n    # Originals\n    hamiltonian ::FullHamiltonian{O}\n    size ::Vector{Int}\n    temperature ::Float64\n\n    # Derivatives\n    momentumgrid ::Array{Vector{Float64}}\n    hfbhamiltonian ::HFBHamiltonian{O}\n    hfbcomputer ::HFBComputer{O}\n    greencollector ::Function\nendHFBSolver contains:The full interacting Hamiltonian,\nSize of the system (i.e. number of k-points in the Brillouin zone)\nTemperature (needed to calculate the ρ and t)And using these properties, BLAH."
},

{
    "location": "guide/hfb/#Using-HFBSolver-to-Find-Self-Consistent-Solution-1",
    "page": "HFB",
    "title": "Using HFBSolver to Find Self-Consistent Solution",
    "category": "section",
    "text": "loop(solver, solution, 100)"
},

{
    "location": "guide/linearizedgap/#",
    "page": "Linearized Gap Equation",
    "title": "Linearized Gap Equation",
    "category": "page",
    "text": ""
},

{
    "location": "guide/linearizedgap/#LinearizedGap-1",
    "page": "Linearized Gap Equation",
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
    "location": "guide/example/#",
    "page": "Example",
    "title": "Example",
    "category": "page",
    "text": ""
},

{
    "location": "internals/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals/#Internal-Documentation-1",
    "page": "Internals",
    "title": "Internal Documentation",
    "category": "section",
    "text": "This page lists all the documented internals of the HartreeFockBogoliubov module and submodules."
},

{
    "location": "internals/#Contents-1",
    "page": "Internals",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]"
},

{
    "location": "internals/#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module.Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]"
},

]}
