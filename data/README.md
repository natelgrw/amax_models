# Table Of AMAX-1 Molecular Desciptors

Each molecule - solvent combination contains two sets of calculated molecular descriptors: one for the target molecule, and one for the solvent.

## Topological Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| BalabanJ | Quantifies molecular complexity based on average distance connectivity and graph branching | RDKit |
| BertzCT | Calculates molecular complexity based on graph connectivity and atomic contributions | RDKit |
| Chi (0-4), Chi_n (1-4) Chi_v (1-4)  | Connectivity indices reflecting molecular topology, branching, and size | RDKit |
| Ipc  | Information content index representing structural complexity | RDKit |
| Kappa (1-3) | Shape indices describing molecular flexibility and overall geometry | RDKit |
| SpMax_Bhv (1-8), SpMin_Bhv (1-8) | Largest and smallest Burden matrix eigenvalues weighted by van der Waals volume | PaDEL Descriptor |
| TopoPSA | Topological polar surface area estimating polar regions accessible on the molecule | PaDEL Descriptor |
| WTPT (1-5) | Weighted path counts of lengths 1 to 5, capturing molecular connectivity | PaDEL Descriptor |
| WPATH | Weighted path descriptor summarizing overall path-based connectivity | PaDEL Descriptor |
| WPOL | Weighted polarity descriptor reflecting polar atom connectivity and distribution | PaDEL Descriptor |
| ETA_Alpha | Measures molecular size and shape using extended topochemical atom indices | PaDEL Descriptor |
| ETA_AlphaP | Modified ETA_Alpha with path corrections | PaDEL Descriptor |
| ETA_dAlpha_A | Differential ETA Alpha descriptor capturing subtle size variations | PaDEL Descriptor |
| ETA_dAlpha_B | Alternative differential ETA Alpha descriptor | PaDEL Descriptor |
| ETA_Epsilon (1-5) | Electronic effect descriptor reflecting polarizability | PaDEL Descriptor |
| ETA_dEpsilon_A | Differential electronic effect descriptor | PaDEL Descriptor |
| ETA_dEpsilon_B | Alternative differential electronic effect descriptor | PaDEL Descriptor |
| ETA_dEpsilon_C | Another variant of differential electronic effect descriptor | PaDEL Descriptor |
| ETA_dEpsilon_D | Additional differential electronic effect descriptor | PaDEL Descriptor |
| ETA_Psi_1 | Molecular shape descriptor based on ETA framework | PaDEL Descriptor |
| ETA_dPsi_A | Differential shape descriptor capturing subtle shape variations | PaDEL Descriptor |
| ETA_dPsi_B | Alternative differential shape descriptor | PaDEL Descriptor |
| ETA_Shape_P | Shape index descriptor quantifying molecular 3D structure | PaDEL Descriptor |
| ETA_Shape_Y | Alternative shape descriptor focusing on molecular topology | PaDEL Descriptor |
| ETA_Shape_X | Additional shape descriptor within ETA framework | PaDEL Descriptor |
| ETA_Beta | Branching descriptor capturing molecular connectivity | PaDEL Descriptor |
| ETA_BetaP | Path-corrected branching descriptor | PaDEL Descriptor |
| ETA_Beta_s | Simplified branching descriptor | PaDEL Descriptor |
| ETA_BetaP_s | Path-corrected simplified branching descriptor | PaDEL Descriptor |
| ETA_Beta_ns | Non-simplified branching descriptor | PaDEL Descriptor |
| ETA_BetaP_ns | Path-corrected non-simplified branching descriptor | PaDEL Descriptor |
| ETA_dBeta | Differential branching descriptor | PaDEL Descriptor |
| ETA_dBetaP | Path-corrected differential branching descriptor | PaDEL Descriptor |
| ETA_Beta_ns_d | Differential non-simplified branching descriptor | PaDEL Descriptor |
| ETA_BetaP_ns_d | Path-corrected differential non-simplified branching descriptor | PaDEL Descriptor |
| ETA_Eta | General molecular topological descriptor | PaDEL Descriptor |
| ETA_EtaP | Path-corrected version of ETA_Eta | PaDEL Descriptor |
| ETA_Eta_R | Radius-corrected ETA descriptor | PaDEL Descriptor |
| ETA_Eta_F | Feature-corrected ETA descriptor | PaDEL Descriptor |
| ETA_EtaP_F | Path and feature-corrected ETA descriptor | PaDEL Descriptor |
| ETA_Eta_L | Length-corrected ETA descriptor | PaDEL Descriptor |
| ETA_EtaP_L | Path and length-corrected ETA descriptor | PaDEL Descriptor |
| ETA_Eta_R_L | Radius and length-corrected ETA descriptor | PaDEL Descriptor |
| ETA_Eta_F_L | Feature and length-corrected ETA descriptor | PaDEL Descriptor |
| ETA_EtaP_F_L | Path, feature, and length-corrected ETA descriptor | PaDEL Descriptor |
| ETA_Eta_B | Branching-corrected ETA descriptor | PaDEL Descriptor |
| ETA_EtaP_B | Path and branching-corrected ETA descriptor | PaDEL Descriptor |
| ETA_Eta_B_RC | Radius-corrected branching ETA descriptor | PaDEL Descriptor |
| ETA_EtaP_B_RC | Path, radius, and branching-corrected ETA descriptor | PaDEL Descriptor |
| topoRadius | Molecular topological radius representing the average distance of atoms from the molecular center | PaDEL Descriptor |
| topoDiameter | Topological diameter indicating the longest shortest path between any two atoms in the molecule | PaDEL Descriptor |
| topoShape | Descriptor quantifying the overall shape of the molecule based on topology | PaDEL Descriptor |
| GGI (1-10) | Topological charge indices capturing charge distribution over the molecular graph | PaDEL Descriptor |
| JGI (1-10) | Mean topological charge indices reflecting local electronic environments | PaDEL Descriptor |
| JGT | Total topological charge index summarizing overall charge distribution | PaDEL Descriptor |

## Electronic Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| MaxAbsEStateIndex | Maximum absolute E-state value in the molecule | RDKit |
| MaxAbsPartialCharge | Maximum absolute atomic partial charge | RDKit |
| MaxEStateIndex | Maximum E-state value in the molecule | RDKit |
| MaxPartialCharge | Highest partial charge in the molecule | RDKit |
| NumValenceElectrons | Total number of valence electrons in the molecule | RDKit |
| NumRadicalElectrons | Total number of unpaired electrons (radicals) | RDKit |
| HallKierAlpha | Atom-type electrotopological descriptor modeling polarity and hybridization | RDKit |
| SpMax_Bhm (1-8), SpMin_Bhm (1-8) | Largest and smallest Burden matrix eigenvalues weighted by atomic mass | PaDEL Desciptor |
| SpMax_Bhe (1-8), SpMin_Bhe (1-8) | Largest and smallest Burden matrix eigenvalues weighted by electronegativity | PaDEL Descriptor |
| SpMax_Bhp (1-8), SpMin_Bhp (1-8) | Largest and smallest Burden matrix eigenvalues weighted by polarizability | PaDEL Descriptor |
| SpMax_Bhi (1-8), SpMin_Bhi (1-8) | Largest and smallest Burden matrix eigenvalues weighted by ionization potential | PaDEL Descriptor |
| SpMax_Bhs (1-8), SpMin_Bhs (1-8) | Largest and smallest Burden matrix eigenvalues weighted by intrinsic state | PaDEL Descriptor |
| SpMax_Dt | Maximum eigenvalue of the molecular graph adjacency matrix indicating molecular branching | PaDEL Descriptor |
| SpDiam_Dt | Spectral diameter from adjacency matrix representing molecular connectivity | PaDEL Descriptor |
| SpAD_Dt | Spectral absolute deviation describing variability in molecular adjacency eigenvalues | PaDEL Descriptor |
| SpMAD_Dt | Spectral mean absolute deviation quantifying spread of adjacency eigenvalues | PaDEL Descriptor |
| EE_Dt | Eccentric connectivity index measuring branching and distance patterns | PaDEL Descriptor |
| VE_Dt (1-3) | Eigenvalue-based descriptors related to molecular connectivity and valence | PaDEL Descriptor |
| VR_Dt (1-3) | Valence-related eigenvalue descriptors reflecting molecular electronic structure | PaDEL Descriptor |

## Surface Area Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| EState_VSA (1-9) | Binned E-state values mapped over atomic surface areas | RDKit |
| LabuteASA | Approximate surface area calculated by Labute’s method | RDKit |
| PEOE_VSA (1-9) | Partial equalization of orbital electronegativities across molecular surface area bins | RDKit |
| SMR_VSA (1-9) | Molar refractivity contributions binned by surface area | RDKit |
| SlogP_VSA (1-12) | SlogP contributions mapped to molecular surface area bins | RDKit |
| TPSA | Topological polar surface area related to hydrogen bonding potential | RDKit |

## Structural Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| NumAliphaticCarbocycles | Aliphatic (non-aromatic) carbon-only rings | RDKit |
| NumAliphaticHeterocycles | Aliphatic rings with at least one heteroatom | RDKit |
| NumAliphaticRings | Number of rings without aromaticity | RDKit |
| NumAromaticCarbocycles | Number of aromatic rings with only carbon atoms | RDKit |
| NumAromaticHeterocycles | Aromatic rings containing heteroatoms | RDKit |
| NumAromaticRings | Number of rings with aromatic character | RDKit |
| NumSaturatedCarbocycles | Fully saturated carbon rings | RDKit |
| NumSaturatedHeterocycles | Fully saturated heteroatom rings | RDKit |
| NumSaturatedRings | Total saturated rings in the molecule | RDKit |
| FractionCSP | Fraction of sp3 hybridized carbons (aliphatic character) | RDKit |
| NumRotatableBonds | Number of freely rotatable single bonds | RDKit |
| NumHAcceptors | Number of hydrogen bond acceptors (e.g. O, N atoms) | RDKit |
| NumHDonors | Number of hydrogen bond donors (e.g. -OH, -NH) | RDKit |
| HybRatio | Ratio of sp³ hybridized atoms to total heavy atoms, indicating molecular saturation and branching | PaDEL Descriptor |
| nAromAtom, nAromBond | Count of aromatic atoms and aromatic bonds in the molecule | PaDEL Descriptor |

## Functional Group Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| fr_Al_COO | Count of aliphatic carboxylic acid groups (-COOH) | RDKit |
| fr_Al_OH | Count of aliphatic alcohol groups (-OH) | RDKit |
| fr_Al_OH_noTert | Count of non-tertiary aliphatic alcohol groups | RDKit |
| fr_ArN | Count of aromatic amines | RDKit |
| fr_Ar_COO | Count of aromatic carboxylic acid groups | RDKit |
| fr_Ar_N | Count of aromatic nitrogen atoms (excluding amines) | RDKit |
| fr_Ar_NH | Count of aromatic amines (-NH- groups) | RDKit |
| fr_Ar_OH | Count of phenol groups (aromatic alcohols) | RDKit |
| fr_COO | Count of carboxylic acid groups (-COOH) | RDKit |
| fr_COO2 | Count of additional carboxylic acid groups (possibly duplicates) | RDKit |
| fr_C_O | Count of carbonyl groups (C=O) | RDKit |
| fr_C_O_noCOO | Count of carbonyls excluding carboxylic acids | RDKit |
| fr_C_S | Count of carbon-sulfur bonds | RDKit |
| fr_HOCCN | Count of hydroxyl connected to carbon connected to nitrogen | RDKit |
| fr_Imine | Count of imine groups (C=NH or C=NR) | RDKit |
| fr_NH0 | Count of nitrogen with zero attached hydrogens (e.g., quaternary) | RDKit |
| fr_NH1 | Count of nitrogen with one attached hydrogen | RDKit |
| fr_NH2 | Count of nitrogen with two attached hydrogens | RDKit |
| fr_N_O | Count of nitrogen-oxygen bonds | RDKit |
| fr_Ndealkylation1 | Count of first type N-dealkylation fragments | RDKit |
| fr_Ndealkylation2 | Count of second type N-dealkylation fragments | RDKit |
| fr_Nhpyrrole | Count of nitrogen atoms in pyrrole rings | RDKit |
| fr_SH | Count of thiol groups (-SH) | RDKit |
| fr_aldehyde | Count of aldehyde groups (-CHO) | RDKit |
| fr_alkyl_carbamate | Count of alkyl carbamate groups | RDKit |
| fr_alkyl_halide | Count of alkyl halide groups | RDKit |
| fr_allylic_oxid | Count of allylic oxidation sites | RDKit |
| fr_amide | Count of amide groups | RDKit |
| fr_amidine | Count of amidine groups | RDKit |
| fr_aniline | Count of aniline groups (aromatic amines) | RDKit |
| fr_aryl_methyl | Count of methyl groups attached to aromatic rings | RDKit |
| fr_azide | Count of azide groups (-N3) | RDKit |
| fr_azo | Count of azo groups (-N=N-) | RDKit |
| fr_barbitur | Count of barbituric acid-like groups | RDKit |
| fr_benzene | Count of benzene rings | RDKit |
| fr_benzodiazepine | Count of benzodiazepine rings | RDKit |
| fr_bicyclic | Count of bicyclic ring systems | RDKit |
| fr_diazo | Count of diazo groups | RDKit |
| fr_dihydropyridine | Count of dihydropyridine rings | RDKit |
| fr_epoxide | Count of epoxide rings | RDKit |
| fr_ester | Count of ester groups | RDKit |
| fr_ether | Count of ether groups | RDKit |
| fr_furan | Count of furan rings | RDKit |
| fr_guanido | Count of guanidine groups | RDKit |
| fr_halogen | Count of halogen atoms (F, Cl, Br, I)  | RDKit   |
| fr_hdrzine | Count of hydrazine groups | RDKit |
| fr_hdrzone | Count of hydrazone groups | RDKit |
| fr_imidazole | Count of imidazole rings | RDKit |
| fr_imid | Count of imide groups | RDKit |
| fr_isocyan | Count of isocyanate groups | RDKit |
| fr_isothiocyan | Count of isothiocyanate groups | RDKit |
| fr_ketone | Count of ketone groups | RDKit |
| fr_ketone_Topliss | Count of ketone groups (Topliss definition) | RDKit |
| fr_lactam | Count of lactam rings | RDKit |
| fr_lactone | Count of lactone rings | RDKit |
| fr_methoxy | Count of methoxy groups (-OCH3) | RDKit |
| fr_morpholine | Count of morpholine rings | RDKit |
| fr_nitrile | Count of nitrile groups (-C≡N) | RDKit |
| fr_nitro | Count of nitro groups (-NO2) | RDKit |
| fr_nitro_arom | Count of aromatic nitro groups | RDKit |
| fr_nitro_arom_nonortho | Count of non-ortho aromatic nitro groups | RDKit |
| fr_nitroso | Count of nitroso groups (-NO) | RDKit |
| fr_oxazole | Count of oxazole rings | RDKit |
| fr_oxime | Count of oxime groups (=NOH) | RDKit |
| fr_para_hydroxylation | Count of para hydroxylation sites | RDKit |
| fr_phenol | Count of phenol groups | RDKit |
| fr_phenol_noOrthoHbond | Count of phenols without ortho hydrogen bonding | RDKit |
| fr_phos_acid | Count of phosphoric acid groups | RDKit |
| fr_phos_ester | Count of phosphoric ester groups | RDKit |
| fr_piperdine | Count of piperidine rings | RDKit |
| fr_piperzine | Count of piperazine rings | RDKit |
| fr_priamide | Count of primary amides | RDKit |
| fr_prisulfonamd | Count of primary sulfonamides | RDKit |
| fr_pyridine | Count of pyridine rings | RDKit |
| fr_quatN | Count of quaternary nitrogen atoms | RDKit |
| fr_sulfide | Count of sulfide groups (thioethers) | RDKit |
| fr_sulfonamd | Count of sulfonamide groups | RDKit |
| fr_sulfone | Count of sulfone groups | RDKit |
| fr_term_acetylene | Count of terminal acetylene groups | RDKi |
| fr_tetrazole | Count of tetrazole rings | RDKit |
| fr_thiazole | Count of thiazole rings | RDKit |
| fr_thiocyan | Count of thiocyanate groups | RDKit |
| fr_thiophene | Count of thiophene rings | RDKit |
| fr_unbrch_alkane | Count of unbranched alkane chains | RDKit |
| fr_urea | Count of urea groups | RDKit |
| nAcid | Counts the number of acidic functional groups in the molecule | PaDEL Descriptor |
| nBase | Counts the number of basic functional groups in the molecule | PaDEL Descriptor |

## Constitutional Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| ExactMolWt | Monoisotopic molecular weight of the compound | RDKit |
| HeavyAtomCount | Number of non-hydrogen atoms in the molecule | RDKit |
| HeavyAtomMolWt | Molecular weight of non-hydrogen atoms only | RDKit |
| MolWt | Molecular weight based on average atomic masses | RDKit |
| NumHeteroatoms | Number of non-carbon and non-hydrogen atoms | RDKit |

## Physiochemical Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| MolLogP | Partition coefficient (logP) estimating hydrophobicity | RDKit |
| MolMR | Molar refractivity related to volume and polarizability | RDKit |
| CrippenLogP | Calculates LogP (octanol-water partition coefficient) via atomic contributions | PaDEL Descriptor |
| CrippenMR | Calculates molar refractivity (a measure of polarizability) via atomic contributions | PaDEL Descriptor |

## 2D Autocorrelation Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| ATS_m (1-8), ATS_v (1-8), ATS_e (1-8), ATS_p (1-8), ATS_i (1-8), ATS_s (1-8) | Broto-Moreau autocorrelation descriptors measuring atomic property correlations at different topological distances | PaDEL Descriptor |
| AATS_m (1-8), AATS_v (1-8), AATS_e (1-8), AATS_p (1-8), AATS_i (1-8), AATS_s (1-8) | Average Broto-Moreau autocorrelation descriptors capturing average correlations of atomic properties over molecular graph distances | PaDEL Desciptor |
| ATSC_m (1-8), ATSC_v (1-8), ATSC_e (1-8), AATS_p (1-8), AATS_i (1-8), AATS_s (1-8) | Centered Broto-Moreau autocorrelation descriptors emphasizing deviations from the mean atomic property values across bonds | PaDEL Descriptor |
| AATSC_m (1-8), AATSC_v (1-8), AATSC_e (1-8), AATSC_p (1-8), AATSC_i (1-8), AATSC_s (1-8) | Average Centered Broto-Moreau autocorrelation descriptors calculating averaged centered correlations of atomic properties by lag | PaDEL Descriptor |
| MATS_m (1-8), MATS_v (1-8), MATS_e (1-8), MATS_p (1-8), MATS_i (1-8), MATS_s (1-8) | Moran autocorrelation descriptors quantifying spatial autocorrelation of atomic properties over molecular topology | PaDEL Descriptor |
| GATS_m (1-8), GATS_v (1-8), GATS_e (1-8), GATS_p (1-8), GATS_i (1-8), GATS_s (1-8) | Geary autocorrelation descriptors measuring molecular graph-based heterogeneity of atomic property distributions | PaDEL Descriptor |

## Fingerprint Descriptors

| Descriptor | Summary | Software Used |
|------------|---------|---------------|
| Morgan_Fingerprint | Circular fingerprints encoding molecular substructures up to a given radius, used for similarity and machine learning | RDKit |
