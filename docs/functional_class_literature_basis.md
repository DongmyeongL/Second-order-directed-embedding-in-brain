# Literature-Guided Functional Classes For Current FCV Plot

This table defines the coarse functional classes used to color nodes in
`geometric_fcv_variants_vs_edge_fcv_all_species.png`.

## Shared Classes

| class_order | class_label | Operational meaning |
|---:|---|---|
| 0 | sensory input | Primary or early sensory-processing nodes. |
| 1 | sensorimotor output | Premotor, descending, locomotor-command, hindbrain/cerebellar, or navigation/visuomotor-output nodes. |
| 2 | integrative relay | Intermediate multimodal relay or association neuropils without a clearly dominant learning/state or motor-output role. |
| 3 | associative / state-dependent | Learning/memory, behavioral-state, neuromodulatory, sleep/arousal, hypothalamic, habenular, or higher associative nodes. |

## C. elegans

Basis:
- WormAtlas individual neuron pages classify neurons by known sensory, interneuron, and motor functions and describe AIB as integrating amphid sensory-neuron information.
- White et al. and subsequent circuit studies describe first-layer amphid interneurons such as AIB, AIY and AIZ as sensory-integration/interneuron nodes.
- Whole-brain activity and circuit studies describe AVA/AVB/AVD/AVE/PVC as locomotor command or premotor interneurons, and RIM/RMG/RIS/ALA-like neurons as state-modulatory or behavioral-state related nodes.

Mapping:
- Sensory input: amphid/cephalic/labial sensory neurons such as ADF, ADL, ADE, AFD, ASE, ASG, ASH, ASI, ASK, AWA, AWB, AWC, BAG, CEP, FLP, IL1, IL2, OLL, OLQ, URX, URY.
- Sensorimotor output: locomotor command interneurons AVA, AVB, AVD, AVE, PVC; head motor/premotor neurons RMD, RME, SMD, SMB, RIV, SAA, URA.
- Associative / state-dependent: AIB, AIM, AIN, AIY, AIZ, ALA, AUA, RIA, RIB, RIC, RIM, RIS, RMG.
- Integrative relay: remaining current-plot interneurons.

## Drosophila

Basis:
- The Insect Brain Name Working Group standardized Drosophila neuropils into optic lobe, antennal lobe/lateral horn, mushroom body, central complex, lateral accessory lobe, superior protocerebrum, and inferior/lateral protocerebrum subdivisions.
- Hemibrain/FlyEM descriptions emphasize the mushroom body for associative learning and memory, the central complex for navigation, and optic-lobe neurons for visual processing.
- Reviews of the central complex and lateral accessory lobe describe these regions as navigation, visuomotor, and premotor/descending interfaces rather than purely associative memory structures.

Mapping:
- Sensory input: AL/LH olfactory stream; ME/LO/AME/AOTU visual/optic stream.
- Sensorimotor output: PB/NO central-complex navigation/visuomotor nodes; LAL/VES/SPS premotor or descending-interface nodes.
- Associative / state-dependent: MB.
- Integrative relay: SIP/SMP/SLP superior protocerebrum; AVLP/PVLP/PLP/WED visual-association protocerebrum; CRE/ICL/IB/ATL/EPA lateral/inferior protocerebrum.

## Zebrafish

Basis:
- Z-Brain defines larval zebrafish anatomical regions and maps whole-brain activity to telencephalon, diencephalon, mesencephalon, rhombencephalon, spinal cord, and ganglia/other.
- Zebrafish optic tectum and pretectal literature identifies TeO and pretectal areas as major visual-processing regions.
- Zebrafish forebrain reviews relate pallium/subpallium to higher associative forebrain functions; hindbrain, reticulospinal, cerebellar, and medullary regions are strongly tied to motor control, escape, and sensorimotor transformations.
- Habenula/interpeduncular/raphe/hypothalamic systems are treated as behavioral-state, neuromodulatory, or valence-related circuits.

Mapping:
- Sensory input: OB/OE/OG olfactory regions; TeO/TL/TS/PT/PrT visual optic-midbrain/pretectal regions.
- Sensorimotor output: Cb, MON, MOS1-5, IO, aRF, imRF, pRF, TG, VR, NX.
- Associative / state-dependent: P/SP telencephalic regions; Hb/Hc/Hi/HR/IPN/Ra state-modulatory forebrain/raphe/habenular system.
- Integrative relay: Th, PO, T, and remaining current-plot relay regions.

## Sources Used

- WormAtlas individual neuron pages, especially AIB: https://www.wormatlas.org/neurons/Individual%20Neurons/AIBmainframe.htm
- White et al., 1986, *The Structure of the Nervous System of the Nematode Caenorhabditis elegans*: https://wormatlas.org/MoW_built0.92/1986%20White%20Southgate%20Thomson%20Brenner%20PhilTransRoySocB_mh.pdf
- Kato et al. / C. elegans whole-brain dynamics context summarized in dynamic-connectome review: https://www.sciencedirect.com/science/article/pii/S0959438821001549
- Ito et al., 2014, *A Systematic Nomenclature for the Insect Brain*: https://www.sciencedirect.com/science/article/pii/S0896627313011781
- Janelia FlyEM hemibrain overview: https://www.janelia.org/project-team/flyem/hemibrain
- Pfeiffer and Homberg, central complex review: https://pubmed.ncbi.nlm.nih.gov/24160424/
- Randlett et al., 2015, Z-Brain atlas: https://www.nature.com/articles/nmeth.3581
- Baier and Scott, 2009/2010, zebrafish optic tectum circuitry review: https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-126
- Cheng et al., zebrafish forebrain review: https://pmc.ncbi.nlm.nih.gov/articles/PMC3895987/
