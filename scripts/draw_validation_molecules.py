"""Draw Fig. 5 from verified molecular graphs rather than a stylized sketch.

PubChem PUG REST structure records (retrieved 2026-10-02):
tetraphenylporphyrin CID 86280046, docetaxel CID 148124,
paclitaxel CID 36314. Stereochemical SMILES are retained for the taxanes.
"""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors


ROOT = Path(__file__).parents[1]
CIDS = [86280046, 148124, 36314]
LABELS = [
    "Tetraphenylporphyrin  C44H30N4 (78 atoms)",
    "Docetaxel  C43H53NO14 (111 atoms)",
    "Paclitaxel  C47H51NO14 (113 atoms)",
]
EXPECTED = ["C44H30N4", "C43H53NO14", "C47H51NO14"]

molecules = [
    Chem.RemoveHs(
        Chem.SDMolSupplier(
            str(ROOT / "data" / "figure5_pubchem_2d" / f"{cid}.sdf"),
            removeHs=False,
        )[0]
    )
    for cid in CIDS
]
assert all(m is not None for m in molecules)
assert [rdMolDescriptors.CalcMolFormula(m) for m in molecules] == EXPECTED
assert sum(atom.GetSymbol() == "N" for atom in molecules[0].GetAtoms()) == 4
assert len(Chem.GetSymmSSSR(molecules[0])) == 9

options = Draw.MolDrawOptions()
options.legendFontSize = 23
options.bondLineWidth = 2.1
options.padding = 0.06
image = Draw.MolsToGridImage(
    molecules,
    molsPerRow=3,
    subImgSize=(850, 680),
    legends=LABELS,
    drawOptions=options,
    useSVG=False,
)
output = ROOT / "figures" / "large_molecule_validation_panel_rdkit.png"
image.save(output)
print(output)
