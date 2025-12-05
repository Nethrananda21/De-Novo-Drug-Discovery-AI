import re
import json
import time
from pathlib import Path
from typing import List, Dict, Tuple

import requests
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey

# Optional rdMolStandardize import (not available in some RDKit builds)
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as std
except Exception:  # pragma: no cover
    try:
        # Older builds sometimes expose it directly under rdkit.Chem
        from rdkit.Chem import rdMolStandardize as std  # type: ignore
    except Exception:
        std = None

ROOT = Path(__file__).resolve().parents[1]
DATA_TXT = ROOT / "Final_Noval_Molecule_Data.txt"
RESULTS_DIR = ROOT / "results"
REPORT_MD = RESULTS_DIR / "mol_001_external_search_report.md"
VARIANTS_CSV = RESULTS_DIR / "mol_001_external_search_variants.csv"

HEADERS = {"User-Agent": "DENOVO/1.0 (+novelty-check; contact: user)"}
TIMEOUT = 20
TOP_N = 25

# Approximate moieties for substructure search
# 1) 6,7-dimethoxyquinazolin-4-amine core
MOIETY_SMILES_RING = "COc1cc2ncnc(N)c2cc1OC"
# 2) cyclohexane-1-carboxylic acid fragment
MOIETY_SMILES_CHEX_ACID = "C1CCC(CC1)C(=O)O"


def read_base_smiles() -> str:
    text = DATA_TXT.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"^SMILES:\s*([^\n\r]+)", text, flags=re.MULTILINE)
    if not m:
        raise RuntimeError("SMILES not found in Final_Noval_Molecule_Data.txt")
    return m.group(1).strip()


def standardize_variants(smiles: str) -> List[Tuple[str, Chem.Mol, str]]:
    variants: List[Tuple[str, Chem.Mol, str]] = []
    m = Chem.MolFromSmiles(smiles)
    if not m:
        raise RuntimeError("Invalid base SMILES, RDKit failed to parse")

    def add(tag: str, mol: Chem.Mol):
        if mol is None:
            return
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        variants.append((tag, Chem.MolFromSmiles(smi), smi))

    add("base", m)

    if std is not None:
        cleanup = std.Cleanup(m)
        add("cleanup", cleanup)

        # Reionization / uncharge / parents
        try:
            reion = std.Reionize(cleanup)
            add("reionized", reion)
        except Exception:
            pass

        try:
            uncharger = std.Uncharger()
            unchg = uncharger.uncharge(cleanup)
            add("uncharged", unchg)
        except Exception:
            pass

        try:
            ch_parent = std.ChargeParent(cleanup)
            add("charge_parent", ch_parent)
        except Exception:
            pass

        try:
            frag_parent = std.FragmentParent(cleanup)
            add("fragment_parent", frag_parent)
        except Exception:
            pass

        # Tautomer enumeration (limited to top 10)
        try:
            te = std.TautomerEnumerator()
            base_for_taut = variants[-1][1] if variants else cleanup
            tautomers = te.Enumerate(base_for_taut)
            for i, t in enumerate(tautomers[:10]):
                add(f"tautomer_{i+1}", t)
        except Exception:
            pass
    else:
        # Fallback variant generation without rdMolStandardize
        # Try deprotonating carboxylic acids to carboxylate
        try:
            rxn = AllChem.ReactionFromSmarts('[CX3:1](=O)[OX2H:2]>>[CX3:1](=O)[O-:2]')
            prods = rxn.RunReactants((m,))
            for i, prod_tuple in enumerate(prods[:3]):
                pm = prod_tuple[0]
                Chem.SanitizeMol(pm)
                add(f"carboxylate_{i+1}", pm)
        except Exception:
            pass

    # Deduplicate by canonical SMILES
    dedup: Dict[str, Tuple[str, Chem.Mol, str]] = {}
    for tag, mol, smi in variants:
        if smi not in dedup:
            dedup[smi] = (tag, mol, smi)
    return list(dedup.values())


def inchikey_info(mol: Chem.Mol) -> Tuple[str, str]:
    try:
        inchi = MolToInchi(mol)
        ik = MolToInchiKey(mol)
        return inchi, ik
    except Exception:
        return "", ""


def q_pubchem_inchikey(ik: str) -> Dict:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{ik}/cids/JSON"
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
    if r.status_code == 404:
        return {"cids": [], "status": 404}
    try:
        js = r.json()
    except Exception:
        return {"cids": [], "status": r.status_code}
    cids = js.get("IdentifierList", {}).get("CID", [])
    return {"cids": cids, "status": r.status_code}


def q_chembl_inchikey(ik: str) -> Dict:
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/molecule?"
        f"molecule_structures__standard_inchi_key={ik}"
    )
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
    if r.status_code != 200:
        return {"count": None, "status": r.status_code}
    # API may return XML/JSON depending on headers; try to parse minimal info
    txt = r.text
    # Simple detection: total_count in XML page_meta
    m = re.search(r"<total_count>(\d*)</total_count>", txt)
    if m:
        cnt = int(m.group(1)) if m.group(1) else 0
        return {"count": cnt, "status": 200}
    # Try JSON
    try:
        js = r.json()
        cnt = js.get("page_meta", {}).get("total_count", 0)
        return {"count": cnt, "status": 200}
    except Exception:
        return {"count": 0, "status": 200}


def q_pubchem_similarity(smiles: str, threshold: int = 95) -> Dict:
    # Fast 2D similarity; returns CIDs if any above threshold
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/"
        f"smiles/{requests.utils.quote(smiles)}/cids/JSON?Threshold={threshold}"
    )
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
    if r.status_code != 200:
        return {"cids": [], "status": r.status_code}
    try:
        js = r.json()
        cids = js.get("IdentifierList", {}).get("CID", [])
        return {"cids": cids, "status": 200}
    except Exception:
        return {"cids": [], "status": r.status_code}


def q_pubchem_substructure(smiles: str, max_records: int = TOP_N) -> Dict:
    # Substructure search by SMILES; may be large, cap via MaxRecords when supported
    base = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/"
        f"smiles/{requests.utils.quote(smiles)}/cids/JSON"
    )
    url = f"{base}?MaxRecords={max_records}" if max_records else base
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
    if r.status_code != 200:
        return {"cids": [], "status": r.status_code}
    try:
        js = r.json()
        cids = js.get("IdentifierList", {}).get("CID", [])
        return {"cids": cids, "status": 200}
    except Exception:
        return {"cids": [], "status": r.status_code}


def fetch_pubchem_props(cids: List[int]) -> Dict[int, Dict]:
    props: Dict[int, Dict] = {}
    if not cids:
        return props
    for i in range(0, len(cids), 25):
        chunk = cids[i:i+25]
        cid_str = ",".join(map(str, chunk))
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid_str}/property/Title,IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"
        )
        r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        if r.status_code != 200:
            continue
        try:
            js = r.json()
            items = js.get("PropertyTable", {}).get("Properties", [])
            for it in items:
                cid = it.get("CID")
                if cid is not None:
                    props[int(cid)] = it
        except Exception:
            continue
    return props


def sim_and_flags(base_smiles: str, target_smiles: str) -> Tuple[float, bool, bool]:
    base = Chem.MolFromSmiles(base_smiles)
    targ = Chem.MolFromSmiles(target_smiles) if target_smiles else None
    if not base or not targ:
        return 0.0, False, False
    fpb = AllChem.GetMorganFingerprintAsBitVect(base, 2, nBits=2048)
    fpt = AllChem.GetMorganFingerprintAsBitVect(targ, 2, nBits=2048)
    try:
        sim = DataStructs.TanimotoSimilarity(fpb, fpt)
    except Exception:
        sim = 0.0
    ring = Chem.MolFromSmiles(MOIETY_SMILES_RING)
    hexacid = Chem.MolFromSmiles(MOIETY_SMILES_CHEX_ACID)
    shares_core = targ.HasSubstructMatch(ring) if ring is not None else False
    has_cyclohex_acid = targ.HasSubstructMatch(hexacid) if hexacid is not None else False
    return float(sim), shares_core, has_cyclohex_acid


def append_neighbor_details(base_smiles: str):
    # Get top-N similarity neighbors (>=90)
    sim = q_pubchem_similarity(base_smiles, 90)
    sim_cids = sim.get("cids", [])[:TOP_N]
    props = fetch_pubchem_props(sim_cids)
    rows = []
    for cid in sim_cids:
        it = props.get(cid, {})
        name = it.get("Title") or it.get("IUPACName") or "—"
        formula = it.get("MolecularFormula") or "—"
        mw = it.get("MolecularWeight")
        smi = it.get("CanonicalSMILES") or ""
        simv, core, hexacid = sim_and_flags(base_smiles, smi)
        note_bits = []
        note_bits.append("shares core" if core else "lacks core")
        note_bits.append("has cyclohexyl-COOH" if hexacid else "lacks cyclohexyl-COOH")
        note = "; ".join(note_bits)
        rows.append({
            "cid": cid,
            "name": name,
            "formula": formula,
            "mw": f"{mw:.2f}" if isinstance(mw, (int, float)) else (mw or "—"),
            "sim": f"{simv:.2f}",
            "note": note,
        })

    # Build Markdown table
    lines = []
    lines.append("\n\n### Similarity Neighbor Details (top %d)" % len(rows))
    lines.append("- Columns: CID, Name, Formula, MW, Sim(Tanimoto), Note")
    lines.append("")
    lines.append("| CID | Name | Formula | MW | Sim | Note |")
    lines.append("| --- | ---- | ------- | --:| ---:| ---- |")
    for r in rows:
        link = f"[CID {r['cid']}](https://pubchem.ncbi.nlm.nih.gov/compound/{r['cid']})"
        lines.append(f"| {link} | {r['name']} | {r['formula']} | {r['mw']} | {r['sim']} | {r['note']} |")

    REPORT_MD.write_text(REPORT_MD.read_text(encoding="utf-8") + "\n" + "\n".join(lines), encoding="utf-8")


def write_csv(rows: List[Dict]):
    import csv
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "variant_tag","smiles","inchi","inchikey",
        "pubchem_status","pubchem_cids","chembl_status","chembl_count",
        "sim95_status","sim95_cids","sim90_status","sim90_cids"
    ]
    with VARIANTS_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def append_report_summary(rows: List[Dict]):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    lines = []
    lines.append("\n\n## Automated Variant Search Summary (" + ts + ")\n")
    total = len(rows)
    pubchem_hits = sum(1 for r in rows if r.get("pubchem_cids"))
    chembl_hits = sum(1 for r in rows if (r.get("chembl_count") or 0) > 0)
    lines.append(f"- Variants checked: {total}")
    lines.append(f"- PubChem exact InChIKey hits: {pubchem_hits}")
    lines.append(f"- ChEMBL exact InChIKey hits: {chembl_hits}")
    sim95 = sum(1 for r in rows if r.get("sim95_cids"))
    sim90 = sum(1 for r in rows if r.get("sim90_cids"))
    lines.append(f"- PubChem similarity hits (>=95): {sim95}")
    lines.append(f"- PubChem similarity hits (>=90): {sim90}")
    lines.append("\nArtifacts\n- CSV: `results/mol_001_external_search_variants.csv`\n")

    REPORT_MD.write_text(REPORT_MD.read_text(encoding="utf-8") + "\n" + "\n".join(lines), encoding="utf-8")


def append_neighbors_section(base_smiles: str):
    # Similarity neighbors (top N from fast similarity >= 90)
    sim = q_pubchem_similarity(base_smiles, 90)
    sim_cids = sim.get("cids", [])[:TOP_N]

    # Substructure neighbors: intersection of ring core and cyclohexyl acid
    ss_ring = q_pubchem_substructure(MOIETY_SMILES_RING, TOP_N * 5)
    ss_hex = q_pubchem_substructure(MOIETY_SMILES_CHEX_ACID, TOP_N * 5)
    set_ring = set(ss_ring.get("cids", []) or [])
    set_hex = set(ss_hex.get("cids", []) or [])
    inter = list(set_ring.intersection(set_hex))[:TOP_N]

    def links(cids: List[int]) -> List[str]:
        return [f"- CID {cid}: https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" for cid in cids]

    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    blk = []
    blk.append(f"\n\n## PubChem Neighbor Snapshot ({ts})\n")
    blk.append(f"- Similarity neighbors (>=90, top {TOP_N}): {len(sim_cids)}")
    blk.extend(links(sim_cids))
    blk.append("")
    blk.append(f"- Substructure neighbors (ring ∩ cyclohexyl-acid, top {TOP_N}): {len(inter)}")
    blk.extend(links(inter))
    blk.append("")

    REPORT_MD.write_text(REPORT_MD.read_text(encoding="utf-8") + "\n" + "\n".join(blk), encoding="utf-8")


def main():
    smiles = read_base_smiles()
    variants = standardize_variants(smiles)

    rows: List[Dict] = []
    for tag, mol, smi in variants:
        try:
            inchi, ik = inchikey_info(mol)
        except Exception:
            inchi, ik = "", ""
        pc = q_pubchem_inchikey(ik) if ik else {"cids": [], "status": None}
        ch = q_chembl_inchikey(ik) if ik else {"count": None, "status": None}
        sim95 = q_pubchem_similarity(smi, 95)
        sim90 = q_pubchem_similarity(smi, 90)
        row = {
            "variant_tag": tag,
            "smiles": smi,
            "inchi": inchi,
            "inchikey": ik,
            "pubchem_status": pc.get("status"),
            "pubchem_cids": ",".join(map(str, pc.get("cids", []))) if pc.get("cids") else "",
            "chembl_status": ch.get("status"),
            "chembl_count": ch.get("count"),
            "sim95_status": sim95.get("status"),
            "sim95_cids": ",".join(map(str, sim95.get("cids", []))) if sim95.get("cids") else "",
            "sim90_status": sim90.get("status"),
            "sim90_cids": ",".join(map(str, sim90.get("cids", []))) if sim90.get("cids") else "",
        }
        rows.append(row)

    write_csv(rows)
    append_report_summary(rows)
    try:
        append_neighbors_section(smiles)
    except Exception:
        # Non-fatal; neighbors are best-effort
        pass
    try:
        append_neighbor_details(smiles)
    except Exception:
        pass
    print(f"Wrote {len(rows)} rows to {VARIANTS_CSV}")
    print(f"Appended summary to {REPORT_MD}")


if __name__ == "__main__":
    main()
