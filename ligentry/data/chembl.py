"""Fetch known active compounds for a target from the ChEMBL REST API (paginated)."""
from __future__ import annotations

import requests

ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def fetch_actives(
    target_chembl_id: str, limit: int = 20, max_nm: float = 100.0
) -> list[tuple[str, str]]:
    """Return up to `limit` unique (chembl_id, smiles) pairs with IC50 <= max_nm."""
    seen: set[str] = set()
    compounds: list[tuple[str, str]] = []
    offset = 0
    page_size = min(1000, max(200, limit))

    while len(compounds) < limit:
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type": "IC50",
            "standard_units": "nM",
            "standard_value__lte": max_nm,
            "limit": page_size,
            "offset": offset,
        }
        resp = requests.get(ACTIVITY_URL, params=params, timeout=60)
        resp.raise_for_status()
        payload = resp.json()
        activities = payload.get("activities", [])
        if not activities:
            break

        for activity in activities:
            mol_id = activity.get("molecule_chembl_id")
            smiles = activity.get("canonical_smiles")
            if not mol_id or not smiles or mol_id in seen:
                continue
            seen.add(mol_id)
            compounds.append((mol_id, smiles))
            if len(compounds) >= limit:
                break

        if not (payload.get("page_meta") or {}).get("next"):
            break
        offset += len(activities)

    return compounds
