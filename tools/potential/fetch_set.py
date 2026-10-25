#!/usr/bin/env python3
"""
Training set of the statistical potential: one chain from each cluster (30%
sequence identity) of high-resolution X-ray protein structures in the PDB.

    fetch_set.py list  training_set.txt          # query RCSB, write "pdbid chain" lines
    fetch_set.py fetch training_set.txt CACHEDIR # download the pdb files (gzipped)

The list is kept in the repository, so that the potential can be rebuilt from
the same structures; the downloaded files are not.
"""
import gzip, json, os, random, sys, time, urllib.request

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"

QUERY = {
    "query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": 1.6}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_entity_polymer_type", "operator": "exact_match", "value": "Protein"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_sample_sequence_length", "operator": "range",
                                                                "value": {"from": 60, "to": 400, "include_lower": True, "include_upper": True}}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "refine.ls_R_factor_R_free", "operator": "less", "value": 0.22}}]},
    "return_type": "polymer_entity",
    "request_options": {"group_by": {"aggregation_method": "sequence_identity", "similarity_cutoff": 30,
                                     "ranking_criteria_type": {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}},
                        "group_by_return_type": "representatives", "paginate": {"start": 0, "rows": 10000},
                        "results_content_type": ["experimental"]}}

NSET = 1500		# chains in the training set
SEED = 1		# the sample of clusters is reproducible


def post(url, data):
    req = urllib.request.Request(url, data=json.dumps(data).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def make_list(fname):
    res = post(SEARCH, QUERY)
    ents = sorted(r["identifier"] for r in res["result_set"])
    print(f"{len(ents)} clusters", file=sys.stderr)
    random.Random(SEED).shuffle(ents)
    out = []
    for k in range(0, len(ents), 200):
        if len(out) >= NSET:
            break
        batch = ents[k:k+200]
        q = '{ polymer_entities(entity_ids: %s) { rcsb_id rcsb_polymer_entity_container_identifiers { auth_asym_ids } } }' % json.dumps(batch)
        d = post(GRAPHQL, {"query": q})
        chains = {e["rcsb_id"]: e["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] for e in d["data"]["polymer_entities"]}
        for e in batch:
            ids = sorted(chains.get(e) or [])
            if ids and len(out) < NSET:
                out.append((e.split("_")[0].lower(), ids[0]))
    with open(fname, "w") as f:
        f.write("# pdbid chain (fetch_set.py: X-ray <= 1.6 A, Rfree < 0.22, 60-400 residues, one per 30% identity cluster)\n")
        for pdb, ch in out:
            f.write(f"{pdb} {ch}\n")
    print(f"{len(out)} chains written to {fname}", file=sys.stderr)


def fetch(fname, cache):
    os.makedirs(cache, exist_ok=True)
    n = bad = 0
    for line in open(fname):
        if line.startswith("#") or not line.strip():
            continue
        pdb = line.split()[0]
        dest = os.path.join(cache, pdb + ".pdb.gz")
        if os.path.exists(dest):
            continue
        try:
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.pdb.gz", timeout=60) as r:
                data = r.read()
            gzip.decompress(data)			# check it is complete
            with open(dest, "wb") as f:
                f.write(data)
            n += 1
        except Exception as e:		# e.g. entries without the pdb format
            bad += 1
            print(f"{pdb}: {e}", file=sys.stderr)
        time.sleep(0.05)
    print(f"downloaded {n}, failed {bad}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "list":
        make_list(sys.argv[2])
    elif len(sys.argv) == 4 and sys.argv[1] == "fetch":
        fetch(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
