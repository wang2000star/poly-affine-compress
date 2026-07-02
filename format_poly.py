"""
Format ANF optimization results as separate transformation and expression files.

Output convention:
  {instance}_opt1_trans.poly — Strategy 1 shared transformation (z = Mx ⊕ b)
  {instance}_opt1_expr.poly  — Strategy 1 expressions (output = g(z))
  {instance}_opt2_trans.poly — Strategy 2 merged transformation
  {instance}_opt2_expr.poly  — Strategy 2 expressions (output = g(z))

z-variables are named z_0, z_1, ... globally.
Strategy 2 merges per-output transformations into a single shared M,b
with deduplication: identical (M_row, b) pairs get the same z-index.
"""
import numpy as np


def merge_transformations(output_results):
    """Merge per-output (g, M, b) into one shared transformation.

    Args:
        output_results: list of (g: SparseANF, M: np.ndarray, b: np.ndarray)

    Returns:
        merged_M: list of rows (list of int)
        merged_b: list of int
        output_g_remapped: list of dict {mask: coeff} with merged z-indices
    """
    rows = []         # (row_tuple, b_val) → z_index
    row_to_z = {}

    output_g_remapped = []
    for g, M_arr, b_arr in output_results:
        M = np.array(M_arr, dtype=np.int64)
        b_flat = np.array(b_arr, dtype=np.int64).flatten()
        m = M.shape[0]

        # Map old z-indices to merged z-indices
        old_to_new = {}
        for j in range(m):
            row_tup = tuple(int(v) for v in M[j])
            bv = int(b_flat[j]) % 2
            key = (row_tup, bv)
            if key not in row_to_z:
                row_to_z[key] = len(rows)
                rows.append(key)
            old_to_new[j] = row_to_z[key]

        # Remap g terms
        remapped = {}
        for mask, coeff in g.terms.items():
            new_mask = 0
            t = mask
            while t:
                lsb = t & -t
                old_j = (lsb.bit_length() - 1)
                t ^= lsb
                new_mask |= (1 << old_to_new.get(old_j, old_j))
            remapped[new_mask] = coeff
        remapped = {k: v for k, v in remapped.items() if v}
        output_g_remapped.append(remapped)

    merged_M = [list(key[0]) for key in rows]
    merged_b = [key[1] for key in rows]
    return merged_M, merged_b, output_g_remapped


def write_transformation(path, inputs, M, b):
    """Write transformation file: z_j = linear_form(inputs)."""
    lines = []
    for j, (M_row, b_val) in enumerate(zip(M, b)):
        terms = []
        for i, v in enumerate(M_row):
            if int(v) % 2:
                terms.append(inputs[i])
        if int(b_val) % 2:
            terms.append("1")
        line = f"z_{j} = {' + '.join(terms) if terms else '0'}"
        lines.append(line)
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def _anf_term_str(mask):
    """Format a single ANF monomial mask as z-variable product string."""
    if mask == 0:
        return "1"
    parts = []
    t = mask
    while t:
        lsb = t & -t
        parts.append(f"z_{(lsb.bit_length() - 1)}")
        t ^= lsb
    return "*".join(parts)


def write_expressions(path, outputs, g_terms_list, verbose=False):
    """Write expression file: output_name = g(z) (compact, one per line)."""
    lines = []
    for name, terms in zip(outputs, g_terms_list):
        if not terms:
            lines.append(f"{name} = 0")
            if verbose:
                print(f"  {name} = 0")
            continue
        monomials = []
        for mask, coeff in sorted(terms.items()):
            if coeff == 0:
                continue
            monomials.append(_anf_term_str(mask))
        line = f"{name} = {' + '.join(monomials)}"
        lines.append(line)
        if verbose:
            print(f"  {name}: T={len(terms)}")
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def anf_dict_to_str(terms, var_names=None):
    """Convert ANF dict to human-readable string (for inline logging)."""
    if not terms:
        return "0"
    parts = []
    for mask, coeff in sorted(terms.items()):
        if coeff == 0:
            continue
        if mask == 0:
            parts.append("1")
        else:
            vars_ = []
            t = mask
            while t:
                lsb = t & -t
                vi = (lsb.bit_length() - 1)
                vars_.append(f"{var_names[vi] if var_names else f'z_{vi}'}")
                t ^= lsb
            parts.append("*".join(vars_))
    return " + ".join(parts)
