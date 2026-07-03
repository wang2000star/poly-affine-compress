"""Gate-level circuit simulator - robust version."""
import re

class GateSim:
    def __init__(self, text):
        self.inputs = []
        self.outputs = []
        self.stmts = []  # (name, rhs_string)
        self._parse(text)

    def _parse(self, text):
        minp = re.search(r'INORDER\s*=\s*([^;]+);', text)
        mout = re.search(r'OUTORDER\s*=\s*([^;]+);', text)
        self.inputs = minp.group(1).strip().split()
        self.outputs = mout.group(1).strip().split()
        for line in text.split('\n'):
            line = line.strip()
            if not line or '=' not in line:
                continue
            if line.startswith('INORDER') or line.startswith('OUTORDER'):
                continue
            parts = line.split('=', 1)
            self.stmts.append((parts[0].strip(), parts[1].strip().rstrip(';')))

    def evaluate(self, input_dict):
        vals = dict(input_dict)
        for name, rhs in self.stmts:
            if name in vals:
                continue
            rhs = rhs.strip()
            if rhs == '0':
                vals[name] = 0; continue
            if rhs == '1':
                vals[name] = 1; continue

            # Build a Python-evaluable expression
            # Replace each token one at a time to avoid number-vs-variable confusion
            py = rhs
            # Replace !var -> (1-val)
            py = re.sub(r'!(\w+)', lambda m: str(1 - vals.get(m.group(1), 0)), py)
            # Replace var -> val (skip numeric tokens)
            py = re.sub(r'\b(\w+)\b', lambda m: str(vals.get(m.group(1), m.group(0))) if not m.group(1).isdigit() else m.group(0), py)
            # Replace operators
            py = py.replace('*', '&').replace('+', '|')
            try:
                vals[name] = eval(py) & 1
            except Exception as e:
                print(f"  EVAL ERROR: {name} = {rhs} -> {py} : {e}")
                vals[name] = 0
        return vals

    def get_outputs(self, vals):
        return {o: vals.get(o, 0) for o in self.outputs}


def probe(sim, description, test_bits_list, output_names):
    print(f"\n--- {description} ---")
    for test_bits in test_bits_list:
        d = dict(zip(sim.inputs, [int(b) for b in test_bits]))
        vals = sim.evaluate(d)
        nz = [o for o in output_names if vals.get(o, 0)]
        print(f"  input={test_bits[:30]}{'...' if len(test_bits)>30 else ''}  nonzero: {nz}")


if __name__ == '__main__':
    # hd09: om_0..om_15 non-trivial
    text = open('/home/wangfz/bool/examples/hd09/hd09.txt').read()
    sim = GateSim(text)
    nz_out = [f'om_{i}' for i in range(16)]

    # Key hypothesis: hd09 is a priority encoder over i0..i15
    # with n97 = NOR(i0..i15) controlling the upper-bit path
    tests = []
    # All zero
    tests.append('0' * 32)
    # Single bit at each position in lower 16
    for k in range(16):
        bits = ['0'] * 32
        bits[15 - k] = '1'  # i_{15-k}=1 (since i0 is bit 0 in our ordering)
        tests.append(''.join(bits))
    # Single bit in upper half (i16 and above)
    for k in range(16, 18):
        bits = ['0'] * 32
        bits[k] = '1'
        tests.append(''.join(bits))
    # HI = all lower 16 zero + some conditions
    # Lower 16 all zero, i16=1, i17=0
    bits = ['0'] * 32
    bits[16] = '1'
    tests.append(''.join(bits))
    # Lower 16 all zero, i16=1, i17=1
    bits = ['0'] * 32
    bits[16] = '1'
    bits[17] = '1'
    tests.append(''.join(bits))
    # i0=1, i16=0, i17=0
    bits = ['0'] * 16 + ['0', '0'] + ['0'] * 14
    bits[0] = '1'
    tests.append(''.join(bits))

    probe(sim, "hd09 priority encoder test", tests, nz_out)

    # hd12: m_26..m_31 non-trivial
    text2 = open('/home/wangfz/bool/examples/hd12/hd12.txt').read()
    sim2 = GateSim(text2)
    nz_out2 = [f'm_{i}' for i in range(26, 32)]
    tests2 = []
    tests2.append('0' * 32)
    for k in range(32):
        bits = ['0'] * 32
        bits[k] = '1'
        tests2.append(''.join(bits))
    probe(sim2, "hd12 single-bit tests", tests2[:10], nz_out2)

    # router: outport0, outport1, outport2
    text3 = open('/home/wangfz/bool/examples/router/router.txt').read()
    sim3 = GateSim(text3)
    nz_out3 = ['outport0', 'outport1', 'outport2']
    tests3 = []
    tests3.append('0' * 60)
    # Single 1 in destx
    for k in range(3):
        bits = ['0'] * 60
        bits[k] = '1'
        tests3.append(''.join(bits))
    # Single 1 in desty
    for k in range(30, 33):
        bits = ['0'] * 60
        bits[k] = '1'
        tests3.append(''.join(bits))
    probe(sim3, "router single-bit tests", tests3, nz_out3)
