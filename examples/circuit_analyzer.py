"""
Circuit analyzer: parse .txt circuit files, build dependency graph,
determine which inputs each output depends on.
"""
import sys, re
sys.path.insert(0, '/home/wangfz/bool')

class CircuitAnalyzer:
    def __init__(self, text):
        self.text = text
        self.inputs = []
        self.outputs = []
        self.nodes = {}  # name -> (type, args)
        self._parse()

    def _parse(self):
        minp = re.search(r'INORDER\s*=\s*([^;]+);', self.text)
        mout = re.search(r'OUTORDER\s*=\s*([^;]+);', self.text)
        self.inputs = minp.group(1).strip().split()
        self.outputs = mout.group(1).strip().split()

        # Parse assignment statements
        for line in self.text.split('\n'):
            line = line.strip()
            if not line or '=' not in line:
                continue
            # Skip INORDER/OUTORDER
            if line.startswith('INORDER') or line.startswith('OUTORDER'):
                continue
            parts = line.split('=', 1)
            lhs = parts[0].strip()
            rhs_str = parts[1].strip()
            self.nodes[lhs] = rhs_str

    def get_deps_for_output(self, output_name):
        """Trace which inputs an output depends on."""
        visited = set()
        deps = set()
        def trace(name):
            if name in visited:
                return
            visited.add(name)
            if name in self.inputs:
                deps.add(name)
                return
            if name not in self.nodes:
                return
            rhs = self.nodes[name]
            # Extract variable names (words that are not operators)
            tokens = re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*', rhs)
            for tok in tokens:
                if tok in ('i', 'o', 'om', 'm', 'n', 'destx', 'desty', 'outport'):
                    continue
                trace(tok)
        trace(output_name)
        return sorted(deps, key=lambda x: self._input_index(x))

    def _input_index(self, name):
        try:
            return self.inputs.index(name)
        except ValueError:
            return 9999

    def analyze(self):
        """Print analysis of the circuit."""
        print(f"  Inputs: {len(self.inputs)}")
        print(f"  Outputs: {len(self.outputs)}")
        print()
        for out in self.outputs:
            deps = self.get_deps_for_output(out)
            print(f"  {out}: depends on {len(deps)} inputs → {', '.join(deps[:5])}{'...' if len(deps) > 5 else ''}")

def main():
    for name, path in [
        ('hd09', '/home/wangfz/bool/examples/hd09/hd09.txt'),
        ('hd11', '/home/wangfz/bool/examples/hd11/hd11.txt'),
        ('hd12', '/home/wangfz/bool/examples/hd12/hd12.txt'),
        ('router', '/home/wangfz/bool/examples/router/router.txt'),
    ]:
        print(f"\n{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")
        text = open(path).read()
        ca = CircuitAnalyzer(text)
        ca.analyze()

if __name__ == '__main__':
    main()
