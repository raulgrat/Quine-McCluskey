class QuineMcCluskey:
    def __init__(self, minterms, dont_cares, num_vars, inputs, outputs, numInputs, numOutputs):
        self.minterms = minterms
        self.dont_cares = dont_cares
        self.num_vars = num_vars
        self.inputs = inputs
        self.outputNames = outputs
        self.numInputs = numInputs
        self.numOutputs = numOutputs
        self.prime_implicants = set()
        self.covering_table = {}

    def get_binary_representation(self, num):
        return format(num, f'0{self.num_vars}b')

    def find_prime_implicants(self):
        # Combine minterms and don't cares for initial processing
        all_terms = set(self.minterms + self.dont_cares)
        groups = {i: [] for i in range(self.num_vars + 1)}
        
        # Group terms by the number of 1s in their binary representation
        for term in all_terms:
            binary = self.get_binary_representation(term)
            groups[binary.count('1')].append(binary)

        combined_any = True
        prime_implicant_candidates = set()

        while combined_any:
            combined_any = False
            new_groups = {i: [] for i in range(self.num_vars + 1)}
            marked = set()

            # Try combining terms in adjacent groups
            for i in range(len(groups) - 1):
                for term1 in groups[i]:
                    for term2 in groups[i + 1]:
                        if self.can_combine(term1, term2):
                            new_term = self.combine_terms(term1, term2)
                            if new_term not in new_groups[new_term.count('1')]:
                                new_groups[new_term.count('1')].append(new_term)
                            marked.update([term1, term2])
                            combined_any = True

            # Add non-marked terms to prime implicant candidates
            for group in groups.values():
                for term in group:
                    if term not in marked:
                        prime_implicant_candidates.add(term)

            groups = new_groups

        # Store the prime implicants
        self.prime_implicants.update(prime_implicant_candidates)
        self.build_covering_table()

    def can_combine(self, term1, term2):
        # Check if exactly one bit differs
        return sum(c1 != c2 for c1, c2 in zip(term1, term2)) == 1

    def combine_terms(self, term1, term2):
        # Combine two terms by replacing differing bit with '-'
        return ''.join(c1 if c1 == c2 else '-' for c1, c2 in zip(term1, term2))

    def build_covering_table(self):
        # Covering table tracks which implicants cover which minterms
        for term in self.prime_implicants:
            self.covering_table[term] = set()

        for minterm in self.minterms:
            binary = self.get_binary_representation(minterm)
            for term in self.prime_implicants:
                if self.covers(term, binary):
                    self.covering_table[term].add(minterm)

        essential_prime_implicants, covered_minterms = self.identify_essential_prime_implicants()
        self.prime_implicants = list(essential_prime_implicants)

        # Cover remaining minterms
        self.cover_remaining_minterms(covered_minterms)

    def covers(self, term, minterm):
        # A term covers a minterm if it matches on all non-dash positions
        return all(t == '-' or t == m for t, m in zip(term, minterm))

    def identify_essential_prime_implicants(self):
        essential_prime_implicants = set()
        covered_minterms = set()

        for minterm in self.minterms:
            covering_implicants = [term for term in self.prime_implicants if minterm in self.covering_table[term]]
            if len(covering_implicants) == 1:
                essential_prime_implicants.add(covering_implicants[0])
                covered_minterms.update(self.covering_table[covering_implicants[0]])

        return essential_prime_implicants, covered_minterms

    def cover_remaining_minterms(self, covered_minterms):
        # Greedily cover the remaining minterms using the best available implicants
        uncovered_minterms = set(self.minterms) - covered_minterms

        while uncovered_minterms:
            best_term = max(self.covering_table, key=lambda term: len(self.covering_table[term].intersection(uncovered_minterms)))
            self.prime_implicants.append(best_term)
            covered_minterms.update(self.covering_table[best_term])
            uncovered_minterms = set(self.minterms) - covered_minterms

    def output_pla(self):
        output = []
        output.append(f".i {self.numInputs}\n")
        output.append(f".o {self.numOutputs}\n")
        output.append(f".ilb {' '.join(self.inputs)}\n")
        output.append(f".ob {' '.join(self.outputNames)}\n")
        output.append(f".p {len(self.prime_implicants)}\n")

        for term in sorted(self.prime_implicants):
            output.append(f"{term} 1\n")

        output.append(".e\n")
        return ''.join(output)

# Parse PLA file
def parse_pla(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    variable_names = []
    output_names = []
    minterms = []
    dont_cares = []

    for line in lines:
        line = line.strip()

        if line.startswith('.i '):
            num_inputs = int(line.split()[1])

        elif line.startswith('.o '):
            num_outputs = int(line.split()[1])

        elif line.startswith('.ilb'):
            variable_names = line.split()[1:]

        elif line.startswith('.ob'):
            output_names = line.split()[1:]

        elif line[0] in '01':
            parts = line.split()
            binary_minterm = parts[0]
            output_value = parts[1]

            if output_value == '1':
                minterms.append(int(binary_minterm, 2))
            elif output_value == '-':
                dont_cares.append(int(binary_minterm, 2))

        elif line.startswith('.e'):
            break

    return minterms, dont_cares, len(variable_names), variable_names, output_names, num_inputs, num_outputs

# Main function
def main(in_pla, out_pla):
    minterms, dont_cares, num_vars, inputs, outputs, numIn, numOut = parse_pla(in_pla)

    if not minterms:
        print("No minterms found in the PLA file.")
        return

    qm = QuineMcCluskey(minterms, dont_cares, num_vars, inputs, outputs, numIn, numOut)
    qm.find_prime_implicants()

    with open(out_pla, 'w') as f:
        f.write(qm.output_pla())

# Entry point
if __name__ == "__main__":
    # Loop through multiple input files
    for i in range(1, 4):
        in_pla = f"inputs/input{i}.pla"
        out_pla = f"outputs/output{i}.pla"

        # Call the main function with each input-output pair
        main(in_pla, out_pla)
