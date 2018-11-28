

    # (1.1) extract the initial mappings between proteins and peptides
    for row_ix, row_values in rnp_df[['Proteins', 'Peptide']].iterrows():

        proteins = row_values['Proteins'].split(",")
        peptide = row_values['Peptide']

        if peptide in pep2pro:
            assert pep2allpro[peptide] == proteins, (
                "The same peptide is observed more than once with different proteins!")
        pep2allpro[peptide] = proteins
        for protein in proteins:
            initial_proteins.add(protein)
            pro2pep[protein].add(peptide)

            if protein.split("|")[0] == "sp":
                protein_level = 1
                top_level_proteins.add(protein)
            elif protein.split("|")[0] == "tr":
                protein_level = 2
            else:
                raise ValueError("Protein does not appear to be either"
                                 "SwissProt(sp) or TrEMBL(tr)")

            pep2pro[peptide][protein_level].add(protein)

    section_blocker = writeSectionHeader(logfile, "Initial file stats")
    logfile.write("# initial peptides: %i\n" % len(pep2pro))
    logfile.write("# initial proteins: %i\n" % len(pro2pep))
    logfile.write("# initial SwissProt proteins: %i\n" % len(top_level_proteins))
    logfile.write("# initial TrEMBL proteins: %i\n" % (
        len(pro2pep)-len(top_level_proteins)))
    logfile.write("%s\n\n" % section_blocker)
            

    # (1.2) find the peptides with only TrEMBL protein matches and
    # 'upgrade' these TrEMBL proteins to being equivalent to SwissProt
    tr_only_peptides = set([x for x in pep2pro.keys() if len(pep2pro[x][1])==0])

    set_upgraded = set()
    for peptide in tr_only_peptides:
        upgraded = pep2pro[peptide][2]
        set_upgraded.update(upgraded)
        top_level_proteins.update(upgraded)
        pep2pro[peptide][2] = pep2pro[peptide][2].difference(set(upgraded))
        pep2pro[peptide][1] = pep2pro[peptide][1].union(set(upgraded))

    section_blocker = writeSectionHeader(
        logfile, "Deciding which TrEMBL proteins to retain:")
    logfile.write("# peptides with only TrEMBL matches: %i\n" % (
        len(tr_only_peptides)))
    logfile.write("# TrEMBL proteins retained as no SwissProt matches for "
                  "peptide: %i\n" % (len(set_upgraded)))
    logfile.write("%s\n\n" % section_blocker)

    # (1.3) Use a parsimonious approach to identifty the minimum number
    # of proteins required to cover all the peptides:
    # Start from the protein(s) with the most peptides and mark these as covered.
    # Continue with remaining proteins in order of peptides per protein
    # until all peptides are covered
    retained_proteins = []
    peptides = copy.deepcopy(set(pep2pro.keys()))
    peptide_counts = {}

    tmppro2pep = copy.deepcopy(pro2pep)
    new_top_level_proteins = copy.deepcopy(top_level_proteins)
    new_pep2pro = collections.defaultdict(set)

    peptide_count = max(map(len, tmppro2pep.values()))

    section_blocker = writeSectionHeader(
        logfile, ("Parsimonious method to identify minimal set of proteins"
                  " to account for all peptides"))

    while True:
        # (1.3.1) If all peptides covered or the maximum peptides per
        # protein = 0, break.
        if len(peptides) == 0 or peptide_count == 0:
            logfile.write("All peptides are now accounted for\n")
            break

        peptide_count -= 1 

        top_proteins = set()
        top_score = 0
        #(1.3.2) Find the proteins with the highest number of peptide matches
        for protein in new_top_level_proteins:
            if len(tmppro2pep[protein]) == top_score:
                top_proteins.add(protein)
            elif len(tmppro2pep[protein]) > top_score:
                top_score = len(tmppro2pep[protein])
                top_proteins = set((protein,))

        logfile.write("%i remaining protein(s) with %i peptides\n" % (
            len(top_proteins), top_score))
        # (1.3.3) Remove the top proteins and the associated peptides
        for top_protein in top_proteins:

            new_top_level_proteins.remove(top_protein)

            retained_proteins.append(top_protein)

            for peptide in pro2pep[top_protein]:
                new_pep2pro[peptide].add(top_protein)
                if peptide in peptides:
                    peptides.remove(peptide)
                for protein in pep2pro[peptide][1]:
                    if protein == top_protein:
                        continue
                    if peptide in tmppro2pep[protein]:
                        tmppro2pep[protein].remove(peptide)

    logfile.write("\n%i proteins retained\n" % len(retained_proteins))
    #logfile.write("\n".join([",".join(map(str, (x, len(tmppro2pep[x]), len(pro2pep[x]))))
    #                         for x in new_top_level_proteins]))
    logfile.write("%i SwissProt proteins retained\n" % len(
        [x for x in retained_proteins if x.split("|")[0] == 'sp']))
    logfile.write("%i TrEMBL proteins retained\n" % len(
        [x for x in retained_proteins if x.split("|")[0] == 'tr']))
    logfile.write("\nNote: If not all SwissProt proteins were retained, this means\n"
                  "these proteins only included peptides which were observed\n"
                  "in other proteins which had a greater number of peptides\n")
    logfile.write("%s\n\n" % section_blocker)

    section_blocker = writeSectionHeader(logfile, "proteins per peptide:")
    counts = collections.Counter([len(x) for x in new_pep2pro.values()])
    sum_counts = sum(counts.values())
    for k, v in counts.items():
        logfile.write("%i peptide(s) (%.2f %%) have %i master protein(s)\n" % (
            v, (100 * v)/sum_counts, k))
    logfile.write("%s\n\n" % section_blocker)

    # Check all the peptides are covered
    assert set(pep2pro.keys()).difference(set(new_pep2pro.keys())) == set()


    # add the top protein and uniprot id annotations
    rnp_df['master_protein(s)'] = [";".join(new_pep2pro[protein]) for protein in rnp_df['Peptide']]
    rnp_df['master_uniprot_id(s)'] = [";".join([protein_id.split("|")[1] for protein_id in new_pep2pro[protein]])
                               for protein in rnp_df['Peptide']]

    # (1.4) Build dictionaries to map from protein id to protein
    # sequence and description using the fasta database
    crap_proteins = set()
    protein2description = {
        entry.title.split(" ")[0]: " ".join(entry.title.split(" ")[1:])
        for entry in FastaIterator(open(args['fasta_db']))}
    protein2seq = {
        entry.title.split(" ")[0]:entry.sequence
        for entry in FastaIterator(open(args['fasta_db']))}
    for entry in FastaIterator(open(args['fasta_crap_db'])):
        protein2seq[entry.title.split(" ")[0]] = entry.sequence
        protein2description[entry.title.split(" ")[0]] = entry.title.split(" ")[0]
        crap_proteins.add(entry.title.split(" ")[0])

    # (1.5) derive further annotations
    protein_lengths = []
    protein_descriptions = []
    crap_protein = []

    position_in_peptide = []
    position_in_protein = []

    motif_13 = []
    motif_15 = []
    motif_17 = []

    for ix, row in rnp_df.iterrows():
        peptide = row['Best localization(s)']
        proteins = row['master_protein(s)'].split(";")
        protein_lengths.append(";".join(map(str, [len(protein2seq[x]) for x in proteins])))
        protein_descriptions.append(";".join([protein2description[x] for x in proteins]))

        # (1.5.1) does peptide match a cRAP protein?
        crap = 0
        for protein in proteins:
            if protein in crap_proteins:
                crap = 1
                break
        crap_protein.append(crap)
			

    new_column_order = [
        "Best localization(s)",
        "RNA",
        "master_protein(s)",
        "master_uniprot_id(s)",
        'protein_description(s)',
        'protein_length(s)',
        'position_in_peptide',
        'position_in_protein(s)', 
        'window_13', 'window_15', 'window_17',
        'crap_protein',
        "Peptide",
        "Proteins"]

    new_column_order.extend([x for x in rnp_df.columns if x not in new_column_order])
    final_rnp_df = rnp_df[new_column_order]

    final_rnp_df.to_csv(args['outfile'], index=False, sep="\t")
    os.chmod(args['outfile'], 0o666)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
