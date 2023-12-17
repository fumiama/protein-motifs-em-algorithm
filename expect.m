function protein_score = expect(score_matrix, amino_acids_rev, protein)
    protein_len = length(protein);
    protein_score = 0;
    for i = 1:protein_len
        protein_score = protein_score + score_matrix(amino_acids_rev(protein(i)), i);
    end
    protein_score = protein_score / protein_len;
end