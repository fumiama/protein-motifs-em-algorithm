function score_matrix = pssm(amino_acids, motifs)
    protein_len = size(motifs, 2);
    amino_len = length(amino_acids);
    score_matrix = zeros(amino_len, protein_len);
    for i = 1:protein_len
        for j = 1:amino_len
            score_matrix(j, i) = sum(motifs(:, i) == amino_acids(j));
        end
    end
    score_matrix = score_matrix ./ sum(score_matrix, 2);
    % define log10(0) = -10
    score_matrix(score_matrix == 0) = 1e-10;
    % take the base-10 logarithm
    score_matrix = log10(score_matrix);
end