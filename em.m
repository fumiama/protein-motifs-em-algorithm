% Parameters
motif_len = 10;
% Read proteins data
data = readcell('proteins.csv', 'Delimiter', ',');
reduced_data = data(2:end, 2:end);
ascii_data = cellfun(@(x) int8(x), reduced_data);
[data_len, protein_len] = size(ascii_data);
% Randomly initialize the starting positions for all the sequences
index = randi([1 protein_len-motif_len], 1, data_len);
% Out-of-loop parameters
motifs = int8(zeros(data_len, motif_len));
amino_acids = unique(ascii_data);
amino_acids_rev = int8(zeros(1, 256));
for i = 1:length(amino_acids)
    amino_acids_rev(amino_acids(i)) = i;
end
expected_socres = zeros(data_len, protein_len-motif_len+1);
loop_count = 0;
while 1 % infinite loop
    loop_count = loop_count+1;
    % Extract the motif from each protein
    for i = 1:data_len
        p = index(i);
        motifs(i, :) = ascii_data(i, p:p+motif_len-1);
    end
    % Calculate PSSM
    score_matrix = pssm(amino_acids, motifs);
    % Expect the score of each position
    expected_socres_len = protein_len-motif_len+1;
    for i = 1:data_len
        for j = 1:expected_socres_len
            expected_socres(i, j) = expect(score_matrix, amino_acids_rev, ascii_data(i, j:j+motif_len-1));
        end
    end
    % Choose the positions that has the highest probability
    [max_expected_socres, max_expected_socres_indices] = max(expected_socres, [], 2);
    % Let it as the updated starting position
    next_index = max_expected_socres_indices';
    % Judge if the solution is stable
    if isequal(index, next_index)
        break;
    end
    index = next_index;
end
% Reached final solution
fprintf('Loop count: %d\n', loop_count);
disp('Starting positions:');
disp(index);
disp('PSSM of the motif:');
disp(score_matrix);
disp('Motifs:');
disp(char(motifs));
