from collections import defaultdict
from openpyxl import Workbook
from Bio import SeqIO


def is_anagram(s1, s2):
    return sorted(s1) == sorted(s2)


def find_anagrams(identical_sequences):
    anagrams = defaultdict(set)
    for sequence, sequence_ids_positions in identical_sequences.items():
        sorted_anagram_seq = ''.join(sorted(sequence))
        anagrams[sorted_anagram_seq].update(sequence_ids_positions)

    anagram_sequences = {k: v for k, v in anagrams.items() if len(v) > 1}
    return anagram_sequences


def save_anagrams_to_excel(anagrams, identical_sequences, output_file, identical_sequences_counts):
    wb_out = Workbook()
    ws_out = wb_out.active
    ws_out.append([
        'Anagram', 'IDs (Avg Position)', 'Count',
        'Original Sequences (count, average position)', 'in X% aa',
        'Most of OS in Anagrams', 'Most of OS Avg Position',
        'Second Most OS', 'Second Most OS Avg Position',
        'Third Most OS', 'Third Most OS Avg Position'
    ])

    for anagram, ids_positions in anagrams.items():
        count = len(set(id for id, pos in ids_positions))

        original_sequences_info = []
        sequence_counts = defaultdict(int)
        sequence_avg_positions = {}

        for seq in identical_sequences:
            if is_anagram(seq, anagram):
                if seq in identical_sequences_counts:
                    av_pos = sum(pos for k2, pos in identical_sequences[seq]) / len(identical_sequences[seq])
                    original_sequences_info.append(f"{seq} ({identical_sequences_counts[seq]:.0f}, {av_pos:.2f}%)")
                    sequence_counts[seq] += identical_sequences_counts[seq]
                    sequence_avg_positions[seq] = av_pos
                else:
                    original_sequences_info.append(f"{seq} (N/A)")

        # Sort the sequences based on counts to find the most, second, and third most common
        sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)

        most_common_seq = sorted_sequences[0][0] if sorted_sequences else 'N/A'
        most_common_seq_avg = sequence_avg_positions.get(most_common_seq, 'N/A')

        second_most_seq = sorted_sequences[1][0] if len(sorted_sequences) > 1 else 'N/A'
        second_most_seq_avg = sequence_avg_positions.get(second_most_seq, 'N/A')

        third_most_seq = sorted_sequences[2][0] if len(sorted_sequences) > 2 else 'N/A'
        third_most_seq_avg = sequence_avg_positions.get(third_most_seq, 'N/A')

        original_seq_with_count = f"{', '.join(original_sequences_info)}"
        avg_position = sum(pos for _id, pos in ids_positions) / count
        formatted_ids_positions = ", ".join([f"{_id} (in {pos:.2f}% aa)" for _id, pos in ids_positions])

        # Append all information to the worksheet
        ws_out.append([
            ' '.join(list(anagram)), formatted_ids_positions, count,
            original_seq_with_count, f'{avg_position:.2f}%',
            most_common_seq, most_common_seq_avg,
            second_most_seq, second_most_seq_avg,
            third_most_seq, third_most_seq_avg
        ])

    wb_out.save(output_file)


def find_sequences(input_file, min_length=9, max_length=9):
    sequences = defaultdict(lambda: defaultdict(list))
    total_sequences = 0
    with open(input_file, 'r', encoding='utf-8') as file:
        for record in SeqIO.parse(file, 'fasta'):
            total_sequences += 1
            sequence = str(record.seq)
            for i in range(len(sequence) - min_length + 1):
                for length in range(min_length, max_length + 1):
                    if i + length <= len(sequence):
                        subseq = sequence[i: i + length]
                        sequences[subseq][record.id].append(i)
    identical_sequences = {k: {(k2, sum(v2) / len(v2)) for k2, v2 in v.items()} for k, v in sequences.items() if
                           len(v) > 1}
    return identical_sequences, total_sequences


def save_sequences_to_excel(identical_sequences, total_sequences, output_file):
    wb = Workbook()
    ws = wb.active
    ws.append(['Sequence', 'IDs (Avg Position)', 'Percentage', 'Count',
               'Average Position'])  # Add the 'Average Position' title

    identical_sequences_counts = {}

    for sequence, ids_positions in identical_sequences.items():
        count = len(set(_id for _id, pos in ids_positions))
        identical_sequences_counts[sequence] = count
        percentage = (count / total_sequences) * 100
        avg_position = sum(pos for _id, pos in ids_positions) / count
        formatted_ids_positions = ', '.join([f"{_id} (in {pos:.2f}% aa)" for _id, pos in ids_positions])
        ws.append([sequence, formatted_ids_positions, f'{percentage:.2f}%', count, f'{avg_position:.2f}%'])
    wb.save(output_file)

    return identical_sequences_counts


if __name__ == '__main__':
    input_file_fasta = 'D:/合成肽/all APD3 AMP.fasta'
    output_file_sequences = 'D://all APD3 AMP_3_se.xlsx'
    identical_sequences, total_sequences = find_sequences(input_file_fasta)
    identical_sequences_counts = save_sequences_to_excel(identical_sequences, total_sequences, output_file_sequences)

    output_file_anagrams = 'D:/all APD3 AMP_9.xlsx'
    anagrams = find_anagrams(identical_sequences)
    save_anagrams_to_excel(anagrams, identical_sequences, output_file_anagrams, identical_sequences_counts)
