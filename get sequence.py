import pandas as pd
import os
import csv
import re
from tqdm import tqdm

def extract_count(original_sequences, most_os):

    sequences = original_sequences.split(', ')
    for seq in sequences:
        if most_os in seq:
            count_match = re.search(r'\((\d+)', seq)
            if count_match:
                return int(count_match.group(1))
    return None

def generate_combinations(target_length, sorted_sequences, unit_length=3):

    N = target_length // unit_length
    if N == 0:
        return

    total_sequences = len(sorted_sequences)
    for i in range(total_sequences - N + 1):
        combo = sorted_sequences[i:i + N]
        concatenated = ''.join([str(seq['Most of OS in Anagrams']) for seq in combo])
        if len(concatenated) == target_length:
            max_rank = max([seq['Global Rank'] for seq in combo])
            sequence_name = f"{unit_length}_{target_length}_{max_rank}_{N}"
            yield f"{sequence_name}: {concatenated}"

def generate_variable_length_sequences(input_excel_file, output_folder, min_length=3, max_length=12, unit_length=3,
                                       count_threshold=0):

    # 读取Excel文件
    try:
        data = pd.read_excel(input_excel_file)
    except FileNotFoundError:
        print(f"not found: {input_excel_file}")
        return
    except Exception as e:
        print(f"error: {e}")
        return


    required_columns = ['Most of OS in Anagrams', 'Count', 'Original Sequences (count, average position)']
    try:
        data = data[required_columns]
    except KeyError as e:
        print(f"未找到列：{e}")
        return


    data['Most of OS Count'] = data.apply(
        lambda row: extract_count(row['Original Sequences (count, average position)'], row['Most of OS in Anagrams']),
        axis=1
    )


    data_filtered = data[data['Count'] > count_threshold]


    data_filtered = data_filtered[data_filtered['Most of OS in Anagrams'].str.len() == unit_length]


    data_sorted = data_filtered.sort_values(
        by=['Count', 'Most of OS Count', 'Most of OS in Anagrams'],
        ascending=[False, False, True]
    ).reset_index(drop=True)


    data_sorted['Global Rank'] = range(1, len(data_sorted) + 1)


    sorted_sequences = data_sorted.to_dict('records')

    print(f"get {len(sorted_sequences)} sequence。")


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    fasta_file_path = os.path.join(output_folder, "Generated_Peptides.fasta")
    total_combinations = 0

    with open(fasta_file_path, 'w', encoding='utf-8') as fasta_file:
        for target_length in range(min_length, max_length + 1, unit_length):
            N = target_length // unit_length
            if N == 0:
                continue
            print(f"combine {target_length} ，need {N} sequence")


            combination_generator = generate_combinations(target_length, sorted_sequences, unit_length)
            count_written = 0

            for combination in tqdm(combination_generator, desc=f"Generating combinations for Length {target_length}"):
                try:

                    parts = combination.split(': ')
                    if len(parts) != 2:
                        raise ValueError(f"error: {combination}")
                    name_part = parts[0]
                    sequence_part = parts[1]
                    unit_length_val, target_length_val, max_rank_val, N_val = name_part.split('_')


                    fasta_header = f">Combination_{total_combinations + 1}_Unit{unit_length_val}_Len{target_length_val}_MaxRank{max_rank_val}_N{N_val}"
                    fasta_file.write(f"{fasta_header}\n{sequence_part}\n")

                    count_written += 1
                    total_combinations += 1
                except Exception as e:
                    print(f"combine error: {combination}. error: {e}")
                    continue

            print(f"get {count_written} length of {target_length} to {fasta_file_path}")

    print(f"all sequence have saved into {fasta_file_path}")
    print(f"get {total_combinations} ")

if __name__ == "__main__":
    input_excel = ''
    output_folder = ''

    generate_variable_length_sequences(
        input_excel,
        output_folder,
        min_length= ,
        max_length= ,
        unit_length= ,
        count_threshold=
    )
