import sys

if len(sys.argv) != 3:
    raise ValueError('Wrong input. Expected usage: correct_sample_file.py input_file_path output_file_path')

sample_file_path = sys.argv[1]
output_filename = sys.argv[2]

with open(sample_file_path, "r") as sample_file:
	with open(output_filename, "w") as corrected_sample_file:
		for i, line in enumerate(sample_file):
			if i >=2: #Skip first two lines
				split_line_list = line.split(" ")
				corrected_line_list = split_line_list
				corrected_line_list[:2] = [split_line_list[1], split_line_list[0]] #Swap IDs
				corrected_line_list[3] = 'NA'
				corrected_line_list[4] = 'NA'
				corrected_line = ' '.join(corrected_line_list)
				corrected_sample_file.write(corrected_line)
			else: #Copy header
				corrected_sample_file.write(line)
		
		
