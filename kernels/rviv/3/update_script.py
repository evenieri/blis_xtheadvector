import re
import sys

# List of specific RVV mnemonics to target
target_mnemonics = [
    'vsetvli', 'vxor', 'vfmacc', 'vfnmsac', 'vfmul', 
    'vfmsub', 'vfmadd', 'vlseg2w', 'vsseg2w', 'vle', 'vse', 
    'vlseg2e', 'vsseg2e', 'vlw', 'vsw'
]

def prepend_th_to_mnemonics(line):
    # Create a regex pattern to match any of the target mnemonics
    pattern = r'\b(' + '|'.join(re.escape(m) for m in target_mnemonics) + r')\b'
    
    # Replace the first occurrence of the target mnemonic with 'th.' prepended
    updated_line = re.sub(pattern, r'th.\1', line, count=1)
    
    return updated_line

def process_file(file_path):
    try:
        # Open the file and process it line by line
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Process each line to prepend 'th.' to only one target mnemonic
        updated_lines = [prepend_th_to_mnemonics(line) for line in lines]

        # Write the updated lines back to the file
        with open(file_path, 'w') as file:
            file.writelines(updated_lines)

        print(f"Processed file: {file_path}")

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_file>")
    else:
        file_path = sys.argv[1]
        process_file(file_path)

