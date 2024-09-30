#!/bin/bash

# Directory to process
DIRECTORY=$(pwd)

# Loop over all .c and .h files in the directory
for file in "$DIRECTORY"/*.c "$DIRECTORY"/*.h "$DIRECTORY"/*.S
do
    if [ -f "$file" ]; then
        echo "Processing $file..."
        python3 update_script.py "$file"
    fi
done

echo "Processing complete."

