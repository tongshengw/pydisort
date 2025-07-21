"""split cdisort.c into cu files"""
import re
import sys
import os


def process_file(input_file_path):
    output_dir = "cu_src"
    os.makedirs(output_dir, exist_ok=True)

    start_pattern = r"/(?!.*end of).*\(\).*"
    end_pattern = r"/.*end of.*\(\).*"
    preprocessor_pattern = r"#.*"

    try:
        with open(input_file_path, "r", encoding="utf-8") as file:
            lines = file.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if re.match(preprocessor_pattern, line):
                print(line)
                i += 1
                continue

            if re.match(start_pattern, line):
                processed_name = process_line_for_filename(line)

                if processed_name:
                    output_file_path = os.path.join(
                        output_dir, f"{processed_name}.cu"
                    )

                    content = []

                    while i < len(lines):
                        current_line = lines[i]
                        content.append(current_line)

                        if re.match(end_pattern, current_line):
                            end_fn_str = process_line_for_filename(
                                current_line
                            )
                            end_fn_name = re.search(r"endof(.*)", end_fn_str)
                            if (
                                end_fn_name
                                and end_fn_name.group(1) != processed_name
                            ):
                                print(
                                    "fn name end start not match",
                                    end_fn_name.group(1),
                                    processed_name,
                                )
                                exit(1)

                            break
                        i += 1

                    try:
                        with open(
                            output_file_path, "w", encoding="utf-8"
                        ) as output_file:
                            output_file.writelines(
                                """#include<configure.h>

#include<memorypool.h>
#include<cdisort.h>

DISPATCH_MACRO
"""
                            )
                            output_file.writelines(content)
                        print(f"{output_file_path}")
                    except Exception as e:
                        print(f"Error writing file {output_file_path}: {e}")

            i += 1

    except FileNotFoundError:
        print(f"Error: File '{input_file_path}' not found.")
    except Exception as e:
        print(f"Error processing file: {e}")


def process_line_for_filename(line):
    chars_to_remove = ["=", "/", "*", " ", "(", ")"]
    processed = line

    for char in chars_to_remove:
        processed = processed.replace(char, "")

    processed = processed.strip()

    processed = re.sub(r"[^\w\-]", "", processed)

    if not processed:
        processed = "unnamed"
    elif len(processed) > 100:
        processed = processed[:100]

    return processed


def main():
    if len(sys.argv) != 2:
        print("usage: <input_file>")
    input_file = sys.argv[1]

    if not input_file:
        print("No file path provided.")
        return

    if not os.path.exists(input_file):
        print(f"File '{input_file}' does not exist.")
        return

    process_file(input_file)


if __name__ == "__main__":
    main()
