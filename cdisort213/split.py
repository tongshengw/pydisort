"""split cdisort.c into cu files"""
import re
import sys
import os


def process_file(input_file_path):
    output_dir = "cu_src"
    os.makedirs(output_dir, exist_ok=True)

    start_pattern = r"/(?!.*end of).*\(\).*"
    end_pattern = r"/.*end of.*\(\).*"
    preprocessor_pattern = r"#.*def.*"
    fprintf_pattern = r".*fprintf.*\(.*\).*"
    fprintf_multiline_pattern = r".*fprintf.*\(.*"
    malloc_pattern = r".*malloc.*"
    calloc_pattern = r".*calloc.*"
    exit_pattern = r".*exit.*"

    try:
        with open(input_file_path, "r", encoding="utf-8") as file:
            lines = file.readlines()

        i = 0
        preprocessor_strings = []
        while i < len(lines):
            line = lines[i].strip()

            if re.match(preprocessor_pattern, line):
                print(line)
                preprocessor_strings.append(line + "\n")
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

                        if re.match(fprintf_pattern, current_line):
                            processed_line = convert_fprintf_to_printf(
                                current_line
                            )
                            content.append(processed_line)
                        elif re.match(fprintf_multiline_pattern, current_line):
                            processed_line = convert_fprintf_partial_line(
                                current_line
                            )
                            content.append(processed_line)
                        elif re.match(malloc_pattern, current_line):
                            processed_line = convert_malloc_to_swappablemalloc(
                                current_line
                            )
                            content.append(processed_line)
                        elif re.match(calloc_pattern, current_line):
                            processed_line = convert_calloc_to_swappablecalloc(
                                current_line
                            )
                            content.append(processed_line)
                        elif re.match(exit_pattern, current_line):
                            processed_line = convert_exit_to_trap(current_line)
                            content.append(processed_line)
                        else:
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
                                [
                                    """// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
"""
                                ]
                                + preprocessor_strings
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


def convert_fprintf_to_printf(code):
    pattern = r"fprintf\s*\(\s*[^,]+\s*,\s*(.*?)\s*\)"

    def replace_fprintf(match):
        # Extract everything after the first comma (format string and args)
        args = match.group(1)
        return f"printf({args})"

    # Use re.sub to replace all fprintf occurrences
    converted_code = re.sub(pattern, replace_fprintf, code)
    return converted_code


def convert_fprintf_partial_line(line):
    pattern = r"fprintf\s*\(\s*([^,]*)\s*,\s*(.*)"

    match = re.search(pattern, line)
    if match:
        remaining_args = match.group(2)

        fprintf_start = line.find("fprintf")

        before_fprintf = line[:fprintf_start]
        converted_line = before_fprintf + "printf(" + remaining_args

        return converted_line

    return line


def convert_malloc_to_swappablemalloc(code):
    return code.replace("malloc", "swappablemalloc")


def convert_calloc_to_swappablecalloc(code):
    return code.replace("calloc", "swappablecalloc")


def convert_exit_to_trap(code):
    return code.replace("exit(1)", "__trap()")


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
