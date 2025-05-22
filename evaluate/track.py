import subprocess
import psutil
import time
import re
import os
import numpy as np


def track_memory_and_runtime(executable_path, argument):
    try:
        start_time = time.time()
        process = subprocess.Popen([executable_path, argument], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pid = process.pid

        # Monitor memory usage
        max_memory = 0
        ps_process = psutil.Process(pid)
        while process.poll() is None:  # While the process is running
            try:
                memory_info = ps_process.memory_info()
                max_memory = max(max_memory, memory_info.rss)  # Track peak memory usage
            except psutil.NoSuchProcess:
                break  # Process has terminated
            time.sleep(0.1)  # Poll every 100ms

        # Wait for the process to finish and capture output
        stdout, stderr = process.communicate()
        end_time = time.time()

        # Calculate runtime
        runtime = end_time - start_time

        # Convert memory usage to MB
        max_memory_mb = max_memory / (1024 * 1024)

        # print(f"Runtime: {runtime:.2f} seconds")
        # print(f"Peak Memory Usage: {max_memory_mb:.2f} MB")
        if stdout:
            output = stdout.decode()
            # print("Program Output:")
            # print(output)

            result = int(re.findall(r"Distinct kmers:\s*(\d+)", output)[0])

            return runtime, max_memory_mb, result

        if stderr:
            print("Program Errors:")
            print(stderr.decode())

    except FileNotFoundError:
        print(f"Error: The executable '{executable_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


def parse_matrix_output(program_output):
    # Split the output into lines
    lines = program_output.strip().split("\n")
    
    # Parse each line into a list of floats
    matrix = []
    for line in lines:
        row = [float(value) for value in line.split()]  # Split by whitespace and convert to float
        matrix.append(row)
    
    # Convert the list of lists into a NumPy array
    matrix_array = np.array(matrix)
    return matrix_array

    
def track_memory_and_runtime_similarity(executable_path, *args):
    try:
        start_time = time.time()
        if args:
            process = subprocess.Popen([executable_path, ''.join(str(arg) for arg in args)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            process = subprocess.Popen([executable_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pid = process.pid

        # Monitor memory usage
        max_memory = 0
        ps_process = psutil.Process(pid)
        while process.poll() is None:  # While the process is running
            try:
                memory_info = ps_process.memory_info()
                max_memory = max(max_memory, memory_info.rss)  # Track peak memory usage
            except psutil.NoSuchProcess:
                break  # Process has terminated
            time.sleep(0.1)  # Poll every 100ms

        # Wait for the process to finish and capture output
        stdout, stderr = process.communicate()
        end_time = time.time()

        # Calculate runtime
        runtime = end_time - start_time

        # Convert memory usage to MB
        max_memory_mb = max_memory / (1024 * 1024)

        # print(f"Runtime: {runtime:.2f} seconds")
        # print(f"Peak Memory Usage: {max_memory_mb:.2f} MB")
        if stdout:
            output = stdout.decode()
            # print("Program Output:")
            # print(output)

            # result = int(re.findall(r"Distinct kmers:\s*(\d+)", output)[0])
            result = parse_matrix_output(output)

            return runtime, max_memory_mb, result

        if stderr:
            print("Program Errors:")
            print(stderr.decode())

    except FileNotFoundError:
        print(f"Error: The executable '{executable_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


