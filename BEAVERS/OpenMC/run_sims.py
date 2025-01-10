import os
import subprocess

def run_in_directories(directories, python_script, system_program):
    """
    Runs a Python script and a system program in each specified directory.

    Parameters:
    directories (list): List of directory paths.
    python_script (str): Path to the Python script to execute.
    system_program (str): Name or path of the system program to execute.
    """
    for directory in directories:
        # Save the current working directory
        original_cwd = os.getcwd()
        try:
            # Change to the target directory
            os.chdir(directory)
            print(f"Changed to directory: {directory}")

            # Run the Python script
            print(f"Running Python script: {python_script}")
            subprocess.run(["python", python_script], check=True)

            # Run the system program
            print(f"Running system program: {system_program}")
            subprocess.run(system_program, shell=True, check=True)

        except Exception as e:
            print(f"An error occurred in directory {directory}: {e}")
        finally:
            # Return to the original directory
            os.chdir(original_cwd)
            print(f"Returned to directory: {original_cwd}")

if __name__ == "__main__":
    # Directories to iterate over
    directories_to_visit = ["Standard", "OTF"]

    # Python script to execute (replace with the actual script path)
    python_script_to_run = "gen_model.py"

    # System program to execute (replace with the actual program name or path)
    system_program_to_run = "openmc > run_log.log"  # Example system program

    run_in_directories(directories_to_visit, python_script_to_run, system_program_to_run)