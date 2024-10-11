import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import subprocess
import threading
import queue
import time

# Constants
RESULTS_DIR = "/mnt/results"
SAMPLE_LIST_PATH = "/mnt/sample_sheets/selected_summary_samples"
TIMEPOINTS = [0.5, 1, 2, 16, 24]
script_version = './metagenomics_summary_report_v1.4.py'
default_directory = "/mnt/summary_reports"
default_directory_2 = "/mnt/sample_sheets"

# Function to list directories in the results folder sorted by creation date
def list_directories():
    # Get the full path of each directory
    dir_paths = [os.path.join(RESULTS_DIR, d) for d in os.listdir(RESULTS_DIR) if os.path.isdir(os.path.join(RESULTS_DIR, d))]
    
    # Sort the directories by their creation time
    dir_paths.sort(key=lambda x: os.path.getctime(x))
    
    # Return the directory names sorted by creation time
    return [os.path.basename(d) for d in dir_paths]

# Function to check file existence for each timepoint
def check_files(sample, q):
    indicators = []
    for timepoint in TIMEPOINTS:
        timepoint_str = str(timepoint).replace('.', '_')
        file_path = os.path.join(
            RESULTS_DIR, sample, f"{timepoint_str}_hours", "centrifuge", "bacterial_centrifuge_report.tsv"
        )
        indicators.append(os.path.isfile(file_path))
    q.put((sample, indicators))

# Function to save selected directories to a file
def save_selected_directories():
    selected_dirs = [d for d in original_sample_list if check_vars[d].get()]
    with open(SAMPLE_LIST_PATH, 'w') as f:
        for dir_name in selected_dirs:
            f.write(f"{dir_name}\n")
    return SAMPLE_LIST_PATH

# Function to load sample list from a file
def load_sample_list():
    file_path = filedialog.askopenfilename(initialdir=default_directory_2, filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    if file_path:
        with open(file_path, 'r') as file:
            sample_dirs = [line.strip() for line in file.readlines()]
        update_sample_list(sample_dirs)
        update_checkboxes(sample_dirs)

# Function to update the sample list and checkboxes
def update_sample_list(sample_dirs):
    global original_sample_list
    original_sample_list = sample_dirs  # Update the original sample list
    # Re-create checkboxes for the new sample list
    for widget in scrollable_frame.winfo_children():
        widget.destroy()
    for dir_name in original_sample_list:
        var = tk.BooleanVar()
        frame = tk.Frame(scrollable_frame)
        tk.Checkbutton(frame, text=dir_name, variable=var).pack(side="left")
        check_vars[dir_name] = var
        # Since we don't have indicator data for imported samples, skip indicators
        frame.pack(anchor='w')

# Function to update checkboxes based on loaded sample list
def update_checkboxes(sample_dirs):
    for dir_name, var in check_vars.items():
        var.set(dir_name in sample_dirs)

# Function to browse output file
def browse_output_file():
    if xlsx_output_format.get():
        file_path = filedialog.asksaveasfilename(initialdir=default_directory, defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
    else:
        file_path = filedialog.asksaveasfilename(initialdir=default_directory, defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
    
    if file_path:
        output_file_path.set(file_path)

# Function to execute the script
def execute_script():
    sample_list_file = save_selected_directories()
    output_file = output_file_path.get()
    if not output_file:
        messagebox.showerror("Error", "Please select an output file.")
        return

    # Validate abundance threshold
    try:
        abundance_threshold = float(abundance_threshold_var.get())
    except ValueError:
        messagebox.showerror("Error", "Abundance threshold must be a number.")
        return
    abundance_threshold = str(abundance_threshold)

    command = [
        "python", script_version,
        "--sample_list", sample_list_file,
        "--results_directory", RESULTS_DIR,
        "--output_directory", output_file,
        "--as_xlsx", str(xlsx_output_format.get()),
        "--abundance_threshold", abundance_threshold,
        "--delimiter", str(delimiter_var.get())
    ]
    print("Parsed command:")
    print(command)
    try:
        subprocess.run(command, check=True)
        messagebox.showinfo("Success", "Script executed successfully!")
    except subprocess.CalledProcessError as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

# Tkinter UI setup
root = tk.Tk()
root.title("Metagenomics Summary Report Generator")

# Output file selection
tk.Label(root, text="Output File:").grid(row=0, column=0, padx=10, pady=5, sticky='e')
output_file_path = tk.StringVar()
tk.Entry(root, textvariable=output_file_path, width=50).grid(row=0, column=1, padx=10, pady=5)
tk.Button(root, text="Browse", command=browse_output_file).grid(row=0, column=2, padx=10, pady=5)

# Sample selection frame with scrollbar
frame = tk.Frame(root)
frame.grid(row=1, column=0, columnspan=3, padx=10, pady=10)

canvas = tk.Canvas(frame, width=600, height=400)
scrollbar = tk.Scrollbar(frame, orient="vertical", command=canvas.yview)
scrollable_frame = tk.Frame(canvas)

scrollable_frame.bind(
    "<Configure>",
    lambda e: canvas.configure(
        scrollregion=canvas.bbox("all")
    )
)

canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
canvas.configure(yscrollcommand=scrollbar.set)

scrollbar.pack(side="right", fill="y")
canvas.pack(side="left", fill="both", expand=True)

# Progress bar
progress = ttk.Progressbar(root, orient="horizontal", length=500, mode="determinate")
progress.grid(row=2, column=0, columnspan=3, pady=5)

check_vars = {}
original_sample_list = []

def load_samples():
    directories = list_directories()
    total_dirs = len(directories)
    global original_sample_list
    original_sample_list = directories  # Update the original sample list
    q = queue.Queue()
    threads = []
    progress['value'] = 0
    root.update_idletasks()
    for i, dir_name in enumerate(directories):
        t = threading.Thread(target=check_files, args=(dir_name, q))
        t.start()
        threads.append(t)
    completed = 0
    while completed < total_dirs:
        sample, indicators = q.get()
        var = tk.BooleanVar()
        frame = tk.Frame(scrollable_frame)
        tk.Checkbutton(frame, text=sample, variable=var).pack(side="left")
        for j, indicator in enumerate(indicators):
            color = "green" if indicator else "red"
            tk.Label(frame, text=f"T{TIMEPOINTS[j]}", bg=color, width=4).pack(side="left", padx=1)
        frame.pack(anchor='w')
        check_vars[sample] = var
        completed += 1
        progress['value'] = (completed / total_dirs) * 100
        root.update_idletasks()
    # Ensure all threads have completed
    for t in threads:
        t.join()

# Start loading samples in a separate thread to keep the UI responsive
threading.Thread(target=load_samples).start()

# Load sample list button
tk.Button(root, text="Load list of sample names", command=load_sample_list).grid(row=3, column=0, columnspan=3, pady=5)

# Relative abundance dropdown
tk.Label(root, text="Relative Abundance threshold:").grid(row=4, column=0, padx=10, pady=5, sticky='e')
abundance_threshold_var = tk.StringVar()
abundance_threshold_var.set("1.0")
abundance_threshold_options = [str(x/10) for x in range(1, 51)]
tk.OptionMenu(root, abundance_threshold_var, *abundance_threshold_options).grid(row=4, column=1, padx=10, pady=5, sticky='w')

# Output format checkbox
xlsx_output_format = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Output as XLSX", variable=xlsx_output_format).grid(row=5, column=0, padx=10, pady=5, sticky='w')

# Delimiter checkbox
delimiter_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Use ';' as delimiter", variable=delimiter_var).grid(row=5, column=1, padx=10, pady=5, sticky='w')

# Execute button
tk.Button(root, text="Generate summary", command=execute_script).grid(row=6, column=0, columnspan=3, pady=10)

# Run the Tkinter event loop
root.mainloop()
