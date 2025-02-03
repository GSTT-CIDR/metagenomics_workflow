from datetime import datetime
import os
import subprocess
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import Tk, Label
from PIL import Image, ImageTk
from tkinter import filedialog

# Function to create the entry fields in a new tab for each sample
def create_tab(index):
    tab = ttk.Frame(notebook)
    notebook.add(tab, text=f'Sample {index}')
    tab_entries = []
    for field in field_names:
        entry_frame = ttk.Frame(tab)
        label = ttk.Label(entry_frame, text=field, width=15, anchor='w')
        if field == "Barcode":
            # Create a dropdown for Barcode with options 1 to 12
            entry = ttk.Combobox(entry_frame, values=[str(x) for x in range(1, 13)], width=47, state="readonly")
            entry.bind('<<ComboboxSelected>>', lambda event, idx=index-1: barcode_selected(event, idx))
        else:
            entry = ttk.Entry(entry_frame, width=50)
        tab_entries.append(entry)
        label.pack(side=tk.LEFT)
        entry.pack(side=tk.LEFT)
        entry_frame.pack(padx=10, pady=5)
    entries.append(tab_entries)

# Function to handle barcode selection
def barcode_selected(event, index):
    selected_value = entries[index][field_names.index("Barcode")].get()
    # Disable the selected value in other barcode dropdowns
    for i, tab_entries in enumerate(entries):
        if i != index:
            barcode_entry = tab_entries[field_names.index("Barcode")]
            barcode_entry['values'] = [v for v in barcode_entry['values'] if v != selected_value]

# Function to choose a file
def choose_file(entry, default_path='/mnt/sample_sheets/'):
    filepath = filedialog.askopenfilename(initialdir=default_path)
    if filepath:
        entry.delete(0, tk.END)
        entry.insert(0, filepath)



# Function to launch the pipeline either from the TSV constructed byt the launcher or a pre existing one.
# Checks if a path has been given for an existing sample sheet


def launch_pipeline():
    # Determine whether the existing sample sheet path has been entered
    entry2_text = entry2.get()
    if entry2_text:
        filename = entry2.get()
        # Counting the number of threads for snakemeake 
        num_lines = sum(1 for line in open(filename, 'r'))
        core_count = int((num_lines - 1) * 2)
        command = f"bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate cmg && for t in {{0.5,1,2,16,24}}; do snakemake --directory /mnt --cores {core_count} -k --config time=$t samples=\"{filename}\" --latency-wait 15 ; done'"
        print(command)
        subprocess.run(command, shell=True)
    else:
        filename = os.path.join("/mnt/sample_sheets/",str(format(datetime.now().strftime("%Y_%m_%d@%H_%M_%S")) + ".tsv"))
        with open(filename, 'w') as f:
            f.write('\t'.join(field_names) + '\n')  # Write headers
            for tab_entries in entries:
                data = [entry.get() for entry in tab_entries]
                f.write('\t'.join(data) + '\n')  # Write data
        messagebox.showinfo("Success", f"Sample sheet saved to {filename} successfully. Starting pipeline...")
        #print(num_samples)
        core_count = int(num_samples.get() * 2)
        command = f"bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate cmg && for t in {{0.5,1,2,16,24}}; do snakemake --directory /mnt --cores {core_count} -k --config time=$t samples=\"{filename}\" --latency-wait 15 ; done'"
        print(command)
        subprocess.run(command, shell=True)

# Initialize the main window
root = tk.Tk()
root.title('CIDR Metagenomics Launcher')

# Add description and image
# Open the image file for header
image = Image.open("/mnt/ref/Template/CIDR_logo.jpeg")
# Resize the image
resized_image = image.resize((432, 136))
# Convert the image for Tkinter
tk_image = ImageTk.PhotoImage(resized_image)
label = Label(root, image=tk_image)
label.pack()

# Text section
description1 = tk.Label(root, text="Launcher for metagenomics workflow. Must be executed from the provided container. If repeating an existing run choose a pre-compiled samplesheet to load. This will overwrite the existing analysis. For a new experiment, fill out the fields below and Launch Pipeline.", wraplength=430)
description1.pack()
description3 = tk.Label(root, text="")
description3.pack()

# Existing sample sheet loader
frame2 = tk.Frame(root)
frame2.pack(padx=10)
label2 = tk.Label(frame2, text="Load existing sample sheet:")
label2.pack(side=tk.TOP)
entry2 = tk.Entry(frame2, width=40)
entry2.pack(side=tk.LEFT, padx=5)
button2 = tk.Button(frame2, text="Choose File", command=lambda: choose_file(entry2))
button2.pack(side=tk.LEFT)

# Main frame
frame1 = tk.Frame(root)
frame1.pack(pady=5)
label1 = tk.Label(frame1, text="Number of samples/barcodes:")
label1.pack(side=tk.LEFT)
# Dropdown for the number of samples
num_samples = tk.IntVar(value=1)
sample_dropdown = ttk.Combobox(frame1, textvariable=num_samples, values=list(range(1, 13)), state="readonly", width=3)
sample_dropdown.pack(padx=10, pady=5)
sample_dropdown.bind('<<ComboboxSelected>>', lambda event: update_tabs())
sample_dropdown.pack(side=tk.RIGHT)

# Notebook for tabs
notebook = ttk.Notebook(root)
notebook.pack(expand=True, fill='both', padx=10, pady=5)

# Entry list for each tab
entries = []
field_names = ["LabID", "Experiment", "SampleID", "Barcode", "SampleType", "PatientID", "Operator"]

# Function to update the number of tabs
def update_tabs():
    while len(entries) < num_samples.get():
        create_tab(len(entries) + 1)
    while len(entries) > num_samples.get():
        removed_entries = entries.pop()
        notebook.forget(notebook.tabs()[-1])
        # Re-enable the barcode value in other tabs
        removed_barcode = removed_entries[field_names.index("Barcode")].get()
        for entry in entries:
            barcode_entry = entry[field_names.index("Barcode")]
            if removed_barcode not in barcode_entry['values']:
                barcode_entry['values'] = list(barcode_entry['values']) + [removed_barcode]

# Create initial tab
create_tab(1)

# Save button
save_button = ttk.Button(root, text="Launch pipeline", command=launch_pipeline)
save_button.pack(pady=10)

# Text section
description2 = tk.Label(root, text="Visit https://github.com/GSTT-CIDR/metagenomics_container for help")
description2.pack()

# Start the main loop
root.mainloop()
