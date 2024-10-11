import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
from datetime import datetime
import subprocess
import csv
import glob
from PIL import Image, ImageTk
import time
from tkinter import Tk, Label
from tkcalendar import DateEntry
import os
from tkinter import PhotoImage

# iBus daemon needs disabling due to performance issues
os.environ['XMODIFIERS'] = ''

# Sample integrity check functions
def check_file_exists(file_path):
    return os.path.isfile(file_path)

def check_string_validity(string):
    valid_strings = ["valid1", "valid2", "valid3"]
    return string in valid_strings

def check_numeric_value(value):
    try:
        num = float(value)
        return 0 <= num <= 100
    except ValueError:
        return False

def get_directories(path):
    """Returns a list of directories in the given path sorted by creation time."""
    try:
        # Retrieve the list of directories with their creation times
        directories = [(name, os.path.getctime(os.path.join(path, name))) for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
        
        # Sort the directories by creation time
        directories.sort(key=lambda x: x[1], reverse=True)
        
        # Return only the directory names in the sorted order
        return [name for name, _ in directories]
    except FileNotFoundError:
        return []

def get_subdirectories(path):
    """Returns a list of subdirectories in the given path."""
    try:
        return [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
    except FileNotFoundError:
        return []

def check_barcode_path_exists(selected_dir, barcode_number):
    """Check if the barcode path exists."""
    pattern = os.path.join(selected_dir, '*', '*', 'fastq_pass', f'barcode{barcode_number:02d}')
    matches = glob.glob(pattern)
    return len(matches) > 0

# Define the path to scan for directories
directory_path = "/data"  # Update this path to your actual directory
directories = get_directories(directory_path)

# Column configuration with the directory dropdown, tooltips, and barcode check
columns = [
    {
        "name": "MinKNOW experiment ID",
        "output_column_name": "Experiment",
        "check": False,
        "width": 30,
        "check_func": None,
        "dropdown": True,
        "values": directories,
        "tooltip": "The experiment ID assigned during setup in MinKNOW (auto-populated from mounted '/data' directory)."
    },
    {
        "name": "MinKNOW sample ID",
        "output_column_name": "SampleID",
        "check": False,
        "width": 30,
        "check_func": None,
        "dropdown": True,
        "values": [],
        "tooltip": "The 'Sample ID' assigned to the sequencing run during setup in MinKNOW."
    },
    {
        "name": "ONT barcode",
        "output_column_name": "Barcode",
        "check": True,
        "width": 5,
        "check_func": None,
        "dropdown": True,
        "values": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12a"],
        "tooltip": "Barcode number (1-12)"
    },
    {
        "name": "Lab/sample ID",
        "output_column_name": "LabID",
        "check": False,
        "width": 30,
        "check_func": None,
        "dropdown": False,
        "values": None,
        "tooltip": "A unique identifier of the lab specimen - Lab accession or Sample ID - can be anonymised if required"
    },
    {
        "name": "Anonymised identifier",
        "output_column_name": "AnonymisedIdentifier",
        "check": False,
        "width": 15,
        "check_func": None,
        "dropdown": False,
        "values": None,
        "tooltip": "Anonymised identifier or hospital number"
    },
    {
        "name": "Collection date (YYYY-MM-DD)",
        "output_column_name": "CollectionDate",
        "check": False,
        "width": 15,
        "check_func": None,
        "dropdown": False,
        "values": None,
        "tooltip": "Date of sample collection"
    },
    {
        "name": "Sample class",
        "output_column_name": "SampleClass",
        "check": False,
        "width": 18,
        "check_func": None,
        "dropdown": True,
        "values": [
            'specimen',
            'positive_control',
            'validation_material',
            'community_standard',
            'negative_control'
        ],
        "tooltip": "Select the category of sample (Additional classes can be added by typing in the field)"
    },
    {
        "name": "Sample type",
        "output_column_name": "SampleType",
        "check": False,
        "width": 8,
        "check_func": None,
        "dropdown": True,
        "values": [
            'other',
            'BAL',
            'SPT',
            'NDL',
            'ETT',
            'NTS',
            'NPA',
            'PFL',
            'TA'
        ],
        "tooltip": "BAL - Bronchoalveolar Lavage<br> SPT - Sputum<br> NDL - Nasal Swab<br> ETT - Endotracheal Tube<br> NTS - Nasopharyngeal Swab<br> NPA - Nasopharyngeal Aspirate<br> PFL - Pleural Fluid<br> TA - Tracheal Aspirate. (Additional classes can be added by typing in the field)"
    },
    {
        "name": "Operator",
        "output_column_name": "Operator",
        "check": False,
        "width": 8,
        "check_func": None,
        "dropdown": False,
        "values": None,
        "tooltip": "User initials"
    },
    {
        "name": "Notes",
        "output_column_name": "Notes",
        "check": False,
        "width": 30,
        "check_func": None,
        "dropdown": False,
        "values": None,
        "tooltip": "Additional notes to be added to results."
    }
]

class Tooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(500, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                         background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

class DataIngestionApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CIDR Metagenomics Launcher")
        self.geometry("2000x700")
        # Taskbar icon
        icon_path = '/mnt/lib/install/CIDR_logo_square_rmg.png'
        icon = PhotoImage(file=icon_path)
        self.iconphoto(True, icon)
        
        self.create_widgets()
    
    def create_widgets(self):
        self.num_rows = tk.IntVar(value=8)
        self.entries = []

        # Create a canvas
        self.canvas = tk.Canvas(self)
        self.canvas.grid(row=0, column=0, sticky="nsew")

        # Add horizontal and vertical scrollbars
        self.h_scrollbar = tk.Scrollbar(self, orient="horizontal", command=self.canvas.xview)
        self.h_scrollbar.grid(row=1, column=0, sticky="ew")
        self.v_scrollbar = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.v_scrollbar.grid(row=0, column=1, sticky="ns")

        # Configure the canvas to use the scrollbars
        self.canvas.configure(xscrollcommand=self.h_scrollbar.set, yscrollcommand=self.v_scrollbar.set)

        # Make the canvas expandable
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # Create a frame inside the canvas
        self.frame = tk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.frame, anchor="nw")

        # Update scroll region when the frame size changes
        self.frame.bind("<Configure>", self.on_frame_configure)

        # Add description and image
        # Open the image file for header
        image = Image.open("/mnt/db/ref/Template/CIDR_logo.jpeg")
        # Resize the image
        resized_image = image.resize((432, 136))
        # Convert the image for Tkinter
        self.tk_image = ImageTk.PhotoImage(resized_image)
        label = Label(self.frame, image=self.tk_image)
        label.grid(row=0, column=0, columnspan=len(columns)*2,  padx=10, pady=10, sticky="ew")

        # Create a label with the text you want to display
        link_text = "Click here for instructions on how to use the launcher"
        url = "https://gstt-cidr.github.io/network_hub/running_metagenomics_workflow/"
        # Use the label as a clickable link
        link = tk.Label(self.frame, text=link_text, fg="blue", cursor="hand2")
        link.grid(row=1, column=0, columnspan=len(columns)*2, padx=10, pady=10, sticky="ew")
        link.bind("<Button-1>", lambda e: self.open_link(url))

        # Existing sample sheet loader
        frame2 = tk.Frame(self.frame)
        frame2.grid(row=2, column=0, columnspan=len(columns)*2, padx=10, pady=10, sticky="ew")

        # Center align label, entry, and button
        label2 = tk.Label(frame2, text="Load existing sample sheet:")
        label2.pack(side=tk.TOP, pady=(0, 5))  # Add some padding for spacing

        # Create a sub-frame to hold the entry and button, and center it
        entry_button_frame = tk.Frame(frame2)
        entry_button_frame.pack(side=tk.TOP, anchor='center')

        self.entry2 = tk.Entry(entry_button_frame, width=40)
        self.entry2.pack(side=tk.LEFT, padx=5)

        button2 = tk.Button(entry_button_frame, text="Choose File", command=lambda: self.choose_file())
        button2.pack(side=tk.LEFT)
        
        # Add a separator
        separator = ttk.Separator(self.frame, orient='horizontal')
        separator.grid(row=3, column=0, columnspan=len(columns)*2, sticky='ew', pady=(10, 20))

        tk.Label(self.frame, text="Select number of rows:").grid(row=4, column=0, padx=10, pady=10, sticky="w")
        self.row_selector = ttk.Combobox(self.frame, textvariable=self.num_rows, values=list(range(1, 25)))
        self.row_selector.grid(row=4, column=1, padx=10, pady=10, sticky="w")
        self.row_selector.bind("<<ComboboxSelected>>", self.update_rows)

        # Add a new text entry for filename suffix
        tk.Label(self.frame, text="Filename Suffix:").grid(row=5, column=0, padx=10, pady=10, sticky="w")
        self.filename_suffix = tk.Entry(self.frame, width=30)
        self.filename_suffix.grid(row=5, column=1, padx=10, pady=10, sticky="w")

        self.entry_frame = tk.Frame(self.frame)
        self.entry_frame.grid(row=6, column=0, columnspan=len(columns)*2, padx=10, pady=10)

        self.submit_button = tk.Button(self.frame, text="Launch Pipeline", command=self.launch_pipeline)
        self.submit_button.grid(row=7, column=6, pady=20)

        self.create_entry_rows(self.num_rows.get())

        # Checkbox for force option
        self.force_var = tk.BooleanVar()
        self.force_check = tk.Checkbutton(self.frame, text="Force overwrite", variable=self.force_var)
        self.force_check.grid(row=7, column=0, pady=10)
        
        # Checkbox for force option
        self.mSCAPE_var = tk.BooleanVar()
        self.mSCAPE_check = tk.Checkbutton(self.frame, text="mSCAPE prompt on completion", variable=self.mSCAPE_var)
        self.mSCAPE_check.grid(row=7, column=1, pady=10)
        
        # Checkbox for sleep option
        #self.disable_ingest_delay_var = tk.BooleanVar()
        #self.disable_ingest_delay_check = tk.Checkbutton(self.frame, text="Disable data ingest sleep", variable=self.disable_ingest_delay_var)
        #self.disable_ingest_delay_check.grid(row=7, column=2, pady=10)
    

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def create_entry_rows(self, num_rows):
        for widget in self.entry_frame.winfo_children():
            widget.destroy()

        self.entries = []

        # Create header row with tooltips
        for col_idx, column in enumerate(columns):
            header = tk.Label(self.entry_frame, text=column["name"])
            header.grid(row=0, column=col_idx*2, padx=5, pady=5, sticky="w")
            Tooltip(header, column["tooltip"])

        for row in range(1, num_rows + 1):
            row_entries = []
            for col_idx, column in enumerate(columns):
                frame = tk.Frame(self.entry_frame)
                frame.grid(row=row, column=col_idx*2, padx=5, pady=5)

                if column["name"] == "Collection Date":
                    entry = DateEntry(frame, width=column["width"])
                    entry.set_date('')  # Set date to an empty string
                    entry.pack(side="left")
                elif column["dropdown"]:
                    entry = ttk.Combobox(frame, values=column["values"], width=column["width"])
                    entry.pack(side="left")
                    if column["name"] == "MinKNOW experiment ID":
                        entry.bind("<<ComboboxSelected>>", lambda e, row=row, col_idx=col_idx: self.update_subdirectory_dropdown(e, row, col_idx))
                    if column["name"] == "Barcode":
                        entry.bind("<<ComboboxSelected>>", lambda e, row=row, col_idx=col_idx: self.check_barcode_path(e, row, col_idx, full_path))
                else:
                    entry = tk.Entry(frame, width=column["width"])
                    entry.pack(side="left")

                if column["check"]:
                    entry.bind("<FocusOut>", lambda e, col_idx=col_idx: self.perform_checks(e, col_idx))
                row_entries.append(entry)

                if column["check"]:
                    indicator = tk.Label(frame, width=2, height=1, bg="yellow")
                    indicator.pack(side="left", padx=5)
                else:
                    indicator = None
                row_entries.append(indicator)

            self.entries.append(row_entries)

    def update_rows(self, event):
        num_rows = self.num_rows.get()
        self.create_entry_rows(num_rows)

    def perform_checks(self, event, col_idx):
        widget = event.widget if event else None
        value = widget.get() if widget else self.entries[0][col_idx * 2].get()
        column = columns[col_idx]

        if column["check"] and column["check_func"] is not None:
            indicator = widget.master.winfo_children()[1] if widget else self.entries[0][col_idx * 2 + 1]
            is_valid = column["check_func"](value)
            if is_valid:
                indicator.config(bg="green")
            else:
                indicator.config(bg="red")
    
    def choose_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("TSV files", "*.tsv"),("CSV files", "*.csv"), ("All files", "*.*")], initialdir='/mnt/sample_sheets/')
        if file_path:
            self.entry2.delete(0, tk.END)
            self.entry2.insert(0, file_path)
            self.load_tsv_data(file_path)

    def load_tsv_data(self, file_path):
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            headers = next(reader)
            data = list(reader)
        
        self.num_rows.set(len(data))
        self.create_entry_rows(len(data))

        for row_idx, row_data in enumerate(data):
            for col_idx, cell_data in enumerate(row_data):
                entry = self.entries[row_idx][col_idx*2]
                # Check if the column is "Collection Date" and set the value accordingly
                if columns[col_idx]["name"] == "Collection Date":
                    try:
                        entry.set_date(datetime.strptime(cell_data, "%d/%m/%Y").date())
                    except ValueError:
                        # Handle the error if the date format is incorrect
                        messagebox.showerror("Error", f"Incorrect date format in row {row_idx+1}, column {col_idx+1}. Expected format: dd-mm-yyyy.")
                else:
                    entry.insert(0, cell_data)
                # Perform integrity check if applicable
                if columns[col_idx]["check"]:
                    self.perform_checks(event=None, col_idx=col_idx)

        self.num_rows.set(len(data))
        self.create_entry_rows(len(data))

        for row_idx, row_data in enumerate(data):
            for col_idx, cell_data in enumerate(row_data):
                entry = self.entries[row_idx][col_idx*2]
                entry.insert(0, cell_data)
                # Perform integrity check if applicable
                if columns[col_idx]["check"]:
                    self.perform_checks(event=None, col_idx=col_idx)

    def any_red_indicators(self):
        for row_entries in self.entries:
            for idx in range(1, len(row_entries), 2):  # Check every indicator
                indicator = row_entries[idx]
                if indicator and indicator.cget("bg") == "red":
                    return True
        return False

    def update_subdirectory_dropdown(self, event, row, col_idx):
        selected_dir = self.entries[row - 1][col_idx * 2].get()  # Get the selected directory
        full_path = os.path.join(directory_path, selected_dir)
        print(f"Selected directory: {selected_dir}")  # Debug print
        print(f"Full path: {full_path}")  # Debug print
        subdirectories = get_subdirectories(full_path)
        print(f"Subdirectories: {subdirectories}")  # Debug print
        
        # Include the selected directory as the first option in subdirectory list
        subdirectories.insert(0, selected_dir)
        
        # Update subdirectory combobox
        subdirectory_combobox = self.entries[row - 1][(col_idx + 1) * 2]
        subdirectory_combobox['values'] = subdirectories
        subdirectory_combobox.set(subdirectories[0])  # Set the default value to the parent directory name

        # Bind the barcode dropdown selection to the check function
        barcode_combobox = self.entries[row - 1][(col_idx + 2) * 2]
        barcode_combobox.bind("<<ComboboxSelected>>", lambda e, row=row, col_idx=col_idx + 2: self.check_barcode_path(e, row, col_idx, full_path))

    def check_barcode_path(self, event, row, col_idx, selected_dir):
        barcode_number = self.entries[row - 1][col_idx * 2].get()  # Get the selected barcode number
        print(f"Checking barcode path for barcode number: {barcode_number}")  # Debug print

        # Perform the check
        path_exists = check_barcode_path_exists(selected_dir, int(barcode_number))
        print(f"Barcode path exists: {path_exists}")  # Debug print

        # Update the indicator
        indicator = self.entries[row - 1][col_idx * 2 + 1]
        if path_exists:
            indicator.config(bg="green")
        else:
            indicator.config(bg="red")

    def proceed_with_action(self):
        # Code to proceed with the action
        pass

    def launch_pipeline(self):
        if self.any_red_indicators():
            proceed = messagebox.askyesno("Warning", "There are red indicators in the entries. Do you want to proceed?")
            if not proceed:
                messagebox.showinfo("Aborted", "The action has been aborted.")
                return

        # Your code to proceed with the action if the user selects yes
        self.proceed_with_action()
        
        # Determine whether the existing sample sheet path has been entered
        entry2_text = self.entry2.get()
        # Force tickboxes and mSCAPE - delay skip is a boolean passed straight to the app
        force_option = '--force' if self.force_var.get() else ''
        mscape_command = 'python mSCAPE_launcher_0.3.py'

        # Added unbuffer to preserve the coloured outputs of snakemake
        filename_suffix = self.filename_suffix.get()
        filename = os.path.join("/mnt/sample_sheets/", str(format(datetime.now().strftime("%Y_%m_%d@%H_%M_%S")) + "_" + filename_suffix + ".tsv"))
        log_path = os.path.join("/mnt/logs/", str(format(datetime.now().strftime("%Y_%m_%d@%H_%M_%S")) + ".txt"))
        with open(filename, 'w') as f:
            f.write('\t'.join(column["output_column_name"] for column in columns) + '\n')  # Write headers
            for row_entries in self.entries:
                data = [entry.get() for entry in row_entries[::2]]  # Skip indicators
                f.write('\t'.join(data) + '\n')  # Write data
        messagebox.showinfo("Success", f"Sample sheet saved to {filename} successfully. Click OK to start.")
        core_count = 16
        command = f"bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate cmg && for t in {{0.5,1,2,16,24}}; do unbuffer snakemake --directory /mnt --cores {core_count} {force_option} -k --config time=$t samples=\"{filename}\" --latency-wait 15 2>&1 | tee -a {log_path} ; done '"
        subprocess.run(command, shell=True)
        
        # Running mSCAPE launcher after the pipeline has completed
        if self.mSCAPE_var.get():
            subprocess.run(mscape_command, shell=True)
            
    def open_link(self, url):
        import webbrowser
        webbrowser.open_new(url)

# Run the application
if __name__ == "__main__":
    app = DataIngestionApp()
    app.mainloop()
