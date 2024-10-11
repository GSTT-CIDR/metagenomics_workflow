import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import yaml
import os
import glob
import keyring
from Crypto.Cipher import AES
from Crypto.Util.Padding import pad
import base64
import gzip
import subprocess
import shutil
from tkinter import PhotoImage

# Path to the YAML configuration file containing terminology dictionaries and settings
config_path = '/mnt/configs/mscape_config.yaml'

# Define the main application class for the mSCAPE Uploader GUI
class MetadataApp:
    def __init__(self, root):
        self.root = root
        self.root.title("mSCAPE Uploader")  # Set the window title
        self.root.geometry("1200x1080")  # Adjust the window size

        # Load and set the taskbar icon
        icon_path = '/mnt/lib/install/CIDR_logo_square_rmg.png'  # Path to the icon image
        icon = PhotoImage(file=icon_path)
        self.root.iconphoto(True, icon)  # Set the window icon

        self.load_config()  # Load configuration settings from the YAML file
        self.create_widgets()  # Create the GUI widgets

    def load_config(self):
        # Load configuration settings from the YAML file
        with open(config_path, 'r') as file:
            self.config = yaml.safe_load(file)

        # Extract configuration parameters
        self.directory_path = self.config['directory_path']  # Directory path for file dialogs
        self.dropdown_options = self.config['dropdown_options']  # Options for dropdown menus
        self.columns = self.config['columns']  # Column definitions for the data table
        self.timepoints = self.config['timepoints']  # Timepoint options for file selection
        self.encryption_columns = self.config.get('encryption_columns', [])  # Columns to be encrypted
        self.initialization_vector = base64.b64decode(self.config.get('initialization_vector', ''))  # IV for encryption
        self.s3_bucket = self.config.get('s3_bucket', '')  # S3 bucket name
        self.aws_profile = self.config.get('aws_profile', '')  # AWS profile name
        self.aws_endpoint = self.config.get('aws_endpoint', '')  # AWS endpoint URL

    def create_widgets(self):
        # Dropdown Section
        self.dropdown_frame = ttk.LabelFrame(self.root, text="Dropdown Options")
        self.dropdown_frame.pack(pady=10, padx=10, fill="x")

        self.dropdown_vars = {}  # Dictionary to hold dropdown variables
        self.linked_vars = {}    # Dictionary to hold linked dropdown variables

        row = 0  # Initialize row counter
        col = 0  # Initialize column counter
        for column in self.columns:
            if column['source'] == 'dropdown' and column.get('show', True):
                if col >= 3:
                    # Move to the next row after 3 columns
                    col = 0
                    row += 1

                # Create a frame for each dropdown menu
                frame = tk.Frame(self.dropdown_frame)
                frame.grid(row=row, column=col, padx=10, pady=5, sticky="w")

                # Create a label for the dropdown menu
                label = tk.Label(frame, text=column['header'])
                label.pack(side=tk.LEFT)

                # Get the key for the dropdown options
                dropdown_key = column['dropdown_key']

                # Create a StringVar to hold the selected value
                dropdown = tk.StringVar(self.root)
                # Set default value to the first option
                dropdown.set(self.dropdown_options[dropdown_key]['values'][0])

                # Create the OptionMenu (dropdown)
                option_menu = tk.OptionMenu(frame, dropdown, *self.dropdown_options[dropdown_key]['values'])
                option_menu.pack(side=tk.LEFT)

                # Store the dropdown variable
                self.dropdown_vars[column['header']] = dropdown

                # Link the dropdown to any linked columns
                if 'linked_values' in self.dropdown_options[dropdown_key]:
                    self.linked_vars[column['header']] = self.dropdown_options[dropdown_key]['linked_values']

                col += 1

            elif column['source'] == 'linked-dropdown' and column.get('show', True):
                # For linked dropdowns, store the linkage information
                linked_to = column['linked_to']
                self.linked_vars[column['header']] = linked_to

        # **Main Table Section**
        self.table_frame = ttk.LabelFrame(self.root, text="Main Table")
        self.table_frame.pack(pady=10, padx=10, expand=True, fill=tk.BOTH)

        # Use a Canvas to enable scrolling
        self.table_canvas = tk.Canvas(self.table_frame)
        self.table_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Vertical scrollbar for the table
        self.scrollbar_y = tk.Scrollbar(self.table_frame, orient=tk.VERTICAL, command=self.table_canvas.yview)
        self.scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)

        # Horizontal scrollbar for the table
        self.scrollbar_x_main = tk.Scrollbar(self.table_frame, orient=tk.HORIZONTAL, command=self.table_canvas.xview)
        self.scrollbar_x_main.pack(side=tk.BOTTOM, fill=tk.X)

        # Configure the canvas to use the scrollbars
        self.table_canvas.configure(xscrollcommand=self.scrollbar_x_main.set, yscrollcommand=self.scrollbar_y.set)

        # Create an inner frame inside the canvas to hold the table widgets
        self.table_inner_frame = tk.Frame(self.table_canvas)
        self.table_canvas_frame = self.table_canvas.create_window((0,0), window=self.table_inner_frame, anchor="nw")

        # Bind events to update the scroll region
        self.table_inner_frame.bind("<Configure>", self.on_table_frame_configure)
        self.table_canvas.bind("<Configure>", self.on_table_canvas_configure)

        # **FASTQ Selection Section**
        self.files_frame = ttk.LabelFrame(self.root, text="FASTQ Selection")
        self.files_frame.pack(pady=10, padx=10, expand=True, fill=tk.BOTH)

        # Use a Canvas to enable scrolling
        self.files_canvas = tk.Canvas(self.files_frame)
        self.files_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Vertical scrollbar for the files section
        self.scrollbar_y_files = tk.Scrollbar(self.files_frame, orient=tk.VERTICAL, command=self.files_canvas.yview)
        self.scrollbar_y_files.pack(side=tk.RIGHT, fill=tk.Y)

        # Horizontal scrollbar for the files section
        self.scrollbar_x_files = tk.Scrollbar(self.files_frame, orient=tk.HORIZONTAL, command=self.files_canvas.xview)
        self.scrollbar_x_files.pack(side=tk.BOTTOM, fill=tk.X)

        # Configure the canvas to use the scrollbars
        self.files_canvas.configure(xscrollcommand=self.scrollbar_x_files.set, yscrollcommand=self.scrollbar_y_files.set)

        # Create an inner frame inside the canvas to hold the file selection widgets
        self.files_inner_frame = tk.Frame(self.files_canvas)
        self.files_canvas_frame = self.files_canvas.create_window((0,0), window=self.files_inner_frame, anchor="nw")

        # Bind events to update the scroll region
        self.files_inner_frame.bind("<Configure>", self.on_files_frame_configure)
        self.files_canvas.bind("<Configure>", self.on_files_canvas_configure)

        # **Buttons Section**
        self.button_frame = ttk.LabelFrame(self.root, text="Actions")
        self.button_frame.pack(pady=10, padx=10, fill="x")

        # Button to load the sample sheet
        self.load_button = tk.Button(self.button_frame, text="1: Load sample sheet", command=self.load_tsv)
        self.load_button.pack(pady=5)

        # Button to update the DataFrame with selections
        self.update_button = tk.Button(self.button_frame, text="Update DataFrame", command=self.update_dataframe)
        self.update_button.pack(pady=5)

        # Button to upload data to S3
        self.upload_button = tk.Button(self.button_frame, text="Upload to S3", command=self.upload_to_s3)
        self.upload_button.pack(pady=5)

        # Button to delete selected rows
        self.delete_button = tk.Button(self.button_frame, text="Delete Selected Rows", command=self.delete_selected_rows)
        self.delete_button.pack(pady=5)

    def on_table_frame_configure(self, event):
        # Update the scroll region when the size of the inner frame changes
        self.table_canvas.configure(scrollregion=self.table_canvas.bbox("all"))

    def on_table_canvas_configure(self, event):
        # Ensure the inner frame width matches the canvas width - still not working
        canvas_width = event.width
        self.table_inner_frame.config(width=canvas_width)
        self.table_canvas.config(scrollregion=self.table_canvas.bbox("all"))

    def on_files_frame_configure(self, event):
        # Update the scroll region when the size of the inner frame changes
        self.files_canvas.configure(scrollregion=self.files_canvas.bbox("all"))

    def on_files_canvas_configure(self, event):
        # Ensure the inner frame width matches the canvas width
        canvas_width = event.width
        self.files_inner_frame.config(width=canvas_width)
        self.files_canvas.config(scrollregion=self.files_canvas.bbox("all"))

    def load_tsv(self):
        # Open a file dialog to select the sample sheet TSV
        file_path = filedialog.askopenfilename(
            initialdir=self.directory_path,
            title="Select TSV file",
            filetypes=(("TSV files", "*.tsv"), ("all files", "*.*"))
        )
        if file_path:
            # Read the sample sheet TSV file into a DataFrame
            self.df = pd.read_csv(file_path, sep='\t')
            self.initialize_dataframe()  # Initialize any missing columns
            self.display_table()  # Display the data in the main table
            self.display_file_checks()  # Display the file selection section

    def initialize_dataframe(self):
        # Initialize the DataFrame with necessary columns
        for column in self.columns:
            if column.get('show', True):
                header = column['header']
                # Add missing columns to the DataFrame
                if header not in self.df.columns:
                    self.df[header] = ""

                if column['source'] == 'linked-tsv':
                    # For linked columns, map values from another column
                    linked_to = column['linked_to']
                    linked_values = column['linked_values']
                    self.df[column['header']] = self.df[linked_to].map(linked_values).fillna("")

    def display_table(self):
        # Clear existing widgets in the table inner frame
        for widget in self.table_inner_frame.winfo_children():
            widget.destroy()

        # Display column headers
        headers = ['Select'] + [column['header'] for column in self.columns if column.get('show', True)]
        for j, header in enumerate(headers):
            label = tk.Label(self.table_inner_frame, text=header, relief=tk.RIDGE)
            label.grid(row=0, column=j, sticky="nsew")

        # Display data rows
        self.selected_vars = []  # List to hold variables for row selection checkboxes
        for i in range(len(self.df)):
            # Create a variable for the checkbox
            selected_var = tk.BooleanVar()
            checkbox = tk.Checkbutton(self.table_inner_frame, variable=selected_var)
            checkbox.grid(row=i+1, column=0, sticky="nsew")
            self.selected_vars.append(selected_var)

            # Display data columns
            for j, column in enumerate([col for col in self.columns if col.get('show', True)]):
                header = column['header']
                # Determine the value based on the source
                if column['source'] == 'tsv':
                    value = self.df.iloc[i][header]
                elif column['source'] == 'dropdown':
                    value = self.dropdown_vars[header].get()
                elif column['source'] == 'linked-dropdown':
                    linked_to = self.linked_vars[header]
                    linked_value = self.dropdown_vars[linked_to].get()
                    value = self.linked_vars[linked_to].get(linked_value, "")
                elif column['source'] == 'generated':
                    value = self.df.iloc[i][header]
                elif column['source'] == 'linked-tsv':
                    value = self.df.iloc[i][header]

                # Create an entry widget to display the value
                entry = tk.Entry(self.table_inner_frame, width=20)
                entry.grid(row=i+1, column=j+1, sticky="nsew")
                entry.insert(tk.END, value)
                # Make non-editable if the source is not 'tsv'
                if column['source'] != 'tsv':
                    entry.configure(state='readonly')

        # Update the scroll region
        self.table_canvas.configure(scrollregion=self.table_canvas.bbox("all"))

    def display_file_checks(self):
        # Clear existing widgets in the files inner frame
        for widget in self.files_inner_frame.winfo_children():
            widget.destroy()

        # Create headers for the file selection table
        headers = ['Sample'] + self.timepoints + ['Selected', 'Upload Status']
        for j, header in enumerate(headers):
            label = tk.Label(self.files_inner_frame, text=header, relief=tk.RIDGE)
            label.grid(row=0, column=j, sticky="nsew")

        self.file_vars = []  # List to hold variables for file selection
        for i in range(len(self.df)):
            # Get the sample ID (using 'LabID' column)
            experiment_id = self.df.iloc[i]['LabID']
            label = tk.Label(self.files_inner_frame, text=experiment_id, relief=tk.RIDGE)
            label.grid(row=i+1, column=0, sticky="nsew")

            row_vars = []
            selected_var = tk.StringVar(value="")
            for j, timepoint in enumerate(self.timepoints):
                # Build the file path pattern to search for
                file_path = f"/mnt/results/{experiment_id}/{timepoint}/microbial/{experiment_id}_{timepoint}_hg38_removed.fastq.gz"
                print(f"Checking for file: {file_path}")  # Debug output
                files_present = glob.glob(file_path)
                print(f"Files found: {files_present}")  # Debug output

                # Create a radio button for each timepoint
                radio = tk.Radiobutton(self.files_inner_frame, variable=selected_var, value=file_path)
                radio.grid(row=i+1, column=j+1, sticky="nsew")
                if files_present:
                    # Enable radio button if file exists
                    radio.config(state='normal')
                else:
                    # Disable radio button if file does not exist
                    radio.config(state='disabled')
                row_vars.append(selected_var)  # Store the StringVar

            # Entry to display the selected file path
            selected_entry = tk.Entry(self.files_inner_frame, textvariable=selected_var, width=30)
            selected_entry.grid(row=i+1, column=len(self.timepoints)+1, sticky="nsew")
            row_vars.append(selected_var)

            # Label to display upload status
            status_var = tk.StringVar(value="Pending")
            status_label = tk.Label(self.files_inner_frame, textvariable=status_var)
            status_label.grid(row=i+1, column=len(self.timepoints)+2, sticky="nsew")
            row_vars.append(status_var)

            self.file_vars.append(row_vars)

        # Update the scroll region
        self.files_canvas.configure(scrollregion=self.files_canvas.bbox("all"))

    def update_dataframe(self):
        if not self.df.empty:
            self.submit_selection()  # Update DataFrame with selected files
            self.extract_fastq_info()  # Extract run_id from FASTQ files
        else:
            messagebox.showwarning("No Data", "Please load a sample sheet first.")

    def submit_selection(self):
        # Update the DataFrame with the selected file paths
        for i, row_vars in enumerate(self.file_vars):
            selected_file = row_vars[0].get()  # Get the selected file path
            if selected_file:
                self.df.at[i, 'SelectedFiles'] = selected_file
            else:
                messagebox.showwarning("No Selection", f"No timepoint selected for sample {self.df.at[i, 'LabID']}")

        self.display_table()  # Refresh the table display

    def extract_fastq_info(self):
        # Extract run_id from the selected FASTQ files
        for i, row_vars in enumerate(self.file_vars):
            selected_file = row_vars[0].get()
            if selected_file:
                run_id = self.get_run_id_from_fastq(selected_file)
                self.df.at[i, 'run_id'] = run_id
            else:
                messagebox.showwarning("No Selection", f"No timepoint selected for sample {self.df.at[i, 'LabID']}")

        self.display_table()  # Refresh the table display

    def get_run_id_from_fastq(self, file_path):
        try:
            # Read the first line of the FASTQ file
            if file_path.endswith('.gz'):
                with gzip.open(file_path, 'rt') as f:
                    first_line = f.readline().strip()
            else:
                with open(file_path, 'rt') as f:
                    first_line = f.readline().strip()
            # Extract the run_id from the first line
            run_id = first_line.split()[0].lstrip('@')
            return run_id
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            return ""

    def get_encryption_key(self, service_name, key_name):
        # Retrieve the encryption key from the keyring
        encoded_key = keyring.get_password(service_name, key_name)
        if encoded_key:
            return base64.b64decode(encoded_key)
        else:
            print("Encryption key not found in keyring")
            return None

    def encrypt_column_data(self):
        # Pre-encryption warning prompt
        proceed = messagebox.askyesno(
            "Encrypt Data",
            "Are you sure you want to encrypt the specified columns? This action cannot be undone."
        )
        if not proceed:
            return

        encryption_key = self.get_encryption_key("your_keyring_service_name", "your_keyring_key_name")
        if not encryption_key:
            return

        for column in self.encryption_columns:
            if column in self.df.columns:
                # Check if any string in the column is longer than 20 characters (suggesting it's unencrypted)
                if any(len(str(value)) > 20 for value in self.df[column]):
                    proceed = messagebox.askyesno(
                        "Has the table already been encrypted?",
                        f"The column '{column}' contains strings longer than 20 characters. "
                        "Do you still want to proceed with encryption? Encrypting twice corrupts data."
                    )
                    if not proceed:
                        return

                iv = self.initialization_vector  # Use the IV from the config file
                encrypted_data = []

                for value in self.df[column]:
                    if pd.isna(value):
                        encrypted_data.append(value)
                    else:
                        # Encrypt the value
                        cipher = AES.new(encryption_key, AES.MODE_CBC, iv)
                        encrypted_bytes = cipher.encrypt(pad(str(value).encode(), AES.block_size))
                        encrypted_value = base64.b64encode(iv + encrypted_bytes).decode('utf-8')
                        encrypted_data.append(encrypted_value)

                self.df[column] = encrypted_data  # Update the column with encrypted data

        self.display_table()  # Refresh the table display
        messagebox.showinfo("Encryption Complete", "The specified columns have been encrypted.")

    def update_dataframe_with_dropdown_values(self):
        # Update the DataFrame with the current values from the dropdown menus
        for column in self.columns:
            if column['source'] == 'dropdown' and column.get('show', True):
                header = column['header']
                self.df[header] = self.dropdown_vars[header].get()
            elif column['source'] == 'linked-dropdown' and column.get('show', True):
                header = column['header']
                linked_to = self.linked_vars[header]
                linked_value = self.dropdown_vars[linked_to].get()
                self.df[header] = self.linked_vars[linked_to].get(linked_value, "")

    def upload_to_s3(self):
        # Update DataFrame with the latest dropdown values
        self.update_dataframe_with_dropdown_values()

        # Create directory for temporary files
        mscape_dir = "/mnt/mscape"
        os.makedirs(mscape_dir, exist_ok=True)

        for i, row in self.df.iterrows():
            # Prepare filenames
            csv_filename = f"mscape.{row['Barcode']}.{row['run_id']}.csv"
            csv_path = os.path.join(mscape_dir, csv_filename)

            # Prepare the row for export
            row_for_export = row[[col['header'] for col in self.columns if col.get('show', True)]]
            row_for_export = row_for_export.rename(
                index={col['header']: col['export_header'] for col in self.columns if col.get('show', True)}
            )
            # Save the row as a CSV file
            row_for_export.to_frame().T.to_csv(csv_path, index=False)

            # Prepare the FASTQ file
            original_fastq_path = row['SelectedFiles']
            fastq_filename = f"mscape.{row['Barcode']}.{row['run_id']}.fastq.gz"
            fastq_path = os.path.join(mscape_dir, fastq_filename)

            # Get the status variable to update the upload status
            status_var = self.file_vars[i][-1]

            try:
                # Copy the FASTQ file to the temporary directory
                shutil.copy(original_fastq_path, fastq_path)

                # Upload CSV file to S3
                subprocess.run([
                    "aws", "--profile", self.aws_profile, "--endpoint", self.aws_endpoint,
                    "s3", "cp", csv_path, f"s3://{self.s3_bucket}/{csv_filename}"
                ], check=True)

                # Upload FASTQ file to S3
                subprocess.run([
                    "aws", "--profile", self.aws_profile, "--endpoint", self.aws_endpoint,
                    "s3", "cp", fastq_path, f"s3://{self.s3_bucket}/{fastq_filename}"
                ], check=True)

                # Remove the copied FASTQ file to save space
                os.remove(fastq_path)

                # Update status to "Complete"
                status_var.set("Complete")
            except subprocess.CalledProcessError as e:
                # Update status to "Failed" on error
                status_var.set("Failed")
                print(f"Error uploading sample {row['Barcode']}: {e}")
            except Exception as e:
                status_var.set("Failed")
                print(f"Error processing sample {row['Barcode']}: {e}")

            self.root.update_idletasks()  # Force update of the GUI

        # Show completion message
        messagebox.showinfo(
            "Upload Complete",
            "All uploads are complete. Please check the terminal outputs and CLIMB to validate uploads."
        )

    def delete_selected_rows(self):
        # Delete the selected rows from the DataFrame
        selected_indices = [i for i, var in enumerate(self.selected_vars) if var.get()]
        if selected_indices:
            self.df.drop(selected_indices, inplace=True)
            self.df.reset_index(drop=True, inplace=True)
            self.display_table()  # Refresh the table display
            self.display_file_checks()  # Refresh the file selection display
        else:
            messagebox.showwarning("No Selection", "No rows selected for deletion.")

if __name__ == "__main__":
    root = tk.Tk()  # Create the main window
    app = MetadataApp(root)  # Create an instance of the application
    root.mainloop()  # Start the Tkinter event loop
