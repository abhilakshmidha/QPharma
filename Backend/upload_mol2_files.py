import os
from tkinter import Tk, filedialog

def upload_mol2_files():
    """Open a dialog to select multiple MOL2 files and return their paths."""
    root = Tk()
    root.withdraw()  # Hide the main window
    file_paths = filedialog.askopenfilenames(
        title="Select MOL2 Files",
        filetypes=[("MOL2 files", "*.mol2")],
        initialdir=os.getcwd()
    )
    return list(file_paths)

