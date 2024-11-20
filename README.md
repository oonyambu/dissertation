# Dissertation Instructions

**A dissertation submitted in partial satisfaction of the requirements for the degree of Doctor of Philosophy in Statistics**  
by **Samuel Onyambu**

---

## Compilation Instructions

To compile the dissertation to a PDF:

1. **Download the folder containing the dissertation files.**
2. **Using a terminal**, navigate to the dissertation folder.
3. Run one of the following commands:

   - **Option A**:  
     ```bash
     make
     ```
   - **Option B**:  
     ```bash
     make BACKEND=biber
     ```
   - **Option C**:  
     ```bash
     make BACKEND=bibtex ENGINE=xelatex
     ```

   ### Notes:
   - The **default backend** is `bibtex`.  
     If you have `biber` installed, you can use `BACKEND=biber` for bibliography management.
   - For the LaTeX engine, you may use any of the following:
     - `xelatex`
     - `pdflatex`
     - `lualatex`

---

## Prerequisites

Ensure the following are properly configured and accessible via your system's PATH:

1. **Rscript**:  
   Verify its availability by running:
   ```bash
   Rscript --version
   ```
   This should return the version of R installed on your system.

2. **PDF Engine:**
  Ensure your selected LaTeX engine is installed and accessible. For example:
  ```bash
  pdflatex --version
  ```
  or
  ```bash
  xelatex --version
  ```
  This should return the version number of the engine.

By following these instructions, you should be able to successfully compile the dissertation into a PDF.
