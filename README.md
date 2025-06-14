# MISTs-Construction-using-MPICH-and-OpenMP

## 📽️ Presentation

The `Deliverable1/Presentation/` folder contains the final presentation materials for this project:

- **MISTs in Bubble-Sort (PDF).pdf** – A PDF version of the presentation slides for quick viewing.  
- **MISTs in Bubble-Sort (PPT).ppt** – The original editable PowerPoint file.

These slides provide a concise overview of the problem, solution, algorithm design, and parallelization strategy used in constructing MISTs for Bubble-Sort Networks.

---

## 📁 Deliverable1 Files

- **Extracted_keypoints.docx**  
  Contains summarized key points and insights extracted from the research paper, serving as a quick reference.

- **Research_Paper.pdf**  
  The original research paper that forms the basis of this project.

---

## 💻 Code Structure

### 🔹 Serial Implementation

**Folder:** `Code/Serial Implementation`

- **File:** `serial_new.cpp`  
  Contains the serial (single-threaded) implementation of the MIST construction algorithm.

- **Compile & Run:**
  ```bash
  g++ -std=c++17 -O2 serial_new.cpp -o mist_bubblesort
  ./serial.out
  ```

## 🔹 Parallel Implementation

**Folder:** `Code/Parallel Implementation`

- **File:** `parallel.cpp`  
  Contains the parallel implementation using both OpenMP (for shared-memory parallelism) and MPICH (for distributed-memory MPI).

### Compile (MPICH + OpenMP)
  ```bash
    mpicxx -fopenmp -O2 parallel.cpp -o parallel
    mpirun -np 2 ./parallel 9
  ```

## ⚙️ Dependencies & Pre-installed Libraries

Before building and running the parallel version, ensure your system has:

- **MPICH**  
  For MPI (Message Passing Interface) support.  
  Check installation with:
  ```bash
  mpicxx -version
  ```
- OpenMP  
Enabled in your C++ compiler (usually via the `-fopenmp` flag in `g++`/`mpicxx`).

- g++ / mpicxx  
C++ compilers that support OpenMP and MPI:  
  ```bash
  g++ --version
  mpicxx --version
  ```

---

## 🚀 Contributing

Contributions, issues and feature requests are welcome! Please take a look at the [Contributing Guidelines](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests.

---

## 📝 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

